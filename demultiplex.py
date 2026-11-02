#!/usr/bin/env python3
"""Species demultiplexing module (mixed_species mode).

Reads the ambient-corrected h5ad (called cells only, post cell_calling +
cleaning or cell_calling_and_cleaning) and assigns each cell to human/mouse/
ambiguous based on the ratio of human- to mouse-prefixed gene counts *in the
cleaned matrix* -- i.e. species assignment runs after ambient removal, not on
raw counts.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent / "src"))
from common import cli  # noqa: E402


def parse_args():
  import argparse
  p = argparse.ArgumentParser(description="Species demultiplexing module (mixed_species)")
  cli.add_base_args(p)                   # --output_dir, --name
  cli.add_stage_args(p, "demultiplexing")  # --filtered_h5ad
  p.add_argument(
    "--quantile", type=float, default=5,
    help="Percentile of majority_ratio below which cells are called ambiguous (default: 5)",
  )
  return p.parse_args()


def assign_species(adata, quantile=5):
  genomes = adata.var_names.str.extract(r"^([^_]+)_")[0]
  unique_genomes = genomes.dropna().unique()

  human_genome = next(g for g in unique_genomes if g.startswith(("hg", "GRCh")))
  mouse_genome = next(g for g in unique_genomes if g.startswith(("mm", "GRCm")))

  is_human = (genomes == human_genome).to_numpy()
  is_mouse = (genomes == mouse_genome).to_numpy()
  human_counts = np.asarray(adata[:, is_human].X.sum(axis=1)).flatten()
  mouse_counts = np.asarray(adata[:, is_mouse].X.sum(axis=1)).flatten()
  total = human_counts + mouse_counts

  human_ratio = human_counts / total
  mouse_ratio = mouse_counts / total

  species_raw = np.where(human_ratio > mouse_ratio, "human", "mouse")
  majority_ratio = np.maximum(human_ratio, mouse_ratio)

  thresh = np.percentile(majority_ratio, quantile)
  species = np.where(majority_ratio > thresh, species_raw, "ambiguous")

  if "pool_id" not in adata.obs:
    raise ValueError("adata.obs is missing 'pool_id' -- can't be joined back after gather/concatenation")

  return pd.DataFrame({
    "pool_id": adata.obs["pool_id"].to_numpy(),
    "cell_id": adata.obs_names,
    "species": species,
    "majority_ratio": majority_ratio,
  })


def main():
  args = parse_args()
  print(f"Full command: {' '.join(sys.argv)}")
  for k in ("output_dir", "name", "filtered_h5ad", "quantile"):
    print(f"  {k}: {getattr(args, k)}")

  out_dir = Path(args.output_dir)
  out_dir.mkdir(parents=True, exist_ok=True)

  print("  reading cleaned h5ad ..")
  adata = sc.read_h5ad(args.filtered_h5ad)

  print("  assigning species (from cleaned counts) ..")
  result = assign_species(adata, quantile=args.quantile)
  counts = result["species"].value_counts()
  print(f"  species counts: {', '.join(f'{k}={v}' for k, v in counts.items())}")

  out_path = out_dir / f"{args.name}_demultiplexed.tsv"
  result.to_csv(out_path, sep="\t", index=False)
  print(f"  wrote: {out_path}")


if __name__ == "__main__":
  main()
