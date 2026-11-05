#!/usr/bin/env python3
"""Gather module: concatenate per-pool cleaned matrices into h5ad per
dataset, joining each cell's species label from demultiplexing.
"""

import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent / "src"))
from common import cli  # noqa: E402


def parse_args():
  import argparse
  # allow_abbrev=False: otherwise argparse treats the unrecognized `--output`
  # (api_version >= 0.8.0's named-outputs flag) as an abbreviation of the
  # declared `--output_dir` and silently overwrites it.
  p = argparse.ArgumentParser(description="Concatenate pools, join species labels (gather module)",
                               allow_abbrev=False)
  cli.add_base_args(p)  # --output_dir, --name
  # Gather flags: one value per upstream producer, not a single path -- the
  # shared schema helper doesn't support nargs, so these are hand-written.
  p.add_argument("--filtered_h5ad", nargs="+", required=True,
                  help="Cleaned h5ad, one per pool")
  p.add_argument("--demultiplexed_tsv", nargs="+", required=True,
                  help="Species assignment TSV (pool_id, cell_id, species, majority_ratio), one per pool")
  # EXPERIMENTAL (api_version >= 0.8.0): the engine also passes
  # `--output <id>=<path>` unconditionally; parse_known_args() tolerates it
  # rather than erroring on an unrecognized flag.
  return p.parse_known_args()[0]


def main():
  args = parse_args()
  print(f"Full command: {' '.join(sys.argv)}")
  print(f"  output_dir: {args.output_dir}")
  print(f"  name: {args.name}")
  print(f"  filtered_h5ad ({len(args.filtered_h5ad)} pools): {args.filtered_h5ad}")
  print(f"  demultiplexed_tsv ({len(args.demultiplexed_tsv)} pools): {args.demultiplexed_tsv}")

  out_dir = Path(args.output_dir)
  out_dir.mkdir(parents=True, exist_ok=True)

  # 10x barcodes aren't globally unique across pools, so rewrite obs_names to
  # "{pool_id}_{barcode}" up front 
  print("  reading + concatenating per-pool h5ads ..")
  adatas = [sc.read_h5ad(p) for p in args.filtered_h5ad]
  for a in adatas:
    if "pool_id" not in a.obs:
      raise ValueError("an input h5ad is missing obs['pool_id']")
    a.obs_names = a.obs["pool_id"].astype(str) + "_" + a.obs_names
  combined = ad.concat(adatas, join="outer")
  print(f"  combined: {combined.n_obs} cells x {combined.n_vars} genes, "
        f"{combined.obs['pool_id'].nunique()} pools")

  print("  reading + concatenating per-pool species assignments ..")
  labels = pd.concat([pd.read_csv(p, sep="\t") for p in args.demultiplexed_tsv], ignore_index=True)
  labels["cell_key"] = labels["pool_id"].astype(str) + "_" + labels["cell_id"]
  labels = labels.set_index("cell_key")

  missing = combined.obs_names.difference(labels.index)
  if len(missing) > 0:
    raise ValueError(f"{len(missing)} cell(s) in the combined matrix have no species assignment")

  joined = labels.loc[combined.obs_names]
  combined.obs["species"] = joined["species"].to_numpy()
  combined.obs["majority_ratio"] = joined["majority_ratio"].to_numpy()

  out_path = out_dir / f"{args.name}.h5ad"
  combined.write_h5ad(out_path)
  print(f"  wrote: {out_path}")


if __name__ == "__main__":
  main()
