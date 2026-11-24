#!/usr/bin/env Rscript
# Concatenate pools, join demultiplexing labels (gather module).
#
# Reads every pool's cleaned h5ad (filtered_h5ad) + demultiplexing TSV
# (demultiplexed_tsv) for one dataset and combines them into a single h5ad:
# rewrites each pool's barcodes to "{pool_id}_{barcode}", cbinds the pools together and joins whichever
# demultiplexing label columns exist 

suppressPackageStartupMessages({
  library(anndataR)
  library(SingleCellExperiment)
  library(data.table)
})

source("src/common/cli.R")

extract_multi_flag <- function(argv, flag) {
  idx <- which(argv == flag)
  if (length(idx) == 0) stop(sprintf("missing required flag %s", flag))
  start <- idx[1] + 1
  end <- start
  while (end <= length(argv) && !startsWith(argv[end], "--")) end <- end + 1
  end <- end - 1
  if (end < start) stop(sprintf("%s requires at least one value", flag))
  list(values = argv[start:end], remaining = argv[-(idx[1]:end)])
}

argv <- strip_output_args()
fh <- extract_multi_flag(argv, "--filtered_h5ad")
dmx <- extract_multi_flag(fh$remaining, "--demultiplexed_tsv")
filtered_h5ad <- fh$values
demultiplexed_tsv <- dmx$values

p <- arg_parser("Concatenate pools, join demultiplexing labels (gather module)")
p <- add_base_args(p)  # --output_dir, --name
args <- parse_args(p, argv = dmx$remaining)

cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("  output_dir: %s\n", args$output_dir))
cat(sprintf("  name: %s\n", args$name))
cat(sprintf("  filtered_h5ad (%d pools): %s\n", length(filtered_h5ad), paste(filtered_h5ad, collapse = ", ")))
cat(sprintf("  demultiplexed_tsv (%d pools): %s\n", length(demultiplexed_tsv), paste(demultiplexed_tsv, collapse = ", ")))

# One label row per GEM barcode, for cells whose ids are GEM-level. A GEM
# called under one probe keeps that probe's labels. A GEM called under several
# is a multiplet (as demultiplex_flex flags it) with no single sample: its
# other label columns become NA.
labels_per_gem <- function(labels) {
  labels <- copy(labels)
  labels[, cell_key := sub("_[^_]+$", "", cell_key)]
  labels[, n_probes := .N, by = cell_key]
  label_cols <- setdiff(colnames(labels), c("pool_id", "cell_id", "cell_key", "n_probes", "multiplet"))
  for (col in label_cols) {
    set(labels, which(labels$n_probes > 1), col, NA)
  }
  if ("multiplet" %in% colnames(labels)) labels[n_probes > 1, multiplet := TRUE]
  labels <- unique(labels, by = "cell_key")
  labels[, cell_id := sub("_[^_]+$", "", cell_id)]
  labels[, n_probes := NULL]
  labels[]
}

main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  cat("  reading + concatenating per-pool h5ads ..\n")
  sces <- lapply(filtered_h5ad, function(path) read_h5ad(path, as = "SingleCellExperiment", x_mapping = "counts"))

  gene_sets <- lapply(sces, rownames)
  if (!all(sapply(gene_sets, identical, gene_sets[[1]]))) {
    stop("per-pool h5ads have different gene sets -- cannot cbind (expected identical features across pools of one dataset)")
  }

  # 10x barcodes aren't globally unique across pools, so rewrite colnames to
  # "{pool_id}_{barcode}" up front
  sces <- lapply(sces, function(sce) {
    if (!"pool_id" %in% colnames(colData(sce))) stop("an input h5ad is missing colData 'pool_id'")
    # anndataR reads X with its own dimnames; left in place they go stale on
    # the rename below and write_h5ad rejects the object.
    for (a in assayNames(sce)) dimnames(assay(sce, a, withDimnames = FALSE)) <- NULL
    colnames(sce) <- paste0(colData(sce)$pool_id, "_", colnames(sce))
    sce
  })
  combined <- do.call(cbind, sces)
  cat(sprintf("  combined: %d cells x %d genes, %d pools\n",
              ncol(combined), nrow(combined), length(unique(colData(combined)$pool_id))))

  cat("  reading + concatenating per-pool demultiplexing assignments ..\n")
  labels <- rbindlist(lapply(demultiplexed_tsv, fread))
  labels[, cell_key := paste0(pool_id, "_", cell_id)]

  missing <- setdiff(colnames(combined), labels$cell_key)
  if (length(missing) > 0) {
    # GEM-merged Flex counts (PREP `merge: gems`): the cells lost their
    # "_{probe_id}" suffix, the labels still carry it.
    cat("  cells do not match the labels by id; matching by GEM barcode ..\n")
    labels <- labels_per_gem(labels)
    missing <- setdiff(colnames(combined), labels$cell_key)
  }
  if (length(missing) > 0) {
    stop(sprintf("%d cell(s) in the combined matrix have no demultiplexing assignment", length(missing)))
  }

  # reorder labels to match combined's column order before attaching
  labels <- labels[match(colnames(combined), labels$cell_key)]
  label_cols <- setdiff(colnames(labels), c("pool_id", "cell_id", "cell_key"))
  for (col in label_cols) {
    colData(combined)[[col]] <- labels[[col]]
  }

  out_path <- file.path(args$output_dir, paste0(args$name, ".h5ad"))
  write_h5ad(combined, out_path, x_mapping = "counts")
  cat(sprintf("  wrote: %s\n", out_path))
}

if (sys.nframe() == 0L) {
  main()
}
