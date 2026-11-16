#!/usr/bin/env Rscript
# No-demultiplexing module (passthrough for non-multiplexed datasets)
#
# For datasets where each pool holds a single sample: there is nothing to
# demultiplex. Every called cell is assigned to that one sample and no cell is
# flagged as multiplet. The output has the same core columns as the other
# demultiplexing modules (pool_id, cell_id, sample_id, multiplet), so
# concatenate_dataset can treat all datasets alike.
#
# Note: multiplet = FALSE here means "not detectable" for now, not "known singlet".
# Do not use this dataset as multiplet ground truth in doublet metrics.
# Possible later: derive singlet ground truth by consensus, i.e. run several
# doublet detection methods and treat cells that all methods call singlet as
# singlets (an approach seen used by another group).

suppressPackageStartupMessages({
  library(anndataR)
  library(SingleCellExperiment)
  library(data.table)
})

# arg parsing
source("src/common/cli.R")
p <- arg_parser("No-demultiplexing module (non-multiplexed datasets)")
p <- add_base_args(p)                    # --output_dir, --name
p <- add_stage_args(p, "demultiplexing") # --rawdata_h5ad, --called_cells_tsv
args <- parse_args(p, argv = strip_output_args())

cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))

assign_sample <- function(sce) {
  if (!"pool_id" %in% colnames(colData(sce))) {
    stop("colData(sce) is missing 'pool_id'")
  }
  pool_id <- colData(sce)$pool_id

  # TODO(sample_id): source still to be agreed with Marusa.
  #   Option A: one-row pool_lookup per pool in metadata(sce)$pool_lookup
  #             (as for dm_hto / dm_flex) -> take its sample_id.
  #   Option B: no pool_lookup -> sample_id = pool_id (one sample per pool).
  # Provisional: option B.
  sample_id <- pool_id

  data.table(
    pool_id = pool_id,
    cell_id = colnames(sce),
    sample_id = sample_id,
    multiplet = FALSE
  )
}

main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  cat("  reading raw h5ad + called cell barcodes ..\n")
  sce <- read_h5ad(args$rawdata_h5ad, as = "SingleCellExperiment", x_mapping = "counts")
  called <- fread(args$called_cells_tsv)$cell_id
  sce <- sce[, colnames(sce) %in% called]

  if (ncol(sce) != length(unique(called))) {
    warning(sprintf("%d called cell(s) not found in the raw h5ad",
                    length(unique(called)) - ncol(sce)))
  }

  cat("  no multiplexing: assigning every called cell to its pool's single sample ..\n")
  result <- assign_sample(sce)
  cat(sprintf("  %d called cells, %d sample(s), 0 multiplets (not detectable)\n",
              nrow(result), uniqueN(result$sample_id)))

  out_path <- file.path(args$output_dir, paste0(args$name, "_demultiplexed.tsv"))
  fwrite(result, out_path, sep = "\t", quote = FALSE)
  cat(sprintf("  wrote: %s\n", out_path))
}

if (sys.nframe() == 0L) {
  main()
}
