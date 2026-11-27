#!/usr/bin/env Rscript
# Flex demultiplexing module
#
# multiplet detection is GEM-barcode based: if the same
# physical GEM was called as a cell for more than one probe_id every cell_id sharing
# that GEM barcode is flagged multiplet

suppressPackageStartupMessages({
  library(anndataR)
  library(SingleCellExperiment)
  library(data.table)
  library(magrittr)
})

# arg parsing
source("src/common/cli.R")
p <- arg_parser("Flex demultiplexing module (flex)")
p <- add_base_args(p)                    # --output_dir, --name
p <- add_stage_args(p, "demultiplexing") # --rawdata_h5ad, --called_cells_tsv
args <- parse_args(p, argv = strip_output_args())

cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))

assign_sample <- function(sce, pool_lookup) {
  if (!"probe_id" %in% colnames(colData(sce))) {
    stop("colData(sce) is missing 'probe_id'")
  }

  out_dt <- data.table(
    pool_id = colData(sce)$pool_id,
    cell_id = colnames(sce),
    probe_id = colData(sce)$probe_id
  ) %>%
  .[, gem_id := sub("_[^_]+$", "", cell_id)] %>%
  merge(pool_lookup, by = c("probe_id", "pool_id"), all.x = TRUE) %>%
  .[, multiplet := .N > 1, by = gem_id] %>%
  .[, c("probe_id", "gem_id") := NULL]

  return(out_dt)
}

main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  cat("  reading raw h5ad + called cell barcodes ..\n")
  sce <- read_h5ad(args$rawdata_h5ad, as = "SingleCellExperiment", x_mapping = "counts")
  called <- fread(args$called_cells_tsv)$cell_id
  sce <- sce[, colnames(sce) %in% called]

  pl_raw <- metadata(sce)$pool_lookup
  if (is.null(pl_raw)) {
    stop("metadata(sce) is missing 'pool_lookup' -- check the anndataR metadata/uns round trip")
  }
  pool_lookup <- as.data.table(pl_raw)

  cat("  assigning samples from probe_id, flagging multi-probe GEM droplets as multiplets ..\n")
  result <- assign_sample(sce, pool_lookup)
  cat(sprintf("  %d / %d called cells flagged as multiplets\n", sum(result$multiplet), nrow(result)))

  out_path <- file.path(args$output_dir, paste0(args$name, "_demultiplexed.tsv"))
  fwrite(result, out_path, sep = "\t", quote = FALSE)
  cat(sprintf("  wrote: %s\n", out_path))
}

if (sys.nframe() == 0L) {
  main()
}
