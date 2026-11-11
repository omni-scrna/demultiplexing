#!/usr/bin/env Rscript
# HTO demultiplexing module
#
# Reads the raw h5ad + called cell barcodes and assigns each called cell to a sample
# using Seurat's HTODemux(). hash.ID is mapped to sample_id via metadata(sce)$pool_lookup (hto_id -> sample_id,
# written by ambient-data-download/run.R). Doublets and
# Negatives (no confident single HTO call) both get NA for sample_id.

suppressPackageStartupMessages({
  library(anndataR)
  library(SingleCellExperiment)
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(magrittr)
})

# arg parsing
source("src/common/cli.R")
p <- arg_parser("HTO demultiplexing module (hto_multiplexed)")
p <- add_base_args(p)                    # --output_dir, --name
p <- add_stage_args(p, "demultiplexing")  # --rawdata_h5ad, --called_cells_tsv
p <- add_argument(p, "--positive_quantile", type = "numeric", default = 0.99,
  help = "Seurat HTODemux() positive.quantile (per-HTO background cutoff)")
args <- parse_args(p, argv = strip_output_args())

cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))

assign_sample <- function(sce, pool_lookup, positive_quantile = 0.99) {
  is_hto <- rowData(sce)$feature_type == "Antibody Capture"
  hto_counts <- as.matrix(counts(sce)[is_hto, , drop = FALSE])
  original_hto_ids <- rownames(hto_counts)

  hto_seu <- CreateSeuratObject(counts = hto_counts, assay = "HTO")
  hto_seu <- NormalizeData(hto_seu, assay = "HTO", normalization.method = "CLR")
  hto_seu <- HTODemux(hto_seu, assay = "HTO", positive.quantile = positive_quantile)

  seu_hto_ids <- rownames(hto_seu[["HTO"]])
  seu_dmux_hto_id <- as.character(hto_seu$hash.ID) # hto assigned to every cell
  
  out_dt = data.table(
    pool_id = colData(sce)$pool_id,
    cell_id = colnames(sce),
    hto_id = original_hto_ids[match(seu_dmux_hto_id, seu_hto_ids)],
    class = tolower(as.character(hto_seu$HTO_classification.global))
  ) %>%
    .[, multiplet := class == 'doublet'] %>%
    merge(pool_lookup, c('hto_id', 'pool_id'), all.x = TRUE)
  
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

  cat("  demultiplexing HTOs (Seurat CLR + HTODemux, pre ambient removal) ..\n")
  result <- assign_sample(sce, pool_lookup, positive_quantile = args$positive_quantile)
  counts_tbl <- result[, .N, by = class]
  cat(sprintf("  classification counts: %s\n",
              paste(sprintf("%s=%d", counts_tbl$class, counts_tbl$N), collapse = ", ")))

  out_path <- file.path(args$output_dir, paste0(args$name, "_demultiplexed.tsv"))
  fwrite(result, out_path, sep = "\t", quote = FALSE)
  cat(sprintf("  wrote: %s\n", out_path))
}

if (sys.nframe() == 0L) {
  main()
}
