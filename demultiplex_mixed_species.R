#!/usr/bin/env Rscript
# Species demultiplexing module.
#
# Reads the raw h5ad + called cell barcodes (post cell_calling)
# and assigns each called cell to a species or declares it a multiplet
#
# classification follows cellranger's barnyard algorithm: take the Nth
# percentile of mouse counts among barcodes where mouse > human (the mouse
# count threshold), and the Nth percentile of human counts among barcodes
# where human > mouse (the human count threshold); a barcode exceeding BOTH
# thresholds is a multiplet (https://kb.10xgenomics.com/s/article/115003517183-How-does-cellranger-count-identify-multiplets)

suppressPackageStartupMessages({
  library(anndataR)
  library(SingleCellExperiment)
  library(Matrix)
  library(data.table)
})

# arg parsing
source("src/common/cli.R")
p <- arg_parser("Species demultiplexing module (mixed_species)")
p <- add_base_args(p)                     # --output_dir, --name
p <- add_stage_args(p, "demultiplexing")  # --rawdata_h5ad, --called_cells_tsv
p <- add_argument(p, "--percentile", type = "numeric", default = 10,
  help = "Percentile of same-species dominant counts used for the mouse/human count thresholds (cellranger default: 10)")
args <- parse_args(p, argv = strip_output_args())

cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))

assign_species <- function(sce, percentile = 10) {
  if (!"species" %in% colnames(rowData(sce))) {
    stop("rowData(sce) is missing 'species")
  }
  is_human <- rowData(sce)$species == "human"
  is_mouse <- rowData(sce)$species == "mouse"

  human_counts <- colSums(counts(sce)[is_human, , drop = FALSE])
  mouse_counts <- colSums(counts(sce)[is_mouse, , drop = FALSE])

  mouse_thresh <- quantile(mouse_counts[mouse_counts > human_counts],
    probs = percentile / 100, na.rm = TRUE)
  human_thresh <- quantile(human_counts[human_counts > mouse_counts],
    probs = percentile / 100, na.rm = TRUE)

  species <- ifelse(human_counts > mouse_counts, "human", "mouse")
  multiplet <- (mouse_counts > mouse_thresh) & (human_counts > human_thresh)

  if (!"pool_id" %in% colnames(colData(sce))) {
    stop("colData(sce) is missing 'pool_id")
  }

  pool_id <- colData(sce)$pool_id
  sample_id <- ifelse(multiplet, NA_character_, paste0(pool_id, "_", species))

  data.table(
    pool_id = pool_id,
    cell_id = colnames(sce),
    species = species,
    sample_id = sample_id,
    multiplet = multiplet
  )
}

main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  cat("  reading raw h5ad + called cell barcodes ..\n")
  sce <- read_h5ad(args$rawdata_h5ad, as = "SingleCellExperiment", x_mapping = "counts")
  called <- fread(args$called_cells_tsv)$cell_id
  sce <- sce[, colnames(sce) %in% called]

  cat("  assigning species (cellranger barnyard classification) ..\n")
  result <- assign_species(sce, percentile = args$percentile)
  counts_tbl <- result[, .N, by = .(species, multiplet)]
  cat(sprintf("  species counts: %s\n",
   paste(sprintf("%s%s=%d", counts_tbl$species,
     ifelse(counts_tbl$multiplet, " (multiplet)", ""), counts_tbl$N),
        collapse = ", ")))

  out_path <- file.path(args$output_dir, paste0(args$name, "_demultiplexed.tsv"))
  fwrite(result, out_path, sep = "\t", quote = FALSE)
  cat(sprintf("  wrote: %s\n", out_path))
}

if (sys.nframe() == 0L) {
  main()
}
