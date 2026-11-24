#!/usr/bin/env Rscript
# concatenate.R on GEM-merged Flex counts. Run from the repo root:
#   Rscript tests/test_concatenate.R
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(anndataR)
  library(data.table)
})
tmp <- tempfile("concat"); dir.create(tmp)
pool <- function(id, cells) {
  m <- matrix(1, nrow = 2, ncol = length(cells), dimnames = list(c("g1", "g2"), cells))
  sce <- SingleCellExperiment(assays = list(counts = as(m, "CsparseMatrix")),
                              colData = DataFrame(pool_id = rep(id, length(cells)), row.names = cells))
  path <- tempfile(id, tmp, ".h5ad")
  write_h5ad(sce, path, x_mapping = "counts")
  path
}
labels <- function(id, cell_id, sample_id, multiplet) {
  path <- file.path(tmp, paste0(id, ".tsv"))
  fwrite(data.table(pool_id = id, cell_id = cell_id, sample_id = sample_id, multiplet = multiplet),
         path, sep = "\t")
  path
}
run <- function(name, h5ads, tsvs) {
  status <- system2("Rscript", c("concatenate.R", "--output_dir", tmp, "--name", name,
                                 "--filtered_h5ad", h5ads, "--demultiplexed_tsv", tsvs),
                    stdout = FALSE, stderr = FALSE)
  if (status != 0) stop("concatenate.R failed")
  colData(read_h5ad(file.path(tmp, paste0(name, ".h5ad")), as = "SingleCellExperiment", x_mapping = "counts"))
}

# DEMUX saw GEM x probe cells: G1 under two probes (a multiplet), G2 under one.
tsv <- labels("p1", c("G1_BC001", "G1_BC002", "G2_BC001"), c("s1", "s2", "s1"), c(TRUE, TRUE, FALSE))

# cells as DEMUX saw them: matched by id, labels untouched
cd <- run("as_is", pool("p1", c("G1_BC001", "G1_BC002", "G2_BC001")), tsv)
stopifnot(identical(as.character(cd$sample_id), c("s1", "s2", "s1")),
          identical(as.logical(cd$multiplet), c(TRUE, TRUE, FALSE)))

# GEM-merged cells: one row per GEM; the multiplet has no sample
cd <- run("gems", pool("p1", c("G1", "G2")), tsv)
stopifnot(identical(rownames(cd), c("p1_G1", "p1_G2")),
          is.na(cd$sample_id[1]), as.character(cd$sample_id[2]) == "s1",
          identical(as.logical(cd$multiplet), c(TRUE, FALSE)))
cat("concatenate.R: all checks passed\n")
