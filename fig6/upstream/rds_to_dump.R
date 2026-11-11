#!/usr/bin/env Rscript
# Dump a Seurat RDS into MatrixMarket + metadata tables, for build_h5ad.py to assemble.
# Usage: Rscript rds_to_dump.R <input.RDS> <dump_dir>
#
# Figure 6 input #1: 3k_dsci_oligo_subset.RDS -> dump -> 3k_dsci_oligo_subset.h5ad
# (the expression channel of the own dsciCAT-seq nuclei; panels A, B, C and the trait-level control).
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
rds_file <- args[1]
dump_dir <- args[2]
if (is.na(rds_file) || is.na(dump_dir)) stop("usage: Rscript rds_to_dump.R <input.RDS> <dump_dir>")

dir.create(dump_dir, showWarnings = FALSE, recursive = TRUE)
obj <- readRDS(rds_file)
cat("[R] object:", class(obj), "| genes x cells:", dim(obj), "\n")
stopifnot(inherits(obj, "Seurat"))

# --- RNA matrix (dgCMatrix, genes x cells) ---
# counts slot = raw / SoupX-corrected counts (non-integer after correction);
# data slot = log1p(counts / column sum * 1e4)  <-- scDRS is called with --flag-raw-count False,
# so it is the data slot that must be exported.
c <- GetAssayData(obj, assay = "RNA", slot = "counts")
d <- GetAssayData(obj, assay = "RNA", slot = "data")
if (is.null(c)) c <- GetAssayData(obj, assay = DefaultAssay(obj), slot = "counts")
if (is.null(d)) d <- c
cat("[R] counts slot max:", max(c), "| data slot max:", max(d), "\n")
writeMM(c, file.path(dump_dir, "matrix.mtx"))
writeMM(d, file.path(dump_dir, "matrix_data.mtx"))
write.table(rownames(c), file.path(dump_dir, "genes.tsv"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(colnames(c), file.path(dump_dir, "barcodes.tsv"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# --- metadata (the columns the figure notebooks need) ---
md <- obj@meta.data
keep <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mito",
          "percent.ribo", "doublet_info", "cell_type")
keep <- intersect(keep, colnames(md))
md2 <- md[, keep, drop = FALSE]
write.table(md2, file.path(dump_dir, "meta.tsv"), sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = TRUE)
cat("[R] metadata columns:", colnames(md2), "\n")

# --- RNA UMAP, if present (panel B draws it; the embedding is never recomputed downstream) ---
umap_key <- NULL
for (r in c("umap", "wnn.umap", "umap.atac")) {
  if (!is.null(obj@reductions[[r]])) { umap_key <- r; break }
}
if (!is.null(umap_key)) {
  um <- as.data.frame(obj@reductions[[umap_key]]@cell.embeddings[, 1:2, drop = FALSE])
  colnames(um) <- c("UMAP_1", "UMAP_2")
  um$barcode <- rownames(um)
  write.table(um, file.path(dump_dir, "umap.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  cat("[R] UMAP taken from reduction:", umap_key, "\n")
} else {
  cat("[R] no UMAP coordinates in this object\n")
}
cat("[R] dumped ->", dump_dir, "\n")
