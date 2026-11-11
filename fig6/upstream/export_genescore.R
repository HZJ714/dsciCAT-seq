#!/usr/bin/env Rscript
# Export the chromatin (ATAC "gene score") channel of the own dsciCAT-seq oligodendrocyte-lineage
# nuclei, as MatrixMarket + metadata, for build_genescore_h5ad.py to assemble into
# genescore_oligo_5977.h5ad — the matrix that panel C's chromatin side is scored from.
#
#   input RDS   3k_dsci_ht_10x_gs_subset.RDS  (Signac object, assay `gene_score`)
#               -> env DSCI_GS_RDS, defaults to $DSCI_FY_ROOT/../3k_dsci_ht_10x_gs_subset.RDS
#   nucleus set lineage_cells_5977.txt        (= obs_names of 3k_dsci_oligo_subset.h5ad, the Fig 2
#               lineage curation; regenerate with the driver script D0_own_lineage_h5ad.sh)
#   matrix      gene_score@data = log1p(ArchR gene score).  Values are on the LOG scale
#               (min 0.0069, max 6.4676), which is what scDRS requires with
#               --flag-raw-count False.  NOTE that this is a different scale from the panel-D
#               cross-chemistry gene-score matrices, which are the un-logged ArchR matrices
#               normalised to 1e4 per nucleus (max 1388.6 / 297.1 / 372.7).  See README.md.
suppressMessages(library(Seurat)); suppressMessages(library(Signac)); suppressMessages(library(Matrix))

D    <- Sys.getenv("DSCI_GS_PROBE_DIR", file.path(Sys.getenv("DSCI_WORK_ROOT", "."),
                   "FY/oligo/fig6_rerun_20261019/dualomics_probe"))
RDS  <- Sys.getenv("DSCI_GS_RDS",
         file.path(Sys.getenv("DSCI_ONLINE_ROOT", "."), "FY/3k_dsci_ht_10x_gs_subset.RDS"))
cells <- readLines(file.path(D, "lineage_cells_5977.txt"))
a <- readRDS(RDS)
idx <- match(cells, colnames(a)); cat("matched", sum(!is.na(idx)), "/", length(cells), "\n")
stopifnot(!any(is.na(idx)))

gs <- a@assays$gene_score
M  <- tryCatch(gs@data, error=function(e) gs[["data"]])
sub <- M[, idx, drop=FALSE]
cat("gene_score subset dim:", paste(dim(sub), collapse=" x "), "\n")
v <- as.vector(sub@x)
cat(sprintf("nonzeros=%d  density=%.4f%%\n", length(v), 100*length(v)/prod(dim(sub))))
cat(sprintf("range min=%.4f max=%.4f | fraction negative=%.4f  (log1p scale => max < 8)\n",
            min(v), max(v), mean(v<0)))
print(quantile(v, c(0,.01,.25,.5,.75,.99,1)))

tmp <- file.path(D, "genescore_5977.mtx")
writeMM(sub, tmp)                       # genes x nuclei, as MatrixMarket always is
system(paste("gzip -f", tmp))
writeLines(rownames(sub), file.path(D, "genescore_features.txt"))
writeLines(cells,         file.path(D, "genescore_barcodes.txt"))
md <- a@meta.data[idx, , drop=FALSE]
write.table(md, gzfile(file.path(D, "genescore_obs.tsv.gz")), sep="\t", quote=FALSE)
cat("obs columns:", paste(colnames(md), collapse=", "), "\n")
cat("Wrote", tmp, ".gz +", file.path(D, "genescore_obs.tsv.gz"), "-> run build_genescore_h5ad.py next\n")
