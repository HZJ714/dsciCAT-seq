# -*- coding: utf-8 -*-
# D1: export the GeneScoreMatrix from Fig 2's ArchR project (all three chemistries share the same score pipeline)
# The recipe is copied verbatim from the (d) block of GitHub/fig2/fig2_out/fig2_up4_archr.R:
#   getMatrixFromProject(useMatrix='GeneScoreMatrix') -> rownames=geneDF$name
#   strip the "sample#" prefix from colnames -> one-to-one with the cell names of Fig 2's light meta
# keep only the OL-lineage nuclei (Oligo / OPC / COP_NFOL) and write the mtx for scDRS scoring.
suppressMessages(library(ArchR)); suppressMessages(library(Matrix))
PRJ <- file.path(Sys.getenv("DSCI_ONLINE_ROOT", "."), "FY/archr-subset-copnfol")
OUT <- file.path(Sys.getenv("DSCI_FIG_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", "."), "figures")), 'routeA_out/D_crosschemistry')
DER <- file.path(OUT, "derived"); TAB <- file.path(OUT, "tables")
say <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n", sep=""); flush.console() }

proj <- loadArchRProject(PRJ, showLogo = FALSE)
say("project loaded | %d cells | matrices: %s", length(proj$cellNames),
    paste(getAvailableMatrices(proj), collapse=", "))

em <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix", threads = 8)
geneDF <- ArchR:::.getFeatureDF(getArrowFiles(proj), "GeneScoreMatrix")
rownames(em) <- geneDF[, "name"]
M <- SummarizedExperiment::assays(em)[["GeneScoreMatrix"]]
cells <- sub("^.*#", "", colnames(M))
say("GeneScoreMatrix: %d genes x %d cells | per-nucleus sum=1e4: %s", nrow(M), ncol(M),
    all(round(Matrix::colSums(M)) == 1e4))
stopifnot(nrow(M) == 24333, ncol(M) == 70643, !anyDuplicated(cells))

meta <- read.csv(file.path(TAB, "fig2_light_meta.csv"))
rownames(meta) <- meta$cell
meta <- meta[cells, , drop = FALSE]
stopifnot(!any(is.na(meta$technique)))
keep <- meta$cell_type_updated %in% c("Oligo", "OPC", "COP/NFOL")
say("OL-lineage nuclei %d / %d | labels: %s", sum(keep), nrow(meta),
    paste(names(table(meta$cell_type_updated[keep])), table(meta$cell_type_updated[keep]), sep="=", collapse=", "))

Ms <- as(M[, keep, drop = FALSE], "CsparseMatrix")
ms <- meta[keep, , drop = FALSE]
for (tk in unique(ms$technique)) {
  idx <- which(ms$technique == tk)
  tag <- gsub("[^A-Za-z0-9]+", "_", tk)
  f <- file.path(DER, paste0("gs_", tag, ".mtx"))
  if (file.exists(paste0(f, ".gz"))) file.remove(paste0(f, ".gz"))
  writeMM(Ms[, idx, drop = FALSE], f)          # must write a path, not a gzfile() connection (fig5 already hit this: a connection yields an empty file)
  system(paste("gzip -f", f))
  f <- paste0(f, ".gz")
  if (file.size(f) < 1000) stop("empty mtx written for ", tk)
  writeLines(rownames(Ms), paste0(sub("\\.mtx\\.gz$", "", f), "_genes.txt"))
  write.csv(data.frame(cell = ms$cell[idx], tech = tk, label = ms$cell_type_updated[idx],
                       orig.ident = ms$orig.ident[idx]),
            file.path(TAB, paste0("cells_", tag, ".csv")), row.names = FALSE)
  say("wrote %s | %d genes x %d cells | nnz=%d", basename(f), nrow(Ms), length(idx),
      sum(Ms@p[idx + 1L] - Ms@p[idx]))
}
saveRDS(ms, file.path(DER, "gs_ol_meta.RDS"))
say("ALL_FINISHED")
