# -*- coding: utf-8 -*-
# D2e: the RNA channel switches to the MEASURED expression layer - Seurat v4 RNA$data (LogNormalize, log1p after CP10k),
#      instead of the ArchR GeneIntegrationMatrix (V1 showed it is a CCA-integrated/imputed product: 48% non-zero vs 6% measured,
#      flattening or even reversing the Plp1/Mbp gradient between OPC and mature OL).
#      dsciCAT takes the FY13/FY14 multiome RDS; HT / 10x reuse fig5's export_tech_expr.R products (the same RNA$data).
suppressMessages({library(Seurat); library(SeuratObject); library(Signac); library(Matrix)})
U   <- Sys.getenv("DSCI_ONLINE_ROOT", ".")
OUT <- file.path(Sys.getenv("DSCI_FIG_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", "."), "figures")), 'routeA_out/D_crosschemistry')
DER <- file.path(OUT,"derived"); TAB <- file.path(OUT,"tables")
say <- function(...) { cat(format(Sys.time(),"%H:%M:%S"), sprintf(...), "\n", sep=""); flush.console() }
ol <- read.csv(file.path(TAB,"cells_dsciCAT_seq.csv"))
parts <- list()
for (s in c("FY13","FY14")) {
  f <- file.path(U,"FY",s,paste0("3k_",s,"_multiome.RDS"))
  say("read %s", basename(f)); o <- readRDS(f)
  want <- ol$cell[startsWith(ol$cell, paste0(s,"_"))]
  hit  <- intersect(want, colnames(o))
  x <- GetAssayData(o, assay="RNA", slot="data")            # v4 slot API (LayerData/$data are unavailable in this environment)
  x <- as(x[, hit, drop=FALSE], "CsparseMatrix")
  say("  %s: object %d nuclei | needed %d | hit %d | RNA data %d genes x %d | range %.3f..%.3f | per-nucleus sum median %.0f",
      s, ncol(o), length(want), length(hit), nrow(x), ncol(x), min(x@x), max(x@x), median(Matrix::colSums(x)))
  p <- file.path(DER, paste0("rnaSeurat_dsciCAT_", s, ".mtx"))
  if (file.exists(paste0(p,".gz"))) file.remove(paste0(p,".gz"))
  writeMM(x, p); system(paste("gzip -f", p))
  writeLines(rownames(x), paste0(p, "_genes.txt"))
  parts[[s]] <- data.frame(cell=hit, genes=nrow(x))
  md <- o@meta.data[hit, , drop=FALSE]
  write.csv(data.frame(cell=hit, nCount_RNA=md$nCount_RNA, nFeature_RNA=md$nFeature_RNA,
                       nCount_ATAC=md$nCount_ATAC, nFeature_ATAC=md$nFeature_ATAC),
            paste0(p, "_meta.csv"), row.names=FALSE)
  rm(x, o, md); invisible(gc())
}
say("ALL_FINISHED")
