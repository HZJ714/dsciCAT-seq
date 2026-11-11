# -*- coding: utf-8 -*-
# D2c: per-nucleus covariates of the three chemistries (RNA depth / ATAC depth / gene_score detection count), unified into one table
#   10x        <- GitHub/fig5/fig5_out/tenx_meta.csv          (WNN Seurat meta)
#   HT         <- GitHub/fig5/fig5_out/HT_harmony_meta.csv    (WNN Seurat meta)
#   dsciCAT    <- FY/FY13,FY14/3k_*_multiome.RDS @meta.data   (the same WNN-family objects)
suppressMessages({library(Seurat); library(SeuratObject); library(Signac)})
U   <- Sys.getenv("DSCI_ONLINE_ROOT", ".")
OUT <- file.path(Sys.getenv("DSCI_FIG_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", "."), "figures")), 'routeA_out/D_crosschemistry')
TAB <- file.path(OUT,"tables")
say <- function(...) { cat(format(Sys.time(),"%H:%M:%S"), sprintf(...), "\n", sep=""); flush.console() }
COLS <- c("nCount_RNA","nFeature_RNA","nCount_ATAC","nFeature_ATAC","nCount_gene_score","nFeature_gene_score")

grab <- function(md, tag){
  k <- intersect(COLS, colnames(md))
  d <- data.frame(cell=rownames(md), tech=tag, check.names=FALSE)
  for (cc in COLS) d[[cc]] <- if (cc %in% colnames(md)) md[[cc]] else NA    # fill missing columns with NA so the sources rbind with one column structure
  say("%s: %d nuclei | valued columns %s", tag, nrow(d), paste(k, collapse=","))
  d[, c("cell","tech",COLS), drop=FALSE]
}
parts <- list()
m10 <- read.csv(file.path(U,"GitHub/fig5/fig5_out/tenx_meta.csv"), row.names=1)
parts[["10x"]] <- grab(m10, "10x")
mht <- read.csv(file.path(U,"GitHub/fig5/fig5_out/HT_harmony_meta.csv"), row.names=1)
parts[["HT"]] <- grab(mht, "HT-scCAT-seq-HT")
for (s in c("FY13","FY14")) {
  f <- file.path(U,"FY",s,paste0("3k_",s,"_multiome.RDS"))
  o <- readRDS(f); parts[[s]] <- grab(o@meta.data, "dsciCAT-seq"); rm(o); invisible(gc())
}
cov <- do.call(rbind, parts)
write.csv(cov, file.path(TAB,"cov_by_cell.csv"), row.names=FALSE)
say("cov_by_cell.csv: %d rows", nrow(cov))

# coverage check: do all 8,948 OL-lineage nuclei have covariates
need <- read.csv(file.path(TAB,"fig2_light_meta.csv"))
need <- need[need$cell_type_updated %in% c("Oligo","OPC","COP/NFOL"),]
hit <- match(need$cell, cov$cell)
say("OL-lineage nuclei %d | covariates hit %d | gaps: %s", nrow(need), sum(!is.na(hit)),
    paste(head(unique(need$technique[is.na(hit)])), collapse=","))
print(table(need$technique[is.na(hit)]))
say("ALL_FINISHED")
