
suppressPackageStartupMessages({library(Seurat)})
DATA <- Sys.getenv("FIGS3_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY"))); OUT <- Sys.getenv('FIGS3_OUT', 'figS3_out')
SRC <- c(FY13 = 'FY13/3k_FY13_subset.RDS', FY14 = 'FY14/3k_FY14_subset.RDS',
         HTp  = 'HT/HT_subset_harmony.RDS', EL14 = 'HT/EL14_subset.RDS')
for (k in names(SRC)) {
  f <- file.path(OUT, paste0('S3_lineage_', k, '.rds'))
  if (file.exists(f) && Sys.getenv('FIGS3_FORCE') != '1') { cat(k, 'cached\n'); next }
  o <- readRDS(file.path(DATA, SRC[[k]])); nall <- ncol(o)
  keep <- as.character(o$cell_type) %in% c('Oligo', 'OPC', 'COP_NFOL', 'COP/NFOL')
  o <- subset(o, cells = colnames(o)[keep]); DefaultAssay(o) <- 'RNA'
  saveRDS(o[['RNA']], f)
  cat(sprintf('%-5s whole brain %6d | lineage %5d (%.2f%%)\n', k, nall, sum(keep), 100 * sum(keep) / nall))
  rm(o); gc(verbose = FALSE)
}
