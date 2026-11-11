
suppressPackageStartupMessages({library(Seurat); library(dplyr)})
DATA  <- Sys.getenv("FIG3_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY")))
LIGHT <- Sys.getenv('FIG3_LIGHT', '')
OUT   <- Sys.getenv('FIG3_OUT', 'figS3_out')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sub <- c(EL14 = 'HT/EL14_subset.RDS', FY13 = 'FY13/3k_FY13_subset.RDS',
         FY14 = 'FY14/3k_FY14_subset.RDS')
copclust <- c(EL14 = NA_integer_, FY13 = 29L, FY14 = 31L)   # the COP/NFOL cluster of each library
libs <- list(); log <- character()
for (nm in names(sub)) {
  t0 <- proc.time()[3]
  o <- readRDS(file.path(DATA, sub[[nm]])); DefaultAssay(o) <- 'RNA'
  ct <- as.character(o$cell_type)
  if (!is.na(copclust[[nm]])) ct[o$RNA_snn_res.1 == copclust[[nm]]] <- 'COP/NFOL'
  d <- GetAssayData(o, assay = 'RNA', slot = 'data')
  cts <- sort(unique(ct))
  avg <- t(sapply(cts, function(k) rowMeans(as.matrix(d[, ct == k, drop = FALSE]))))
  libs[[nm]] <- list(n = length(ct), ct = ct, cells = colnames(o),
                     rna = as.data.frame(Embeddings(o, 'umap')), avg = avg,
                     genes = rownames(o), hvg = VariableFeatures(o),
                     counts = table(ct)[cts])
  log <- c(log, sprintf('%s cells=%d types=%d genes=%d hvg=%d loaded %.1f min', nm, ncol(o),
                        length(cts), length(libs[[nm]]$genes), length(libs[[nm]]$hvg),
                        (proc.time()[3] - t0) / 60))
  cat(tail(log, 1), '\n'); rm(o, d); gc(verbose = FALSE)
}
common <- Reduce(intersect, lapply(libs, function(s) s$genes))
hv <- intersect(unique(Reduce(union, lapply(libs, function(s) s$hvg))), common)
log <- c(log, sprintf('common genes %d | shared HVG used %d', length(common), length(hv)))
cat(tail(log, 1), '\n')
for (nm in names(libs)) {
  s <- libs[[nm]]
  s$avg <- s$avg[, hv, drop = FALSE]
  write.csv(round(s$avg, 5), file.path(OUT, paste0('S3_avgexpr_', nm, '.csv')))
  write.csv(as.data.frame(s$counts), file.path(OUT, paste0('S3_counts_', nm, '.csv')))
  s$genes <- s$hvg <- NULL; libs[[nm]] <- s
}
# ATAC embedding: Figure 3's cache already holds it (same cells, matched by barcode)
if (!is.na(LIGHT) && nzchar(LIGHT) && file.exists(LIGHT)) {
  S3 <- readRDS(LIGHT)
  for (nm in names(libs)) {
    idx <- match(libs[[nm]]$cells, S3$libs[[nm]]$cells)
    libs[[nm]]$atac <- S3$libs[[nm]]$atac[idx, , drop = FALSE]
    log <- c(log, sprintf('%s ATAC from Figure 3 cache: %d / %d cells', nm, sum(!is.na(idx)),
                          libs[[nm]]$n))
  }
} else {
  ht <- readRDS(file.path(DATA, '3k_dsci_ht_10x_gs_subset.RDS'))
  E <- Embeddings(ht, 'umap.atac')
  for (nm in names(libs)) {
    idx <- match(libs[[nm]]$cells, colnames(ht))
    libs[[nm]]$atac <- as.data.frame(E[idx, , drop = FALSE])
    log <- c(log, sprintf('%s ATAC from joint object: %d / %d cells', nm, sum(!is.na(idx)),
                          libs[[nm]]$n))
  }
  rm(ht, E); gc(verbose = FALSE)
}
cat(paste(log, collapse = '\n'), '\n')
types <- sort(unique(unlist(lapply(libs, function(s) s$ct))))
cols <- unique(c(as.vector(pals::glasbey(32)), as.vector(pals::alphabet(26))))[1:length(types)]
names(cols) <- types
saveRDS(list(libs = libs, cols = cols, types = types, hv = hv), file.path(OUT, 'S3_light.RDS'))
writeLines(log, file.path(OUT, 'S3_summary.txt'))
cat('wrote', file.path(OUT, 'S3_light.RDS'), '\n')
