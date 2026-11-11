#!/usr/bin/env Rscript
# ============================================================================
# s3c_cop_umap_export.R —— 50k / 100k COP-NFOL 再聚类对象的 UMAP 坐标 + "同参数"证据
# ============================================================================
# 作者要的新版式：最左放 UMAP（50k 与 100k 的 COP/NFOL 再聚类，标出 sub-0/sub-1），
# 用来说"同一套参数下高通量才分得出第二群"。所以这里除了坐标，必须把图里写
# "identical parameters" 的**依据**原样打出来核对，不能凭记忆写。
suppressMessages({library(Seurat)})
P <- Sys.getenv("DSCI_SEURAT_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "pi/plots-seurat")); OUT <- file.path(Sys.getenv("DSCI_FIG_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "figures")), 'u13_cop_4k/s3c')
objs <- list(fy13_50k = 'fy13_cop_nfol_subset_reclustered.RDS',
             fy14_100k = 'fy14_cop_nfol_subset_reclustered.RDS')
rows <- list()
for (nm in names(objs)) {
  o <- readRDS(file.path(P, objs[[nm]])); DefaultAssay(o) <- 'RNA'
  cat('\n===', nm, ' 核数 =', ncol(o), ' 降维 =', paste(Reductions(o), collapse = ','), '\n')
  cat('   meta.data 列:', paste(colnames(o@meta.data), collapse = ', '), '\n')
  cl <- intersect(c('cluster_id', 'seurat_clusters', 'RNA_snn_res.0.5'), colnames(o@meta.data))
  for (k in cl) cat('   ', k, ':', paste(names(table(o[[k, drop = TRUE]])), table(o[[k, drop = TRUE]]),
                                         sep = '=', collapse = ' '), '\n')
  for (g in intersect(c('snn', 'RNA_snn'), names(o[[]]))) {
    m <- tryCatch(o[[g]]@graph@misc, error = function(e) NULL)
    if (!is.null(m)) cat('   ', g, '参数: dims =', paste(m$dims, collapse = ','),
                         '| algorithm =', m$algorithm, '| resolution =', m$resolution,
                         '| random.seed =', m$random.seed, '| method =', m$method, '\n')
  }
  ja <- tryCatch(SeuratObject::ToolAssays(o[['RNA']]), error = function(e) list())
  for (w in c('RunPCA', 'FindNeighbors', 'FindClusters', 'RunUMAP'))
    if (length(ja[[w]])) cat('   ', w, ':', paste(names(ja[[w]]), unlist(ja[[w]]), sep = ' = ',
                                                  collapse = '; '), '\n')
  if ('umap' %in% Reductions(o)) {
    e <- as.data.frame(Embeddings(o, 'umap')); e$cell <- rownames(e); e$lib <- nm
    e$cluster <- if (length(cl)) as.character(o[[cl[1], drop = TRUE]]) else NA
    e$nCount_RNA <- o$nCount_RNA
    rows[[nm]] <- e
    cat('   umap x [', round(min(e$umap_1), 1), ',', round(max(e$umap_1), 1), '] y [',
        round(min(e$umap_2), 1), ',', round(max(e$umap_2), 1), ']  簇数 =',
        length(unique(na.omit(e$cluster))), '\n')
  }
}
if (length(rows)) {
  d <- do.call(rbind, rows)
  write.csv(d, file.path(OUT, 's3c_cop_umap.csv'), row.names = FALSE)
  cat('\n导出 s3c_cop_umap.csv 行数', nrow(d), '\n')
  print(table(lib = d$lib, cluster = d$cluster))
}
