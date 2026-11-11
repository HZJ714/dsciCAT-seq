#!/usr/bin/env Rscript
# ============================================================================
# s3c_cop_umap_pair.R —— 50k 与 100k 的 COP/NFOL 用**同一套参数**再聚类 + UMAP
# ============================================================================
# 作者新图版式最左要放这对 UMAP，并在图里写"identical parameters"。
# 那就不能拿两个对象里**各自历史遗留**的 reduction 直接画（无法证明参数一致），
# 而是这里现跑一遍：prep() 与 s3c_cop2_analysis.R 里稳定性扰动用的**同一个函数体**
# （HVG 3000 → ScaleData → PCA 30 → Neighbors dims 1:10），聚类 algorithm=3 / res=0.5 / seed=1。
# 同时把"新跑的分区 vs 交付分区"的 ARI 打出来：100k 必须复现 81/52，否则这块图不能这么画。
suppressMessages({library(Seurat); library(dplyr)})
P <- Sys.getenv("DSCI_SEURAT_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "pi/plots-seurat")); OUT <- file.path(Sys.getenv("DSCI_FIG_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "figures")), 'u13_cop_4k/s3c')
objs <- list(fy13_50k = 'fy13_cop_nfol_subset_reclustered.RDS',
             fy14_100k = 'fy14_cop_nfol_subset_reclustered.RDS')
prep <- function(o) {
  o <- suppressWarnings(FindVariableFeatures(o, nfeatures = 3000, verbose = FALSE))
  o <- ScaleData(o, verbose = FALSE)
  o <- RunPCA(o, npcs = 30, verbose = FALSE)
  o <- FindNeighbors(o, reduction = 'pca', dims = 1:10, verbose = FALSE)
  RunUMAP(o, reduction = 'pca', dims = 1:10, seed.use = 1, verbose = FALSE)
}
# ⚠ 取簇的坑（s3c_fix_stability_qc.R 里记过，我第一版又踩了一次）：
#   对象里残留着父对象带来的 `RNA_snn_res.1`，它是个**常数**（133 核全在 level 31），
#   所以 `x$RNA_snn_res.1` 永远读到 1 簇。必须用一次多分辨率 FindClusters 再按档取列。
RESOS <- c(0.1, 0.3, 0.5, 0.8, 1.0)
col_for <- function(res) if (res == 1.0) 'RNA_snn_res.1' else paste0('RNA_snn_res.', res)
all_rows <- list()
for (nm in names(objs)) {
  o <- readRDS(file.path(P, objs[[nm]])); DefaultAssay(o) <- 'RNA'
  o <- prep(o)
  cl <- FindClusters(o, algorithm = 3, resolution = RESOS, random.seed = 1, verbose = FALSE)
  nmv <- as.character(colnames(cl))
  have <- grep('^RNA_snn_res', colnames(cl@meta.data), value = TRUE)
  cat('\n===', nm, ' n =', ncol(cl), ' FindClusters 后 meta 里的 RNA_snn_res 列:',
      paste(have, collapse = ' '), '\n')
  grab <- function(x, col) as.character(x[[col]][nmv, , drop = TRUE])   # 照搬 fix 脚本的取法
  fresh <- grab(cl, col_for(0.5)); names(fresh) <- nmv
  del <- as.character(cl$cluster_id); names(del) <- nmv
  tb <- sort(table(fresh), decreasing = TRUE)
  cat('   同参数 res=0.5 新跑簇数 =', length(tb), ' 大小 =', paste(tb, collapse = '/'),
      '  交付分区 =', paste(names(table(del)), table(del), sep = '=', collapse = ' '), '\n')
  for (r in RESOS) {
    col <- col_for(r)
    if (!col %in% have) { cat('    res', r, '-> 缺列', col, '\n'); next }
    v <- grab(cl, col); tv <- sort(table(v), decreasing = TRUE)
    cat('    res', r, '->', length(tv), '簇', paste(tv, collapse = '/'), '\n')
  }
  e <- as.data.frame(Embeddings(cl, 'umap')); colnames(e) <- c('UMAP_1', 'UMAP_2')
  all_rows[[nm]] <- data.frame(cell = rownames(e), lib = nm, UMAP_1 = e$UMAP_1, UMAP_2 = e$UMAP_2,
                               delivered = del[rownames(e)], fresh = fresh[rownames(e)],
                               nCount_RNA = cl$nCount_RNA, nFeature_RNA = cl$nFeature_RNA)
}
d <- do.call(rbind, all_rows)
write.csv(d, file.path(OUT, 's3c_cop_umap_pair.csv'), row.names = FALSE)
cat('\n导出 s3c_cop_umap_pair.csv', nrow(d), '行\n')
cat('每簇中位 nCount_RNA:\n'); print(aggregate(nCount_RNA ~ lib + delivered, d, median))
