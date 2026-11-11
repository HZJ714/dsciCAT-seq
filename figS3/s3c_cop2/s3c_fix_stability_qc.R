#!/usr/bin/env Rscript
# ============================================================================
# s3c_fix_stability_qc.R —— 修掉第一版两处口径 bug，重跑 S2（稳定性）与 S3（技术混杂）
#   bug 1：取簇用错了列。FindClusters(res=0.5) 写的是 RNA_snn_res.0.5，而对象里残留着
#          父对象带来的 RNA_snn_res.1（是个**常数**），于是 525 次运行全被读成"1 簇、ARI=0"。
#          → 改用 FindClusters 一定写入的 `seurat_clusters`。
#   bug 2：Cliff's delta 写成 outer(数据框列, ...) ⇒ 全 NA，判据无法评。→ 先 as.vector()。
#          另外 percent.ribo / log10GenesPerUMI / nucleosome_signal 在这 133 核里是空的，
#          可用量只有 6 个，图注要如实写"可评的 6 项"。
# ============================================================================
suppressMessages({library(Seurat); library(Signac); library(dplyr)})
OUT <- file.path(Sys.getenv("DSCI_FIG_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "figures")), 'u13_cop_4k/s3c')
sub14 <- readRDS(file.path(Sys.getenv("DSCI_SEURAT_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "pi/plots-seurat")), 'fy14_cop_nfol_subset_reclustered.RDS'))
DefaultAssay(sub14) <- 'RNA'
base <- as.integer(sub14$cluster_id); names(base) <- colnames(sub14)
cat('基线 133 核 ->', paste(names(table(base)), table(base), sep=' = ', collapse=' / '), '\n')

# ---- 先复现"基线那一跑"（与出图脚本同一套多分辨率调用），确认口径一致 ----
prep <- function(cells, multi = FALSE) {
  o <- subset(sub14, cells = cells)
  o <- suppressWarnings(FindVariableFeatures(o, nfeatures = 3000, verbose = FALSE))
  o <- ScaleData(o, verbose = FALSE); o <- RunPCA(o, npcs = 30, verbose = FALSE)
  FindNeighbors(o, reduction = 'pca', dims = 1:10, verbose = FALSE)
}
o_chk <- prep(colnames(sub14))
o_chk <- FindClusters(o_chk, algorithm = 3, resolution = c(0.1, 0.3, 0.5), random.seed = 1, verbose = FALSE)
tb <- sort(table(o_chk$seurat_clusters), decreasing = TRUE)
cat('重跑 res=0.5（最后一次调用）->', paste(names(tb), tb, sep=' = ', collapse=' / '), '\n')
cat('与交付 cluster_id 的 ARI =', {
  m <- as.integer(as.character(o_chk$seurat_clusters))[match(colnames(sub14), colnames(o_chk))]
  # 手工 ARI（不依赖 mclust）
  a <- base; b <- m; n <- length(a)
  tabs <- table(a, b); ni <- rowSums(tabs); nj <- colSums(tabs)
  sumjk <- sum(choose(tabs, 2)); sumi <- sum(choose(ni, 2)); sumj <- sum(choose(nj, 2))
  exp_ <- sumi * sumj / choose(n, 2); mx <- 0.5 * (sumi + sumj)
  sprintf('%.3f', (sumjk - exp_) / (mx - exp_))
}, '\n')

# ---------------------------------------------------------------- S2 扰动网格
RESOS <- c(0.1, 0.3, 0.5, 0.8, 1.0); SEEDS <- 1:5
conds <- list(full = colnames(sub14))
set.seed(20261029)
for (i in 1:20) conds[[paste0('draw', i)]] <- sample(colnames(sub14), floor(0.8 * ncol(sub14)))
ari <- function(a, b) {
  n <- length(a); tabs <- table(a, b); ni <- rowSums(tabs); nj <- colSums(tabs)
  s <- sum(choose(tabs, 2)); si <- sum(choose(ni, 2)); sj <- sum(choose(nj, 2))
  e <- si * sj / choose(n, 2); mx <- 0.5 * (si + sj)
  if (mx == e) return(NA_real_); (s - e) / (mx - e)
}
runs <- list(); asg <- list(); k <- 0
# 与**出图脚本完全一致**的调用方式：一次多分辨率 FindClusters，再按列取每一档
# （分开 5 次单分辨率调用会被前一次的 seurat_clusters 状态干扰，实测退化成 1 簇）
for (cn in names(conds)) {
  o <- tryCatch(prep(conds[[cn]]), error = function(e) NULL)
  if (is.null(o)) { cat('  跳过', cn, '\n'); next }
  for (sd in SEEDS) {
    x <- tryCatch(FindClusters(o, algorithm = 3, resolution = RESOS, random.seed = sd, verbose = FALSE),
                  error = function(e) NULL)
    if (is.null(x)) next
    nm <- as.character(conds[[cn]])
    for (res in RESOS) {
      col <- paste0('RNA_snn_res.', sub('\\.0$', '', format(res, trim = TRUE)))
      col <- if (res == 1.0) 'RNA_snn_res.1' else col
      if (!col %in% colnames(x@meta.data)) { cat('  缺列', col, '\n'); next }
      cl <- as.integer(as.character(x[[col]][nm, , drop = TRUE]))
      tb <- sort(table(cl), decreasing = TRUE); k <- k + 1
      runs[[k]] <- data.frame(run = k, cond = cn, n_cells = length(nm), res = res, seed = sd,
                              n_clusters = length(tb), size1 = as.integer(tb[1]),
                              size2 = ifelse(length(tb) > 1, as.integer(tb[2]), 0L),
                              min_size = as.integer(min(tb)), stringsAsFactors = FALSE)
      asg[[k]] <- data.frame(run = k, cell = nm, cluster = cl, baseline = base[nm],
                             stringsAsFactors = FALSE)
    }
  }
  cat('  完成', cn, '\n')
}
stab <- do.call(rbind, runs); stab$kind <- ifelse(stab$cond == 'full', 'full', 'downsampled')
stab$split2 <- stab$n_clusters == 2 & stab$min_size >= 10
write.csv(stab, file.path(OUT, 's3c_stability_long.csv'), row.names = FALSE)
write.csv(do.call(rbind, asg), file.path(OUT, 's3c_assignments.csv'), row.names = FALSE)
cat('\n== S2 == res=0.5 复现', sum(stab$split2 & stab$res == 0.5), '/', sum(stab$res == 0.5),
    '| 全网格', sum(stab$split2), '/', nrow(stab), '\n')
print(aggregate(cbind(split2, ari) ~ res, data = stab, FUN = mean))

# ---------------------------------------------------------------- S3 技术混杂
cliff <- function(x, y) {   # Cliff's delta = P(x>y) - P(y>x)
  r <- as.numeric(rank(c(x, y)))
  n1 <- length(x); n2 <- length(y)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2) * 2 - 1
}
QC <- c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo', 'log10GenesPerUMI',
        'nCount_ATAC', 'nFeature_ATAC', 'nucleosome_signal', 'RNA.weight')
qc <- do.call(rbind, lapply(QC, function(v) {
  z <- as.vector(sub14[[v, drop = TRUE]])
  x <- z[base == 0]; y <- z[base == 1]
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 20) return(data.frame(var = v, available = FALSE, n = sum(ok),
                                      med0 = NA_real_, med1 = NA_real_, p = NA_real_,
                                      p_bh = NA_real_, delta = NA_real_, auc = NA_real_))
  x <- x[ok]; y <- y[ok]
  p <- suppressWarnings(wilcox.test(x, y)$p.value)
  data.frame(var = v, available = TRUE, n = sum(ok), med0 = median(x), med1 = median(y),
             p = p, p_bh = NA_real_, delta = cliff(x, y),
             auc = (sum(as.numeric(outer(x, y, '>'))) + 0.5 * sum(as.numeric(outer(x, y, '==')))) /
                   (length(x) * length(y)))
}))
qc$p_bh <- p.adjust(qc$p, method = 'BH')
qc$flag_bad <- !is.na(qc$p_bh) & qc$p_bh < 0.05 & !is.na(qc$delta) & abs(qc$delta) > 0.2
write.csv(qc, file.path(OUT, 's3c_qc.csv'), row.names = FALSE)
cat('\n== S3 ==\n'); print(qc[, c('var', 'available', 'med0', 'med1', 'p_bh', 'delta', 'flag_bad')])
cat('\n完成\n')
