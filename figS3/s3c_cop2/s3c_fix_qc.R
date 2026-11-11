#!/usr/bin/env Rscript
# S3（技术混杂）重算：第一版 `sub14[[v]]` 返回的是 data.frame，用 [logical] 取行取错了
# （警告 "longer object length is not a multiple..."、delta 全 NA）⇒ 改用 FetchData 拿数值向量。
suppressMessages({library(Seurat); library(Signac)})
OUT <- file.path(Sys.getenv("DSCI_FIG_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "figures")), 'u13_cop_4k/s3c')
sub14 <- readRDS(file.path(Sys.getenv("DSCI_SEURAT_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "pi/plots-seurat")), 'fy14_cop_nfol_subset_reclustered.RDS'))
base <- as.integer(sub14$cluster_id)
cliff <- function(x, y) {
  r <- rank(c(x, y)); n1 <- length(x); n2 <- length(y)
  2 * ((sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2)) - 1
}
QC <- c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo', 'log10GenesPerUMI',
        'nCount_ATAC', 'nFeature_ATAC', 'nucleosome_signal', 'RNA.weight')
rows <- lapply(QC, function(v) {
  z <- tryCatch(FetchData(sub14, v)[, 1], error = function(e) NULL)
  if (is.null(z) || all(is.na(z))) return(data.frame(var = v, available = FALSE, n = 0L,
        med0 = NA_real_, med1 = NA_real_, p = NA_real_, p_bh = NA_real_, delta = NA_real_, auc = NA_real_))
  # 两侧各自去 NA（原来 !is.na(x) & !is.na(y) 长度 81 vs 52 会被循环补全 ⇒ 统计量全错）
  x <- z[base == 0]; y <- z[base == 1]
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 10 || length(y) < 10) return(data.frame(var = v, available = FALSE, n = length(x) + length(y),
        med0 = NA_real_, med1 = NA_real_, p = NA_real_, p_bh = NA_real_, delta = NA_real_, auc = NA_real_))
  p <- suppressWarnings(wilcox.test(x, y)$p.value)
  auc <- (sum(as.numeric(outer(x, y, '>'))) + 0.5 * sum(as.numeric(outer(x, y, '==')))) / (length(x) * length(y))
  data.frame(var = v, available = TRUE, n = length(x) + length(y), med0 = median(x), med1 = median(y),
             p = p, p_bh = NA_real_, delta = cliff(x, y), auc = auc)
})
qc <- do.call(rbind, rows); qc$p_bh <- p.adjust(qc$p, method = 'BH')
qc$flag_bad <- !is.na(qc$p_bh) & qc$p_bh < 0.05 & !is.na(qc$delta) & abs(qc$delta) > 0.2
write.csv(qc, file.path(OUT, 's3c_qc.csv'), row.names = FALSE)
print(qc[, c('var', 'available', 'n', 'med0', 'med1', 'p', 'p_bh', 'delta', 'auc', 'flag_bad')], digits = 3)
cat('\n可评项 =', sum(qc$available), '| 判为技术差异(flag_bad) =', sum(qc$flag_bad, na.rm = TRUE),
    '| S3 ', if (sum(qc$flag_bad, na.rm = TRUE) == 0) 'PASS' else 'FAIL', '\n')
