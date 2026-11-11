#!/usr/bin/env Rscript
# 导出 133 核的 raw counts + 元信息，供 Python 端做「深度匹配重检」（cell 2b）
suppressMessages({library(Seurat); library(Signac); library(Matrix)})
sub <- readRDS(file.path(Sys.getenv("DSCI_SEURAT_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "pi/plots-seurat")), 'fy14_cop_nfol_subset_reclustered.RDS'))
DefaultAssay(sub) <- 'RNA'
cnt <- GetAssayData(sub, assay = 'RNA', slot = 'counts')
keep <- rowSums(cnt > 0) >= 10                       # 至少 10 核表达，避免全零基因
cat('基因', nrow(cnt), '-> 保留', sum(keep), '\n')
m <- as.matrix(cnt[keep, , drop = FALSE])
df <- data.frame(gene = rownames(m), m, check.names = FALSE)
write.csv(df, file.path(Sys.getenv("DSCI_FIG_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "figures")), 'u13_cop_4k/s3c/s3c_counts.csv'), row.names = FALSE)
meta <- data.frame(cell = colnames(sub), cluster = as.integer(sub$cluster_id),
                   nCount_RNA = sub$nCount_RNA, nFeature_RNA = sub$nFeature_RNA,
                   percent.mt = sub$percent.mt)
write.csv(meta, file.path(Sys.getenv("DSCI_FIG_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "figures")), 'u13_cop_4k/s3c/s3c_meta.csv'), row.names = FALSE)
cat('已导出\n')
