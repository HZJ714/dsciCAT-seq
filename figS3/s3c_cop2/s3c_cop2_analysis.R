#!/usr/bin/env Rscript
# ============================================================================
# s3c_cop2_analysis.R —— 100k COP/NFOL「两簇」的四格证据（判据冻结见
#   <work>/figures/S3C_COP2_FROZEN_CRITERIA_20261029.md，先写规则后看结果）
# 输出（全部 CSV/JSON，不画图；画图由 s3c_cop2_figure.py 负责）:
#   s3c_marker_table.csv  独立 marker × 各参照簇
#   s3c_de.csv            两簇无先验差异（Wilcoxon + AUC + 均值差）
#   s3c_stability_long.csv 扰动网格逐次结果（ARI 在 Python 里用 sklearn 算）
#   s3c_assignments.csv   每次运行的细胞级分配（算 ARI 用）
#   s3c_qc.csv            9 个技术量 × 两簇（p / BH / Cliff's delta）
#   s3c_summary.json      S1/S2/S3 三条判据的通过与否（含负结果）
# ============================================================================
suppressMessages({library(Seurat); library(Signac); library(dplyr)})
set.seed(1)
OUT <- file.path(Sys.getenv("DSCI_FIG_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "figures")), 'u13_cop_4k/s3c'); dir.create(OUT, showWarnings = FALSE)
P   <- Sys.getenv("DSCI_SEURAT_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "pi/plots-seurat")); FYP <- Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY"))

# 判据里用到的基因（S1 必须排除，避免自证）
CALL_GENES <- c('Bmp4', 'Gpr17', 'C1ql1', 'Ptprz1')
# 独立 marker（预先指定，不看结果挑；Ccl4aat/Htraf1/Asf1 在矩阵里没有，会自动落 NA 并被剔除）
INDEP_GENES <- c('Pdgfra', 'Cspg4', 'Sox6', 'Myc', 'Pdgfrb', 'Dcn', 'Gpc5', 'Mobp', 'Mal',
                 'Sox10', 'Olig1', 'Cdk1', 'Mki67', 'Top2a', 'Fn1', 'Sox11', 'Pdgfc',
                 'Trf', 'Cnp', 'Mog', 'Dcx')
QC_VARS <- c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo', 'log10GenesPerUMI',
             'nCount_ATAC', 'nFeature_ATAC', 'nucleosome_signal', 'RNA.weight')

sub14 <- readRDS(file.path(P, 'fy14_cop_nfol_subset_reclustered.RDS'))   # 100k, 133 核
DefaultAssay(sub14) <- 'RNA'
sub13 <- readRDS(file.path(P, 'fy13_cop_nfol_subset_reclustered.RDS'))   # 50k, 60 核
DefaultAssay(sub13) <- 'RNA'
opc14 <- readRDS(file.path(FYP, 'FY14/3k_FY14_subset.RDS'))              # 同库全脑（参照）
DefaultAssay(opc14) <- 'RNA'

base <- sub14$cluster_id
names(base) <- colnames(sub14)
cat('基线：133 核 ->', paste(names(table(base)), table(base), sep=' = ', collapse=' / '), '\n')

# ---------------------------------------------------------------- 1. marker 表
tab_of <- function(obj, cells, label) {
  d <- GetAssayData(obj, assay = 'RNA', slot = 'data')
  gs <- intersect(c(INDEP_GENES, CALL_GENES), rownames(d))
  m <- as.matrix(d[gs, cells, drop = FALSE])
  data.frame(group = label, gene = rownames(m),
             mean = rowMeans(m), pct = 100 * rowMeans(m > 0), n = length(cells),
             independent = rownames(m) %in% INDEP_GENES, stringsAsFactors = FALSE)
}
oligo_ref <- colnames(opc14)[opc14$cell_type %in% c('Oligo', 'OPC')]
mk <- rbind(
  tab_of(sub14, names(base)[base == 0], '100k COP/NFOL sub-0'),
  tab_of(sub14, names(base)[base == 1], '100k COP/NFOL sub-1'),
  tab_of(sub13, colnames(sub13),        '50k COP/NFOL (single)'),
  tab_of(opc14, colnames(opc14)[opc14$cell_type == 'OPC'],   '100k OPC (same lib)'),
  tab_of(opc14, colnames(opc14)[opc14$cell_type == 'Oligo'], '100k Oligo (same lib)'))
mk$gene <- factor(mk$gene, levels = c(INDEP_GENES, CALL_GENES))
write.csv(mk, file.path(OUT, 's3c_marker_table.csv'), row.names = FALSE)

# ---------------------------------------------------------------- 2. 无先验 DE
Idents(sub14) <- 'cluster_id'
# Seurat 4.3.0：检验名是 'wilcox'（写成 'wilcoxon' 会报 Unknown test），参数名是 min.pct
de <- tryCatch(FindMarkers(sub14, ident.1 = 0, ident.2 = 1, assay = 'RNA',
                           test.use = 'wilcox', logfc.threshold = 0, min.pct = 0),
               error = function(e) {cat('FindMarkers 失败:', conditionMessage(e), '\n'); NULL})
if (is.null(de) || nrow(de) == 0) stop('DE 算不出来，直接终止（不硬凑结论）')
# AUC = |U| / (n1*n2)（与 Wilcoxon 同一条秩统计量，直接换成效应量）
d14 <- GetAssayData(sub14, assay = 'RNA', slot = 'data')
auc_of <- function(g) {
  x <- as.numeric(d14[g, ]); a <- x[base == 0]; b <- x[base == 1]
  U <- sum(as.numeric(outer(a, b, '>'))) + 0.5 * sum(as.numeric(outer(a, b, '==')))
  abs(U - length(a) * length(b) / 2) / (length(a) * length(b)) + 0.5
}
de$mean_diff <- sapply(rownames(de), function(g) mean(d14[g, base == 0]) - mean(d14[g, base == 1]))
de$AUC <- sapply(rownames(de), auc_of)
de$gene <- rownames(de)
de$in_call_set <- de$gene %in% CALL_GENES
de$qualifies <- (!de$in_call_set) & abs(de$mean_diff) >= 0.5 & de$p_val_adj < 0.05
de <- de[order(-abs(de$mean_diff)), ]
write.csv(de, file.path(OUT, 's3c_de.csv'), row.names = FALSE)

# ---------------------------------------------------------------- 3. 稳定性扰动
# 每个条件只做一次 HVG/Scale/PCA/Neighbors（21 次），FindClusters 只跑图上的社区划分（21×25）
RESOS <- c(0.1, 0.3, 0.5, 0.8, 1.0); SEEDS <- 1:5
prep <- function(cells) {
  o <- subset(sub14, cells = cells)
  o <- suppressWarnings(FindVariableFeatures(o, nfeatures = 3000, verbose = FALSE))
  o <- ScaleData(o, verbose = FALSE)
  o <- RunPCA(o, npcs = 30, verbose = FALSE)
  o <- FindNeighbors(o, reduction = 'pca', dims = 1:10, verbose = FALSE)
  o
}
conds <- list(full = colnames(sub14))
set.seed(20261029)
for (i in 1:20) conds[[paste0('draw', i)]] <- sample(colnames(sub14), floor(0.8 * ncol(sub14)))
runs <- list(); asg <- list(); k <- 0
for (cn in names(conds)) {
  o <- tryCatch(prep(conds[[cn]]), error = function(e) NULL)
  if (is.null(o)) { cat('  跳过', cn, '\n'); next }
  for (res in RESOS) for (sd in SEEDS) {
    a <- tryCatch({
      x <- FindClusters(o, algorithm = 3, resolution = res, random.seed = sd, verbose = FALSE)
      setNames(as.character(x$RNA_snn_res.1), colnames(x))
    }, error = function(e) NULL)
    if (is.null(a)) next
    tb <- sort(table(a), decreasing = TRUE); k <- k + 1
    runs[[k]] <- data.frame(run = k, cond = cn, n_cells = length(a), res = res, seed = sd,
                            n_clusters = length(tb), size1 = as.integer(tb[1]),
                            size2 = ifelse(length(tb) > 1, as.integer(tb[2]), 0L),
                            min_size = as.integer(min(tb)), stringsAsFactors = FALSE)
    asg[[k]] <- data.frame(run = k, run_id = sprintf('%s_r%s_s%d', cn, res, sd),
                           cell = names(a), cluster = as.integer(as.character(a)),
                           baseline = as.integer(base[names(a)]), stringsAsFactors = FALSE)
  }
  cat('  完成', cn, '(', length(conds[[cn]]), '核 )\n')
}
stab <- do.call(rbind, runs)
write.csv(stab, file.path(OUT, 's3c_stability_long.csv'), row.names = FALSE)
write.csv(do.call(rbind, asg), file.path(OUT, 's3c_assignments.csv'), row.names = FALSE)

# ---------------------------------------------------------------- 4. 技术混杂
cliff <- function(x, y) {
  m <- mean(outer(x, y, '>')) - mean(outer(x, y, '<'))
  m
}
qc <- do.call(rbind, lapply(QC_VARS, function(v) {
  x <- sub14[[v, drop = TRUE]][base == 0]; y <- sub14[[v, drop = TRUE]][base == 1]
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 20) return(data.frame(var = v, available = FALSE, n0 = sum(!is.na(x)), n1 = sum(!is.na(y)),
                                      med0 = NA_real_, med1 = NA_real_, p = NA_real_, p_bh = NA_real_, delta = NA_real_))
  wt <- suppressWarnings(wilcox.test(x[ok], y[ok])$p.value)
  data.frame(var = v, available = TRUE, n0 = sum(ok & (base == 0)), n1 = sum(ok & (base == 1)),
             med0 = median(x[ok]), med1 = median(y[ok]), p = wt, p_bh = NA_real_,
             delta = cliff(x[ok], y[ok]))
}))
qc$p_bh <- p.adjust(qc$p, method = 'BH')
qc$flag_bad <- (!is.na(qc$p) & qc$p_bh < 0.05) & (!is.na(qc$delta) & abs(qc$delta) > 0.2)
write.csv(qc, file.path(OUT, 's3c_qc.csv'), row.names = FALSE)

# ---------------------------------------------------------------- 5. 判据结论
S1 <- list(n_qual = sum(de$qualifies),
           genes = head(de$gene[de$qualifies], 20),
           pass = sum(de$qualifies) >= 3)
cat('\n== S1 独立基因 ==', ifelse(S1$pass, 'PASS', 'FAIL'), ' 合格基因数 =', S1$n_qual, '\n')
print(head(de[de$qualifies & !is.na(de$qualifies), c('gene', 'mean_diff', 'p_val_adj', 'AUC')], 10))
cat('\n== S3 技术混杂 ==\n'); print(qc[, c('var', 'available', 'p_bh', 'delta', 'flag_bad')])
writeLines('{ "note": "S2 的复现率与 ARI 由 s3c_score_stability.py 计算后回填" }',
           file.path(OUT, 's3c_summary_partial.json'))
cat('\n输出目录:', OUT, '\n')
