#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# fig1c_v4.R —— Figure 1 panel C 重画（Codex #237：缺坐标/基因组版本/归一化口径 + 基因名错）
# ---------------------------------------------------------------------------------------------
# 为什么重画而不是修补：交付版 C 里只有 ATAC 覆盖行与基因模型是矢量，
#   「300 single cells」三行 = 3 张 300x300 px 位图（175 dpi），RNA 三行 = 1 张 1149x379 px 位图，
#   且 RNA 行的来源/归一化在任何脚本里都查不到（track plots*.ipynb 只留了 bedgraph->bw 的中间步骤），
#   正文也从未引用 Figure 1C。=> 位图块 + 不可复现的轨道 = 审稿人一问就崩，故整块重画。
# 本版做的事：
#   1) 坐标轴：CoveragePlot 自带的 x 轴不再 blank 掉，刻度落在 5-kb 整倍数，标签 kb（万位带千分位）；
#   2) 基因组版本：人 hg38 / 鼠 mm10 印在细胞系标题里，区间写在 x 轴标题；
#   3) 归一化口径：三个覆盖行**共用一个 y 量程**（原图是各自 auto-scale，行间不可比），量程写在左列；
#      tile 行显式 limits 0-2（Signac 把每细胞片段归一到 10,000）；
#   4) 基因名由 EnsDb 给出：人 MIR17HG（HGNC 全大写），鼠 Hes7（MGI 首字母大写，不是 Hex7）。
# 数据：<raw> 下 per-species peak 矩阵 + fragments（与 D/E/F 同源）。
# 输出：replot/F1C_tracks_v4.pdf —— 尺寸 = fig1_layout 的 C 槽位 11.85 x 6.53 cm
suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(Matrix); library(readr); library(magrittr)
  library(ggplot2); library(patchwork); library(GenomicRanges); library(scales)
  library(EnsDb.Hsapiens.v86); library(EnsDb.Mmusculus.v79)
})

DATA <- Sys.getenv("FIG1_DATA", "")
OUT  <- file.path(Sys.getenv("DSCI_FIG_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "figures")), 'replot/F1C_tracks_v4.pdf')
PNG  <- file.path(Sys.getenv("DSCI_FIG_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "figures")), 'replot/F1C_tracks_v4.png')
PW <- 11.85; PH <- 6.53                       # = fig1_layout.C_SLOT 的 w/h
NCOL <- 300                                   # tile 每方法细胞数（沿用原图 300）
COL  <- c('#E15759', '#59A14F', '#4E79A7')    # 与原图/图例一致：50k 红 / 100k 绿 / HT 蓝
LBL  <- c('dsciCAT-seq-50k', 'dsciCAT-seq-100k', 'HT-scCAT-seq')
FS   <- 4.6                                   # 行内刻度字号（本室下限 4pt）
PT   <- function(n) n * .pt                     # ggplot2 3.4+ 的 size 单位是 mm，* .pt 才回到 pt

SP <- list(
  list(id = 'HEK293T', build = 'hg38', reg = 'chr13-91340000-91370000', ens = 'hs',
       chrom = 'chr13', lo = 91340000, hi = 91370000,
       frag = c('HG/HG12/sub/human_fragments.tsv.gz', 'HG/HG14/sub/human_fragments.tsv.gz',
                'scCAT-seq-HT/IL9/ATAC-sub/human_fragments.tsv.gz'),
       mtx  = c('HG/HG12/sub/ATAC/Human', 'HG/HG14/sub/ATAC/Human',
                'scCAT-seq-HT/IL9/ATAC-sub/Human')),
  list(id = 'NIH3T3', build = 'mm10', reg = 'chr11-69090000-69130000', ens = 'mm',
       chrom = 'chr11', lo = 69090000, hi = 69130000,
       frag = c('HG/HG12/sub/mouse_fragments.tsv.gz', 'HG/HG14/sub/mouse_fragments.tsv.gz',
                'scCAT-seq-HT/IL9/ATAC-sub/mouse_fragments.tsv.gz'),
       mtx  = c('HG/HG12/sub/ATAC/Mouse', 'HG/HG14/sub/ATAC/Mouse',
                'scCAT-seq-HT/IL9/ATAC-sub/Mouse'))
)

ann_for <- function(which) {
  a <- if (which == 'hs') GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
       else GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevels(a) <- paste0('chr', seqlevels(a))
  genome(a) <- if (which == 'hs') 'hg38' else 'mm10'
  a
}

CACHEO <- file.path(Sys.getenv("DSCI_FIG_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "figures")), 'scripts/cache/fig1c_objs')
dir.create(CACHEO, showWarnings = FALSE, recursive = TRUE)
load_one <- function(id, frags, mtx, which) {
  cf <- file.path(CACHEO, paste0(id, '_', which, '.rds'))
  if (file.exists(cf)) return(readRDS(cf))
  f <- as.character(read_tsv(file.path(mtx, 'features.tsv'), col_names = FALSE,
                             show_col_types = FALSE)[[1]])
  b <- as.character(read_tsv(file.path(mtx, 'barcodes.tsv'), col_names = FALSE,
                             show_col_types = FALSE)[[1]])
  X <- readMM(file.path(mtx, 'matrix.mtx'))
  if (nrow(X) == length(b)) X <- t(X)
  pk  <- GenomicRanges::GRanges(f, sep = c(':', '-'))
  sel <- which(width(pk) > 100 & width(pk) < 5000)
  ca <- CreateChromatinAssay(X[sel, , drop = FALSE] %>% set_rownames(f[sel]) %>% set_colnames(b),
                             sep = c(':', '-'), fragments = frags, min.cells = 10)
  o <- CreateSeuratObject(ca, assay = 'ATAC'); o$orig.ident <- id
  Annotation(o) <- ann_for(which)
  DefaultAssay(o) <- 'ATAC'; Idents(o) <- 'orig.ident'
  saveRDS(o, cf); o
}

theme_row <- function(x_axis = FALSE, y_axis = TRUE)
  theme_classic(base_size = FS, base_family = 'Arial') +
  theme(axis.text.x  = if (x_axis) element_text(size = FS, colour = 'black',
                                                margin = margin(t = 0.3)) else element_blank(),
        axis.ticks.x = if (x_axis) element_line(colour = 'black', linewidth = 0.15) else element_blank(),
        axis.title.x = element_blank(),
        axis.text.y  = if (y_axis) element_text(size = FS, colour = 'black') else element_blank(),
        axis.ticks.y = if (y_axis) element_line(colour = 'black', linewidth = 0.15) else element_blank(),
        axis.title.y = element_blank(), axis.line = element_line(colour = 'black', linewidth = 0.2),
        strip.text   = element_blank(), plot.margin = margin(0.6, 1.2, 0.6, 0.6),
        legend.position = 'none', panel.background = element_rect(fill = 'white', colour = NA),
        panel.grid = element_blank())

set_area_fill <- function(x, colour) {           # CoveragePlot 会自己装 fill scale，只能改 geom 参数
  for (i in seq_along(x$layers))
    if (inherits(x$layers[[i]]$geom, 'GeomArea')) x$layers[[i]]$aes_params$fill <- colour
  x
}
ymax_of <- function(x) {
  for (i in seq_along(x$layers)) {
    if (!inherits(x$layers[[i]]$geom, 'GeomArea')) next
    d <- tryCatch(ggplot2::layer_data(x, i), error = function(e) NULL)
    if (is.null(d) || !is.numeric(d$y)) next
    return(max(d$y, na.rm = TRUE))
  }
  NA
}
txt_plot <- function(...) {                       # 一个只装文字的占位 plot（左列行标签用）
  ggplot() + annotate('text', 0, 0.5, ...) + xlim(0, 1) + ylim(0, 1) +
    theme_void(base_size = FS) + theme(plot.margin = margin(0, 0, 0, 0))
}

build_species <- function(s) {
  objs <- lapply(1:3, function(i) load_one(LBL[i], file.path(DATA, s$frag[i]),
                                          file.path(DATA, s$mtx[i]), s$ens))
  cat(sprintf('  %s: cells %s\n', s$id, paste(sapply(objs, ncol), collapse = ' / ')))
  cov <- vector('list', 3); genes <- NULL
  for (i in 1:3) {
    p <- CoveragePlot(object = objs[[i]], region = s$reg, annotation = (i == 3), peaks = FALSE,
                      extend.upstream = 0, extend.downstream = 0, group.by = 'orig.ident')
    if (inherits(p, 'patchwork')) { cov[[i]] <- p[[1]]; if (i == 3) genes <- p[[2]] }
    else cov[[i]] <- set_area_fill(p, COL[i])
    cov[[i]] <- set_area_fill(cov[[i]], COL[i])
  }
  rng <- signif(max(sapply(cov, ymax_of), na.rm = TRUE) * 1.02, 3)
  brk <- seq(ceiling(s$lo / 5000) * 5000, s$hi, by = 5000)
  lab <- function(b) ifelse(b %% 10000 == 0, formatC(b / 1000, format = 'd', big.mark = ','), '')
  # 方法名直接印在每行左上角（省掉图例，行-方法一一对应）
  for (i in 1:3)
    cov[[i]] <- cov[[i]] + ylim(0, rng) + theme_row(x_axis = FALSE, y_axis = (i == 1)) +
      annotate('text', x = s$lo + 0.012 * (s$hi - s$lo), y = rng * 0.88, label = LBL[i],
               hjust = 0, vjust = 0.5, size = PT(FS * 0.95), colour = COL[i])
  genes <- genes + theme_row(x_axis = TRUE, y_axis = FALSE) +
    scale_x_continuous(breaks = brk, labels = lab, limits = c(s$lo, s$hi)) +
    theme(axis.text.x = element_text(size = FS, colour = 'black', angle = 0, hjust = 0.5))
  til <- lapply(1:3, function(i)
    TilePlot(object = objs[[i]], region = s$reg, tile.cells = NCOL, order.by = 'top',
             group.by = 'orig.ident') +
      scale_fill_gradientn(colours = c('white', 'blue', 'darkblue'), limits = c(0, 2),
                           breaks = c(0, 1, 2)) + theme_row(x_axis = FALSE, y_axis = FALSE))
  stack <- wrap_plots(cov[[1]], cov[[2]], cov[[3]], til[[1]], til[[2]], til[[3]], genes,
                      ncol = 1, heights = c(1.0, 1.0, 1.0, 0.72, 0.72, 0.72, 0.95))
  list(stack = stack, rng = rng, ncell = sapply(objs, ncol), s = s)
}

if (Sys.getenv('FIG1C_BUILD_ONLY', '') != '') {                 # 只建对象+缓存，便于快速迭代绘图
  for (sp in SP) for (i in 1:3) {
    o <- load_one(LBL[i], file.path(DATA, sp$frag[i]), file.path(DATA, sp$mtx[i]), sp$ens)
    cat(sprintf('  cached %-14s %-2s cells=%d peaks=%d annot=%d\n', LBL[i], sp$ens, ncol(o), nrow(o), length(Annotation(o)))
    ) }
  q(status = 0)
}
cat('[C] human\n'); H <- build_species(SP[[1]])
cat('[C] mouse\n'); M <- build_species(SP[[2]])

# 左列行标签（与 7 行同高）：覆盖行给共用量程，tile 行给细胞数与色标
lab_stack <- function(x) wrap_plots(
  txt_plot(label = sprintf('ATAC\ncoverage\n(0-%s)', formatC(x$rng, format = 'f', digits = 0)),
           size = PT(FS * 0.95), lineheight = 0.85),
  txt_plot(), txt_plot(),
  txt_plot(label = sprintf('single cells\n(%d / method)\ncolour 0-2', NCOL),
           size = PT(FS * 0.95), lineheight = 0.85),
  txt_plot(), txt_plot(),
  txt_plot(label = 'gene models', size = PT(FS * 0.95), lineheight = 0.85),
  ncol = 1, heights = c(1.0, 1.0, 1.0, 0.72, 0.72, 0.72, 0.95)) +
  theme(plot.margin = margin(0, 0, 0, 0))

title_of <- function(x) sprintf('%s  (%s)   %s:%s-%s', x$s$id, x$s$build, x$s$chrom,
                                formatC(x$s$lo, format = 'd', big.mark = ','),
                                formatC(x$s$hi, format = 'd', big.mark = ','))
half <- function(x) x$stack + plot_annotation(title = title_of(x),
  theme = theme(plot.title = element_text(family = 'Arial', size = 5.6, hjust = 0,
                                          margin = margin(b = 0.3))))

fig <- lab_stack(H) | half(H) | half(M)
fig <- fig + plot_layout(widths = c(0.80, 5.4, 5.4))
ggsave(OUT, fig, width = PW, height = PH, units = 'cm', dpi = 600, limitsize = FALSE)
ggsave(PNG, fig, width = PW, height = PH, units = 'cm', dpi = 500, limitsize = FALSE)
cat('wrote', OUT, '\n')
cat(sprintf('shared coverage range: human 0-%.0f  mouse 0-%.0f\n', H$rng, M$rng))
