#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
replot_v4.py — Fig1_v3 四张源图重绘 (B 污染率 / D RNA QC 箱线 / E ATAC 小提琴 / F 相关热图)
================================================================================
v3.1 (2026-10-12, 按反馈改版, 直接覆盖 Fig1_v3.pdf):
  新版式: A=Workflow_combined(17x7) | B=污染率(竖排两图, 方形子图, 无 RNA/ATAC 标题) + C=tracks |
          下排 = 参考图布局: 左 D=RNA QC 2x2 箱线(上下两行), 右列 E=ATAC 小提琴原图(上,
          plots-cellline 原文件直接嵌入, 不重绘) + F=相关热图(下); 行高对齐:
          D 上行 = E 小提琴行高(4.46), D 下行 = F 热图方块行高(3.00)
  与 replot_DEFG.py (v2) 的差别:
  1) 原 D 拆出为 B: **去掉 RNA/ATAC 面板标题** (用户要求); 子图改**正方形** (2.46x2.46) 并整体右移
     (左边留出 ylabel 与面板字母 B 的间距); 轴标签拆两行防溢出
  2) 原 E -> D: 8.45x8.84, 两行行高分别对齐右侧小提琴/热图 (add_axes 精确定位)
  3) 原 F -> E: 不再重绘, 直接嵌入 plots-cellline 原图 (compose 阶段处理)
  4) 原 G -> F: 8.35x4.25 (右列宽, 与小提琴同宽); 方块 3.00 cm
  字号规范不变: 刻度 6.5 / 图例 6.0 / 轴标题 7.0 / 面板内标题 8.0 pt
  方法名缩写表 ABBREV 同 v2 (图注需给出全称)

用法: python replot_v3.py
输出: <work>/figures/replot/{B,D,E,F}_v4_*.pdf
缓存: 复用 scripts/cache/corr_{HEK293T,NIH3T3}.csv (v2 已算好, 无需重算)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.io import mmread
from scipy.cluster.hierarchy import linkage, leaves_list

HZJ    = os.environ.get('FIG1_DATA', '')
HZJ_HG = f'{HZJ}/HG'
OUT    = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'replot')
CACHE  = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'scripts/cache')
os.makedirs(OUT, exist_ok=True)

# ---- 槽位尺寸 (cm) —— 唯一真源 scripts/fig1_layout.py，与 compose_fig1_v3.py 共用 --------
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig1_layout as FL
SLOT = FL.SLOT

FS_TICK, FS_LEG, FS_LBL, FS_TITLE = 6.5, 6.0, 7.0, 8.0

ABBREV = {'dsciCAT-seq-50k': 'dsci50k', 'dsciCAT-seq-100k': 'dsci100k',
          'HT-scCAT-seq': 'HT-scCAT', 'ISSAAC-seq(Droplet)': 'ISSAAC',
          'Paired-seq': 'Paired', 'sciCAR-seq': 'sciCAR',
          'SUM-seq': 'SUM', '10xMultiome': '10x'}

# 字体规范: 全图 6-8pt (最低 4pt); [2026-10-13] 改用真 Arial (figkit 从 <work>/fonts 注册)
import sys, os as _os
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from figkit import use_matplotlib_arial as _arial
_FAM = _arial()          # 真 Arial (无则 Liberation Sans) + pdf.fonttype=42
plt.rcParams.update({'savefig.dpi': 600})
# v4: 作者要求轴标签用正式的 log10 而不是 lg。真 Arial 没有 U+2080/2081 下标数字，
#     而 mathtext 的 \log 命令会强制嵌入 DejaVuSans（混字体）。
#     解法：mathtext 的 rm/it/bf/sf/tt/cal 全部指向 Arial，标签写 'log$_{10}$'
#     —— 'log' 是普通文本(Arial)，只有下标 '{10}' 走 mathtext，实测整页只嵌 ArialMT。
for _k in ('rm', 'it', 'bf', 'sf', 'tt', 'cal'):
    matplotlib.rcParams['mathtext.' + _k] = _FAM
matplotlib.rcParams['mathtext.default'] = 'regular'
LOG = 'log$_{10}$'                       # 统一从这里取，避免各处写错
# 只有「矩阵型」公开数据集（无 FASTQ，不能下采样到 50,000 reads/cell 的共同深度）——
#   正文 §2.7: "All FASTQ-based libraries were down-sampled to a common depth of 50,000 read
#   pairs per cell per modality ... whereas the matrix-only methods were analyzed at the depth
#   of their released matrices"。数据路径也印证：sciCAR/ SUM-seq corrected_matrix/ Paired_seq_Cell_Mix_matrix。
MATRIX_ONLY = {'sciCAR-seq', 'Paired-seq', 'SUM-seq'}
# sci-CAR 的 NIH3T3 ATAC 列不可比：GSE117089 的 ATAC 矩阵是对 **人 hg19 峰值集** 定量的
#   （见 SCICAR_MOUSE_ATAC_FORENSICS_20261028.md：189,603 个 hg19 峰、含 hs37d5、无鼠峰；
#    711 个鼠核里 75 个为 0，中位数 6 片段）=> 画成灰底斜纹 + ‡，不参与"谁更强"的读法。
NOT_COMPARABLE = {'sciCAR-seq'}
# 但只给**小鼠侧**加灰标：GSE117089 的 ATAC 矩阵是按人源 hg19 峰集定量的，所以 HEK293T 那一列
#   仍是真实的（只是矩阵型、非等深度，用 * 标）；只有 NIH3T3 的 sci-CAR ATAC 列不可比。
NOT_COMPARABLE_ATAC = {'NIH3T3'}

# ================================================================================
# 数据: 26 个 metrics (与 replot_DEFG.py 完全一致)
# ================================================================================
metrics_location = [
    ['ISSAAC-seq(Droplet)', 'HEK293T', 'RNA',  f'{HZJ}/ISAAC/sub/RNA-sub/human_metrics.csv'],
    ['ISSAAC-seq(Droplet)', 'HEK293T', 'ATAC', f'{HZJ}/ISAAC/sub/ATAC-sub/human_metrics.csv'],
    ['ISSAAC-seq(Droplet)', 'NIH3T3',  'RNA',  f'{HZJ}/ISAAC/sub/RNA-sub/mouse_metrics.csv'],
    ['ISSAAC-seq(Droplet)', 'NIH3T3',  'ATAC', f'{HZJ}/ISAAC/sub/ATAC-sub/mouse_metrics.csv'],
    ['10xMultiome', 'NIH3T3', 'RNA',  f'{HZJ}/10x/RNA-sub0.74/metrics.csv'],
    ['10xMultiome', 'NIH3T3', 'ATAC', f'{HZJ}/10x/ATAC_sub0.62/ATAC-sub/metrics.csv'],
    ['HT-scCAT-seq', 'HEK293T', 'RNA',  f'{HZJ}/scCAT-seq-HT/IL9/RNA-sub/Human/metrics.csv'],
    ['HT-scCAT-seq', 'HEK293T', 'ATAC', f'{HZJ}/scCAT-seq-HT/IL9/ATAC-sub/Human/metrics.csv'],
    ['HT-scCAT-seq', 'NIH3T3',  'RNA',  f'{HZJ}/scCAT-seq-HT/IL9/RNA-sub/Mouse/metrics.csv'],
    ['HT-scCAT-seq', 'NIH3T3',  'ATAC', f'{HZJ}/scCAT-seq-HT/IL9/ATAC-sub/Mouse/metrics.csv'],
    ['dsciCAT-seq-50k', 'HEK293T', 'ATAC', f'{HZJ_HG}/HG12/sub/ATAC/Human/metrics.csv'],
    ['dsciCAT-seq-50k', 'NIH3T3',  'RNA',  f'{HZJ_HG}/HG12/sub/RNA/mouse_metrics.csv'],
    ['dsciCAT-seq-50k', 'NIH3T3',  'ATAC', f'{HZJ_HG}/HG12/sub/ATAC/Mouse/metrics.csv'],
    ['dsciCAT-seq-50k', 'HEK293T', 'RNA',  f'{HZJ_HG}/HG12/sub/RNA/human_metrics.csv'],
    ['dsciCAT-seq-100k', 'HEK293T', 'ATAC', f'{HZJ_HG}/HG14/sub/ATAC/Human/metrics.csv'],
    ['dsciCAT-seq-100k', 'NIH3T3',  'RNA',  f'{HZJ_HG}/HG14/sub/RNA/mouse_metrics.csv'],
    ['dsciCAT-seq-100k', 'NIH3T3',  'ATAC', f'{HZJ_HG}/HG14/sub/ATAC/Mouse/metrics.csv'],
    ['dsciCAT-seq-100k', 'HEK293T', 'RNA',  f'{HZJ_HG}/HG14/sub/RNA/human_metrics.csv'],
    ['Paired-seq', 'HEK293T', 'RNA',  f'{HZJ}/paired-seq/Paired_seq_Cell_Mix_matrix/HEK293T/RNA/metrics.csv'],
    ['Paired-seq', 'HEK293T', 'ATAC', f'{HZJ}/paired-seq/Paired_seq_Cell_Mix_matrix/HEK293T/ATAC/metrics.csv'],
    ['sciCAR-seq', 'HEK293T', 'RNA',  f'{HZJ}/sciCAR/HEK293T/RNA/metrics.csv'],
    ['sciCAR-seq', 'HEK293T', 'ATAC', f'{HZJ}/sciCAR/HEK293T/ATAC/metrics.csv'],
    ['sciCAR-seq', 'NIH3T3',  'RNA',  f'{HZJ}/sciCAR/NIH3T3/RNA/metrics.csv'],
    ['sciCAR-seq', 'NIH3T3',  'ATAC', f'{HZJ}/sciCAR/NIH3T3/ATAC/metrics.csv'],
    ['SUM-seq', 'NIH3T3', 'RNA',  f'{HZJ}/SUM-seq/RNA/3t3filtered/corrected_matrix/metrics.csv'],
    ['SUM-seq', 'NIH3T3', 'ATAC', f'{HZJ}/SUM-seq/ATAC/3t3filtered/metrics.csv'],
]
dfs = []
for t, c, m, q in metrics_location:
    tmp = pd.read_csv(q); tmp['Method'] = t; tmp['Cell_line'] = c; tmp['Modality'] = m
    dfs.append(tmp)
cell_lines = pd.concat(dfs, ignore_index=True)
cell_lines['logC'] = np.log10(cell_lines['nCounts'] + 1)
cell_lines['logF'] = np.log10(cell_lines['nFeatures'] + 1)
print('[data] cell_lines =', cell_lines.shape)

# ================================================================================
# 面板 B: HG14 RNA / ATAC 污染率竖排两图, 3.82 x 7.05 cm
#   (v3.2: **无 RNA/ATAC 标题**; 子图正方形 2.64 并放大; 左边距 1.00 防标签遮挡;
#    上下间距加大 (gap 0.85); 允许比 C 行更高 (7.95-15.00))
#   模态可由轴标签区分 (GRCh38 RNA reads / GRCh38 ATAC fragments)
# ================================================================================
df_HG14 = pd.read_csv(f'{HZJ_HG}/HG14/sub/ISSAAC_joint_metrics_sub.csv')
sp_order, sp_color = ['Human', 'Mouse', 'Mix'], {'Human': 'red', 'Mouse': 'blue', 'Mix': 'grey'}
W, H = SLOT['B']
fig, axs = plt.subplots(2, 1, figsize=(W / 2.54, H / 2.54))
spec = [('species_rna', 'count_hs_rna', 'count_mm_rna',
         'GRCh38 RNA reads\n(\u00d710\u2074)', 'mm10 RNA reads\n(\u00d710\u2074)',
         (-5000, 100000), [0, 20000, 40000, 60000, 80000, 100000], ['0', '2', '4', '6', '8', '10']),
        ('species_atac', 'count_hs_atac', 'count_mm_atac',
         'GRCh38 ATAC fragments\n(\u00d710\u2074)', 'mm10 ATAC fragments\n(\u00d710\u2074)',
         (-10000, 150000), [0, 50000, 100000, 150000], ['0', '5', '10', '15'])]
for ax, (spc, xc, yc, xl, yl, lim, tk, tl) in zip(axs, spec):
    cnt = df_HG14[spc].value_counts(); rat = (cnt / len(df_HG14) * 100).round(2)
    for sp in sp_order:
        sub = df_HG14.query(f'{spc} == @sp')
        ax.scatter(sub[xc], sub[yc], c=sp_color[sp], ec=None, s=1.2, alpha=0.18,
                   label=f'{sp}: {cnt.get(sp,0)}\n({rat.get(sp,0)}%)')   # 两行: 方形子图内放得下
    ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_xticks(tk); ax.set_yticks(tk)
    ax.set_xticklabels(tl); ax.set_yticklabels(tl)
    ax.set_xlabel(xl, fontsize=6.5, labelpad=1.5)   # B 面板窄, 轴标签用 6.5pt 防溢出
    ax.set_ylabel(yl, fontsize=6.5, labelpad=1.0)
    ax.tick_params(labelsize=FS_TICK, length=1.6, width=0.4, pad=1.2)
    for s_ in ax.spines.values(): s_.set_linewidth(0.4)
    leg = ax.legend(loc='upper right', fontsize=FS_LEG, frameon=False, markerscale=2.6,
                    handletextpad=0.2, labelspacing=0.28, borderpad=0.1)
    leg.set_zorder(4)
fig.subplots_adjust(**FL.b_adj(W, H))     # 方形子图边距由版式模块算（行高变了不用重手调分数）
fig.savefig(f'{OUT}/B_v4_HG14_RNA_ATAC_contamination.pdf', transparent=True)
plt.close(fig); print('[B] saved (%.2fx%.2f, 方形子图 %.2f, 无 RNA/ATAC 标题)' % (W, H, FL.B_SIDE))

# ================================================================================
# 面板 D: RNA 捕获 lg(UMIs) / lg(detected genes) 2x2 箱线, 8.45 x 10.30 cm
#   v3.3 (2026-10-20 作者意见): 上行 = 下行 = H_ROW 3.50（原 4.04 / 3.00 不等高，
#   从 add_axes 的箱体高度真改，不是只换外框）；细胞系标记从上行左上角移到**下行左下角**，
#   与 E 同一约定 (transAxes (0.02,0.02) va=bottom，加粗 8pt)。几何全取 fig1_layout。
# ================================================================================
W, H = SLOT['D']
fig = plt.figure(figsize=(W / 2.54, H / 2.54))
LX, RX, AXW = FL.D_LX, FL.D_RX, FL.D_AW
AXPOS = FL.D_AXPOS
ax00, ax01, ax10, ax11 = [fig.add_axes([x / W, y / H, w / W, hh / H])
                          for x, y, w, hh in AXPOS]
rna_order_HEK293T = ['dsciCAT-seq-50k', 'dsciCAT-seq-100k', 'HT-scCAT-seq',
                     'ISSAAC-seq(Droplet)', 'Paired-seq', 'sciCAR-seq']
rna_pal_HEK293T   = ['#AF2337', '#AF2337', '#2F3C28', '#3A7CA5', '#EF476F', '#5F4BB6']
rna_order_NIH3T3  = ['dsciCAT-seq-50k', 'dsciCAT-seq-100k', 'SUM-seq', 'HT-scCAT-seq',
                     'ISSAAC-seq(Droplet)', '10xMultiome', 'sciCAR-seq']
rna_pal_NIH3T3    = ['#AF2337', '#AF2337', '#FF6B35', '#2F3C28', '#3A7CA5', '#ECC342', '#5F4BB6']
kw = dict(width=0.5, saturation=1, fliersize=0.4,
          boxprops=dict(linewidth=0.4, edgecolor='black'), medianprops=dict(linewidth=0.5),
          whiskerprops=dict(linewidth=0.4), capprops=dict(linewidth=0.4),
          flierprops=dict(marker='o', markeredgecolor='black', markersize=0.5,
                          alpha=0.5, markeredgewidth=0.4))
def star(order_):
    """v4: 矩阵型（非等深度）方法名后加 *，与图例符号对应"""
    return [m + ('*' if m in MATRIX_ONLY else '') for m in order_]


def lighten(hexc, f=0.55):
    """向白插值 f，用于把“不参与比较”那一列的提琴体淡下去（斜纹要看得清，底色就不能深）"""
    c = [int(hexc.lstrip('#')[i:i + 2], 16) for i in (0, 2, 4)]
    return '#%02X%02X%02X' % tuple(round(v + (255 - v) * f) for v in c)


def darken(hexc, f=0.35):
    """向黑插值 f —— hatchsame 样式用：斜纹/描边换成**同色系深色**，色相不变，颜色键不破。"""
    c = [int(hexc.lstrip('#')[i:i + 2], 16) for i in (0, 2, 4)]
    return '#%02X%02X%02X' % tuple(round(v * (1 - f)) for v in c)

def d_plot(a, cl, order_, pal, ycol):
    sns.boxplot(data=cell_lines.query(f'Modality == "RNA" & Cell_line == "{cl}"'),
                x='Method', y=ycol, ax=a, palette=pal, order=order_, **kw)
    # v4.1（作者 2026-10-28：“有网格感觉不好看”）：**去掉 D 的横向网格**。
    #   原 v4 加网格是为了跨列读数，但箱线本来就有 y 刻度 + 上下两行共用同一把尺，网格属多余。
    a.tick_params(axis='both', labelsize=FS_TICK, length=1.6, width=0.4, pad=1.2)
    for s_ in a.spines.values(): s_.set_linewidth(0.4)
    return a


d_plot(ax00, 'HEK293T', rna_order_HEK293T, rna_pal_HEK293T, 'logC')
ax00.set_ylim(top=FL.D_YLIM_TOP)   # 上行不再放标题，仍留 0.5 头位防外点贴框
ax00.set_ylabel(LOG + '(UMIs)', fontsize=FS_LBL, labelpad=1.0)
ax00.set_yticks([0, 1, 2, 3, 4, 5])
ax00.set_xlabel(''); ax00.set_xticks([])
# v3.3: 原上行左上角的 'HEK293T' 已下移到 ax10 左下角（与 E 统一）
d_plot(ax10, 'HEK293T', rna_order_HEK293T, rna_pal_HEK293T, 'logF')
ax10.set_ylim(top=FL.D_YLIM_BOT)
ax10.set_ylabel(LOG + '(detected genes)', fontsize=FS_LBL, labelpad=1.0)
ax10.set_yticks([0, 1, 2, 3, 4])
ax10.set_xlabel('')
ax10.text(0.02, 0.02, 'HEK293T', transform=ax10.transAxes, fontsize=FS_TITLE,
          fontweight='bold', va='bottom', ha='left')      # v3.3 与 E 同约定：左下角
ax10.set_xticklabels(star(rna_order_HEK293T), rotation=90, fontsize=6.0)   # v4: 全称 + 矩阵型 *
for _t in ax10.get_xticklabels():       # 锚点修正: 全称较长, 从轴下缘向下排, 不捅出图底
    _t.set_rotation_mode('anchor'); _t.set_ha('right')
d_plot(ax01, 'NIH3T3', rna_order_NIH3T3, rna_pal_NIH3T3, 'logC')
ax01.set_ylim(top=5.5)
ax01.set_ylabel(''); ax01.set_yticks([0, 1, 2, 3, 4, 5])
ax01.tick_params(axis='y', which='both', labelleft=False, left=False)
ax01.set_xlabel(''); ax01.set_xticks([])
# v3.3: 原上行左上角的 'NIH3T3' 已下移到 ax11 左下角（与 E 统一）
d_plot(ax11, 'NIH3T3', rna_order_NIH3T3, rna_pal_NIH3T3, 'logF')
ax11.set_ylim(top=FL.D_YLIM_BOT)
ax11.set_ylabel(''); ax11.set_yticks([0, 1, 2, 3, 4])
ax11.tick_params(axis='y', which='both', labelleft=False, left=False)
ax11.set_xlabel('')
ax11.text(0.02, 0.02, 'NIH3T3', transform=ax11.transAxes, fontsize=FS_TITLE,
          fontweight='bold', va='bottom', ha='left')       # v3.3 与 E 同约定：左下角
ax11.set_xticklabels(star(rna_order_NIH3T3), rotation=90, fontsize=6.0)     # v4: 全称 + 矩阵型 *
for _t in ax11.get_xticklabels():
    _t.set_rotation_mode('anchor'); _t.set_ha('right')
fig.savefig(f'{OUT}/D_v4_RNA_umis_genes.pdf', transparent=True)
plt.close(fig); print('[D] saved (%.2fx%.2f, 上下行等高 %.2f, 细胞系标记在下行左下角)' % (W, H, FL.H_ROW))

# ================================================================================
# 面板 E: ATAC lg(unique fragments) 小提琴 1x2 —— v3.3 改为**自绘**
#   原来直接嵌 plots-cellline 的 22x12 原图（缩放 0.437），高度被原图白边钉死，无法与 D 上行对齐。
#   数据 / 方法顺序 / 配色与 notebook `cellline_QC_v2_no_scifi_dsciATAC_no_cellbar.ipynb` 的
#   ATAC 段逐项一致（同一 26 份 metrics.csv，logC = lg(nCounts+1)，已实测 HEK293T 最大 5.504 /
#   NIH3T3 最大 5.586，与原图尾巴一致）。上下零边距 => 绘图框整高 = H_ROW = D 上行。
# ================================================================================
atac_order_HEK293T, atac_order_NIH3T3 = rna_order_HEK293T, rna_order_NIH3T3
atac_pal_HEK293T, atac_pal_NIH3T3 = rna_pal_HEK293T, rna_pal_NIH3T3   # notebook 里两模态 order/palette 相同
W, H = SLOT['E']
e_fig = plt.figure(figsize=(W / 2.54, H / 2.54))
e_aw = (W - FL.E_LBL_W - FL.E_RGT - FL.E_GAP) / 2
e_kw = dict(width=0.60, inner='box', linewidth=0.4, linecolor='black',
            density_norm='area')     # seaborn 0.13: 旧 scale='area' 现在是 density_norm，写死才与原图同口径
# ↑ **width 0.8→0.60 的原因**（作者：“图小了之后提琴难看”）：交付版 E 是把 19.12x10.20 cm 的
#   原图缩到 4.46 cm 高，提琴最大宽实测 0.42-0.46 cm；新版绘图框只有 3.50 cm 高，
#   沿用 width=0.8 提琴最大宽仍为 0.48 cm —— **同样宽、矮了 18%，视觉上就胖成一团**。
#   按高度比例 3.50/4.46=0.78 缩到 width≈0.62（实测 0.60 时最大宽 0.36 cm）就恢复原来的瘦长观感。
#   历史坑（2026-09/10 会话）：`bw_method` 是 Scott 带宽的**倍数因子**，0.5 会更碎不是更平滑；
#   本数据用默认 scott 已经平滑，实测 bw_method=1.5/2.0 反而把提琴撑成圆棒，故**不动 bw**。
for _i, (_cl, _order, _pal) in enumerate([('HEK293T', atac_order_HEK293T, atac_pal_HEK293T),
                                          ('NIH3T3', atac_order_NIH3T3, atac_pal_NIH3T3)]):
    ea = e_fig.add_axes([(FL.E_LBL_W + _i * (e_aw + FL.E_GAP)) / W, FL.E_BOT / H,
                         e_aw / W, FL.H_ROW / H])
    sns.violinplot(data=cell_lines.query('Modality == "ATAC" & Cell_line == @_cl'),
                   x='Method', y='logC', ax=ea, palette=_pal, order=_order, **e_kw)
    # 0.13 的 violinplot 不再收 boxprops/medianprops/flierprops，事后统一成 0.3 细线（与原图一致）
    from matplotlib.collections import LineCollection as _LC
    for _lc in ea.findobj(_LC):
        _lc.set_linewidth(0.3)
    ea.set_ylim(FL.E_YLIM_BOT, FL.E_YLIM_TOP)   # 数据最大 5.586 -> 顶 5.9；底 -0.30 让 cut=2 的 KDE 下尾留在框内
    ea.set_yticks([0, 1, 2, 3, 4, 5])
    ea.set_xlabel(''); ea.set_xticks([])
    ea.tick_params(axis='both', labelsize=FS_TICK, length=1.6, width=0.4, pad=1.2)
    ea.text(0.02, 0.02, _cl, transform=ea.transAxes, fontsize=FS_TITLE,
            fontweight='bold', va='bottom', ha='left')      # 左下角：D 现在同此约定
    for s_ in ea.spines.values(): s_.set_linewidth(0.4)
    # ---- v4.1 (Codex #222 + sci-CAR 取证 + 作者 2026-10-28“sciCAR 的小提琴怎么这么奇怪”):
    #   v4 旧做法是在整列上盖一个**顶天立地的白纱+斜纹矩形**，视觉上变成“一个灰色长条”，
    #   把提琴形状完全遮住，而且贴着右脊柱、看着像坐标轴的一部分 —— 就是作者看到的“奇怪”。
    #   现在改成**只给提琴本体**上斜纹 + 淡色：seaborn 0.13 的 violinplot 每列就是一个
    #   PolyCollection（实测 ax.collections 按类别顺序一列一个、每个 1 条 path），
    #   所以直接 set_facecolor/set_hatch 即可，不需要覆盖层（旧注释说“collections 为空”
    #   是在另一个旧环境里量的，本环境 0.13.2+mpl3.9 不是这样）。
    # ---- v4.2：不可比列的**样式开关**（作者 2026-10-28：“还是太突兀、颜色不对”）
    #   背景：颜色是跨 D/E/F 的方法键（sciCAR = #5F4BB6 紫），把这一列改淡紫/灰色**会破掉颜色编码**，
    #   这就是“颜色不对”的根。所以候选样式都围绕“**不动色相，只降视觉权重**”。
    #   FIG1E_SCICAR_STYLE = plain | fade | hatchsame | pale | omit（默认 fade）
    #   取证结论见 SCICAR_MOUSE_ATAC_FORENSICS_20261028.md：数据没用错，是 GSE117089 的 ATAC
    #   矩阵在**仅人源 hg19 峰集**（189,603 峰、含 hs37d5、无鼠峰）上定量，鼠核无峰可数
    #   → 这个数度量的是参考集属性，不是 sci-CAR 性能，所以不能与 dsciCAT/10x 并排读。
    # ---- v4.3（作者 2026-10-28：“颜色还是不对，干脆用回最初版 fig1_v3 的，数据没弄错就直接用原版”）：
    #   默认改为 **plain** —— 提琴完全不动（与 E_v3 一致），只保留顶端一个灰色 ‡ 作为图注锚点。
    #   理由：取证已排除“用错数据”（见 SCICAR_MOUSE_ATAC_FORENSICS_20261028.md），剩下的只是
    #   “不能当性能读”，这属于**图注层面**要说清的事（取证报告 §4 当初的建议也是“图不动、图注如实标”）；
    #   而颜色是跨 D/E/F 的方法键，改色相会破颜色编码，淡色/斜纹又都抢戏。
    #   一个 ‡ 都不想要：FIG1E_SCICAR_STYLE=none（与 E_v3 逐像素一致）。
    #   要回到上一版淡色：FIG1E_SCICAR_STYLE=fade；其余 hatchsame / pale / omit 见 §11 表。
    if 'sciCAR-seq' in _order and _cl in NOT_COMPARABLE_ATAC:
        _j = _order.index('sciCAR-seq')
        _st = os.environ.get('FIG1E_SCICAR_STYLE', 'plain')
        _pc = ea.collections[_j]
        _top = _pc.get_paths()[0].vertices[:, 1].max()
        matplotlib.rcParams['hatch.linewidth'] = 0.45
        if _st == 'none':
            pass                                        # 连 ‡ 也不画（下面 _top 用不到）
        elif _st == 'plain':
            pass                                        # 只放 ‡，提琴完全不动
        elif _st == 'fade':
            _pc.set_alpha(0.38)                         # 保持原色相，只淡下去
            # 内部 inner='box' 画了 3 条 Line2D（ whisker / IQR 盒 / 白色中位线），
            #  whisker 延伸到 3.4、还是全黑，会从淡提琴里“戳”出来 → 只把它一起淡掉，
            #  盒体与中位线不动（否则白中位线一淡就看不见，盒变实心黑块）。
            _ls = [l for l in ea.lines if abs(float(np.mean(l.get_xdata())) - _j) < 1e-6]
            if _ls:
                max(_ls, key=lambda l: np.ptp(l.get_ydata())).set_alpha(0.45)
        elif _st == 'hatchsame':
            _pc.set_hatch('///'); _pc.set_edgecolor(darken(_pal[_j], 0.35))
            _pc.set_linewidth(0.5)                      # 斜纹用**同色系深色**，不改色相
        elif _st == 'pale':
            _pc.set_facecolor(lighten(_pal[_j], 0.62))  # v4.1 旧样式（作者否了：颜色不对）
            _pc.set_edgecolor('#7A7A7A'); _pc.set_hatch('///'); _pc.set_linewidth(0.5)
        elif _st == 'omit':
            _pc.set_visible(False)                       # 不画提琴（“不可数”=没有分布），只留一个中位刻线
            _med = np.median(cell_lines.query('Modality == "ATAC" & Cell_line == @_cl & Method == "sciCAR-seq"').logC)
            ea.plot([_j - 0.19, _j + 0.19], [_med, _med], color=_pal[_j], lw=1.0, alpha=0.9, zorder=6)
            _top = _med
        if _st != 'none':
            ea.text(_j, _top + 0.18, '\u2021', ha='center', va='bottom',
                    fontsize=FS_TICK + 1.0, color='#4D4D4D', zorder=7)
            # ---- v4.4（作者 2026-10-29：“如果真的有问题，图里也标一下，不能只在图注写”）：
            #   ‡ 单放在图里是无定义符号，投稿惯例必须就地解释。右半绘图框底部那条
            #   y < 2.2 的空带（除最右 sciCAR 列外，其余方法最小值 >= 2.64）正好放两行小字，
            #   字号沿用 C 面板的 note 口径（3.2 pt）。字体只用 ASCII，避免 mathtext 回退到非 Arial。
            if os.environ.get('FIG1E_SCICAR_NOTE', 'on') != 'off':
                ea.text(0.345, 0.028,
                        '\u2021 not comparable: the GSE117089\n'
                        'ATAC matrix is quantified on a\n'
                        'human-only hg19 peak set',
                        transform=ea.transAxes, ha='left', va='bottom', fontsize=3.0,
                        color='#4D4D4D', linespacing=1.30, zorder=7)
        print('[E] %s: sciCAR ATAC 列样式 = %s' % (_cl, _st))
    if _i == 0:
        ea.set_ylabel(LOG + '(unique fragments)', fontsize=FS_LBL, labelpad=1.0)
        ea.tick_params(axis='y', which='both', labelright=False, right=False)
    else:
        # v4.1（作者：“右边加了个多余的轴，看着和左边差不多，不能左右统一坐标轴么”）：
        #   两半本来就是**同一 ylim、同一组刻度**（FL.E_YLIM_BOT/TOP），v4 又在右半右脊柱上
        #   复制了一套一模一样的数字 → 冗余。现在右半不画任何 y 刻度，共用左轴读数（同 D 右列）。
        ea.set_ylabel('')
        ea.tick_params(axis='y', which='both', labelleft=False, left=False,
                       labelright=False, right=False)
e_fig.savefig(f'{OUT}/E_v4_ATAC_fragments.pdf', transparent=True)
plt.close(e_fig)
print('[E] saved (%.2fx%.2f, 自绘小提琴, 绘图框 %.2f = D 上行, 标记左下角)' % (W, H, FL.H_ROW))

# ================================================================================
# (历史) 面板 E v3.1/v3.2：直接嵌 plots-cellline/ATAC_fragments_no_scifi_dsciATAC_nocellbar.pdf
#   (19.12x10.20, 缩放 0.437 进 8.35x4.46 槽位)。**v3.3 已改为上面的自绘块** —— 原图白边
#   把高度钉死在 4.46，无法做到「E 总高 = D 上行」；数据口径未变（同一批 metrics.csv）。
# ================================================================================

# ================================================================================
# 面板 F: RNA similarity 方法间相关热图 x2 (HEK293T + NIH3T3), 8.35 x 4.25 cm
#   (原 G 内容; v3.1 右列宽 = 小提琴同宽 8.35, 高 4.35 = 标题0.40+方块3.00+标签0.95)
# ================================================================================
methods_cfg = {
    'dsciCAT-seq-50k':   (f'{HZJ_HG}/HG12/sub/RNA/{{sp}}', False),
    'dsciCAT-seq-100k':  (f'{HZJ_HG}/HG14/sub/RNA/{{sp}}', False),
    'HT-scCAT-seq':      (f'{HZJ}/scCAT-seq-HT/IL9/RNA-sub/{{sp}}', False),
    'ISSAAC-seq(Droplet)': (f'{HZJ}/ISAAC/sub/RNA-sub/{{sp}}', False),
    'sciCAR-seq':        (f'{HZJ}/sciCAR/{{line}}/RNA', True),
    'SUM-seq':           (f'{HZJ}/SUM-seq/RNA/3t3filtered/corrected_matrix', False),
    'Paired-seq':        (f'{HZJ}/paired-seq/Paired_seq_Cell_Mix_matrix/HEK293T/RNA', False),
    # 2026-10-21 作者补漏：10xMultiome 在 D/E 里有（rna_order_NIH3T3 第 6 项、色 #ECC342），
    # 但 F 的 methods_cfg 从来没登记它，所以 NIH3T3 侧一直只有 6×6，与 D/E 对不上。
    # 矩阵 = cellranger-arc 的 outs/Solo.out/Gene/filtered（32285 基因 × 6603 细胞），与 D/E 用的
    # RNA-sub0.74/metrics.csv 属同一子集（barcodes 6603 == metrics 行数 6603，已核）；
    # mtx 是基因×细胞，故 genes_x_cells=True（同 sciCAR-seq）。features.tsv 第 0 列是 ENSMUSG，
    # 与其它小鼠侧方法同规范 -> 加入后公共基因集 29621 不变，原 6×6 数值一字不改。
    '10xMultiome':       (f'{HZJ}/10x/RNA-sub0.74/outs/Solo.out/Gene/filtered', True),
}
method_set = {'HEK293T': ['dsciCAT-seq-50k', 'dsciCAT-seq-100k', 'HT-scCAT-seq',
                          'ISSAAC-seq(Droplet)', 'Paired-seq', 'sciCAR-seq'],
              'NIH3T3':  ['dsciCAT-seq-50k', 'dsciCAT-seq-100k', 'SUM-seq', 'HT-scCAT-seq',
                          'ISSAAC-seq(Droplet)', '10xMultiome', 'sciCAR-seq']}   # 7 种，与 D/E 的 rna_order_NIH3T3 同序同集合

def load_gene_sums(prefix, genes_x_cells=False):
    mtx = mmread(os.path.join(prefix, 'matrix.mtx')).tocsr()
    var = pd.read_csv(os.path.join(prefix, 'features.tsv'), sep='\t', header=None, index_col=0)
    return pd.Series(np.asarray(mtx.sum(axis=1 if genes_x_cells else 0)).ravel(), index=var.index)

def get_methods(cell_line):
    sp = 'Human' if cell_line == 'HEK293T' else 'Mouse'
    return [(m, methods_cfg[m][0].format(sp=sp, line=cell_line), methods_cfg[m][1])
            for m in method_set[cell_line]]

corr_r_store = {}
for cell_line in ['HEK293T', 'NIH3T3']:
    cpath = f'{CACHE}/corr_{cell_line}.csv'
    if os.path.exists(cpath):
        cached = pd.read_csv(cpath, index_col=0)
        # 缓存口径自检：形状与行列名必须与当前 method_set 完全一致才用。
        # （2026-10-21 补 NIH3T3 的 10xMultiome 时，旧 6×6 缓存就在同一个路径上，
        #   不自检就会被静默复用、新列白补。备份留作 corr_NIH3T3.csv.bak_6x6_20261021。）
        if (cached.shape == (len(method_set[cell_line]),) * 2
                and list(cached.index) == method_set[cell_line]
                and list(cached.columns) == method_set[cell_line]):
            corr_r_store[cell_line] = cached
            continue
        print(f'[F] {cell_line}: 缓存口径不符 {cached.shape} vs 当前 {len(method_set[cell_line])} 法 -> 重算')
    per = [load_gene_sums(p, g).rename(m) for m, p, g in get_methods(cell_line)]
    common = sorted(set.intersection(*map(set, [s.index for s in per])))
    counts = np.log2(pd.concat([s[common] for s in per], axis=1) + 1)
    normalized = np.log2(counts / counts.sum(axis=0) * 1_000_000 + 1)
    c = normalized.corr(); c.to_csv(cpath)
    corr_r_store[cell_line] = c
    print(f'[F] {cell_line}: common genes {len(common)}')

W, H = SLOT['F']
SIDE, X0 = FL.F_SIDE, FL.F_X0
LBLH, TITH, CBAND = FL.F_LBLH, FL.F_TITH, FL.F_CBAND
SW_D, SW_W, BAR_H = FL.F_SW_D, FL.F_SW_W, FL.F_BAR_H
SW_H, SW_GAP = FL.F_SWH, FL.F_SWD2            # 列色块条（贴方块下缘，旋转全称排在它下面）
Y0 = LBLH                                     # 方块下缘（自图底）
F_TITLE_DROP = 0.09                           # 标题下移量（cm）—— 见下面 g_fig.text 处的实测注释
CB_Y = Y0 + SIDE + FL.F_BAR_GAP               # 横向渐变条下缘
# 方法 -> 颜色：直接取 D/E 图例定义过的那套 palette（作者要求"色块是 de 图例定义过的颜色"）
M2C = dict(zip(rna_order_HEK293T, rna_pal_HEK293T))
M2C.update(dict(zip(rna_order_NIH3T3, rna_pal_NIH3T3)))
g_fig = plt.figure(figsize=(W / 2.54, H / 2.54))
norm = matplotlib.colors.Normalize(vmin=0.2, vmax=1.0)
sm = matplotlib.cm.ScalarMappable(cmap='OrRd', norm=norm); sm.set_array([])
LAB_W = 1.22                                  # 'correlation' 文字宽（6pt 实测 ≈1.16）
for i, cell_line in enumerate(['HEK293T', 'NIH3T3']):
    c = corr_r_store[cell_line]
    order_ = leaves_list(linkage(c.values, method='average'))
    names = [method_set[cell_line][j] for j in order_]
    ah = g_fig.add_axes([X0[i] / W, Y0 / H, SIDE / W, SIDE / H])
    sns.heatmap(c.iloc[order_, order_], ax=ah, cmap='OrRd', vmin=0.2, vmax=1, cbar=False,
                xticklabels=star(names), yticklabels=False, linewidths=0)   # 底部旋转 90 全称；矩阵型加 *
    # 旋转全称整体下移，给列色块条让位（pad 单位 pt：0.20 cm = 5.67 pt）
    ah.tick_params(axis='x', labelsize=6.0, length=0, pad=1.0 + (SW_GAP + SW_H) * 28.346,
                   rotation=90, bottom=False)
    ah.tick_params(axis='y', length=0, left=False, labelleft=False)
    for s_ in ah.spines.values(): s_.set_linewidth(0.4)
    # 每行右侧的方法色块（颜色 = D/E 图例；行序跟随聚类后的 names）
    sax = g_fig.add_axes([(X0[i] + SIDE + SW_D) / W, Y0 / H, SW_W / W, SIDE / H])
    sax.imshow(np.arange(len(names))[:, None], aspect='auto',
               cmap=matplotlib.colors.ListedColormap([M2C[n] for n in names]),
               vmin=-0.5, vmax=len(names) - 0.5)
    sax.set_xticks([]); sax.set_yticks([])
    for s_ in sax.spines.values(): s_.set_visible(False)
    # 列下方色块条：与列一一对应，紧接其下就是该方法的全称
    kax = g_fig.add_axes([X0[i] / W, (Y0 - SW_GAP - SW_H) / H, SIDE / W, SW_H / H])
    kax.imshow(np.arange(len(names))[None, :], aspect='auto',
               cmap=matplotlib.colors.ListedColormap([M2C[n] for n in names]),
               vmin=-0.5, vmax=len(names) - 0.5)
    kax.set_xticks([]); kax.set_yticks([])
    for s_ in kax.spines.values(): s_.set_visible(False)
    # 色条移到热图**上方**（横条）：左 'correlation'，刻度只标两端 0.2 / 1.0
    cbx = g_fig.add_axes([(X0[i] + LAB_W) / W, CB_Y / H, (SIDE - LAB_W) / W, BAR_H / H])
    cb = plt.colorbar(sm, cax=cbx, orientation='horizontal')
    cb.outline.set_linewidth(0.3)                     # 参考样式：细边（默认 1pt 太重）
    cbx.tick_params(axis='x', labelsize=FS_TICK, length=1.4, width=0.4, pad=1.2,
                    top=True, bottom=False, labeltop=True, labelbottom=False)
    cbx.set_xticks([0.2, 1.0])
    for s_ in cbx.spines.values(): s_.set_linewidth(0.4)
    cbx.text(-0.02, 0.5, 'correlation', transform=cbx.transAxes, ha='right', va='center',
             fontsize=FS_LEG)
    # 标题下移 F_TITLE_DROP = 0.09 cm（2026-10-21 作者指出“HEK293T / NIH3T3 上面被截掉了一点”）。
    # 原来的 +0.07 只把**基线**放到距顶 0.19 cm，而 8pt Arial-Bold 大写高 0.202 cm，
    # 于是墨顶 = -0.0094 cm，被面板自己的上页边切平（1200 dpi 实测第 0 行就有墨）。
    # 下移后：上留白 0.078 cm、下留白（到 0.2/1.0 刻度墨顶 0.3579）0.075 cm，上下基本均衡。
    g_fig.text(X0[i] / W, (H - TITH + 0.07 - F_TITLE_DROP) / H, cell_line, fontsize=FS_TITLE,
               fontweight='bold', va='baseline', ha='left')
g_fig.savefig(f'{OUT}/F_v4_RNA_similarity_heatmap.pdf', transparent=True)
plt.close(g_fig); print('[F] saved (%.2fx%.2f, 方块 %.2f 底边对齐 D 下行, 色条上移 + 行/列色块)' % (W, H, SIDE))
