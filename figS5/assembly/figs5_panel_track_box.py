#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
⚠ 2026-10-29（Codex #92 前半，实测纠错）：
  旧标题写 "ATAC insertions per cell state (RPM-normalized, 200bp bins)" —— **RPM 是错的**。
  上游 `grn_mef2c_nrxn3_track.py`：rpkm = 每底 counts × 1e6 / 该态 chr12 总插入；
  画图时 y = 200 底求和 × (1e6/200) ⇒ y ∝ 每底平均插入数 × 1e12/total_chr。
  即：**每行按各自的 chr12 文库总量归一，再各自 autoscale**（实测各行 y 最大 293,123–622,101，
  交付版刻度 500000/500000/200000/200000/200000/200000）——不是 RPM（RPM 量级是个位数到百）。
  同处第二处纠错：框内 "MWU p<1e-22" 与正文 "6.2 × 10⁻⁸" 不是同一个比较：
  实测（`scripts/cache/nrxn3_expr.npz`，单侧 MWU）MOL3 vs 其余五态**合并** = 6.15e-80（=正文数），
  而 1.1e-22 是**最弱的逐对**比较（vs COP/NFOL，n=203）。图内只留合并口径的 6.2e-80，
  两个口径同时写进图注。

figs5_panel_track_box.py — FigS5_noc2 的 C 面板: Nrxn3 locus ATAC track (boxplot 版)
================================================================================
内容 = 原图 <work>/pi/plots-scenic+/grn_mef2c_nrxn3_track_box.pdf
       (由 scripts-scenic+/grn_mef2c_nrxn3_track.py --style box 产出)。
原图画布 47.07 x 27.02 cm, 放进 18cm 版面须缩 0.382 倍 -> 5.5-10.5pt 文字只剩
2.1-4.0pt, 不满足投稿 6-8pt。本脚本用 **同一份数据、同一套几何** 在槽位尺寸重画,
字号全部 6-8pt (Liberation Sans = Arial 度量, pdf.fonttype=42), 并修掉原图两处排版问题:
  1) 中间各行的 x 刻度被下一行遮住, 只剩右端 "96" 碎片 -> 只在最后一行下方画一次刻度;
  2) 底部注释 (enhancer 位置说明 / 距离标注 / MEF2 motif 文字 / 图例) 在 18cm 上互相
     压字 -> 分成 4 个文字行 (E 标签 / 位置说明 x2 / 距离 x2) + 图例移到表达图下方
     竖排 4 条, 'MEF2 motif hit (…)' 由文字改为图例条目 (信息不丢)。
数据依赖: ArchR GroupCoverages h5 (12 个, ~2s, 解码后缓存) + Nrxn3 表达缓存
scripts/cache/nrxn3_expr.npz —— **不需要 scenicplus / SCENIC+ 环境**。
用法: python3 figs5_panel_track_box.py [宽cm 高cm 输出pdf]
"""
import os, sys, importlib.util
import numpy as np

sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib                                    # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.patches import FancyArrowPatch       # noqa: E402
from matplotlib.lines import Line2D                  # noqa: E402
from figkit import use_matplotlib_liberation         # noqa: E402

SRC = os.path.join(os.environ.get('DSCI_SCENIC_SCRIPTS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pi/scripts-scenic+')), 'grn_mef2c_nrxn3_track.py')
CACHE = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'scripts/cache')
EXPR_NPZ = f'{CACHE}/nrxn3_expr.npz'                  # 由 scplus_obj_final.pkl 预提取
TRACK_NPZ = f'{CACHE}/nrxn3_tracks.npz'               # ArchR h5 解码缓存
OUT = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'replot/S5c_nrxn3_track_box.pdf')
W_CM, H_CM = 18.0, 11.2

# 状态名两行折法 (轴外左侧, 宽仅 0.75cm)
STATE_LBL = {'COP_NFOL': 'COP\nNFOL', 'MFOL': 'MFOL', 'MOL2': 'MOL2',
             'MOL3': 'MOL3', 'MOL3_6': 'MOL3\n6', 'MOL6': 'MOL6'}

# ---- 版面 (cm, y 自顶向下) ----
XL, XR = 1.62, 13.75          # track 绘图区
VXL, VXR = 14.45, 17.85       # 右侧表达图
T0, TRH, TGP = 0.92, 1.12, 0.12
TBOT = T0 + 6 * TRH + 5 * TGP              # 8.24  最后一行底边
GM_Y, GM_H = 8.62, 1.20       # gene model 行
E_LBL_Y = 8.00                # E1..E4 标签 (压在最后一行 track 底部, 刻度行上方)
NOTE_ROW = {'E1': 10.05, 'E3': 10.05, 'E2': 10.32, 'E4': 10.32}
DIST_ROW = {'E1': 10.60, 'E2': 10.60, 'E3': 10.60, 'E4': 10.88}
LEG_Y = 8.75                  # 图例 (表达图下方)

_spec = importlib.util.spec_from_file_location('nrxmod', SRC)
mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mod)                           # 只取常量/函数, 不跑 __main__

CELLTYPES, COLORS, ENH = mod.CELLTYPES, mod.COLORS, mod.ENHANCERS
WIN_S, WIN_E, TSS = mod.WIN_S, mod.WIN_E, mod.TSS
EXONS, GENE_START, GENE_END = mod.NRXN3_EXONS, mod.GENE_START, mod.GENE_END
LOC_NOTE = {'E1': 'intergenic · TSS −14.8 kb',
            'E2': "NM_172544 5'UTR edge (alt. first exon)",
            'E3': 'intron (exon6-7)', 'E4': 'intron (exon7-8)'}
DISTS = {88_707_581: '-14.8 kb (TSS)', 88_952_720: '+229 kb (intron)',
         89_284_141: '+561 kb (intron)', 89_376_068: '+653 kb (intron)'}


# ---------------------------------------------------------------- data
def get_tracks():
    """6 态 chr12:88.55-89.60Mb 的 RPM coverage (h5 解码一次后缓存到 npz)。"""
    if os.path.exists(TRACK_NPZ):
        z = np.load(TRACK_NPZ)
        return {ct: z[ct].astype(float) for ct in CELLTYPES}
    tr = mod.load_tracks()
    np.savez_compressed(TRACK_NPZ, **{k: v.astype('float32') for k, v in tr.items()})
    return tr


def get_expr():
    z = np.load(EXPR_NPZ)
    return {k: z[k].astype(float) for k in z.files if k in CELLTYPES}


# ---------------------------------------------------------------- figure
def main(w_cm=W_CM, h_cm=H_CM, out=OUT):
    use_matplotlib_liberation()
    plt.rcParams.update({'axes.linewidth': 0.5, 'xtick.major.width': 0.4,
                         'ytick.major.width': 0.4, 'xtick.major.size': 1.6,
                         'ytick.major.size': 1.6, 'font.size': 6.5})
    tracks, expr, hits = get_tracks(), get_expr(), mod.load_motif_hits()
    fig = plt.figure(figsize=(w_cm / 2.54, h_cm / 2.54))

    def axrect(x0, ytop, w, h):
        return [x0 / w_cm, 1 - (ytop + h) / h_cm, w / w_cm, h / h_cm]

    def dx(bp):                                   # 基因组坐标 -> cm
        return XL + (bp - WIN_S) / (WIN_E - WIN_S) * (XR - XL)

    def fy(ycm):                                  # cm(自顶) -> 图分数
        return 1 - ycm / h_cm

    rows = [T0 + i * (TRH + TGP) for i in range(6)]
    TBOT = rows[-1] + TRH

    # ---- 6 行 coverage ----
    BIN = 200
    n_bins = (WIN_E - WIN_S) // BIN
    for i, ct in enumerate(CELLTYPES):
        a = fig.add_axes(axrect(XL, rows[i], XR - XL, TRH))
        y = tracks[ct][:n_bins * BIN].reshape(n_bins, BIN).sum(axis=1) * (1e6 / BIN)
        xs = WIN_S + BIN * (np.arange(n_bins) + 0.5)
        c = COLORS[ct]
        a.fill_between(xs, y, color=c, alpha=0.75, linewidth=0)
        a.plot(xs, y, color=c, lw=0.5)
        a.set_ylim(0, max(y.max() * 1.25, 0.02))
        a.set_xlim(WIN_S, WIN_E)
        a.tick_params(axis='y', labelsize=6, pad=1.0)
        a.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        for s, e, lab, nes in ENH:
            a.axvspan(s, e, color='#d62728', alpha=0.30, lw=0, zorder=0)
        for sp in ('right', 'top'):
            a.spines[sp].set_visible(False)
        # 状态名: 轴外左侧 (7pt, 必要时两行)
        a.text(-0.075, 0.5, STATE_LBL[ct], transform=a.transAxes, ha='right',
               va='center', fontsize=7, linespacing=1.0)
    axes = fig.axes

    # ---- x 刻度: 只画一次 (最后一行下方) ----
    axx = axes[-1]
    axx.spines['bottom'].set_visible(True)
    axx.tick_params(axis='x', which='both', bottom=True, labelbottom=True,
                    labelsize=6.5, pad=1.0)
    axx.set_xticks([88_600_000 + k * 200_000 for k in range(6)])
    axx.set_xticklabels([f'{(88.60 + 0.20 * k) / 10:.2f}' for k in range(6)])
    fig.text(XR / w_cm + 0.004, fy(TBOT + 0.30), '1e7', fontsize=6.5, ha='left',
             va='top')

    # ---- gene model ----
    ax = fig.add_axes(axrect(XL, GM_Y, XR - XL, GM_H))
    ax.axis('off'); ax.set_xlim(WIN_S, WIN_E); ax.set_ylim(-0.4, 1.8)
    (px0, _), (px1, _) = ax.transData.transform([(WIN_S, 0), (WIN_E, 0)])
    min_w = 6 / ((px1 - px0) / (WIN_E - WIN_S))
    ax.plot([max(GENE_START, WIN_S), min(GENE_END, WIN_E)], [0.35, 0.35],
            color='#b0b0b0', lw=1.0, solid_capstyle='butt', zorder=1)
    ax.add_patch(FancyArrowPatch((TSS, 0.35), (min(TSS + 8000, WIN_E), 0.35),
                                 arrowstyle='-|>', mutation_scale=9,
                                 color='#333333', lw=1.0))
    ax.text(TSS, 1.02, f'TSS {TSS:,}', fontsize=6, ha='center', color='#333333')
    ntx = len(EXONS)
    for ti, (tname, ex) in enumerate(EXONS):
        yc = 0.35 + 0.27 * (1 - 2 * ti)          # 两条转录本行 (上: NM_001198587)
        for s, e, cds in ex:
            if e < WIN_S or s > WIN_E:
                continue
            s2, e2 = max(s, WIN_S), min(e, WIN_E)
            if e2 - s2 < min_w:
                mid = (s2 + e2) / 2
                s2, e2 = mid - min_w / 2, mid + min_w / 2
            ax.plot([s2, e2], [yc, yc], color='#222222' if cds else '#cccccc',
                    lw=5.0, solid_capstyle='butt', zorder=3)
        ax.text(WIN_S + 400, yc, tname, fontsize=6, color='#666666', va='center')
    for s, e, lab, nes in ENH:                    # enhancer 顶栏 (红色短横)
        ax.plot([s, e], [1.72, 1.72], color='#d62728', lw=3)
    if hits:                                      # MEF2 motif hits (单独一行)
        for s, e, lab, nes in ENH:
            if lab not in hits:
                continue
            for k, (_, r) in enumerate(hits[lab].iterrows()):
                hs, he = max(int(r['start']), WIN_S), min(int(r['end']), WIN_E)
                ax.plot([hs, he], [1.38, 1.38],
                        color='#8b0000' if k == 0 else '#d62728',
                        lw=3.5 if k == 0 else 1.2, alpha=1 if k == 0 else 0.8,
                        solid_capstyle='butt', zorder=4)

    # ---- 底部注释 (E 标签 / 位置说明 / 距离) ----
    for s, e, lab, nes in ENH:
        mid = (s + e) / 2
        shift = {'E4': 55_000}.get(lab, 0)        # E4 右移, 与 E3 标签分开
        fig.text(dx(mid + shift) / w_cm, fy(E_LBL_Y), f'{lab}  NES={nes}', fontsize=6,
                 ha='center', va='center', color='#d62728', zorder=6,
                 bbox=dict(boxstyle='square,pad=0.12', fc='white', ec='none', alpha=0.85))
        fig.text(dx(mid) / w_cm, fy(NOTE_ROW[lab]), LOC_NOTE[lab], fontsize=6,
                 ha='center', va='center', color='#8b0000')
        fig.text(dx(s) / w_cm, fy(DIST_ROW[lab]), DISTS[s], fontsize=6,
                 ha='center', va='center', color='#d62728')

    # ---- 右侧表达图 (boxplot + jitter, 与原 _box.pdf 同画法) ----
    axv = fig.add_axes(axrect(VXL, rows[0], VXR - VXL, TBOT - rows[0]))
    lin = {ct: 2.0 ** expr[ct] - 1.0 for ct in CELLTYPES}
    datasets = [lin[ct] for ct in CELLTYPES]
    bp = axv.boxplot(datasets, positions=[5.5, 4.5, 3.5, 2.5, 1.5, 0.5], vert=False,
                     widths=0.55, patch_artist=True, showfliers=False,
                     medianprops=dict(color='#111111', lw=1.2, zorder=5),
                     whiskerprops=dict(color='#555555', lw=0.7, zorder=3),
                     capprops=dict(color='#555555', lw=0.7, zorder=3),
                     boxprops=dict(lw=0.7))
    for patch, ct in zip(bp['boxes'], CELLTYPES):
        patch.set_facecolor(COLORS[ct])
        patch.set_alpha(0.90 if ct == 'MOL3' else 0.60)
        patch.set_edgecolor('#1b7a24' if ct == 'MOL3' else '#444444')
        patch.set_linewidth(1.0 if ct == 'MOL3' else 0.7)
    rng = np.random.default_rng(42)               # 固定 seed, 与原图一致
    for yy, v, ct in zip([5.5, 4.5, 3.5, 2.5, 1.5, 0.5], datasets, CELLTYPES):
        v = np.asarray(v)
        if len(v) > 300:
            v = v[rng.choice(len(v), 300, replace=False)]
        axv.scatter(v + rng.uniform(-0.13, 0.13, len(v)), np.full(len(v), yy),
                    s=1.6, color=COLORS[ct], alpha=0.30, linewidths=0, zorder=2)
    hi = max(c.get_xdata()[0] for c in bp['caps']) * 1.06
    axv.set_xlim(0, hi); axv.set_ylim(-0.2, 6.2); axv.set_yticks([])
    axv.set_xticks([t for t in (0, 5, 10, 15, 20) if t <= hi])
    axv.tick_params(axis='x', labelsize=6, pad=1.0)
    axv.set_xlabel('Nrxn3 expr. (2ˣ−1)', fontsize=6.5, labelpad=1.5)
    axv.set_title('Nrxn3 expression\n(boxplot + jitter)', fontsize=7, loc='left',
                  pad=2.0, linespacing=1.05)
    axv.text(hi * 0.40, 5.05, 'MOL3-specific\nmedian ≈7.0\nvs ≈3.1–3.8\n≈2×; MWU p = 6.2e-80',
             fontsize=6, color='#1b7a24', ha='left', va='center', linespacing=1.15,
             bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='#1b7a24', lw=0.6))
    for sp in ('right', 'top'):
        axv.spines[sp].set_visible(False)

    # ---- 图例 (表达图下方, 竖排 4 条) ----
    handles = [Line2D([0], [0], color='#222222', lw=5.0, label='CDS'),
               Line2D([0], [0], color='#cccccc', lw=5.0, label='UTR'),
               Line2D([0], [0], color='#b0b0b0', lw=1.0, label='intron line'),
               Line2D([0], [0], color='#8b0000', lw=4, label='MEF2 motif hit')]
    axl = fig.add_axes(axrect(VXL, LEG_Y, VXR - VXL, 1.6)); axl.axis('off')
    axl.legend(handles=handles, loc='upper left', frameon=False, fontsize=6,
               handlelength=1.1, columnspacing=1.0, handletextpad=0.4,
               labelspacing=0.45, borderaxespad=0)

    # ---- 标题 (3 行 6.5pt, 左对齐) ----
    fig.text(0.004, fy(0.14),
             f'Nrxn3 locus (chr12:{WIN_S:,}-{WIN_E:,})  |  ATAC insertion density per '
             '200-bp bin; each row normalized to its own library size, drawn on its own y-scale\n'
             'Mef2c eRegulon → Nrxn3 enhancers: 4 regions (2000hvr)  ·  MEF2 motif '
             'metacluster_179.4 (JASPAR MA0052.4), NES/hits shown\n'
             'Right: Nrxn3 expression (boxplot + jitter, linearized 2ˣ−1) per state '
             '— MOL3-specific (median ≈7.0 vs ≈3.1–3.8, ≈2×)',
             fontsize=6.5, ha='left', va='top', linespacing=1.45)

    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out)                              # 不用 tight -> 画布恰为槽位尺寸
    plt.close(fig)
    print('saved', out)
    return out


if __name__ == '__main__':
    a = sys.argv[1:]
    main(float(a[0]) if len(a) > 0 else W_CM,
         float(a[1]) if len(a) > 1 else H_CM,
         a[2] if len(a) > 2 else OUT)
