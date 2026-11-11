import os
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""s3c_cop_figure_v9.py —— 100k COP/NFOL 验证图的新三格版式（作者定）
=========================================================================================
作者看后口径：
  * 得先看到"两簇长什么样" ⇒ 最左放 **50k 与 100k 的 COP/NFOL 再聚类 UMAP**，
    标出 sub-0/sub-1，并写清是**同一套参数**（实测：50k 60 核 → 1 簇；100k 133 核 → 81/52，
    与交付分区完全一致，见 s3c_cop_umap_pair.R 的输出与下面的 ARI 断言）。
  * 中间放 marker 气泡图，**去掉 Oligo 那列**（作者要求；Oligo 的均值改由图注给）。
  * 右边放稳定性 J。
  * 原 K（技术混杂 6 行 lollipop）作者说没信息量 ⇒ 撤图，**但两条数不能丢**：
    中位 nCount_RNA 2,682 vs 4,292 直接标在 UMAP 那格的注释里（它就是"只作描述"的理由）。
  * 原 I（深度匹配重检）作者问"这是啥" ⇒ 撤图，结论压成图注一句：306 → 251 个基因仍合格、
    方向 100% 不变。
判据条（FROZEN RULE…）仍然进图注（compact 口径）。
输出 replot/S3c/FigS3c_v9{,_HIJK}.pdf，18.00 × 7.75 cm，供 merge 到 S3 页底部。
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'scripts'))
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, Normalize
import figkit as K
K.use_matplotlib_arial()      # ⚠ 少了这一行 matplotlib 会静默用 DejaVu Sans（宽 7~10%、高 1.7em），
                              #   轴标/题注全变成非 Arial —— 本轮就先踩过一次，靠 inkcheck 盒高
                              #   从 6.1 pt 涨到 9.3 pt 才暴露。

D = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'u13_cop_4k/s3c')
OUTDIR = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'replot/S3c')
os.makedirs(OUTDIR, exist_ok=True)
W, H = 18.00, 7.75
RED, GREY, INK = (0.835, 0.122, 0.149), (0.55, 0.55, 0.55), (0.13, 0.13, 0.13)
FAINT, BLUE, LGTEAL, PALE = (0.80, 0.80, 0.80), (0.20, 0.42, 0.70), (0.13, 0.47, 0.31), (0.66, 0.77, 0.90)
FS_N = 5.5
LETTER_Y, PLOT_TOP, PLOT_H = 0.36, 1.44, 4.73   # 顶部 1.44 cm 让给 J 的图例（作者：别写文字，用图例）
SUF = '_HIJK' if '--letters=HIJK' in sys.argv else ''
LET = dict(zip('HIJ', 'HIJ')) if SUF else dict(zip('HIJ', 'HIJ'))   # 独立编号 H/I/J；并进 S3 时同名

# 三格位置：x = 格最左沿, gut = 行名沟槽, w = 绘图区宽
CELLS = {'H': dict(x=0.42, gut=0.00, w=5.55),
         'I': dict(x=6.35, gut=1.45, w=2.95),
         'J': dict(x=11.45, gut=1.15, w=4.90)}
RIGHT = {c: CELLS[c]['x'] + CELLS[c]['gut'] + CELLS[c]['w'] for c in CELLS}
assert RIGHT['H'] < CELLS['I']['x'] and RIGHT['I'] < CELLS['J']['x'] and RIGHT['J'] < W - 0.30, RIGHT

um = pd.read_csv(f'{D}/s3c_cop_umap_pair.csv')
mk = pd.read_csv(f'{D}/s3c_marker_table.csv')
stab = pd.read_csv(f'{D}/s3c_stability_scored.csv')
S2 = json.load(open(f'{D}/s3c_S2.json'))
qc = pd.read_csv(f'{D}/s3c_qc.csv')

# ---- 上版式的前提断言：同参数新跑的分区必须与交付分区一致，否则 UMAP 不能这么画 -------
from sklearn.metrics import adjusted_rand_score
for lib in um.lib.unique():
    d = um[um.lib == lib]
    a = adjusted_rand_score(d.delivered.astype(str), d.fresh.astype(str))
    print('  %s: %d 核, 交付 %d 簇 / 同参数新跑 %d 簇, ARI=%.3f'
          % (lib, len(d), d.delivered.nunique(), d.fresh.nunique(), a))
    if a < 0.999:
        raise RuntimeError('同参数再聚类没复现交付分区（%s ARI=%.3f）——UMAP 格不能这么画' % (lib, a))

fig = plt.figure(figsize=(W / 2.54, H / 2.54))


def header(c, title):
    g = CELLS[c]
    fig.text(g['x'] / W, (H - LETTER_Y) / H, c, fontsize=8.0, family=K.FAM, color=INK,
             ha='left', va='baseline')
    x1 = g['x'] + K.text_width(c, 8.0) + 0.20
    nxt = {'H': CELLS['I']['x'], 'I': CELLS['J']['x'], 'J': W - 0.40}[c]
    w = K.text_width(title, FS_N + 0.7)
    if x1 + w > nxt - 0.08:
        raise RuntimeError('标题越格 %s: 右沿 %.2f，隔壁在 %.2f（"%s"）' % (c, x1 + w, nxt, title))
    fig.text(x1 / W, (H - LETTER_Y) / H, title, fontsize=FS_N + 0.7, family=K.FAM, color=INK,
             ha='left', va='baseline')


def ax_cell(c, y0=PLOT_TOP, h=PLOT_H, x=None, w=None):
    g = CELLS[c]
    x = g['x'] + g['gut'] if x is None else x
    w = g['w'] if w is None else w
    return fig.add_axes([x / W, (H - y0 - h) / H, w / W, h / H])


def clean(ax, left=True, bottom=True):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.spines['left'].set_visible(left)
    ax.spines['bottom'].set_visible(bottom)


# ============================================================================ H 再聚类 UMAP 对
header('H', 're-clustering, identical parameters')
LIM = (min(um.UMAP_1) - 0.5, max(um.UMAP_1) + 0.5, min(um.UMAP_2) - 0.5, max(um.UMAP_2) + 0.5)
AW, GAPX = 2.30, 0.45
titles = {'fy13_50k': '50k COP/NFOL', 'fy14_100k': '100k COP/NFOL'}
for k, lib in enumerate(['fy13_50k', 'fy14_100k']):
    d = um[um.lib == lib]
    ax = ax_cell('H', x=CELLS['H']['x'] + k * (AW + GAPX), w=AW)
    ax.set_xlim(LIM[0], LIM[1]); ax.set_ylim(LIM[2], LIM[3])
    ax.set_aspect('equal'); ax.axis('off')
    for cl, col in zip(sorted(d.delivered.unique()), [BLUE, RED]):
        s = d[d.delivered == cl]
        ax.scatter(s.UMAP_1, s.UMAP_2, s=2.6, color=[col], edgecolors='none', zorder=3)
        if len(d.delivered.unique()) > 1:            # 两簇才标名字，单簇标在下方
            cx, cy = s.UMAP_1.mean(), s.UMAP_2.mean()
            ax.text(cx, cy, 'sub-%d' % cl, fontsize=FS_N, color='white', zorder=4,
                    ha='center', va='center', bbox=dict(fc=col, ec='none', pad=0.6))
    fig.text((CELLS['H']['x'] + k * (AW + GAPX) + AW / 2) / W, (H - 0.72) / H, titles[lib],
             fontsize=FS_N + 0.7, color=INK, ha='center', va='baseline')
    fig.text((CELLS['H']['x'] + k * (AW + GAPX) + AW / 2) / W, (H - PLOT_TOP - PLOT_H - 0.22) / H,
             '%d nuclei \u2192 %d cluster%s' % (len(d), d.delivered.nunique(),
                                               '' if d.delivered.nunique() == 1 else 's'),
             fontsize=FS_N, color=INK, ha='center', va='top')
med = qc[qc['var'].isin(['nCount_RNA', 'nFeature_RNA'])].set_index('var')
n0, n1 = med.loc['nCount_RNA', ['med0', 'med1']]
nf0, nf1 = med.loc['nFeature_RNA', ['med0', 'med1']]
# 作者：核心参数可以写在 UMAP 下面，但原来那两行 5.5 pt 灰字沉在页底、太小看不清；
# 中位 UMI 那行是"数"不是"参数"，按"图里不写"移进图注。
note = '3,000 HVG \u00b7 PCA dims 1\u201310 \u00b7 Leiden res 0.5'
FS_P = FS_N + 0.7
if K.text_width(note, FS_P) > CELLS['H']['w'] - 0.10:
    raise RuntimeError('H 参数行越格: %.2f cm' % K.text_width(note, FS_P))
fig.text((CELLS['H']['x'] + CELLS['H']['w'] / 2) / W, (H - PLOT_TOP - PLOT_H - 0.60) / H, note,
         fontsize=FS_P, color=GREY, ha='center', va='top')

# ============================================================================ I 独立 marker 气泡
axI = ax_cell('I')
GRP = ['100k COP/NFOL sub-0', '100k COP/NFOL sub-1', '50k COP/NFOL (single)', '100k OPC (same lib)']
SHORT = ['sub-0', 'sub-1', '50k', 'OPC']
ind = mk[mk.independent]
Pm = ind.pivot(index='gene', columns='group', values='pct')[GRP]
Mm = ind.pivot(index='gene', columns='group', values='mean')[GRP]
Z = Mm.sub(Mm.mean(1), axis=0).div(Mm.std(1) + 1e-9, axis=0)   # 4 组重算 z（去掉 Oligo 之后）
nrm = TwoSlopeNorm(vmin=-2.2, vcenter=0.0, vmax=2.2)
gl0 = sorted(ind.gene.unique())
LOWEXP = 0.10
# 作者：基因给得多了点。原来 21 行按字母序排，7 行"四组全测不到"的空行夹在中间，
# 看起来既长又随机。改成：**能测到的按表达量降序在前，测不到的收成底部一小块（灰名）**，
# 并只留 4 行有论证作用的阴性对照（增殖 Cdk1/Mki67/Top2a、基质 Dcn）——
# 它们证明这两簇不是增殖细胞或成纤维/基质污染的混入；剩下 3 行（Fn1/Myc/Pdgfrb）进图注。
# 要全表：NEG_KEEP = []。
NEG_KEEP = ['Cdk1', 'Mki67', 'Top2a', 'Dcn']
det = sorted([g for g in gl0 if Mm.loc[g].max() >= LOWEXP], key=lambda g: -Mm.loc[g].max())
neg = [g for g in gl0 if Mm.loc[g].max() < LOWEXP]
NEG_DROP = [g for g in neg if g not in NEG_KEEP]
gl = det + [g for g in neg if g in NEG_KEEP]
faint_rows = set(g for g in gl if Mm.loc[g].max() < LOWEXP)
for i, g in enumerate(gl):
    y = len(gl) - 1 - i
    for j in range(len(GRP)):
        s = 2.0 + 34.0 * (Pm.loc[g, GRP[j]] / 100.0) ** 1.2     # 最大墨径 0.239 cm < 行距 0.263
        c = FAINT if g in faint_rows else plt.cm.coolwarm(
            nrm(float(np.clip(Z.loc[g, GRP[j]], -2.2, 2.2))))
        axI.scatter([j], [y], s=s, color=[c], edgecolors=(0.25, 0.25, 0.25), linewidths=0.2, zorder=3)
axI.set_xlim(-0.62, len(GRP) - 0.38); axI.set_ylim(-0.55, len(gl) - 0.45)
axI.set_xticks(range(len(GRP)))
NS = [int(mk[mk.group == g].n.iloc[0]) for g in GRP]
axI.set_xticklabels(['%s\nn=%s' % (SHORT[k], f'{NS[k]:,}') for k in range(len(GRP))],
                    fontsize=FS_N, color=INK)
for t in axI.get_xticklabels():
    t.set_linespacing(1.55)
axI.tick_params(axis='x', length=0, pad=2.2)
axI.set_yticks([len(gl) - 1 - i for i in range(len(gl))])
axI.set_yticklabels(gl, fontsize=FS_N)
for t, g in zip(axI.get_yticklabels(), gl):
    if g in faint_rows:
        t.set_color(GREY)        # ⚠ 不能用 FAINT(0.80)：204 灰阶在 final_audit 的 <200 墨量阈下
                                 #   等于**看不见**（29 个"隐藏字符"全是这 7 行基因名）。
                                 #   淡色只用在点上，标签保持中灰可读。
axI.tick_params(axis='y', length=0, pad=1.4)
axI.axhline(len(gl) - len(det) - 0.5, color=(0.78, 0.78, 0.78), lw=0.5, ls=(0, (2, 2)), zorder=2)
clean(axI, left=False, bottom=False)
header('I', 'independent markers')

# ============================================================================ J 分裂稳定性
axJ = ax_cell('J')
g = stab.groupby('res').agg(split=('split2', 'mean'), ari=('ari', 'mean'), n=('split2', 'size'))
full = stab[stab.kind == 'full'].groupby('res').split2.mean() * 100
draw = stab[stab.kind == 'downsampled'].groupby('res').split2.mean() * 100
xs = np.arange(len(g))
axJ.bar(xs - 0.16, full.values, 0.30, color=BLUE, zorder=3)
axJ.bar(xs + 0.16, draw.values, 0.30, color=PALE, zorder=3)
axJ.plot(xs, g.ari.values * 100, marker='D', ms=3.0, color=LGTEAL, ls='none', zorder=4)
axJ.axhline(70, color=RED, lw=0.6, ls=(0, (2, 2)), zorder=2)
axJ.set_xticks(xs); axJ.set_xticklabels(['%.1f' % r for r in g.index], fontsize=FS_N)
axJ.set_ylim(0, 118); axJ.set_xlim(-0.55, len(g) - 0.45)
axJ.set_xlabel('sub-clustering resolution', fontsize=FS_N, color=INK, labelpad=1.8)
axJ.set_ylabel('split recovered (%) / ARI \u00d7100', fontsize=FS_N, color=INK, labelpad=1.4)
axJ.tick_params(axis='both', labelsize=FS_N, length=2, pad=1.4, colors=INK)
for s in ('top', 'right'):
    axJ.spines[s].set_visible(False)
nfull = int(stab[stab.kind == 'full'].groupby('res').size().iloc[0])
ndraw = int(stab[stab.kind == 'downsampled'].groupby('res').size().iloc[0])
print('  J: full %d 次/档, downsampled %d 次/档（写进图注）' % (nfull, ndraw))
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
HS = [Patch(fc=BLUE, lw=0), Patch(fc=PALE, lw=0),
      Line2D([], [], marker='D', ms=3.0, color=LGTEAL, ls='none'),
      Line2D([], [], color=RED, lw=0.6, ls=(0, (2, 2)))]
LB = ['whole library', '80% downsamples', 'ARI \u00d7100', 'frozen rule, 70%']
leg = axJ.legend(HS, LB, loc='lower left', bbox_to_anchor=(0.0, 1.005), ncol=2, frameon=False,
                 fontsize=FS_N, handlelength=1.0, handletextpad=0.45, columnspacing=1.2,
                 labelspacing=0.45, borderaxespad=0)
fig.canvas.draw()
bb = leg.get_window_extent(fig.canvas.get_renderer())
x0c, y0c = bb.x0 / fig.dpi * 2.54, (fig.get_size_inches()[1] - bb.y1 / fig.dpi) * 2.54
x1c, y1c = bb.x1 / fig.dpi * 2.54, (fig.get_size_inches()[1] - bb.y0 / fig.dpi) * 2.54
if x1c > W - 0.05 or x0c < CELLS['J']['x']:
    raise RuntimeError('J 图例越格: %.2f–%.2f cm（格 %s–%.2f）' % (x0c, x1c, CELLS['J']['x'], W - 0.05))
if y0c < LETTER_Y + 0.16:
    raise RuntimeError('J 图例顶到标题行: 上沿 %.2f cm' % y0c)
print('  J 图例 %.2f–%.2f × %.2f–%.2f cm' % (x0c, x1c, y0c, y1c))
header('J', 'split stability')

OUT = f'{OUTDIR}/FigS3c_v10{SUF}.pdf'
fig.savefig(OUT, format='pdf', bbox_inches=None, pad_inches=0)
plt.close(fig)
print('OK ->', OUT)
print('H 中位 UMI %s vs %s（进图注）| I 行 %d = 可测 %d + 阴性 %d（另 %s 进图注）| J 网格 %d 次'
      % (f'{n0:,.0f}', f'{n1:,.0f}', len(gl), len(det), len(faint_rows), '/'.join(NEG_DROP), int(g.n.sum())))
