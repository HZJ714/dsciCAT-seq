import os
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""s3c_cop2_figure.py —— 方案 A：独立小图 Figure S3c（18.00 × 8.30 cm，四格一排）
=========================================================================================
四个格子对应**在看数据之前冻结**的三条判据（S3C_COP2_FROZEN_CRITERIA_20261029.md）：
  A 独立 marker 点阵   —— 21 个没参与 COP/NFOL 判定的基因 × 5 组（同库 OPC/Oligo 作参照）
  B 无先验差异 + 深度匹配敏感性 —— 两簇差异基因在"把两组一起下采样到同一深度"前后的变化
  C 分裂稳定性         —— 分辨率 × 种子 × 下采样网格上的复现率与 ARI
  D 技术混杂           —— 9 项里可评的 6 项，Cliff's delta + BH（中位数两栏数字）
底部一段给出**判据结论**（S1/S2/S3 逐条 PASS/FAIL），不粉饰：三条必须同时过才保留"两亚群"。
作者口径：图面文字尽量少、长解释进图注 ⇒ 图内只留读图必需的小图例，没有通栏说明带。

⚠ 三条踩过的坑（都改在这里）：
  1) merge 之后列名是 `mean_diff_a/mean_diff_m`，不是 `diff_a/diff_m`（上一版 KeyError 的来源）；
  2) 图面数字一律从 CSV 现算。上一版手写 "146 of 31,263 → 36"：146 是 Seurat FindMarkers 的
     p_val_adj（Bonferroni 口径），冻结判据写的是 BH-FDR ⇒ 真口径 260；"→36" 那对数来自
     一个「秩与 p 错配」的 BH bug，修好后是 306 → 251（8,081 基因宇宙）；
  3) 版式：v1 把 3 行说明写在每格下面，结果说明文字横向压进隔壁格、标题压进隔壁字母；
     现在标题用 K.text_width 实测宽度并 assert 不越格，图例改放进各坐标轴自己的留白里。
"""
import json, os, sys, textwrap
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figkit as K
K.use_matplotlib_arial()
from matplotlib import pyplot as plt
from matplotlib.colors import TwoSlopeNorm

D = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'u13_cop_4k/s3c')
OUTDIR = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'replot/S3c')
os.makedirs(OUTDIR, exist_ok=True)
W = 18.00
RED, GREY, INK = (0.835, 0.122, 0.149), (0.55, 0.55, 0.55), (0.13, 0.13, 0.13)
FAINT, BLUE, LGTEAL, PALE = (0.80, 0.80, 0.80), (0.20, 0.42, 0.70), (0.13, 0.47, 0.31), (0.66, 0.77, 0.90)
FS_N = 5.5
CALL = ['Bmp4', 'C1ql1', 'Gpr17', 'Ptprz1']

# 两种版式（作者二选一）：
#   完整版 18.00 × 9.75 cm —— 26 行（21 独立 + 4 判据基因灰行）+ 底部判据条
#   compact 18.00 × 7.55 cm —— 去掉 4 行灰色判据基因（判据基因的值进图注）+ 判据条进图注
COMPACT = '--compact' in sys.argv
# 面板字母可换：并进 Figure S3 时 A–D 会和母版 A–G 撞号，必须 --letters=HIJK。
# 只改**画出来的字**和图内自引用（"Depth matching (B)"），内部键名仍是 A–D。
LET = {}
for a in sys.argv[1:]:
    if a.startswith('--letters='):
        v = a.split('=', 1)[1].strip()
        assert len(v) == 4 and v.isalpha(), '--letters 要 4 个字母，收到 %r' % v
        LET = dict(zip('ABCD', v.upper()))


def L(c):
    return LET.get(c, c)


SUF = ('_' + ''.join(LET.values())) if LET else ''
LETTER_Y, PLOT_TOP = 0.36, 1.02
PLOT_H = 5.15 if COMPACT else 6.45     # 行距下限 0.245 cm（5.5 pt 上下不贴字）
RULE_Y, VERD_Y = (0, 0) if COMPACT else (8.42, 8.60)
H = 7.75 if COMPACT else 9.75
# x = 格最左沿, gut = 行名沟槽, w = 绘图区宽；右沿 = x+gut+w
CELLS = {'A': dict(x=0.42, gut=0.66, w=3.30),     # 右沿 4.38
         'B': dict(x=4.62, gut=0.86, w=2.30),     # 右沿 7.78
         'C': dict(x=8.00, gut=0.90, w=2.75),     # 右沿 11.65（原来 3.00 宽会越过 D 的左沿 11.80，
                                                  #   于是 C 的 rule-70 虚线一直画进 D 的行名沟槽）
         'D': dict(x=11.95, gut=1.34, w=4.20)}    # 右沿 17.49
RIGHT_EDGE = {c: CELLS[c]['x'] + CELLS[c]['gut'] + CELLS[c]['w'] for c in CELLS}

mk = pd.read_csv(f'{D}/s3c_marker_table.csv')
de = pd.read_csv(f'{D}/s3c_de.csv')
da = pd.read_csv(f'{D}/s3c_de_asis.csv')
dm = pd.read_csv(f'{D}/s3c_de_matched.csv')
stab = pd.read_csv(f'{D}/s3c_stability_scored.csv')
S2 = json.load(open(f'{D}/s3c_S2.json'))
qc = pd.read_csv(f'{D}/s3c_qc.csv')
meta = pd.read_csv(f'{D}/s3c_meta.csv')


def bh(p):
    p = np.asarray(p, float); n = p.size
    o = np.argsort(p, kind='stable'); ps = p[o]
    qs = np.minimum.accumulate((ps * n / np.arange(1, n + 1))[::-1])[::-1]
    q = np.empty(n); q[o] = np.clip(qs, 0, 1); return q


bh_seurat = bh(de.p_val.values)
n_s1_bh = int(((~de.in_call_set) & (de.mean_diff.abs() >= 0.5) & (bh_seurat < 0.05)).sum())
n_s1_bonf = int(de.qualifies.sum())
n_a, n_m, n_gene = int(da.qualifies.sum()), int(dm.qualifies.sum()), len(da)
TARGET = int(min(np.median(meta.loc[meta.cluster == 0, 'nCount_RNA']),
                 np.median(meta.loc[meta.cluster == 1, 'nCount_RNA'])))
P2, P5 = S2['primary_frozen_grid'], S2['secondary_res_matched']
mid = stab[stab.res.isin([0.3, 0.5])]
n_qc_ok, n_qc_bad = int(qc.available.sum()), int(qc.flag_bad.sum())
bad_vars = ', '.join(qc[qc.flag_bad]['var'])

fig = plt.figure(figsize=(W / 2.54, H / 2.54))


def ax_cell(c, y0=PLOT_TOP, h=PLOT_H):
    g = CELLS[c]
    return fig.add_axes([(g['x'] + g['gut']) / W, (H - y0 - h) / H, g['w'] / W, h / H])


def header(c, title):
    """面板字母 8 pt Arial 常规（house rule）+ 短标题；实测宽度，越格就报错而不是叠字"""
    g = CELLS[c]
    fig.text(g['x'] / W, (H - LETTER_Y) / H, L(c), fontsize=8.0, family=K.FAM, color=INK,
             ha='left', va='baseline')
    x1 = g['x'] + K.text_width(L(c), 8.0) + 0.20
    nxt = {'A': CELLS['B']['x'], 'B': CELLS['C']['x'], 'C': CELLS['D']['x'], 'D': W - 0.40}[c]
    if x1 + K.text_width(title, FS_N + 0.7) > nxt - 0.08:
        raise RuntimeError('标题越格 %s: 右沿 %.2f cm，隔壁在 %.2f cm（标题 "%s"）'
                           % (c, x1 + K.text_width(title, FS_N + 0.7), nxt, title))
    fig.text(x1 / W, (H - LETTER_Y) / H, title, fontsize=FS_N + 0.7, family=K.FAM,
             color=INK, ha='left', va='baseline')


def clean(ax, left=True, bottom=True):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.spines['left'].set_visible(left)
    ax.spines['bottom'].set_visible(bottom)


# ============================================================================ A 独立 marker
axA = ax_cell('A')
ind, call = mk[mk.independent], mk[~mk.independent]
gl, cl = sorted(ind.gene.unique()), sorted(call.gene.unique())
GRP = ['100k COP/NFOL sub-0', '100k COP/NFOL sub-1', '50k COP/NFOL (single)',
       '100k OPC (same lib)', '100k Oligo (same lib)']
NS = [int(mk[mk.group == g].n.iloc[0]) for g in GRP]
Pm = ind.pivot(index='gene', columns='group', values='pct')[GRP]
Mm = ind.pivot(index='gene', columns='group', values='mean')[GRP]
PC, MC = call.pivot(index='gene', columns='group', values='pct')[GRP], call.pivot(index='gene', columns='group', values='mean')[GRP]
Z = Mm.sub(Mm.mean(1), axis=0).div(Mm.std(1) + 1e-9, axis=0)
nrm = TwoSlopeNorm(vmin=-2.2, vcenter=0.0, vmax=2.2)
LOWEXP = 0.10
rows = gl + ([] if COMPACT else [None] + cl)
IG = len(gl)
for i, g in enumerate(rows):
    y = len(rows) - 1 - i - (0.45 if i == IG else 0)
    if g is None:
        axA.plot([-0.55, 4.55], [y + 0.5, y + 0.5], color=(0.82, 0.82, 0.82), lw=0.4, zorder=1)
        continue
    is_call = i > IG
    Pi, Mi = (PC if is_call else Pm).loc[g], (MC if is_call else Mm).loc[g]
    faint = (not is_call) and Mi.max() < LOWEXP
    for jj in range(5):
        s = 2.0 + 30.0 * (Pi[GRP[jj]] / 100.0) ** 1.2   # 最大墨径 0.22 cm < 行距 0.245 cm
        c = GREY if is_call else (FAINT if faint else
                                  plt.cm.coolwarm(nrm(float(np.clip(Z.loc[g, GRP[jj]], -2.2, 2.2)))))
        axA.scatter([jj], [y], s=s, color=[c], edgecolors=(0.25, 0.25, 0.25), linewidths=0.2, zorder=3)
axA.set_xlim(-0.62, 4.62); axA.set_ylim(-0.55, len(rows) - 0.45)
axA.set_xticks(range(5))
SHORT = ['sub-0', 'sub-1', '50k', 'OPC', 'Oligo']
axA.set_xticklabels(['%s\nn=%d' % (SHORT[k], NS[k]) for k in range(5)], fontsize=FS_N, color=INK)
for t in axA.get_xticklabels():
    t.set_linespacing(1.55)          # 两行刻度标签默认 1.2 行距会把 'Oligo' 和 'n=3614' 的墨贴上
axA.tick_params(axis='x', length=0, pad=2.2)
axA.set_yticks([len(rows) - 1 - i - (0.45 if i == IG else 0) for i, g in enumerate(rows) if g])
axA.set_yticklabels([g for g in rows if g], fontsize=FS_N)
for t, g in zip(axA.get_yticklabels(), [g for g in rows if g]):
    is_call = g in cl
    Mi = (MC if is_call else Mm).loc[g]
    t.set_color(GREY if (is_call or Mi.max() < LOWEXP) else INK)
axA.tick_params(axis='y', length=0, pad=1.4)
clean(axA, left=False, bottom=False)
header('A', 'independent markers')

# ============================================================================ B DE + 深度匹配
axB = ax_cell('B')
j = da.merge(dm, on='gene', suffixes=('_a', '_m'))
q = j[j.qualifies_a & ~j.gene.isin(CALL)].copy()
q['am'] = q[['mean_diff_a', 'mean_diff_m']].abs().mean(1)
q = q.sort_values('am', ascending=False).head(11)
yy = np.arange(len(q))[::-1]
axB.set_ylim(-0.85, len(q) + 2.6)
for k, (_, r) in enumerate(q.iterrows()):
    axB.plot([r.mean_diff_a, r.mean_diff_m], [yy[k], yy[k]], color=(0.72, 0.72, 0.72), lw=0.7, zorder=2)
    axB.scatter([r.mean_diff_a], [yy[k]], s=11, facecolors='none', edgecolors=BLUE, linewidths=0.6, zorder=3)
    axB.scatter([r.mean_diff_m], [yy[k]], s=11, color=RED, zorder=4)
axB.axvline(0, color=INK, lw=0.5, zorder=1)
axB.set_yticks(yy); axB.set_yticklabels(q.gene.values, fontsize=FS_N, color=INK)
axB.set_xlim(-2.9, 2.7)
axB.set_xlabel('mean difference,\nlog-normalised (sub-0 \u2212 sub-1)', fontsize=FS_N, color=INK, labelpad=1.8,
               linespacing=1.55)      # 默认 1.2 行距下两行紧墨盒重叠 0.36 pt（inkcheck 报出来的）
axB.tick_params(axis='both', labelsize=FS_N, length=2, pad=1.4, colors=INK)
yb = len(q) + 1.9
axB.scatter([-2.70, -2.70], [yb, yb - 1.15], s=11, facecolors=['none', RED],
            edgecolors=[BLUE, RED], linewidths=0.6)
axB.text(-2.30, yb, 'as analysed', fontsize=FS_N, color=INK, va='center')
axB.text(-2.30, yb - 1.15, 'depth matched', fontsize=FS_N, color=INK, va='center')
clean(axB)
header('B', 'depth-matched re-test')

# ============================================================================ C 稳定性
axC = ax_cell('C')
res = sorted(stab.res.unique())
bw = 0.33
axC.set_ylim(0, 150)
for i, r in enumerate(res):
    sf, sd = (stab.res == r) & (stab.kind == 'full'), (stab.res == r) & (stab.kind == 'downsampled')
    axC.bar(i - bw / 2, 100 * stab.split2[sf].mean(), bw, color=BLUE, edgecolor='none')
    axC.bar(i + bw / 2, 100 * stab.split2[sd].mean(), bw, color=PALE, edgecolor='none')
    axC.scatter([i - bw / 2, i + bw / 2], [100 * stab.ari[sf].mean(), 100 * stab.ari[sd].mean()],
                s=8, marker='D', color=LGTEAL, zorder=5, linewidths=0)
axC.axhline(70, color=RED, lw=0.6, ls=(0, (2, 1.5)), zorder=1)
axC.text(-0.56, 71.5, 'rule 70', fontsize=FS_N, color=RED, ha='left', va='bottom')
for k, txt in enumerate(['bars: whole library (dark) /', '80% downsamples (pale)',
                         'diamonds: ARI vs Figure 3']):
    axC.text(-0.56, 146 - k * 11.0, txt, fontsize=FS_N, color=INK, va='center', ha='left')
axC.set_xticks(range(len(res))); axC.set_xticklabels(['%.1f' % r for r in res], fontsize=FS_N, color=INK)
axC.set_xlim(-0.62, 4.62)
axC.set_ylabel('split recovered / ARI (%)', fontsize=FS_N, color=INK, labelpad=1.5)
axC.set_xlabel('sub-clustering resolution', fontsize=FS_N, color=INK, labelpad=1.8)
axC.tick_params(axis='both', labelsize=FS_N, length=2, pad=1.4, colors=INK)
axC.set_yticks([0, 20, 40, 60, 80, 100])
clean(axC)
header('C', 'split stability')

# ============================================================================ D 技术混杂
axD = ax_cell('D')
qq = qc[qc.available].sort_values('p_bh')
yy = np.arange(len(qq))[::-1]
axD.axvspan(-0.2, 0.2, color=(0.94, 0.94, 0.94), zorder=0)
axD.axvline(0, color=INK, lw=0.5, zorder=1)
for k, (_, r) in enumerate(qq.iterrows()):
    col = RED if r.flag_bad else INK
    axD.plot([0, r.delta], [yy[k], yy[k]], color=col, lw=0.9, zorder=2, solid_capstyle='round')
    axD.scatter([r.delta], [yy[k]], s=15, color=[col], zorder=3, edgecolors='none')
    axD.text(0.80, yy[k], '{:,}'.format(int(round(r.med0))) if r.med0 >= 1 else '%.3f' % r.med0,
             fontsize=FS_N, color=col, va='center', ha='right')
    axD.text(1.72, yy[k], '{:,}'.format(int(round(r.med1))) if r.med1 >= 1 else '%.3f' % r.med1,
             fontsize=FS_N, color=col, va='center', ha='right')
axD.text(0.80, len(qq) + 0.05, 'median\nsub-0', fontsize=FS_N, color=GREY, va='bottom', ha='right',
         linespacing=1.55)
axD.text(1.72, len(qq) + 0.05, 'sub-1', fontsize=FS_N, color=GREY, va='bottom', ha='right')
axD.text(0.0, -0.72, '|delta| \u2264 0.2', fontsize=FS_N, color=GREY, va='top', ha='center')
axD.set_yticks(yy); axD.set_yticklabels(qq['var'], fontsize=FS_N, color=INK)
axD.set_xlim(-0.55, 1.80); axD.set_ylim(-1.30, len(qq) + 1.60)
axD.set_xticks([-0.4, -0.2, 0.0])
axD.set_xlabel("Cliff's delta (sub-0 vs sub-1)", fontsize=FS_N, color=INK, labelpad=1.8)
axD.tick_params(axis='both', labelsize=FS_N, length=2, pad=1.4, colors=INK)
clean(axD)
header('D', 'technical confounders')

if not COMPACT:
    # --- 判据条（compact 版里这段进图注，不画图）    # ============================================================================ 判据条
    VERD = ('FROZEN RULE, ALL THREE MUST HOLD TO NAME TWO SUBSTATES.   '
            'S1 independent genes >=3: PASS (%d of %s tested).   '
            'S2 recovery >=70%% and ARI >=0.70 over the frozen grid: FAIL (%.0f%% / ARI %.2f) \u2014 the split is '
            'recovered in 100%% of runs with ARI %.2f at res 0.3\u20130.5, and the grid also counts res 0.1/0.8/1.0, '
            'where a two-cluster answer is not expected.   '
            'S3 no technical metric both BH<0.05 and |delta|>0.2: FAIL (%s).   '
            'Depth matching (%s) keeps %d of the %d genes, all in the same direction.   '
            'VERDICT: two spatially separable sub-clusters within one library and one animal, running along a '
            'depth/maturation gradient; reported descriptively, not as two substates.'
            % (n_s1_bh, f'{len(de):,}', P2['pct'], P2['ari'], P5['ari'], bad_vars, L('B'), n_m, n_a))
    fig.text(0.42 / W, (H - RULE_Y) / H, '\n'.join(textwrap.wrap(VERD, 172)),
             fontsize=FS_N, color=INK, va='top', ha='left', linespacing=1.55)
    axv = fig.add_axes([0.42 / W, (H - (RULE_Y - 0.18)) / H, (W - 0.84) / W, 0.0001])
    axv.axhline(0.5, color=(0.82, 0.82, 0.82), lw=0.4); axv.set_axis_off()

if COMPACT:
    print('compact 版：判据条与灰色判据基因行都进图注，图面只留四格')

OUT = f'{OUTDIR}/FigS3c_{"compact" if COMPACT else "full"}{SUF}.pdf'
fig.savefig(OUT, format='pdf', bbox_inches=None, pad_inches=0)
plt.close(fig)
print('OK ->', OUT)
print('S1 BH=%d (Bonferroni 口径 %d) | 深度匹配 %d → %d / %d 基因 | S2 网格 %.1f%% ARI %.3f | '
      'S2 res0.3-0.5 %.0f%% ARI %.3f (%d 次) | S3 可评 %d 项挂 %d 项 (%s)'
      % (n_s1_bh, n_s1_bonf, n_a, n_m, n_gene, P2['pct'], P2['ari'],
         100 * mid.split2.mean(), mid.ari.mean(), len(mid), n_qc_ok, n_qc_bad, bad_vars))
