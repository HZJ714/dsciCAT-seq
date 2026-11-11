#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p0_8c_fig6c_c1.py —— Fig6C 定案 C1 版：38 个 *cis*-window 风险基因的落点(Ca) + 程序分解(Cb)

依据 FIG6C_FIG7D_NARRATIVE_DESIGN_20261116.md §1.3 / §5 P0-a；取代 v25(top12 分数偶联) 的主图位。

回答 A/B/D 都没答的问题：B 里富集起来的那部分遗传风险，落在成熟 OL 的哪些分子程序上？

防循环（本格要害，必须进图注）
  行 = Table 2 的 38 个 *cis*-window 风险基因（外部输入，人源符号，38/38 同源映射进 RNA 与 ATAC 矩阵）；
  选择与排序都不用 scDRS 分数（排序 = 功能类分块，类内按亚型间 log2 差降序）。
  Ca 纯描述（每亚型平均表达/可及性，行内 z），**不出现逐基因 p**；
  Cb 只在类水平检验（38 基因内置换 competitive test），hypothesis-generating。

P0-a 实测（p0_8c_c1_precheck.py，本脚本据此定版式）
  L2「OL 亚型特异表达」二值口径：可检出 34/38=89.5%（几乎全实心）、log2FC>=1 仅 3/38、严格 marker 1/38
  ⇒ 按 §3 双边规则，二值层全部不占列，改连续色热图。
  top 亚型分布 MFOL 17 / MOL2 9 / COP-NFOL 7 / MOL6 4 / MOL3 1 ⇒ 只有 4/38 在富集亚型达峰。

诚实要点（本版本的核心修订）
  类**均值**被单基因绑架：lipid/energy 均值 +0.34 全部来自 PRR5L（+1.97），中位仅 +0.06；
  myelin structure 均值 −0.59 全部来自 BCAS1（−1.91），中位 +0.04。
  ⇒ 主统计量用**中位数**（实心条），均值只画空心方块，单基因一律叠小点，读者自己看得见是谁在撑。

版式 17.00 × 5.50 cm；字号两档 8.0/6.0；左上角留给拼页 panel letter；色标不压数据、不出页。
"""
import os, sys, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

CM = 1 / 2.54
WORK = os.environ.get('DSCI_WORK_ROOT', '')
BASE = WORK + "/figures/routeA_out/P0_8C"
TAB, PAN = BASE + "/tables", BASE + "/panels"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arial_style import use_arial, use_two_tier
use_arial(); use_two_tier()

SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]      # 成熟轴，与 Fig6A 同序
ENA, ENB = "MOL3/6", "MOL6"                                        # Fig6A group_q<0.05 的富集列
N_ENA, N_ENB, N_COP = 991, 2574, 203
ANCHORS = ["Mbp", "Plp1", "Mog", "Mag", "Cnp", "Mobp", "Sox10", "Olig2"]

CLS = {                                                            # 人工归类、非穷举（图注声明）
    "lipid / energy":     ["APOE", "PLIN3", "EFR3B", "PHGDH", "PRR5L", "ATP1A1"],
    "myelin structure":   ["CNP", "BCAS1", "PLLP"],
    "cytoskeleton / MT":  ["CAMSAP2", "TUBB4A", "SCHIP1", "YPEL2", "CEP97"],
    "vesicle / membrane": ["SEC11C", "PHLDB1", "ENDOD1", "CLIC4", "TMCC3"],
    "immune / adhesion":  ["LITAF", "TRIM59", "SKAP2", "AIG1", "MCAM", "IGSF8", "CD9", "ADAMTS4", "SPON1"],
    "lineage / TF":       ["OLIG1", "ZNF365", "CREB5", "ZBTB38"],
    "other signalling":   ["DDC", "TMEM151A", "SEZ6L2", "KCTD13", "SBF1", "RRAS2"],
}
CCOL = {"lipid / energy": "#3f7fb0", "myelin structure": "#6f9c3f", "cytoskeleton / MT": "#b58a3c",
        "vesicle / membrane": "#8e6fa8", "immune / adhesion": "#c1544f", "lineage / TF": "#4f8f8a",
        "other signalling": "#8f8f8f"}
g2c = {g: c for c, gs in CLS.items() for g in gs}
assert len(g2c) == 38

# ---------------- 数据 ----------------
S = pd.read_csv(TAB + "/p08c_c1_rna_atac_by_subtype.csv")
M = pd.read_csv(TAB + "/p08c_c1_gene38_map.csv")
S["cls"], S["mouse"] = S.gene.map(g2c), S.gene.map(dict(zip(M.gene, M.mouse)))
S["dist_kb"] = S.gene.map(dict(zip(M.gene, M.dist_kb)))
L = pd.read_csv(TAB + "/p08c_subtype_stats_long.csv")
A = L[L.gene.isin(ANCHORS)].pivot_table(index="gene", columns="subtype", values=["rna_mean", "atac_mean"])

def pdelta(vals, atac=False):                                     # 富集池(细胞数加权) − COP/NFOL
    pre = "atac_" if atac else "rna_"
    return ((vals[pre + ENA] * N_ENA + vals[pre + ENB] * N_ENB) / (N_ENA + N_ENB)) - vals[pre + SUB[0]]
S["d_rna"] = pdelta(S, False); S["d_atac"] = pdelta(S, True)
anc_r = np.array([pdelta({("rna_" + s): float(A.loc[g, ("rna_mean", s)]) for s in SUB}) for g in ANCHORS])
anc_a = np.array([pdelta({("atac_" + s): float(A.loc[g, ("atac_mean", s)]) for s in SUB}, True) for g in ANCHORS])

# 类水平 competitive test（38 基因内置换，保持类大小；中位数为主统计量）
rng = np.random.default_rng(7); NP = 20000
lab = S.cls.values
rows_c = []
for c, gs in CLS.items():
    m = lab == c
    for ch, v in (("rna", S.d_rna.values), ("atac", S.d_atac.values)):
        null = np.array([v[rng.permutation(len(lab))][m] for _ in range(NP)])
        med, mean = np.median(v[m]), np.mean(v[m])
        rows_c.append(dict(cls=c, ch=ch, n=int(m.sum()), med=med, mean=mean,
                           p_med=max((np.median(null, 1) >= med).mean(), 1 / (NP + 1)),
                           p_mean=max((null.mean(1) >= mean).mean(), 1 / (NP + 1))))
C = pd.DataFrame(rows_c).pivot(index="cls", columns="ch", values=["med", "mean", "p_med", "p_mean", "n"])
C.columns = ["%s_%s" % (a, b) for a, b in C.columns]
C = C.reset_index().sort_values("med_rna", ascending=False).reset_index(drop=True)
C["n"] = [len(CLS[k]) for k in C.cls]
C.to_csv(TAB + "/p08c_c1_class_effect.csv", index=False)
print(C.round(4).to_string(index=False))
print("锚基因参照（同口径）：RNA 中位 %+0.3f 均值 %+0.3f | ATAC 中位 %+0.3f 均值 %+0.3f"
      % (np.median(anc_r), anc_r.mean(), np.median(anc_a), anc_a.mean()))

# ---------------- 版式预算（全部显式，不再靠估） ----------------
PW, PH = 17.00, 6.00      # 38 行×两叠 ⇒ 6pt 行标签需要 0.236 cm 行距，5.50 cm 排不下（回执 §11 记偏差）
fig = plt.figure(figsize=(PW * CM, PH * CM))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, PW); ax.set_ylim(0, PH); ax.axis("off")
for a in fig.axes[1:]:
    a.axis("off")

YA, YT = 5.74, 5.42                     # 8pt 主标题基线 / 6pt 通道组名基线
BODY_T, BODY_B = 5.34, 0.86             # 热图体（下方 0.86 cm 给横排轴标，不再用旋转标签）
NR = 19
RH = (BODY_T - BODY_B) / NR
LBLW, RIBW, DISW, COLW, CHGAP = 1.10, 0.14, 0.40, 0.215, 0.07
STKW = LBLW + RIBW + DISW + 12 * COLW + CHGAP + 0.06
X_S1, X_S2 = 0.30, 0.30 + STKW + 0.30
X_CBAR = X_S2 + STKW + 0.10
X_CB0 = X_CBAR + 0.68                                   # Cb 类名左边界
X_CB1 = 16.72                                          # Cb 右边界
print("Ca 两叠 x=%.2f/%.2f 宽 %.2f | 色标 %.2f | Cb %.2f–%.2f (宽 %.2f)"
      % (X_S1, X_S2, STKW, X_CBAR, X_CB0, X_CB1, X_CB1 - X_CB0))
assert X_CB1 <= PW - 0.20

CMAP = plt.get_cmap("coolwarm")
ax.text(1.55, YA, "Risk-set landing across OL subtypes", fontsize=8.0, va="baseline", weight="bold")

def stack_x(k):
    x0 = X_S1 if k == 0 else X_S2
    xr = x0 + LBLW + RIBW + DISW
    return x0, xr, xr + 6 * COLW + CHGAP

def draw_stack(k, rows):
    x0, x_rna, x_atac = stack_x(k)
    x_cls, x_dist = x0 + LBLW, x0 + LBLW + RIBW
    if k == 0:
        ax.text(x0, YT, "(a)", fontsize=6.0, ha="left", va="baseline", color="#333333")
    ax.text(x_rna + 3 * COLW, YT, "RNA", fontsize=6.0, ha="center", va="baseline", color="#333333")
    ax.text(x_atac + 3 * COLW, YT, "ATAC", fontsize=6.0, ha="center", va="baseline", color="#333333")
    ax.text(x_dist + DISW / 2, YT, "\u00b1kb", fontsize=6.0, ha="center", va="baseline", color="#666666")
    y, prev = BODY_T, None
    for i, (_, r) in enumerate(rows.iterrows()):
        y -= RH
        ax.text(x0 + LBLW - 0.05, y + RH / 2, r.gene, fontsize=6.0, ha="right", va="center")
        ax.add_patch(Rectangle((x_cls, y), RIBW, RH, fc=CCOL[r.cls], ec="none"))
        d = float(r.dist_kb); c0 = x_dist + DISW / 2
        ax.plot([c0, c0 + np.sign(d) * min(abs(d), 500) / 500 * (DISW / 2 - 0.03)],
                [y + RH / 2] * 2, lw=0.5, color="#b5502a" if d > 0 else "#3f7fb0", zorder=3)
        ax.plot(c0, y + RH / 2, ms=0.8, color="#444444", zorder=4)
        for xg, pre, pct in ((x_rna, "rna_", "pct_"), (x_atac, "atac_", None)):
            vv = [float(r[pre + s]) for s in SUB]
            mu, sd = np.nanmean(vv), np.nanstd(vv)
            for j, s in enumerate(SUB):
                ax.add_patch(Rectangle((xg + j * COLW, y), COLW - 0.012, RH - 0.012,
                                       fc=CMAP((np.clip((vv[j] - mu) / (sd + 1e-12), -2, 2) + 2) / 4), ec="none"))
                if pct and float(r[pct + s]) < 10:
                    ax.plot(xg + j * COLW + COLW / 2, y + RH / 2, marker=".", ms=1.0, color="#555555", zorder=5)
        if prev is not None and r.cls != prev:
            ax.plot([x_cls, x_atac + 6 * COLW], [y + RH, y + RH], lw=0.35, color="#b0b0b0", zorder=6)
        prev = r.cls
    for xg in (x_rna, x_atac):                         # 富集列顶托（与 Fig6A group_q<0.05 同口径）
        for s in (ENA, ENB):
            ax.add_patch(Rectangle((xg + SUB.index(s) * COLW - 0.005, BODY_T), COLW, 0.05, fc="#222222", ec="none"))
    for xg in (x_rna, x_atac):                         # 分通道托线 + 亚型分界小刻度
        ax.plot([xg, xg + 6 * COLW], [BODY_B - 0.08, BODY_B - 0.08], lw=0.4, color="#999999")
        for j in range(1, 6):
            ax.plot([xg + j * COLW], [BODY_B - 0.08], marker="|", ms=1.5, color="#999999", zorder=3)
    ax.text(x_rna, BODY_B - 0.40, "COP/NFOL", fontsize=6.0, ha="left", va="top", color="#333333")
    ax.text(x_atac + 6 * COLW, BODY_B - 0.40, "MOL6", fontsize=6.0, ha="right", va="top", color="#333333")
    ya_ = BODY_B - 0.49
    ax.annotate("", xy=(x_atac + 6 * COLW - 0.58, ya_), xytext=(x_rna + 1.06, ya_),
                arrowprops=dict(arrowstyle="-|>", lw=0.5, color="#666666", shrinkA=0, shrinkB=0))

def stack_of(cs):
    return pd.concat([S[S.cls == c].sort_values("d_rna", ascending=False) for c in cs], ignore_index=True)
h1 = stack_of(["lipid / energy", "myelin structure", "cytoskeleton / MT", "vesicle / membrane"])
h2 = stack_of(["immune / adhesion", "lineage / TF", "other signalling"])
assert len(h1) == 19 and len(h2) == 19, (len(h1), len(h2))
draw_stack(0, h1); draw_stack(1, h2)

# ---------------- 共享色标 ----------------
gh = np.linspace(-2, 2, 64)[:, None]
ax.imshow(np.hstack([gh, gh]), extent=(X_CBAR, X_CBAR + 0.22, BODY_B + 0.86, BODY_T - 0.86),
          cmap=CMAP, vmin=-2, vmax=2, aspect="auto", zorder=2)
CB_T, CB_B = BODY_T - 0.80, BODY_B + 0.80
ax.text(X_CBAR + 0.10, CB_T + 0.06, "row z", fontsize=6.0, ha="center", va="baseline", color="#333333")
ax.text(X_CBAR + 0.24, CB_T, "+2", fontsize=6.0, ha="left", va="top", color="#333333")
ax.text(X_CBAR + 0.24, CB_B, "−2", fontsize=6.0, ha="left", va="bottom", color="#333333")
ax.plot(X_CBAR + 0.05, BODY_B + 0.52, marker=".", ms=1.2, color="#555555")
ax.text(X_CBAR + 0.13, BODY_B + 0.52, "pct<10%", fontsize=6.0, ha="left", va="center", color="#555555")

# ---------------- Cb：类水平效应量 ----------------
ax.text(X_CB0 + 0.02, YA, "(b) Programme decomposition", fontsize=8.0, va="baseline", weight="bold")
ax.text(X_CB0, YT, "\u0394 enriched (MOL3/6+MOL6) \u2212 COP/NFOL, log2", fontsize=6.0, va="baseline", color="#333333")
ax.text(X_CB0, YT - 0.20, "bar = class median; \u25a1 = mean; dots = genes", fontsize=6.0, va="baseline", color="#666666")
LBL2 = 2.22
x0 = X_CB0 + LBL2
x1 = X_CB1 - 0.44
LO, HI = -0.62, 0.62          # 类比较用压缩轴；越界单基因以三角标出，另在注脚给真值
sc = (x1 - x0) / (HI - LO)
def px(v):
    return x0 + (v - LO) * sc
ax.plot([px(0), px(0)], [1.62, 5.06], lw=0.6, color="#333333", zorder=3)
BH = (4.98 - 1.68) / (len(C) + 1)
y = 5.06
for _, r in C.iterrows():
    y -= BH
    yc = y + BH / 2
    ax.text(x0 - 0.06, yc, "%s (%d)" % (r.cls, r.n), fontsize=6.0, ha="right", va="center", color=CCOL[r.cls])
    for ch, col, off in (("rna", CCOL[r.cls], +BH * 0.17), ("atac", "#8a8a8a", -BH * 0.17)):
        med, mean = float(r["med_" + ch]), float(r["mean_" + ch])
        ax.plot([px(0), px(med)], [yc + off] * 2, lw=2.0, color=col, solid_capstyle="butt", zorder=4)
        ax.plot(px(mean), yc + off, marker="s", ms=1.6, mfc="none", mec=col, mew=0.5, zorder=5)
        gv = S.loc[S.cls == r.cls, "d_" + ch].values
        for t, gvv in enumerate(gv):
            yy = yc + off + (t % 3 - 1) * BH * 0.08
            if LO < gvv < HI:
                ax.plot(px(gvv), yy, ms=0.8, mfc="none", mec=col if ch == "rna" else "#b0b0b0", mew=0.35, zorder=6)
            else:                                   # 越界基因：贴边三角，方向即符号
                ax.plot(px(HI if gvv > 0 else LO), yy, ms=1.3, marker=">" if gvv > 0 else "<",
                        mfc="none", mec=col if ch == "rna" else "#b0b0b0", mew=0.4, zorder=6)
        if r["p_med_" + ch] < 0.05:
            ax.text(px(med) + 0.10, yc + off, "P<0.05", fontsize=6.0, va="center", ha="left", color=col)
ax.plot([px(0), px(0)], [y, 5.06], lw=0.6, color="#333333", zorder=3)
y -= BH
ax.text(x0 - 0.06, y + BH / 2, "myelin anchors* (8)", fontsize=6.0, ha="right", va="center", color="#333333")
for v_, col, off in ((np.median(anc_r), "#222222", +BH * 0.17), (np.median(anc_a), "#999999", -BH * 0.17)):
    ax.plot([px(0), px(v_)], [y + BH / 2 + off] * 2, lw=2.0, color=col, solid_capstyle="butt", zorder=4)
for tv in (-0.5, 0.0, 0.5):
    ax.text(px(tv), 1.56, "%+0.1f" % tv if tv else "0", fontsize=6.0, ha="center", va="top", color="#333333")
ax.text(X_CB0, 1.24, "* 8 canonical OL/myelin genes, same contrast, not in the risk set",
        fontsize=6.0, ha="left", va="top", color="#666666")
ax.text(X_CB0, 0.98, "off-scale: PRR5L +1.97, BCAS1 \u22121.91 (they drive the class means)",
        fontsize=6.0, ha="left", va="top", color="#666666")
ax.text(X_CB0, 0.72, "permutation within the 38 genes, 20,000 draws; hypothesis-generating",
        fontsize=6.0, ha="left", va="top", color="#666666")
ax.text(X_CB0, 0.46, "neither selection nor ordering of these rows uses the scDRS score",
        fontsize=6.0, ha="left", va="top", color="#666666")

OUT_P = PAN + "/fig6v27_C_riskset_landing.pdf"
fig.savefig(OUT_P, transparent=False)
import fitz
d = fitz.open(OUT_P); p = d[0]
print("[saved] %s\n尺寸 %.2f × %.2f cm" % (OUT_P, p.rect.width * 2.54 / 72, p.rect.height * 2.54 / 72))
