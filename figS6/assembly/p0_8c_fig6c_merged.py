#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p0_8c_fig6c_merged.py —— Fig6C 合并版（Ca+Cb 合一）：38 个 *cis*-window 风险基因，两列热图 + 右侧程序汇总
用法：python3 p0_8c_fig6c_merged.py [delta|cluster]
  delta   行序 = 显示量 Δ 降序（bootstrap 秩相关 0.98；相邻行 Δ 差中位 0.025 ⇒ 又稳又成块）
  cluster 行序 = 表达谱 ward 聚类模块（作者要求先看效果；实测 ARI 中位 0.56–0.71，峰亚型复现仅 3/38）

设计要点
  * 两列一列 RNA 一列 ATAC，显示的就是被检验的那个量：Δ = 富集(MOL3/6+MOL6，细胞数加权) − COP/NFOL（log2）。
    —— 修掉旧版"Ca 显示行内 z、Cb 显示 Δ"两套口径不一致的毛病（Fig4B 的"显示量=被检量"性质）。
  * 38 行全在体里，只标 8–12 个名（|合成 z| 最大者，且被标行至少隔 3 行）；完整名单 = Table 2。
  * 符号不可分辨的格（|Δ| ≤ 2×bootstrap SE）压成灰格：RNA 8/38、ATAC 3/38。
  * 越界格（|Δ| > 1.0）画三角，真值进注脚。
  * 右侧程序汇总以**中位数**为主视觉（均值只画空心方块，单基因一律叠点）：
    lipid/energy 均值 +0.34 全靠 PRR5L(+1.97)、myelin structure 均值 −0.59 全靠 BCAS1(−1.91)。
  * 参照行 myelin anchors 用 7 个（去掉 Cnp——它本身在 Table 2 的 38 行里，避免重复计数）。
  * 防循环：行来自 Table 2 的 ±500 kb 窗口指派，不用 scDRS 分数；行序用显示量 Δ（数据驱动，但不是分数）。

版式带（10.72 × 5.50，自上而下）：标题 8pt / 通道名 6pt / 热图体 38 行 / 注脚 5 行 6pt。
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, fcluster

CM = 1 / 2.54
WORK = os.environ.get('DSCI_WORK_ROOT', '')
BASE = WORK + "/figures/routeA_out/P0_8C"
TAB, PAN = BASE + "/tables", BASE + "/panels"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arial_style import use_arial, use_two_tier
use_arial(); use_two_tier()

MODE = (sys.argv[1] if len(sys.argv) > 1 else "delta").lower()
SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
ENA, ENB = "MOL3/6", "MOL6"
N_ENA, N_ENB, N_COP = 991, 2574, 203
LO, HI = -1.0, 1.0
CLS = {"lipid / energy": ["APOE", "PLIN3", "EFR3B", "PHGDH", "PRR5L", "ATP1A1"],
       "myelin structure": ["CNP", "BCAS1", "PLLP"],
       "cytoskeleton / MT": ["CAMSAP2", "TUBB4A", "SCHIP1", "YPEL2", "CEP97"],
       "vesicle / membrane": ["SEC11C", "PHLDB1", "ENDOD1", "CLIC4", "TMCC3"],
       "immune / adhesion": ["LITAF", "TRIM59", "SKAP2", "AIG1", "MCAM", "IGSF8", "CD9", "ADAMTS4", "SPON1"],
       "lineage / TF": ["OLIG1", "ZNF365", "CREB5", "ZBTB38"],
       "other signalling": ["DDC", "TMEM151A", "SEZ6L2", "KCTD13", "SBF1", "RRAS2"]}
g2c = {g: c for c, gs in CLS.items() for g in gs}
CCOL = {"lipid / energy": "#3f7fb0", "myelin structure": "#6f9c3f", "cytoskeleton / MT": "#b58a3c",
        "vesicle / membrane": "#8e6fa8", "immune / adhesion": "#c1544f", "lineage / TF": "#4f8f8a",
        "other signalling": "#8f8f8f"}
MCOL = ["#2f4b6e", "#7a5c9e", "#3f8f6f", "#b5651d"]
ANCH = ["Mbp", "Plp1", "Mog", "Mag", "Mobp", "Sox10", "Olig2"]

# ---------------------------------------------------------------- 数据
S = pd.read_csv(TAB + "/p08c_c1_rna_atac_by_subtype.csv")
S["cls"] = S.gene.map(g2c)
d_r0 = S.rna_delta_enrich_vs_cop.values.astype(float)
d_a0 = S.atac_delta_enrich_vs_cop.values.astype(float)
comp0 = (d_r0 / d_r0.std() + d_a0 / d_a0.std()) / 2.0

CACHE = TAB + "/p08c_c2_celllevel38.npz"
if os.path.exists(CACHE):
    zz = np.load(CACHE, allow_pickle=True)
    XG, AG, key = zz["XG"], zz["AG"], zz["key"]
    rows = {s: np.flatnonzero(key == s) for s in SUB}
    rng = np.random.default_rng(20261116)
    bs_r, bs_a = [], []
    for _ in range(500):
        idx = {s: rng.choice(rows[s], len(rows[s]), replace=True) for s in SUB}
        for M, acc in ((XG, bs_r), (AG, bs_a)):
            p = np.column_stack([M[idx[s]].mean(0) for s in SUB])
            acc.append(((p[:, 4] * N_ENA + p[:, 5] * N_ENB) / (N_ENA + N_ENB)) - p[:, 0])
    SE_R, SE_A = np.array(bs_r).std(0), np.array(bs_a).std(0)
    ok_r0, ok_a0 = np.abs(d_r0) > 2 * SE_R, np.abs(d_a0) > 2 * SE_A
    print("bootstrap SE 中位 RNA %.3f ATAC %.3f | 符号可分辨 RNA %d/38 ATAC %d/38"
          % (np.median(SE_R), np.median(SE_A), ok_r0.sum(), ok_a0.sum()))
else:
    ok_r0 = ok_a0 = np.ones(len(S), bool)
    print("[warn] 无细胞级缓存，全部按可分辨处理")

L = pd.read_csv(TAB + "/p08c_subtype_stats_long.csv")
Al = L[L.gene.isin(ANCH)].pivot_table(index="gene", columns="subtype", values=["rna_mean", "atac_mean"])
def pdelta(v):
    return ((v["MOL3/6"] * N_ENA + v["MOL6"] * N_ENB) / (N_ENA + N_ENB)) - v["COP/NFOL"]
anc_r = np.array([pdelta(Al.loc[g, "rna_mean"]) for g in ANCH])
anc_a = np.array([pdelta(Al.loc[g, "atac_mean"]) for g in ANCH])

# ---------------------------------------------------------------- 版式常量
PW, PH = 10.72, 5.50
YA, YT = 5.23, 4.95                                   # 8pt 标题基线 / 6pt 通道名基线
BODY_T, BODY_B = 4.88, 1.62                            # 热图体
NR = len(S)
LBLW, RIBW, COLW, CGAP = 1.10, 0.13, 0.86, 0.06
MODW = 0.11 if MODE == "cluster" else 0.0
x_lbl = 0.30
x_rib = x_lbl + LBLW
x_mod = x_rib + RIBW
x_rna = x_mod + MODW + 0.06
x_atac = x_rna + COLW + CGAP
X_CBAR = x_atac + COLW + 0.16
X_CB0 = X_CBAR + 0.86
RH = (BODY_T - BODY_B) / NR

# ---------------------------------------------------------------- 行序 + 标签挑选
def pick_labels(cands, comp, gap_cm=0.26):
    need = max(1, int(np.ceil(gap_cm / RH)))
    keep = []
    for i in sorted(cands, key=lambda x: -abs(comp[x])):
        if all(abs(i - j) >= need for j in keep):
            keep.append(i)
    return sorted(keep)


if MODE == "cluster":
    P = S[["rna_" + s for s in SUB]].values.astype(float)
    Pa = S[["atac_" + s for s in SUB]].values.astype(float)
    Zj = np.column_stack([(P - P.mean(1, keepdims=True)) / np.maximum(P.std(1, keepdims=True), 1e-9),
                          (Pa - Pa.mean(1, keepdims=True)) / np.maximum(Pa.std(1, keepdims=True), 1e-9)])
    mod0 = fcluster(linkage(Zj, "ward"), 4, "maxclust")
    order = np.lexsort((-d_r0, mod0))
    cand = [i for m in sorted(set(mod0)) for i in order[mod0[order] == m][:3]]
else:
    mod0 = None
    order = np.argsort(-comp0)
    cand = list(np.argsort(-np.abs(comp0)))          # 全部 38 行按 |合成 d| 排，贪心间隔取名
S = S.iloc[order].reset_index(drop=True)
mod = None if mod0 is None else mod0[order]
d_r, d_a, comp = d_r0[order], d_a0[order], comp0[order]
ok_r, ok_a = ok_r0[order], ok_a0[order]
lab_rows = pick_labels(cand, comp)
n_lab = len(lab_rows)

# ---------------------------------------------------------------- 绘制
fig = plt.figure(figsize=(PW / 2.54, PH / 2.54))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, PW); ax.set_ylim(0, PH); ax.axis("off")
CMAP = LinearSegmentedColormap.from_list("by", ["#2166ac", "#f7f7f7", "#b2182b"])

ax.text(0.75, YA, "Risk-set landing, 38 genes", fontsize=8.0, ha="left", va="baseline",
        color="#000000", fontweight="bold")

def cell_y(i):
    return BODY_T - (i + 1) * RH

for d, ok, xx, nm in ((d_r, ok_r, x_rna, "RNA"), (d_a, ok_a, x_atac, "ATAC")):
    for i in range(NR):
        v, y, h = d[i], cell_y(i), RH - 0.006
        if ok[i]:
            t = max(min(v, HI), LO)
            ax.add_patch(Rectangle((xx, y), COLW, h, fc=CMAP((t - LO) / (HI - LO)), ec="none"))
            if abs(v) > HI:
                ax.plot(xx + (COLW - 0.08 if v > 0 else 0.08), y + h / 2,
                        marker=">" if v > 0 else "<", ms=1.6, color="#ffffff", mew=0.4, zorder=6)
        else:
            ax.add_patch(Rectangle((xx, y), COLW, h, fc="#e6e6e6", ec="none"))
    ax.text(xx + COLW / 2, YT, nm, fontsize=6.0, ha="center", va="baseline", color="#333333")
    ax.plot([xx, xx + COLW, xx + COLW, xx, xx], [BODY_T, BODY_T, BODY_B, BODY_B, BODY_T], lw=0.4, color="#333333")

for i in range(NR):
    y, h = cell_y(i), RH - 0.006
    ax.add_patch(Rectangle((x_rib, y), RIBW, h, fc=CCOL[S.cls[i]], ec="none"))
    if mod is not None:
        ax.add_patch(Rectangle((x_mod, y), MODW, h, fc=MCOL[mod[i] - 1], ec="none"))
    if i in lab_rows:
        ax.text(x_lbl + LBLW - 0.04, y + h / 2, S.gene[i], fontsize=6.0, ha="right", va="center", color="#111111")

# 色标
CB_T, CB_B = BODY_T, BODY_B + 1.05
gh = np.linspace(HI, LO, 64)[:, None]
ax.imshow(np.hstack([gh, gh]), extent=(X_CBAR, X_CBAR + 0.20, CB_B, CB_T), cmap=CMAP, aspect="auto", zorder=2)
ax.add_patch(Rectangle((X_CBAR, CB_B), 0.20, CB_T - CB_B, fc="none", ec="#333333", lw=0.4))
for tv in (1.0, 0.0, -1.0):
    yy = CB_B + (tv - LO) / (HI - LO) * (CB_T - CB_B)
    ax.plot([X_CBAR + 0.20, X_CBAR + 0.24], [yy] * 2, lw=0.4, color="#333333")
    ax.text(X_CBAR + 0.27, yy, "%+.1f" % tv if tv else "0", fontsize=6.0, ha="left", va="center", color="#333333")
ax.plot(X_CBAR + 0.10, CB_B - 0.18, marker="s", ms=2.2, mfc="#e6e6e6", mec="#999999", mew=0.3)
ax.text(X_CBAR + 0.22, CB_B - 0.18, "|d|<2SE", fontsize=6.0, ha="left", va="center", color="#666666")

# 右侧程序汇总
C = []
for c, gs in CLS.items():
    m = (S.cls == c).values
    C.append(dict(cls=c, n=int(m.sum()), med_r=float(np.median(d_r[m])), med_a=float(np.median(d_a[m])),
                  mean_r=float(d_r[m].mean()), mean_a=float(d_a[m].mean())))
C = pd.DataFrame(C).sort_values("med_r", ascending=False).reset_index(drop=True)
B_T, B_B = BODY_T, BODY_B + 0.26
BH = (B_T - B_B) / (len(C) + 2)
bx0 = X_CB0 + 2.26
bx1 = PW - 0.56
sc = (bx1 - bx0) / (HI - LO)
def px(v):
    return bx0 + (min(max(v, LO), HI) - LO) * sc
ax.text(X_CB0, YT, "Programme mix", fontsize=6.0, ha="left", va="baseline", color="#333333")
ax.plot([px(0), px(0)], [B_B - 0.06, B_T], lw=0.6, color="#333333", zorder=3)
y = B_T
for _, r in C.iterrows():
    yc = y - BH / 2
    ax.text(X_CB0, yc, "%s (%d)" % (r.cls, r.n), fontsize=6.0, ha="left", va="center", color=CCOL[r.cls])
    for med, mean, col, off in ((r.med_r, r.mean_r, CCOL[r.cls], +BH * 0.17),
                                (r.med_a, r.mean_a, "#8a8a8a", -BH * 0.17)):
        ax.plot([px(0), px(med)], [yc + off] * 2, lw=2.0, color=col, solid_capstyle="butt", zorder=4)
        ax.plot(px(mean), yc + off, marker="s", ms=1.5, mfc="none", mec=col, mew=0.4, zorder=5)
    m = (S.cls == r.cls).values
    for vals, col, off in ((d_r[m], CCOL[r.cls], +BH * 0.17), (d_a[m], "#8a8a8a", -BH * 0.17)):
        for v in vals:
            if abs(v) <= HI:
                ax.plot(px(v), yc + off, marker=".", ms=1.0, color=col, zorder=6)
            else:
                ax.plot(bx1 + 0.10 if v > 0 else bx0 - 0.10, yc + off, marker=">" if v > 0 else "<", ms=1.5,
                        mfc="none", mec=col, mew=0.4, zorder=6)
    y -= BH
ax.text(X_CB0, B_B + BH / 2, "myelin anchors* (7)", fontsize=6.0, ha="left", va="center", color="#333333")
for med, col, off in ((np.median(anc_r), "#333333", +BH * 0.17), (np.median(anc_a), "#8a8a8a", -BH * 0.17)):
    ax.plot([px(0), px(med)], [B_B + BH / 2 + off] * 2, lw=2.0, color=col, solid_capstyle="butt", zorder=4)
for v, col, off in ((anc_r, "#333333", +BH * 0.17), (anc_a, "#8a8a8a", -BH * 0.17)):
    for u in v:
        ax.plot(px(u), B_B + BH / 2 + off, marker=".", ms=1.0, color=col, zorder=6)
ax.plot([bx0, bx1], [B_B - 0.10, B_B - 0.10], lw=0.4, color="#999999")
for tv in (-1.0, 0.0, 1.0):
    ax.text(px(tv), B_B - 0.15, "%+.1f" % tv if tv else "0", fontsize=6.0, ha="center", va="top", color="#333333")

# 注脚（通栏 5 行，全部左对齐 x=0.30）
CL = pd.read_csv(TAB + "/p08c_c1_class_effect.csv")
pmin = min(CL.p_med_rna.min(), CL.p_med_atac.min())
notes = ["\u0394 = enriched (MOL3/6+MOL6, cell-weighted) \u2212 COP/NFOL, log2; blue = lower in enriched",
         "all 38 Table 2 genes shown; %d named (largest |\u0394|, \u22653 rows apart)" % n_lab,
         "one scale for both channels: ATAC shifts are ~2\u00d7 smaller than RNA (p95 0.31 vs 0.76)",
         "off-scale: PRR5L +1.97, BCAS1 \u22121.91, CAMSAP2 \u22120.92, YPEL2 +0.90, RRAS2 \u22120.80",
         "class permutation within the 38: min P = %.3f, none survives BH; order = \u0394" % pmin]
if MODE == "cluster":
    notes[1] = "rows = 4 expression modules (ward, joint RNA+ATAC z); %d named" % n_lab
    notes[4] = "class permutation min P = %.3f, none survives BH; order = modules" % pmin
assert len(notes) == 5, len(notes)
for k, t in enumerate(notes):
    ax.text(x_lbl, 1.44 - 0.26 * k, t, fontsize=6.0, ha="left", va="top", color="#666666")

out = PAN + "/fig6v29_C_merged_%s.pdf" % MODE
fig.savefig(out, dpi=600)
print("[saved] %s | 体 x=%.2f-%.2f 行距 %.3f | 色标 %.2f | 右块 %.2f-%.2f | 标签 %d"
      % (MODE, x_rna, x_atac + COLW, RH, X_CBAR, X_CB0, PW - 0.30, n_lab))
import fitz
d = fitz.open(out); p = d[0]
print("尺寸 %.2f × %.2f cm" % (p.rect.width * 2.54 / 72, p.rect.height * 2.54 / 72)); d.close()
