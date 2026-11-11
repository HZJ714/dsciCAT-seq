#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p0_8c_fig6c_top12.py —— Fig6C v25 主图版：top 12 被检风险基因 + 4 个髓鞘结构锚基因，RNA|ATAC 双热图

为什么是 12 不是 20（行距预算，不是挑肥拣瘦）：
  主图第二行槽位高 5.50 cm（与 A/B 齐平，拼页后页高仍是 11.31，不动 v23 的版面），
  正文净高 3.85 cm ⇒ 每行 0.24 cm；6 pt 行标签字高 ~0.21 cm ⇒ 只有"每行都标名"这一种排法，
  能标全的行数上限就是 16 行。20 行会被挤成隔行标，反而更像 cherry-picking。
  ⇒ 12 风险基因 + 4 髓鞘锚 = 16 行，全部带名字。全 223 行的透明版留附图（v24）。

统计口径（P0_8C_RECEIPT §8.7 之后）：
  行序 = **消融（LOO）r 降序**，即把该基因从 MS 基因集里删掉、同一 scdrs compute-score 参数
  重打分之后，再与该基因表达求 Pearson。交付口径 r（含自己）只作浅/深灰长条叠在同一轴上，
  用来让读者看见"长条被机械项抬高、短条才是真值"。红点 = ATAC 通道 r（外部量，无自我 inclusion）。
  注意：LOO 后没有任何基因过 BH FDR（最小 FDR≈0.79），黑色短条只是**名义 MC P<0.05**，
  图注必须这么写，不能写成"显著"。

4 个锚基因（PLP1/MBP/MOG/MAG）不在 MS 基因集内 ⇒ 它们的 r 天然无需消融，是全图唯一干净的读数。

口径：主池 all6_pool（与 Fig6A 同源同序）；高亮列 = Fig6A group_q<0.05（MOL3/6、MOL6）。
纪律：字号两档 8.0/6.0；文字不压数据；色标不得压热图（v24 首版踩过，见 §8.8）。
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

SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
NTOP, ANCH = 12, ["Plp1", "Mbp", "Mog", "Mag"]
# 功能标签：为读者定向用的策展注释，不进任何统计（图注里要写"curated for orientation"）
TAG = {"Prr5l": "?", "Edil3": "ECM", "Arsg": "lipid", "Kctd13": "K+ chan.",
       "Vmp1": "ER/aut.", "Litaf": "TNF", "Jak2": "JAK-STAT", "Rps6ka5": "kinase",
       "Gatm": "creatine", "Abca2": "lipid", "Smarcc2": "chrom.", "Abca1": "lipid",
       "Plp1": "myelin", "Mbp": "myelin", "Mog": "myelin", "Mag": "myelin"}

L = pd.read_csv(TAB + "/p08c_subtype_stats_long.csv")
K = pd.read_csv(TAB + "/p08c_genelevel_all6_pool.csv")
O = pd.read_csv(TAB + "/p08c_genelevel_LOOexact.csv")
D = pd.read_csv(TAB + "/p08c_atac_genelevel_all6_pool.csv")
AN = pd.read_csv(TAB + "/p08c_anchor_genelevel.csv")
meta = json.load(open(TAB + "/p08c_meta.json"))
ENRICH = meta["enriched_by_group_q"]

M = K.merge(O[["gene", "r_LOO", "mc_p_LOO", "fdr_LOO"]], on="gene").merge(
    D[["gene", "r_atac", "fdr_atac"]], on="gene", how="left")
M = M.sort_values("r_LOO", ascending=False).head(NTOP).reset_index(drop=True)
A4 = AN.set_index("gene").reindex(ANCH).reset_index()
nT, nA = len(M), len(A4)
genes = list(M.gene) + list(A4.gene)
print("[行序] LOO r 降序 top%d：%s" % (nT, ", ".join(M.gene_human.astype(str))))
print("      锚块：%s（不在 MS 基因集内 ⇒ r 无需消融）" % ", ".join(A4.gene_human.astype(str)))

# ---------------------------------------------------------------- 版式（cm）
FW, FH = 10.72, 5.50                       # 与 A/B 同高，拼页后页高不变
BODY_B, BODY_H = 1.00, 3.80                # 16 行 ⇒ 0.238 cm/行；下沿 1.00 cm 才够 45° 亚型名（实测需 0.96）
XLAB, WLAB = 0.06, 2.08
XR, WR = 2.44, 1.65   # 右移 0.10：功能标签实测最宽到 2.33 cm，与统计栏左沿要留缝
XH1, WH = 4.40, 2.55
XH2 = 7.15
XCB, WCB = 9.88, 0.26                       # 色标左沿 9.88 > ax2 右沿 9.70，不压图
YA, YB = 5.22, 4.88                         # 8 pt 标题 / 6 pt 副行（正文顶 4.80；5.40 会把标题上伸量顶出页边）

def axrec(x, w):
    return [x / FW, BODY_B / FH, w / FW, BODY_H / FH]
fig = plt.figure(figsize=(FW * CM, FH * CM))
axL = fig.add_axes(axrec(XLAB, WLAB)); axL.axis("off")
axR = fig.add_axes(axrec(XR, WR))
ax1 = fig.add_axes(axrec(XH1, WH)); ax2 = fig.add_axes(axrec(XH2, WH))
YMAX = nT + nA - 0.5
for a in (axR, ax1, ax2):
    a.set_ylim(YMAX, -0.5); a.tick_params(length=0)
def cx(x, w):
    return (x + w / 2) / FW
def cy(y):
    return y / FH
def ry(r):
    return 1.0 - (r + 0.5) / (nT + nA)

P = L.pivot_table(index="gene", columns="subtype", values="rna_mean").reindex(genes)[SUB]
Q = L.pivot_table(index="gene", columns="subtype", values="atac_mean").reindex(genes)[SUB]
def rowz(Dm):
    return ((Dm.T - Dm.mean(1)) / (Dm.std(1, ddof=0).replace(0, np.nan))).T
V = 2.0
for ax, Z in ((ax1, rowz(P)), (ax2, rowz(Q))):
    ax.imshow(Z.values, cmap="coolwarm", vmin=-V, vmax=V, aspect="auto", zorder=1)
    ax.set_xlim(-0.5, 5.5); ax.set_yticks([])
    ax.set_xticks(range(6)); ax.set_xticklabels(SUB, fontsize=6.0, rotation=45, ha="right")
    for j, s in enumerate(SUB):
        if s in ENRICH:
            ax.add_patch(Rectangle((j - 0.5, -0.5), 1, nT + nA, fill=False, ec="0.15", lw=0.8, zorder=4))
    for sp_ in ax.spines.values():
        sp_.set_linewidth(0.5)
for i in range(nA):
    for ax in (ax1, ax2):
        ax.add_patch(Rectangle((-0.5, nT + i - 0.5), 12, 1, color="0.92", zorder=0.2))
    axR.add_patch(Rectangle((-10, nT + i - 0.5), 40, 1, color="0.92", zorder=0.2))
for ax in (axR, ax1, ax2):
    ax.axhline(nT - 0.5, color="0.15", lw=0.8, zorder=6)

# ---------------------------------------------------------------- 统计栏
r_test = np.array(list(M.r) + list(A4.r), float)
r_loo = np.array(list(M.r_LOO) + list(A4.r), float)      # 锚块：分数里没有它 ⇒ 交付 r 即消融 r
r_atac = np.array(list(M.r_atac) + list(A4.r_atac), float)
sig_test = np.array(list(M.fdr_bh < 0.05) + [False] * nA)
sig_loo = np.array(list(M.mc_p_LOO < 0.05) + list(A4.mc_p.values < 0.05))
y = np.arange(nT + nA)
axR.barh(y, r_test, height=0.90, color=np.where(sig_test, "0.42", "0.80"), edgecolor="none", zorder=2)
axR.barh(y, r_loo, height=0.50, color=np.where(sig_loo, "0.05", "0.35"), edgecolor="none", zorder=3)
axR.plot(r_atac, y, ls="none", marker="o", ms=1.7, color="#B03030", zorder=5)
axR.axvline(0, color="0.2", lw=0.6, zorder=6)
axR.set_xlim(-0.06, 0.14); axR.set_ylim(YMAX, -0.5)
axR.set_xticks([0]); axR.set_xticklabels(["0"], fontsize=6.0)
axR.set_yticks([]); axR.tick_params(length=1.4, pad=1.0)
for sp_ in axR.spines.values():
    sp_.set_linewidth(0.5)

# ---------------------------------------------------------------- 文字
fig.text(cx(XR, WR), cy(YA), "r vs MS score", fontsize=8.0, ha="center", va="baseline")
fig.text(cx(XH1, WH), cy(YA), "RNA  log2(CP10k+1)", fontsize=8.0, ha="center", va="baseline")
fig.text(cx(XH2, WH), cy(YA), "ATAC  gene activity", fontsize=8.0, ha="center", va="baseline")
fig.text(cx(XR, WR), cy(YB), "grey = as tested, dark = LOO", fontsize=6.0, ha="center", va="baseline")
# 用绝对 cm 坐标：符号右对齐 0.75、标签左对齐 0.85 ⇒ RPS6KA5(0.65 cm) 起点 0.10 不再被页左边裁掉
def fy(r):
    return (BODY_B + (1.0 - (r + 0.5) / (nT + nA)) * BODY_H) / FH
for i, g in enumerate(genes):
    sym = str(M.gene_human.values[i]) if i < nT else str(A4.gene_human.values[i - nT])
    fig.text(1.28 / FW, fy(i), sym, fontsize=6.0, style="italic", ha="right", va="center", color="0.1")
    fig.text(1.36 / FW, fy(i), TAG.get(g, ""), fontsize=6.0, ha="left", va="center", color="0.35")

axc = fig.add_axes([XCB / FW, BODY_B / FH + BODY_H / FH * 0.28,
                    WCB / FW, BODY_H / FH * 0.44])
cb = fig.colorbar(matplotlib.cm.ScalarMappable(matplotlib.colors.Normalize(-V, V), "coolwarm"), cax=axc)
cb.set_ticks([-V, 0, V]); cb.ax.tick_params(labelsize=6.0, length=1.4, pad=1.0)
cb.outline.set_linewidth(0.5)
fig.text((XCB + WCB / 2) / FW, (BODY_B + BODY_H * 0.72 + 0.10) / FH, "z",
         fontsize=6.0, ha="center", va="baseline")

out = PAN + "/fig6v25_C_dualheatmap_top12.pdf"
fig.savefig(out)
print("[写出] %s  %.2f x %.2f cm | 行 %d（top%d + 锚 %d）| 高亮 %s"
      % (os.path.basename(out), FW, FH, nT + nA, nT, nA, "+".join(ENRICH)))
print(" 风险行：中位 r 交付 %.4f → LOO %.4f | 名义 P<0.05 %d/%d | FDR<0.05 %d/%d | ATAC 中位 %.4f"
      % (M.r.median(), M.r_LOO.median(), int((M.mc_p_LOO < 0.05).sum()), nT,
         int((M.fdr_LOO < 0.05).sum()), nT, M.r_atac.median()))
print(" 锚块：中位 r %.4f | 名义 P<0.05 %d/4 | ATAC 中位 %.4f" % (A4.r.median(), int((A4.mc_p < 0.05).sum()), A4.r_atac.median()))
