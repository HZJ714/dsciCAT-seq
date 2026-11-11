#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
P0-8c S6C mockup — Fig6C 两口径稳健性对照（v4 角色互换后）

规格 FIG6_7_SYNTHESIS_RECEIPT_20261112.md / FIG7_DRUG_HEATMAP_RECEIPT_20261112.md §8.3 P0-8c v4：
  S6C = Fig6C 两口径稳健性 —— 富集池（MOL3/6+MOL6）名单/排序 vs 主图全池
  top20 交集 13/20；top10 沿成熟轴升 8/10（全池）vs 3/10（富集池）
  = 定义域差异可视证据 + 亚型混杂可见性；数据 = P0-8c 副产品，零新计算。

设计：同 6 列、同 row-wise z（z 只依赖亚型均值，与检验池无关）——
      两池唯一的差别是"检验定义域选出的基因名单/排序"。
      左组 = 全池 top10（主图口径，沿轴升），右组 = 富集池 top10（对照，方向与轴脱耦）。
      每组 RNA + ATAC 两块，与 Fig6C 同型，读者一眼可对比。
纪律：措辞 complementary（禁 upgrade/driver）；r 为 marginal Pearson；两档字号 8/6。
版式：cm 精算（上一版分数布局 LEFT+4*(w+gap)=1.50 越页，旋转列名被裁，已改）。
"""
import os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.insert(0, HERE)
from arial_style import use_arial, use_two_tier, CM
use_arial(); use_two_tier()

WORK = os.environ.get('DSCI_WORK_ROOT', '')
TAB = WORK + "/figures/routeA_out/P0_8C/tables"
PAN = WORK + "/figures/routeA_out/P0_8C/panels"
SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
N_TOP = 10

L = pd.read_csv(TAB + "/p08c_subtype_stats_long.csv")
meta = json.load(open(TAB + "/p08c_meta.json"))
subm = L.pivot_table(index="gene", columns="subtype", values="rna_mean")
atacm = L.pivot_table(index="gene", columns="subtype", values="atac_mean")

def rowz(D):
    return ((D.T - D.mean(1)) / (D.std(1, ddof=0).replace(0, np.nan))).T
zr, za = rowz(subm), rowz(atacm)

STAR = lambda p: "***" if p < 0.001 else ("**" if p < 0.01 else "*")

def pool_top10(pool):
    G = pd.read_csv(TAB + "/p08c_genelevel_%s.csv" % pool)
    return G.sort_values(["mc_p", "r"], ascending=[True, False]).head(N_TOP).reset_index(drop=True)

A6, EP = pool_top10("all6_pool"), pool_top10("enriched_pool")
axis_up = lambda df: int(sum(zr.loc[g]["MOL6"] > zr.loc[g]["COP/NFOL"] for g in df.gene))

# ---------------------------------------------------------------- 版式（cm 精算）
LAB_W, HM_W, GAP_IN, GAP_GRP = 1.35, 1.45, 0.10, 1.40
TOP_CM, BOT_CM = 1.18, 1.42          # 上：组标题(8pt)+双行通道名(6pt)；下：45° 列名 + 注释行          # 上留标题、下留 45° 列名 + 一行注释
CB_X, CB_W = 0.25, 0.42
FW = LAB_W + HM_W + GAP_IN + HM_W + GAP_GRP + HM_W + GAP_IN + HM_W + CB_X + CB_W + 0.15
FH = TOP_CM + (FH_HM := 3.80) + BOT_CM
fig = plt.figure(figsize=(FW * CM, FH * CM))

def ax_at(x_cm):
    return fig.add_axes([x_cm / FW, (BOT_CM) / FH, HM_W / FW, FH_HM / FH])

x0 = LAB_W
axl, axl2 = ax_at(x0), ax_at(x0 + HM_W + GAP_IN)
xg = x0 + HM_W + GAP_IN + HM_W + GAP_GRP
axr, axr2 = ax_at(xg), ax_at(xg + HM_W + GAP_IN)
axcb = fig.add_axes([xg + HM_W + GAP_IN + HM_W + CB_X, BOT_CM + 0.55, CB_W, FH_HM - 1.10])

def draw(ax, df, chan):
    genes = list(df.gene)
    Z = (zr if chan == "RNA" else za).reindex(genes)[SUB].values
    ax.imshow(Z, cmap="coolwarm", vmin=-2, vmax=2, aspect="auto", zorder=1)
    ax.set_xticks(range(6)); ax.set_xticklabels(SUB, fontsize=6.0, rotation=45, ha="right")
    ax.set_yticks(range(len(genes))); ax.set_yticklabels([])
    for j in (4, 5):                                   # 富集亚型列高亮 MOL3/6、MOL6
        ax.add_patch(Rectangle((j - 0.5, -0.5), 1, len(genes), fill=False, edgecolor="0.2", lw=0.7))
    for i in range(len(genes) + 1):
        ax.axhline(i - 0.5, color="w", lw=0.3, zorder=2)
    ax.tick_params(length=0, pad=1.6)
    for sp_ in ax.spines.values():
        sp_.set_linewidth(0.5)
    ax.set_xlim(-0.5, 5.5); ax.set_ylim(len(genes) - 0.5, -0.5)

draw(axl, A6, "RNA"); draw(axl2, A6, "ATAC")
draw(axr, EP, "RNA"); draw(axr2, EP, "ATAC")

for ax, df in ((axl, A6), (axr, EP)):                   # 基因名 + MC P 星标（斜体=富集基因）
    labels = ["%s%s" % (h, STAR(p)) for h, p in zip(df.gene_human, df.mc_p)]
    ax.set_yticklabels(labels, fontsize=6.0)
    for t in ax.get_yticklabels():
        t.set_style("italic"); t.set_color("0.15")

# 标题：组名 8pt 居中跨两块；通道名 6pt 双行、只在左组写一次（§8.3 并列面板只写一次）
for ax, txt in ((axl, "all 6 subtypes  (main)"), (axr, "MOL3/6 + MOL6  (enriched pool)")):
    x0 = ax.get_position().x0 * FW
    w2 = (HM_W * 2 + GAP_IN) / 2.0
    fig.text((x0 + w2) / FW, (BOT_CM + FH_HM + 0.72) / FH, txt,
             ha="center", va="bottom", fontsize=8.0, color="0.1")
for ax, txt in ((axl, "RNA\nlog2(CP10k+1)"), (axl2, "ATAC\ngene activity")):
    cx = (ax.get_position().x0 + ax.get_position().width / 2) * FW
    fig.text(cx / FW, (BOT_CM + FH_HM + 0.06) / FH, txt,
             ha="center", va="bottom", fontsize=6.0, color="0.25", linespacing=1.15)

cb = fig.colorbar(matplotlib.cm.ScalarMappable(matplotlib.colors.Normalize(-2, 2), "coolwarm"),
                  cax=axcb)
cb.set_ticks([-2, 0, 2]); cb.ax.tick_params(labelsize=6.0, length=1.4, pad=1.6)
cb.outline.set_linewidth(0.5)
cb.set_label("row-wise z", fontsize=6.0, labelpad=2)

fig.text(0.005, 0.020,
         "top20 overlap %d/20  ·  top10 along maturation axis %d/10 vs %d/10  ·  "
         "RNA×ATAC ρ median 0.71 / 0.69"
         % (meta["top20_two_pool_overlap"], axis_up(A6), axis_up(EP)),
         fontsize=6.0, color="0.2", ha="left", va="bottom")

out = PAN + "/figS6v13_C_robustness_mockup.pdf"
fig.savefig(out); plt.close(fig)
print("写出 %s | 页面 %.2f x %.2f cm" % (os.path.basename(out), FW, FH))
print("all6   top10:", " ".join(A6.gene_human), "| axis-up", axis_up(A6))
print("enrich top10:", " ".join(EP.gene_human), "| axis-up", axis_up(EP))
