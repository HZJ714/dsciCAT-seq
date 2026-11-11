#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
P0-8c Fig6C 双热图 mockup（先不拼进 Fig6 页面，单独出片看效果）

规格 §8.3 P0-8c v3：
  行 = 富集基因（按 MC P 升序，同分按 r 降序）取 top N + 底部锚块 8 行（髓鞘经典基因，C∩E=∅）
  列 = 6 亚型成熟轴 COP/NFOL → MOL6（与 6A 同序同源同归一化）
  左 = RNA log2(CP10k+1) 亚型均值，右 = ATAC gene activity 同名单均值，均跨亚型缩放（row-wise z）
  富集基因名旁给 MC P 星标；富集亚型列高亮
纪律：措辞用 concordance / decomposition，禁 driver。
"""
import os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

CM = 1 / 2.54          # inch per cm（figsize 用英寸）
WORK = os.environ.get('DSCI_WORK_ROOT', '')
BASE = WORK + "/figures/routeA_out/P0_8C"
TAB, PAN = BASE + "/tables", BASE + "/panels"
SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
N_TOP = 10
POOL = os.environ.get("HEAT_POOL", "all6_pool")          # enriched_pool | all6_pool
sys_path = os.path.dirname(os.path.abspath(__file__))

import sys
sys.path.insert(0, sys_path)
from arial_style import use_arial, use_two_tier
use_arial(); use_two_tier()

L = pd.read_csv(TAB + "/p08c_subtype_stats_long.csv")
G = pd.read_csv(TAB + "/p08c_genelevel_%s.csv" % POOL)
A = pd.read_csv(TAB + "/p08c_anchor_block.csv")
meta = json.load(open(TAB + "/p08c_meta.json"))
ENRICH = meta["enriched_by_group_q"]                     # 高亮列 = Fig6A group_q<0.05 的亚型

sig = G.sort_values(["mc_p", "r"], ascending=[True, False]).head(N_TOP)
genes = list(sig.gene) + list(A.gene)
STAR = lambda p: "***" if p < 0.001 else ("**" if p < 0.01 else "*")
labels = ["%s%s" % (h, STAR(p)) for h, p in zip(sig.gene_human, sig.mc_p)] + \
         [str(h) for h in A.gene_human]
is_anchor = [False] * N_TOP + [True] * len(A)

P = L.pivot_table(index="gene", columns="subtype", values="rna_mean").reindex(genes)[SUB]
Q = L.pivot_table(index="gene", columns="subtype", values="atac_mean").reindex(genes)[SUB]
def rowz(D):
    return ((D.T - D.mean(1)) / (D.std(1, ddof=0).replace(0, np.nan))).T
Zr, Za = rowz(P).values, rowz(Q).values

V = 2.0
fig = plt.figure(figsize=(10.72 * CM, 5.50 * CM))
LM, TOP, BOT = 0.175, 0.085, 0.178   # 页高 5.30→5.50：18 行行距 ≥6.5pt 且 45° 列名不裁（§9 加呼吸）   # BOT↑：45° 列名竖直净空；TOP 留标题（§9 加呼吸不缩字）
gap = 0.055
w1, w2 = 0.30, 0.30
ax1 = fig.add_axes([LM, BOT, w1, 1 - TOP - BOT]); ax2 = fig.add_axes([LM + w1 + gap, BOT, w2, 1 - TOP - BOT])
for ax, Z, ttl in ((ax1, Zr, "RNA  log2(CP10k+1)"), (ax2, Za, "ATAC  gene activity")):
    ax.imshow(Z, cmap="coolwarm", vmin=-V, vmax=V, aspect="auto", zorder=1)
    ax.set_xticks(range(6)); ax.set_xticklabels(SUB, fontsize=6.0, rotation=45, ha="right")
    ax.set_yticks(range(len(genes))); ax.set_yticklabels([])
    for k, j in enumerate([SUB.index(s) for s in ENRICH]):          # 富集亚型列高亮
        ax.add_patch(Rectangle((j - 0.5, -0.5), 1, len(genes), fill=False, edgecolor="0.2", lw=0.9))
    for i in range(len(genes) + 1):
        ax.axhline(i - 0.5, color="w", lw=0.35, zorder=2)
    ax.set_title(ttl, fontsize=8.0, pad=3)
    ax.tick_params(length=0)
    for sp_ in ax.spines.values():
        sp_.set_linewidth(0.5)
    ax.set_xlim(-0.5, 5.5); ax.set_ylim(len(genes) - 0.5, -0.5)
ax1.set_yticklabels(labels, fontsize=6.0)
for t, s, p in zip(ax1.get_yticklabels(), is_anchor, list(sig.mc_p) + [np.nan] * len(A)):
    t.set_style("italic" if not s else "normal")
    t.set_color("0.15" if not s else "0.45")
for i in range(len(genes)):                                          # 锚块行灰底（§8.3：底部两行解释性小字删掉，交图注）
    if is_anchor[i]:
        for ax in (ax1, ax2):
            ax.add_patch(Rectangle((-0.5, i - 0.5), 12, 1, color="0.92", zorder=0.2))
sep = N_TOP - 0.5
for ax in (ax1, ax2):
    ax.axhline(sep, color="0.2", lw=0.9, zorder=3)
cb = fig.colorbar(matplotlib.cm.ScalarMappable(matplotlib.colors.Normalize(-V, V), "coolwarm"),
                  ax=[ax1, ax2], fraction=0.035, pad=0.03)
cb.set_ticks([-V, 0, V]); cb.ax.tick_params(labelsize=6.0, length=1.6, pad=1.0); cb.outline.set_linewidth(0.5)
cb.set_label("scaled across subtypes (row-wise z)", fontsize=6.0, labelpad=2)
out = PAN + "/fig6v23_C_dualheatmap_%s.pdf" % ("enrichpool" if POOL == "enriched_pool" else "all6pool")
fig.savefig(out); plt.close(fig)
print("写出 %s | 检验池 %s | 行 %d（富集 %d + 锚 %d）| 高亮列 %s"
      % (os.path.basename(out), POOL, len(genes), N_TOP, len(A), "+".join(ENRICH)))
print(sig[["gene_human", "mc_p", "mc_z", "fdr_bh", "r", "mean_pct"]].round(4).to_string(index=False))
