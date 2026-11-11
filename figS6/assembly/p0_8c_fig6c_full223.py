#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p0_8c_fig6c_full223.py - Fig6C v24: all 223 tested genes + the 8 anchor block, RNA|ATAC dual heatmap
                        (following the Fig 4B precedent: the heatmap body carries all genes; only a few names are drawn)

Three differences from the old v23 (top-10 only), all from the P0_8C_RECEIPT §8.4/§8.6/§8.7 audit:
 1. rows = ALL 223 tested genes, ordered by the test statistic r descending, so the displayed and tested quantities are one and the same
    (fix 8.4b: the old heatmap drew subtype means while the test measured per-cell co-variation, and only 10 rows were picked).
 2. each statistics bar stacks three layers: the light/dark grey long bar = the delivered-convention r (as tested); the dark short bar = r after the gene is
     REMOVED FROM ITS GENE SET and the score recomputed with identical scdrs parameters (the c2_logo ablation product, LOO);
     the red dot = the ATAC-channel r (the scores are constructed from RNA; ATAC is an external quantity, so no self-inclusion by construction).
     The face itself shows "long bars +0.04, short bars collapsed to zero" - which is §8.7's verdict: a mechanical share of 101%.
 3. labels go only to genes with MS/immune-literature names (greedy + MIN_GAP=11 rows against collisions, 6 pt italic);
    the 8 myelin/OL anchor genes get a whole grey block, not per-gene labels (the list goes to the caption). The anchors sit outside the MS gene set,
    so their r reads directly without ablation - the figure's only clean gene-level evidence (RNA median +0.021, 1/8 nominal P<0.05).

Convention: the main pool all6_pool (same source and order as Fig 6A); highlighted columns = Fig6A group_q<0.05 (MOL3/6, MOL6).
Discipline: two tiers 8.0/6.0; text never on data; wording concordance / decomposition, "driver" banned.
Row-pitch budget: body 5.60 cm / 231 rows = 0.0242 cm/row; a 6 pt label is ~0.254 cm tall, so >= 11 rows of spacing are needed.
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
L = pd.read_csv(TAB + "/p08c_subtype_stats_long.csv")
K = pd.read_csv(TAB + "/p08c_genelevel_all6_pool.csv")        # the 223 tested (delivered convention)
O = pd.read_csv(TAB + "/p08c_genelevel_LOOexact.csv")         # the ablation convention (exact LOO)
D = pd.read_csv(TAB + "/p08c_atac_genelevel_all6_pool.csv")   # ATAC-channel per-gene
AN = pd.read_csv(TAB + "/p08c_anchor_genelevel.csv")          # the 8 anchors
meta = json.load(open(TAB + "/p08c_meta.json"))
ENRICH = meta["enriched_by_group_q"]

M = K.merge(O[["gene", "r_LOO", "mc_p_LOO"]], on="gene").merge(
    D[["gene", "r_atac", "fdr_atac"]], on="gene", how="left")
M = M.sort_values("r", ascending=False).reset_index(drop=True)
genes = list(M.gene) + list(AN.gene)
nT, nA = len(M), len(AN)

# ---------------------------------------------------------------- label placement (greedy + minimum row pitch)
order = {g: i for i, g in enumerate(M.gene)}
PRI = ["Jak2", "Tmem132a", "Stat3", "Nfkb1", "Bach2", "Ptpn2", "Ubash3b", "Apoe",
       "Smarca4", "Atxn2", "Hdac7", "Abca1", "Arsg", "Vmp1", "Litaf", "Jak1", "Klf3", "Prr5l"]
MIN_GAP = 11
lab_row, lab_txt = [], []
for g in PRI:
    r = order.get(g)
    if r is None:
        continue
    if all(abs(r - x) >= MIN_GAP for x in lab_row):
        lab_row.append(r); lab_txt.append(str(M.gene_human.values[r]))
print("[labels] candidates %d -> placed %d (MIN_GAP=%d, positions %s): %s"
      % (len(PRI), len(lab_txt), MIN_GAP, sorted(lab_row), ", ".join(lab_txt)))

# ---------------------------------------------------------------- layout (cm)
FW, FH = 10.72, 7.75   # page height +0.15: 45° subtype names need 1.0 cm clearance, or "COP/NFOL"'s bottom gets clipped
BODY_B, BODY_H = 1.00, 5.60
XLAB, WLAB = 0.06, 2.08
XR, WR = 2.34, 1.75
XH1, WH = 4.40, 2.55
XH2 = 7.15
XCB, WCB = 9.88, 0.26   # v24's first colour bar sat at 10.00, pressing on the ATAC heatmap's right edge (10.25);
                     # narrowing the heatmap and moving the bar out put its left edge at 9.88 > ax2's 9.70, fully clear of the plot
YTOP = BODY_B + BODY_H                      # 6.46
YA, YB, YC = 7.45, 7.10, 6.78   # 8pt title / 6pt sub-line / 6pt third line (body top 6.60, 0.18 clearance)                         # 8pt title baseline / 6pt sub-line baseline
def axrec(x, w):                     # cm -> figure fraction (FW/FH are already in cm; do not multiply by CM again)
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
def ry(r):                            # row index -> label-axis y (row 0 at the top, matching the heatmap)
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
    for j, s in enumerate(SUB):                      # highlight the enriched-subtype columns (Fig6A group_q<0.05)
        if s in ENRICH:
            ax.add_patch(Rectangle((j - 0.5, -0.5), 1, nT + nA, fill=False, ec="0.15", lw=0.8, zorder=4))
    for sp_ in ax.spines.values():
        sp_.set_linewidth(0.5)
for i in range(nA):                                  # anchor block grey background
    for ax in (ax1, ax2):
        ax.add_patch(Rectangle((-0.5, nT + i - 0.5), 12, 1, color="0.92", zorder=0.2))
    axR.add_patch(Rectangle((-10, nT + i - 0.5), 40, 1, color="0.92", zorder=0.2))
for ax in (axR, ax1, ax2):
    ax.axhline(nT - 0.5, color="0.15", lw=0.8, zorder=6)

# ---------------------------------------------------------------- the statistics bar (three layers on one axis)
r_test = np.array(list(M.r) + list(AN.r), float)
r_loo = np.array(list(M.r_LOO) + [np.nan] * nA, float)
r_atac = np.array(list(M.r_atac) + list(AN.r_atac), float)
sig_test = np.array(list(M.fdr_bh < 0.05) + [False] * nA)
sig_loo = np.array(list(M.mc_p_LOO < 0.05) + list(AN.mc_p.values < 0.05))
y = np.arange(nT + nA)
axR.barh(y, r_test, height=0.92, color=np.where(sig_test, "0.42", "0.80"), edgecolor="none", zorder=2)
axR.barh(y, r_loo, height=0.52, color=np.where(sig_loo, "0.05", "0.35"), edgecolor="none", zorder=3)
for i in range(nA):                                  # anchors: never in the scores, so the delivered r IS the "ablated" r
    axR.barh(nT + i, AN.r.values[i], height=0.52, color="0.35", edgecolor="none", zorder=3)
axR.plot(r_atac, y, ls="none", marker="o", ms=1.4, color="#B03030", zorder=5)
axR.axvline(0, color="0.2", lw=0.6, zorder=6)
axR.set_xlim(-0.06, 0.14); axR.set_ylim(YMAX, -0.5)
axR.set_xticks([0]); axR.set_xticklabels(["0"], fontsize=6.0)   # a 0.10 tick would fight the 45° subtype names
axR.set_yticks([]); axR.tick_params(length=1.4, pad=1.0)
for sp_ in axR.spines.values():
    sp_.set_linewidth(0.5)

# ---------------------------------------------------------------- text (all outside the data area)
fig.text(cx(XR, WR), cy(YA), "r vs MS score", fontsize=8.0, ha="center", va="baseline")
fig.text(cx(XH1, WH), cy(YA), "RNA  log2(CP10k+1)", fontsize=8.0, ha="center", va="baseline")
fig.text(cx(XH2, WH), cy(YA), "ATAC  gene activity", fontsize=8.0, ha="center", va="baseline")
fig.text(cx(XR, WR), cy(YC), "+0.040  \u2192  \u22120.000 (LOO)", fontsize=6.0, ha="center", va="baseline")
# the panel's own notes all move to the statistics-bar title band: the top-left is reserved for the compose-time letter C (v24's first version was pressed by C)
fig.text(cx(XR, WR), cy(YB), "all 223 tested, r desc", fontsize=6.0, ha="center", va="baseline")
axL.text(0.02, ry(nT + nA / 2 - 0.5), "8 anchors", fontsize=6.0, transform=axL.transAxes,
         ha="left", va="center", color="0.2")
for r, t in zip(lab_row, lab_txt):
    axL.text(0.985, ry(r), t, fontsize=6.0, style="italic", ha="right",
             va="center", transform=axL.transAxes, color="0.1")

axc = fig.add_axes([XCB / FW, BODY_B / FH + BODY_H / FH * 0.30,
                    WCB / FW, BODY_H / FH * 0.40])
cb = fig.colorbar(matplotlib.cm.ScalarMappable(matplotlib.colors.Normalize(-V, V), "coolwarm"), cax=axc)
cb.set_ticks([-V, 0, V]); cb.ax.tick_params(labelsize=6.0, length=1.4, pad=1.0)
cb.outline.set_linewidth(0.5)
fig.text((XCB + WCB / 2) / FW, (BODY_B + BODY_H * 0.70 + 0.10) / FH, "z",
         fontsize=6.0, ha="center", va="baseline")

out = PAN + "/fig6v24_C_dualheatmap_full223.pdf"
fig.savefig(out)   # no tight: keep the 10.72 x 7.60 cm slot; plt.close(fig)
print("[written] %s | rows %d (tested %d + anchors %d) | highlighted %s" % (os.path.basename(out), nT + nA, nT, nA, "+".join(ENRICH)))
print(" median r: delivered %.4f | LOO %.4f | ATAC %.4f | anchors delivered %.4f"
      % (M.r.median(), M.r_LOO.median(), M.r_atac.median(), AN.r.median()))
print(" significance counts: FDR(delivered)<0.05 %d | MC P(LOO)<0.05 %d | FDR(ATAC)<0.05 %d | anchors MC P<0.05 %d/8"
      % (int((M.fdr_bh < 0.05).sum()), int((M.mc_p_LOO < 0.05).sum()),
         int((M.fdr_atac < 0.05).sum()), int((AN.mc_p < 0.05).sum())))
