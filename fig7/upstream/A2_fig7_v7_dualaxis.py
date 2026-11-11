# -*- coding: utf-8 -*-
"""Fig 7c (or S7): the two mirror cell-identity axes.
x = d1 (the OL maturation axis, mature OLIGO - precursor OPC)
y = d2A (AST - IMM) / d2B (NEU - GLIA)
Same cells, same 74 traits, same algorithm - only the cell-identity axis changes, so the
"""
import os
import sys

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
sys.path.insert(0, WORK + "/scripts")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import arial_style
plt = arial_style.use_arial()
from scipy.stats import spearmanr

BASE = WORK + "/figures/routeA_out"
OUT = os.path.join(BASE, "FIG7_v1", "pdf")
m = pd.read_csv(os.path.join(BASE, "A2_GSE277430_v1/tables/trait_sweep_dual_axis.csv"))

LAB = {"PASS_Multiple_sclerosis": ("MS", "#c1121f"),
       "UKB_460K.blood_PLATELET_COUNT": ("platelets", "#e07b39"),
       "PASS_Schizophrenia_Pardinas2018": ("schizophrenia", "#1d3557"),
       "PASS_Intelligence_SavageJansen2018": ("intelligence", "#1d3557"),
       "PASS_SWB": ("well-being", "#1d3557"),
       "PASS_MDD_Howard2019": ("depression", "#1d3557"),
       "PASS_BIP_Stahl2019": ("bipolar", "#1d3557"),
       "PASS_SleepDuration_Dashti2019": ("sleep", "#1d3557"),
       "UKB_460K.repro_MENARCHE_AGE": ("menarche", "#8a8a8a"),
       "PASS_Alzheimers_Jansen2019": ("Alzheimer", "#6f6f6f"),
       "PASS_Rheumatoid_Arthritis": ("RA", "#6f6f6f")}


def decollide(rows, gap):
    ys = np.array([r[2] for r in rows], float)
    order = np.argsort(ys)
    for _ in range(400):
        moved = 0.0
        for a, b in zip(order[:-1], order[1:]):
            g = ys[b] - ys[a]
            if g < gap:
                p = (gap - g) / 2.0 + 1e-4
                ys[a] -= p; ys[b] += p; moved = max(moved, p)
        if moved < 1e-4:
            break
    return ys


def panel(ax, ycol, rcol, ylab, tag):
    x, y = m.d1_OLmat.values, m[ycol].values
    rho, prho = spearmanr(x, y)
    lo_x, hi_x = x.min() - 0.5, x.max() + 0.9
    lo_y, hi_y = y.min() - 0.4, y.max() + 0.5
    ax.axvline(0, color="0.8", lw=0.5)
    ax.axhline(0, color="0.8", lw=0.5)
    # MS guide lines: MS sits high on d1 and low on the d2 axes
    ms = m[m.trait == "PASS_Multiple_sclerosis"].iloc[0]
    scz = m[m.trait == "PASS_Schizophrenia_Pardinas2018"].iloc[0]
    ax.plot([ms.d1_OLmat, ms.d1_OLmat], [lo_y, ms[ycol]], color="#c1121f", lw=0.6, ls=":")
    ax.plot([lo_x, ms.d1_OLmat], [ms[ycol], ms[ycol]], color="#c1121f", lw=0.6, ls=":")
    ax.plot([scz.d1_OLmat, scz.d1_OLmat], [lo_y, scz[ycol]], color="#1d3557", lw=0.6, ls=":")
    ax.plot([lo_x, scz.d1_OLmat], [scz[ycol], scz[ycol]], color="#1d3557", lw=0.6, ls=":")
    for _, r in m.iterrows():
        is_ms = r.trait == "PASS_Multiple_sclerosis"
        is_scz = r.trait == "PASS_Schizophrenia_Pardinas2018"
        ax.scatter(r.d1_OLmat, r[ycol], s=58 if (is_ms or is_scz) else 12,
                   color=LAB.get(r.trait, (None, "#bcbcbc"))[1],
                   zorder=6 if is_scz else (7 if is_ms else 3),
                   edgecolors="#333333" if (is_ms or is_scz) else "none", linewidths=0.6)
    mid = (lo_x + hi_x) / 2.0
    rows = []
    for _, r in m.iterrows():
        if r.trait in LAB:
            right = r.d1_OLmat < mid
            rows.append([r.d1_OLmat + (0.25 if right else -0.25), LAB[r.trait][0], r[ycol],
                         LAB[r.trait][1], int(r["r1"]), int(r[rcol]), right])
    ys = decollide(rows, (hi_y - lo_y) / 14.5)
    for (px, lab, y0, c, r1, r2, right), yl in zip(rows, ys):
        if abs(yl - y0) > (hi_y - lo_y) / 90.0:
            ax.plot([px + (-0.22 if right else 0.22), px + (-0.04 if right else 0.04)],
                    [y0, yl], color=c, lw=0.35, alpha=0.55)
        ax.text(px, yl, "%s (%d / %d)" % (lab, r1, r2), fontsize=4.7, color=c, va="center",
                ha="left" if right else "right",
                fontweight="bold" if lab in ("MS", "schizophrenia") else "normal")
    ax.set_xlim(lo_x, hi_x); ax.set_ylim(lo_y, hi_y)
    ax.set_xlabel("\u03941 OL-maturation axis\n(mature OLIGO \u2212 precursor OPC)", fontsize=5.3)
    ax.set_ylabel(ylab, fontsize=5.3)
    ax.text(0.025, 0.045, "Spearman \u03c1 = %.2f\nP = %.1e" % (rho, prho), transform=ax.transAxes,
            fontsize=5.2, va="bottom")
    ax.text(0.975, 0.035, tag, transform=ax.transAxes, fontsize=5.0, va="bottom", ha="right",
            color="#c1121f")
    arial_style.open_spines(ax)


fig = arial_style.cm2fig(16.6, 6.6)
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)
panel(ax1, "d2_NEU_GL", "rB", "\u03942B neuron axis\n(median NEU \u2212 median GLIA)",
      "rank shown as (\u03941 / \u03942B)")
panel(ax2, "d2_AST_IM", "rA", "\u03942A astrocyte\u2013immune axis\n(median AST \u2212 median IMM)",
      "rank shown as (\u03941 / \u03942A)")
fig.suptitle("Fig 7c  Same cells, same 74 traits, different cell-identity axis: "
             "MS peaks only on the OL-maturation axis,\nneuropsychiatric traits peak on the neuron / "
             "astrocyte-immune axes (median rank of 11 neuro-cognitive traits: 61 on \u03941 vs 11 on \u03942B)",
             fontsize=6.3)
fig.tight_layout(rect=[0, 0, 1, 0.935])
fig.savefig(os.path.join(OUT, "FIG7_v7_dualaxis.pdf"))
fig.savefig(os.path.join(OUT, "FIG7_v7_dualaxis.png"), dpi=300)
print("saved FIG7_v7_dualaxis")
