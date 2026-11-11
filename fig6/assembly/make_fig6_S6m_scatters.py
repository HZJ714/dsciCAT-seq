# -*- coding: utf-8 -*-
"""Fig S6m (supplement): per-trait scatters across library chemistries = the v17 main-figure D top-row triplet, moved here wholesale.

  - syntax identical to C: points = the 74 traits, x = chromatin-channel Δ1, y = expression-channel Δ1
    Δ1 = median(scDRS norm_score | Oligo) - median(| OPC); labels = the Fig 2 33-class cell_type_updated
  - three cells = three independent experiments, three library chemistries; all three share one scale (else the visual comparison is void)
  - the three-line provenance small print moves here from the main-figure D (the main figure keeps the forest plot only)
  - every figure number reads tables/*.csv; the script never recomputes
"""
import os, sys

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # arial_style.py, check_text_overlap.py

import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from arial_style import use_arial

ARIAL = os.path.join(os.environ.get('DSCI_FONTS', WORK + "/fonts"), 'ARIAL.TTF')
mpl.font_manager.fontManager.addfont(ARIAL)
mpl.rcParams["font.family"] = mpl.font_manager.FontProperties(fname=ARIAL).get_name()
plt.rcParams["pdf.fonttype"] = 42
use_arial({"xtick.labelsize": 5.8, "ytick.labelsize": 5.8})
mpl.rcParams.update({"font.size": 6.0, "axes.linewidth": 0.55, "xtick.major.width": 0.5,
                    "ytick.major.width": 0.5, "xtick.major.size": 2.2, "ytick.major.size": 2.2,
                    "axes.spines.top": False, "axes.spines.right": False})
CM = 2.54
RED, BLU, MUT, GRID = "#D51F26", "#272E6A", "0.42", "0.85"
WB = dict(fc="white", ec="none", alpha=0.9, pad=0.25)
TABD = os.environ.get("TABD", os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "fig6_out"))
FIG = os.environ.get("OUTFIG", os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "fig6_out"))   # the supplementary S6m panel is drawn from the same two tables
os.makedirs(FIG, exist_ok=True)
DL = pd.read_csv(TABD + "/fig6D_crosschem_ladder.csv")
DT = pd.read_csv(TABD + "/fig6D_delta1_74traits_per_tech.csv")
ORDER = ["dsciCAT-seq (own)", "10x Multiome", "HT-scCAT-seq"]
DL["tech"] = pd.Categorical(DL.tech, ORDER, ordered=True); DL = DL.sort_values("tech")
DT["tech"] = pd.Categorical(DT.tech, ORDER, ordered=True)
NT = int(DL.n_traits.iloc[0])

fig = plt.figure(figsize=(14.6 / CM, 6.4 / CM))
gs = fig.add_gridspec(1, 3, wspace=0.30, left=0.055, right=0.985, top=0.855, bottom=0.20)
sub = [DT[DT.tech == t] for t in ORDER]
lo = min(min(s.d1_atac.min(), s.d1_rna.min()) for s in sub)
hi = max(max(s.d1_atac.max(), s.d1_rna.max()) for s in sub)
pad = 0.12 * (hi - lo); LO, HI = lo - pad, hi + pad

for k, (t, T) in enumerate(zip(ORDER, sub)):
    ax = fig.add_subplot(gs[0, k])
    r = DL[DL.tech == t].iloc[0]
    ax.grid(True, lw=0.32, color=GRID, zorder=0); ax.set_axisbelow(True)
    ax.axhline(0, color=MUT, lw=0.5, alpha=0.6, zorder=1); ax.axvline(0, color=MUT, lw=0.5, alpha=0.6, zorder=1)
    ax.plot([LO, HI], [LO, HI], ls="--", lw=0.6, color=MUT, zorder=1)
    ax.scatter(T.d1_atac, T.d1_rna, s=11, facecolor="0.76", edgecolor="none", zorder=2)
    rd = []
    for key, col, nm in [("Multiple_sclerosis", RED, "MS"), ("Schizophrenia", BLU, "SCZ"), ("HEIGHTz", "0.25", "HEIGHT")]:
        q = T[T.trait.str.contains(key, regex=False)].iloc[0]
        ax.scatter([q.d1_atac], [q.d1_rna], s=30, facecolor=col, edgecolor="white", lw=0.4, zorder=4)
        rd.append("%s #%d\u2192#%d" % (nm, int(q.rank_rna_adj), int(q.rank_atac_adj)))
    ax.text(0.03, 0.975, "\u03c1 %.3f \u2192 %.3f\nperm P %.4f\n95%% CI %.2f\u2013%.2f\n%s"
            % (r.rho_raw, r.rho_adj, r.perm_P_adj, r.ci_lo_adj, r.ci_hi_adj, "   ".join(rd)),
            transform=ax.transAxes, va="top", ha="left", fontsize=5.5, linespacing=1.45, bbox=WB, zorder=6)
    ax.set_xlim(LO, HI); ax.set_ylim(LO, HI)
    ax.set_xlabel("chromatin \u03941 (Oligo \u2212 OPC)", fontsize=6.0)
    ax.set_ylabel("expression \u03941 (Oligo \u2212 OPC)" if k == 0 else "", fontsize=6.0)
    ax.set_title("%s\nn=%s Oligo / %s OPC" % (t, format(int(r.n_oligo), ","), format(int(r.n_opc), ",")),
                 fontsize=6.2, pad=3, linespacing=1.3)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontsize(5.6)

fig.text(0.055, 0.955, "Fig S6m  Cross-library-chemistry replication, per-trait scatter \u00b7 same \u03941 criterion as Fig 6C \u00b7 points = %d traits per panel" % NT,
         fontsize=6.8, fontweight="bold", va="center")
fig.savefig(FIG + "/fig6S6m_crosschem_scatters.pdf", bbox_inches="tight", pad_inches=0.05)
if os.environ.get("PNGOUT"):
    fig.savefig(os.environ["PNGOUT"] + "/fig6S6m_crosschem_scatters.png", dpi=600, bbox_inches="tight", pad_inches=0.05)
print("fig6S6m_crosschem_scatters ok")
