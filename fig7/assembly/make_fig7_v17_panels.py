#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_fig7_v17_panels.py -- the four Figure 7 panels of the FINAL (v17) delivery, one PDF each.

Delivered sheet (assembled outside this file): A | B | C on row 1, D on row 2, inside a
17.00 x 14.34 cm two-column page.  This script draws panels only -- no letters, no composing, no QA
gates -- because those are journal assembly, not content (house spec SS1, SS9).

Typography: the house two-tier scale, 8.0 pt (titles, axis titles, row heads) and 6.0 pt (everything
else), Arial, black on the face except the trait colours.  Panels are drawn at their FINAL width; no
delivered PDF is rescaled (SS9 forbids vector scaling).

Run from the fig7/ directory.  Environment:
  DSCI_WORK_ROOT  analysis tree (default: two levels above the repo) -- needed only for the Arial
                  TTF; every table this script reads is vendored under data/ or tables_v17/
                  (see data/MANIFEST.csv for path -> producing script -> panel)
  FIG7_OUT        where the panel PDFs go (default fig7_out/)
  FIG7_TABLES     vendored small tables for D (default tables_v17/)
  FIG7_DATA       vendored inputs for A/B/C (default data/)

Vendored inputs (data/, mirrored from the analysis tree, byte-identical, all <= 75 rows):
  data/FIG7_v11/tables/fig7A_donors_delta1.csv                upstream/D9_fig7_assembly.py
  data/FIG7_v11/tables/fig7A_between_donor_P.csv              upstream/D9_fig7_assembly.py
  data/FIG7_v1/tables/seven4_confound_ALL_cohorts.csv         upstream/A8c_confound_unified.py
  data/A2_GSE277430_v1/tables/trait_sweep_dual_axis.csv        scripts/routeA/A2_GSE277430_delta2_mirror.py
  data/A3_GSE118257_v1/tables/trait_sweep_delta2_GSE118257_74.csv  scripts/routeA/A3_GSE118257_delta2_full74.py
  data/A4_GSE301908_v1/tables/trait_sweep_delta2_GSE301908_74.csv  scripts/routeA/A4_GSE301908_delta2_full74.py
  data/A7_GSE231585_rna_v1/tables/trait_sweep_delta2_GSE231585_74.csv  upstream/A7_GSE231585_delta2_full74.py
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.stats import mannwhitneyu, spearmanr, gaussian_kde

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # arial_style.py, the font policy
from arial_style import use_arial, use_two_tier, FAMILY, open_spines   # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
OUTD = os.environ.get("FIG7_OUT", os.path.join(REPO, "fig7_out"))
TABD = os.environ.get("FIG7_TABLES", os.path.join(REPO, "tables_v17"))
os.makedirs(OUTD, exist_ok=True)
_ROOT = os.environ.get("DSCI_WORK_ROOT", os.path.abspath(os.path.join(REPO, "..", "..")))
DATA = os.environ.get("FIG7_DATA", os.path.join(REPO, "data"))   # vendored A/B/C inputs (G1, 2026-11-18)

matplotlib.rcParams.update({
    # NOTE (2026-11-18, author: "why do A and C have frames?"): the house style loaded by use_arial()
    # below sets axes.spines.top/right = True, and it is applied AFTER this dict, so it wins. Panel A
    # and the four panel-C tiles are statistical panels and must be open-frame; they now call
    # open_spines() explicitly instead of inheriting the UMAP-era default.
    "axes.linewidth": 0.55, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.pad": 1.6, "xtick.major.size": 2.2, "ytick.major.size": 2.2})
use_arial({"savefig.bbox": "tight", "savefig.pad_inches": 0.04, "savefig.dpi": 600, "figure.dpi": 120})
use_two_tier()                                # registers the 8.0 / 6.0 tier names used below
S8, S6 = 8.0, 6.0                       # the only two sizes on the face
CM = 2.54
RED, BLU, ORG, GRN, MUT, GRID = "#D51F26", "#272E6A", "#F47D2B", "#0B8457", "0.42", "0.82"
WB = dict(fc="white", ec="none", alpha=0.92, pad=0.25)
MS = "PASS_Multiple_sclerosis"
SCZ = "PASS_Schizophrenia_Pardinas2018"
ADt = "PASS_Alzheimers_Jansen2019"
COH = ["A2", "A3", "A4", "A7"]
GSE = {"A2": "277430", "A3": "118257", "A4": "301908", "A7": "231585"}
STUDIES = ["A2 GSE277430", "A3 GSE118257", "A4 GSE301908", "A7 GSE231585", "M GSE250589"]


def savefig(fig, name, tight=True):
    """Write the panel PDF and report the size MEASURED from the file (tight bbox != figsize).
    A/B/C are trimmed to their ink; D is laid out on an exact canvas (its own coordinates are the
    layout), so it is saved with bbox_inches=None at exactly FW x FH."""
    import fitz
    fig.savefig(os.path.join(OUTD, name + ".pdf"),
                bbox_inches=(None if not tight else "tight"),
                pad_inches=(0 if not tight else 0.04))
    plt.close(fig)
    p = fitz.open(os.path.join(OUTD, name + ".pdf"))[0]
    print("%-28s %.2f x %.2f cm" % (name + ".pdf", p.rect.width * CM / 72, p.rect.height * CM / 72))


# ---------------------------------------------------------------- shared geometry helpers
def twcm(f, a, st, fs=6.0):
    """rendered width of a string, in cm (never trust a char-count estimate)"""
    t = a.text(0, -9999, st, fontsize=fs, clip_on=False)
    bb = t.get_window_extent(renderer=f.canvas.get_renderer()); t.remove()
    return (bb.x1 - bb.x0) / f.dpi * CM


def dxcm(a, d_cm):
    """cm along x -> that axes' data units"""
    p0 = a.transData.transform((0.0, 0.0)); p1 = a.transData.transform((1.0, 0.0))
    return d_cm / 2.54 * 72.0 / (p1[0] - p0[0])


def rot_anchor(a, x_data, st, cos_half):
    """rotated tick labels read DOWN-LEFT from their anchor; shift the anchor right by half the
    horizontal projection so the label body lands centred under its own column"""
    return x_data + dxcm(a, twcm(a.figure, a, st) * cos_half)


def mwu_p(x, y):
    try:
        return float(mannwhitneyu(x, y).pvalue)
    except ValueError:
        return float("nan")


def robust_z(v):
    med = np.nanmedian(v); mad = 1.4826 * np.nanmedian(np.abs(v - med))
    return (v - med) / (mad if mad > 1e-12 else 1.0)


# ---------------------------------------------------------------- inputs (all delivered tables)
R = DATA                                          # vendored sweep tables (mirrors the routeA_out layout)
T11 = os.path.join(DATA, "FIG7_v11", "tables")
T7 = os.path.join(DATA, "FIG7_v1", "tables")
for _p in [T11 + "/fig7A_donors_delta1.csv", T11 + "/fig7A_between_donor_P.csv",
           T7 + "/seven4_confound_ALL_cohorts.csv",
           R + "/A2_GSE277430_v1/tables/trait_sweep_dual_axis.csv",
           R + "/A3_GSE118257_v1/tables/trait_sweep_delta2_GSE118257_74.csv",
           R + "/A4_GSE301908_v1/tables/trait_sweep_delta2_GSE301908_74.csv",
           R + "/A7_GSE231585_rna_v1/tables/trait_sweep_delta2_GSE231585_74.csv"]:
    if not os.path.exists(_p):
        raise SystemExit("vendored input missing: %s (see data/MANIFEST.csv)" % _p)
A = pd.read_csv(T11 + "/fig7A_donors_delta1.csv")
AP = pd.read_csv(T11 + "/fig7A_between_donor_P.csv")
CONF = pd.read_csv(T7 + "/seven4_confound_ALL_cohorts.csv")

# Confounding of the raw ranking by regulatory-domain size, recomputed here from the same table Panel B
# is drawn from (columns A2..A7 are the cohort-level d1, lnp = ln n_peak_raw).  These are the four rho
# values quoted in the caption; Fig S6/S7 store the identical numbers in their cross-concordance ledger
# D3_crossconcordance_domainadjust.csv, so the main figure and the supplement cannot disagree.
from scipy.stats import spearmanr as _spearmanr
RHO_LNPEAK = {c: round(_spearmanr(CONF[c].dropna(), CONF.loc[CONF[c].notna(), "lnp"]).statistic, 3) for c in COH}
assert [RHO_LNPEAK[c] for c in COH] == [-0.838, -0.691, -0.778, -0.873], RHO_LNPEAK


# %% PANEL A
# ================================================================ A  donor-level delta1
PVALS = os.environ.get("FIG7A_PVALS", "bare")     # "bare" = 0.091 | "P=" = P=0.091 | "off" = not on face
_mg = A[A.study == "M GSE250589"]
HEADS = [st.split(" ")[1] for st in STUDIES[:4]] + ["Early", "Peak", "Late"]
TPS = ["Early", "Peak", "Late"]
COLG = [A[A.study == st] for st in STUDIES[:4]] + [_mg[_mg.id.str.startswith(tp)] for tp in TPS]
# Canvas defaults are the geometry that reproduces the delivered sheet: the panels are saved with a
# tight ink bbox, so the requested width feeds back into the ink width, and the composer needs the row-1
# ink to leave a 0.20-0.40 cm inter-panel gap at the fixed 17.00 cm page. Do not "tidy" these numbers.
W_A, H_A = float(os.environ.get("W_A", 6.36)), float(os.environ.get("H_A", 6.30))
fig, ax = plt.subplots(figsize=(W_A / CM, H_A / CM))
open_spines(ax)                    # statistical panel: left + bottom only (author 2026-11-18)
fig.subplots_adjust(left=0.175, right=0.985, top=0.845, bottom=0.305)
ax.axvspan(-0.5, 3.5, color="#EAF1FA", lw=0, zorder=0)
ax.axvspan(3.5, 6.5, color="#FDF3E4", lw=0, zorder=0)
ax.axvline(3.5, color="0.62", ls="--", lw=0.7, zorder=1)
ax.grid(True, axis="y", lw=0.32, color=GRID, zorder=1); ax.set_axisbelow(True)
for x0, x1, cc in ((-0.42, 3.42, "#1F3E6E"), (3.58, 6.42, "#8A4B12")):
    ax.plot([x0, x1], [3.14, 3.14], color=cc, lw=0.8, zorder=6)
    for xe in (x0, x1):
        ax.plot([xe, xe], [3.14, 3.03], color=cc, lw=0.8, zorder=6)
ax.text(1.50, 3.22, "Human MS (brain)", fontsize=S6, ha="center", va="bottom", color="#1F3E6E", zorder=6)
ax.text(5.00, 3.22, "Mouse EAE (spinal cord)", fontsize=S6, ha="center", va="bottom", color="#8A4B12", zorder=6)
ax.text(5.00, 3.00, "GSE250589", fontsize=S6, ha="center", va="top", color="#8A4B12", zorder=6)
rng = np.random.default_rng(11)
PP = {r.study: float(r.P) for _, r in AP.iterrows()}
for i, g in enumerate(COLG):
    for cond, sub in g.groupby("cond"):
        x = i + rng.uniform(-0.16, 0.16, len(sub))
        ax.scatter(x, sub.d1, s=13, marker="o", color={"MS": RED, "EAE": ORG, "Ctrl": BLU}.get(cond, MUT),
                   alpha=0.9, edgecolor="white", lw=0.35, zorder=4)
    for cond in ["Ctrl", "MS", "EAE"]:
        s_ = g[g.cond == cond].d1.values
        if len(s_) == 0:
            continue
        xm = i + (0.22 if cond == "Ctrl" else -0.22)
        ax.plot([xm - 0.13, xm + 0.13], [np.median(s_)] * 2,
                color={"MS": RED, "EAE": ORG, "Ctrl": BLU}[cond], lw=1.6, zorder=5)
    pv = PP[STUDIES[i]] if i < 4 else mwu_p(g[g.cond == "EAE"].d1.values, g[g.cond == "Ctrl"].d1.values)
    if PVALS != "off":
        ax.text(i, 2.62, ("P=%.3f" if PVALS == "P=" else "%.3f") % pv, fontsize=S6, va="top",
                ha="center", color="0.15", zorder=6)
ax.axhline(0, color=MUT, ls=":", lw=0.7, zorder=2)
ax.set_xlim(-0.5, 6.5); ax.set_ylim(-0.45, 3.42)
ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
ax.set_xticks(range(len(HEADS))); ax.set_xticklabels([""] * len(HEADS))
# column heads: full accession / stage name, 45 deg (house spec SS14.2 step 2: two-line stacking no
# longer clears the 0.75 cm column pitch at 6 pt -- "301908" and "231585" ended 0.03 cm apart).
for i, hd in enumerate(HEADS):
    ax.text(rot_anchor(ax, i, hd, 0.3536), -0.95, hd, fontsize=S6, ha="right", va="top",
            rotation=45, rotation_mode="anchor")
ax.set_ylabel("mature OL \u2212 OPC\n(median score)", fontsize=S8, linespacing=1.35)
handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markeredgecolor=c,
                  markersize=5, label=c) for c in ["#D51F26", "#F47D2B", "#272E6A"]]
handles[0].set_label("MS"); handles[1].set_label("EAE"); handles[2].set_label("Ctrl")
# legend moved to the top-right of the panel (same convention as panel C): the rotated 45 deg column
# heads need the whole band below the axes.
fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.995, 1.005), ncol=3, fontsize=S6,
           frameon=False, handletextpad=0.2, columnspacing=1.15)
savefig(fig, "fig7v17_A_donors")

# %% PANEL B
# ================================================================ B  74 x 4 rank heatmap (portrait, deconfounded ranks)
B_ = CONF.copy()
B_["mean_adj"] = B_[["rank_adj_%s" % c for c in COH]].mean(axis=1)
order = B_.sort_values("mean_adj", ascending=True).trait.tolist()
M = B_.set_index("trait").loc[order]
assert len(M) == 74
Z = np.array([M["rank_adj_%s" % c].values for c in COH], dtype=float).T
imsi = int(np.where(M.index == MS)[0][0])
W_B, H_B = float(os.environ.get("W_B", 4.00)), float(os.environ.get("H_B", 5.98))
fig, ax = plt.subplots(figsize=(W_B / CM, H_B / CM))
fig.subplots_adjust(left=0.20, right=0.830, top=0.975, bottom=0.255)
im = ax.imshow(Z, cmap="viridis_r", aspect="auto", vmin=1, vmax=74)
ax.set_xlim(-0.5, 3.5); ax.set_ylim(73.5, -0.5)
assert abs(ax.get_ylim()[0] - 73.5) < 1e-6 and len(M) == 74, "74 rows must all be drawn"
ax.set_xticks(range(4))
ax.set_xticklabels(["GSE%s" % GSE[c] for c in COH], fontsize=S6, rotation=90, ha="right")
ax.tick_params(axis="x", pad=1.5)
for _t in ax.get_xticklabels():
    _t.set_fontweight("normal")
ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)
# short forms as in the delivered v16 face and panel C (house spec SS11.3: full names in caption)
for t_, nm_, col_ in [(MS, "MS", RED), (ADt, "AD", GRN), (SCZ, "SCZ", BLU)]:
    ry = list(M.index).index(t_)
    ax.text(-0.018, ry, nm_, transform=ax.get_yaxis_transform(), ha="right", va="center",
            fontsize=S6, color=col_, zorder=6)
    ax.plot([-0.004, 0.028], [ry, ry], transform=ax.get_yaxis_transform(), color=col_, lw=0.9,
            clip_on=False, zorder=6)
cax = fig.add_axes([0.868, 0.295, 0.022, 0.45])
cb = fig.colorbar(im, cax=cax, orientation="vertical", ticks=[1, 20, 40, 60, 74])
cb.outline.set_linewidth(0.4)
cb.ax.tick_params(labelsize=S6, length=1.6, pad=1.5)
fig.text(0.868, 0.79, "rank", fontsize=S6, ha="left", va="bottom")
savefig(fig, "fig7v17_B_adjranks")

# %% PANEL C
# ================================================================ C  4 datasets x 3 axes, 2x2 violins
d2a = pd.read_csv(R + "/A2_GSE277430_v1/tables/trait_sweep_dual_axis.csv")
d3t = pd.read_csv(R + "/A3_GSE118257_v1/tables/trait_sweep_delta2_GSE118257_74.csv")
d4t = pd.read_csv(R + "/A4_GSE301908_v1/tables/trait_sweep_delta2_GSE301908_74.csv")
d7t = pd.read_csv(R + "/A7_GSE231585_rna_v1/tables/trait_sweep_delta2_GSE231585_74.csv")
TILES = [("GSE277430", d2a, "d1_OLmat", "d2_AST_IM", "d2_NEU_GL", "r1", "rA", "rB", None),
         ("GSE118257", d3t, "d1_OLIGO_minus_OPC", "dA_AST_minus_IMM", "dB_NEU_minus_GLIA", "r1", "rA", "rB", None),
         ("GSE301908", d4t, "d1", "dA", "dB", "r1", "rA", "rB", None),
         # `note` = the AXIS INDEX whose pool is thin, not a dataset-level flag: in GSE231585 only
         # the neuron pool is small (717 of 62,911 nuclei; 4-133 per sample; 2 of 15 samples drop out
         # under the >=10-nuclei rule, hence nB=13 vs n1=nA=14).  v13 printed the mark in the tile
         # title, which reads as a caveat about the whole cohort; it now hangs on the axis it qualifies.
         ("GSE231585", d7t, "d1", "dA", "dB", "r1", "rA", "rB", 2)]
# Axis labels of the three violin columns.  These are the **cell-type contrasts**, written with the same
# symbols the manuscript abbreviation table uses and the same strings as Fig S7 E/G and Table S7, so a
# reader can map one axis across the three places without a dictionary:
#   OL−OPC        mature oligodendrocyte vs OPC            (the old face text "OL maturation")
#   AST−(MG+LYM)  astrocyte vs the pooled immune fractions (the old face text "astro − immune")
#   NEU−non-NEU   neuron vs every non-neuronal nucleus     (the old face text "neuron − glia", which
#                 was wrong on both sides: the comparator is all non-neuronal nuclei, not glia only)
# Why the change (2026-11-18, author): panel D of the same figure already prints AST / MG / LYM / NEU,
# while C spelled out astro / immune / neuron / glia - one figure, two naming systems.
AXL = ["OL\u2212OPC", "AST\u2212(MG+LYM)", "NEU\u2212non-NEU"]
AXCOL = ["#4E7A3A", "#A6651E", "#5B5688"]


def violin_clipped(ax, v, x, half, col, ylo, yhi):
    from scipy.stats import gaussian_kde
    g0, g1 = max(v.min(), ylo + 0.06), min(v.max(), yhi - 0.06)
    if g1 <= g0:
        return
    grid = np.linspace(g0, g1, 220)
    den = gaussian_kde(v, bw_method=0.35)(grid)
    den = den / den.max() * half
    ax.fill_betweenx(grid, x - den, x + den, facecolor=col, edgecolor=col, lw=0.5, alpha=0.30, zorder=3)


def robust_z(v):
    med = np.nanmedian(v); mad = 1.4826 * np.nanmedian(np.abs(v - med))
    return (v - med) / (mad if mad > 1e-12 else 1.0)


def separate_rank_labels(ax, fig, maxit=18, step=1.4, cap=11.0):
    """DEAD CODE kept for provenance: the fixed MS-up / SCZ-down stagger could not clear a corner
    collision between '#71' and '#4' without pushing a label past the 0.35 cm anchor limit
    (house spec SS14.8).  Replaced by place_rank_labels(), which searches bounded positions."""
    return


def place_rank_labels(fig, axes, rings=(7.4, 9.9)):
    """Place every '#rank' label on a small ring around its own dot.

    Why not the delivered v13 offsets: v13 tiles were 7.88 cm wide (column pitch ~0.9 cm) so a label
    0.42 cm to the right of its dot never reached the next column.  In the 17 cm two-column figure the
    tiles are 5.82 cm (pitch ~0.65 cm) and a fixed right/left offset collides with the neighbour's dot
    or label.  Here each label tries 8 bearings on a 0.26 cm ring, then the same 8 on a 0.35 cm ring
    (house spec SS14.8 cap), and takes the first position whose MEASURED box touches neither another
    label, nor another dot, nor leaves the figure.  Placement is global over the four tiles, so
    cross-tile collisions are seen too."""
    rend = fig.canvas.get_renderer()
    fig_bb = fig.get_window_extent(renderer=rend)
    items = []
    for ax in axes:
        dots = []
        for coll in ax.collections:
            for x0, y0 in coll.get_offsets():
                dots.append(ax.transData.transform((float(x0), float(y0))))
        for a in ax.texts:
            if a.get_text().startswith("#"):
                items.append((a, ax, dots))
    from matplotlib.transforms import Bbox

    def infl(b, m=0.8):
        return Bbox.from_bounds(b.x0 - m, b.y0 - m, b.width + 2 * m, b.height + 2 * m)
    placed = []
    for a, ax, dots in items:
        a.set_ha("center"); a.set_va("center")
        px, py = ax.transData.transform((a.xy[0], a.xy[1]))
        chosen = None
        for rad in rings:
            cand = [(rad, 0.0), (rad * .707, rad * .707), (0.0, rad), (-rad * .707, rad * .707),
                    (-rad, 0.0), (-rad * .707, -rad * .707), (0.0, -rad), (rad * .707, -rad * .707)]
            for dx, dy in cand:
                a.xyann = (dx, dy)
                bb = a.get_window_extent(renderer=rend)
                if (bb.x0 < fig_bb.x0 or bb.x1 > fig_bb.x1 or bb.y0 < fig_bb.y0 or bb.y1 > fig_bb.y1):
                    continue
                if any(infl(bb).overlaps(q) for q in placed):
                    continue
                if any(fitz_dot(infl(bb), qx, qy) for qx, qy in dots
                       if not (abs(qx - px) < 1e-6 and abs(qy - py) < 1e-6)):
                    continue
                chosen = (dx, dy, bb)
                break
            if chosen:
                break
        if chosen is None:
            a.xyann = (rings[0], 0.0)
            chosen = (rings[0], 0.0, a.get_window_extent(renderer=rend))
        placed.append(infl(chosen[2], 0.0))


def fitz_dot(bb, x, y, r=3.4):
    """True if the (r pt) dot at (x, y) touches box bb."""
    return bb.x0 - r <= x <= bb.x1 + r and bb.y0 - r <= y <= bb.y1 + r


def place_axis_marks(fig, marks,
                     try_offsets=((0.0, 1.7), (0.0, 3.0), (0.22, 1.7), (-0.22, 1.7), (0.0, 0.6))):
    """Hang the thin-pool marker DIRECTLY ABOVE the violin it qualifies.

    Author 2026-11-17, two rejected positions: (a) v13/v16 appended it to the tile title -> reads as a
    caveat about the whole cohort; (b) appending it to the rotated axis label put it in the empty corner
    between the two right-hand tiles, which reads as "both of these cohorts are suspect".  Straight over
    the column, just outside the top spine, is the only spot that is unambiguously inside THIS tile and
    over THIS axis.  Placement is measured, not assumed: the box must clear every label in the panel,
    the tile title, this tile's dots and the figure edge; the first free candidate wins."""
    from matplotlib.transforms import Bbox
    rend = fig.canvas.get_renderer()
    fig_bb = fig.get_window_extent(renderer=rend)
    boxes = []
    for ax in fig.axes:
        for t in list(ax.texts) + [ax.title, ax.xaxis.label, ax.yaxis.label]:
            if t.get_text().strip():
                boxes.append(t.get_window_extent(renderer=rend))
    for lg in fig.legends:
        for t in lg.get_texts():
            boxes.append(t.get_window_extent(renderer=rend))
    out = []
    for ax, k in marks:
        dots = [ax.transData.transform((float(x0), float(y0)))
                for coll in ax.collections for x0, y0 in coll.get_offsets()]
        top = ax.get_ylim()[1]
        chosen = None
        for dx, dy in try_offsets:
            t = ax.annotate("*", (k + dx, top), xytext=(0.0, dy), textcoords="offset points",
                            fontsize=S6, ha="center", va="bottom", color="black",
                            clip_on=False, zorder=8)
            bb = t.get_window_extent(renderer=rend)
            if bb.x0 < fig_bb.x0 or bb.x1 > fig_bb.x1 or bb.y1 > fig_bb.y1:
                t.remove(); continue
            if any(bb.overlaps(b) for b in boxes):
                t.remove(); continue
            if any(fitz_dot(bb, x, y) for x, y in dots):
                t.remove(); continue
            boxes.append(bb); chosen = (t, k + dx, dy); break
        if chosen is None:
            raise SystemExit("axis mark has no free position over column %d" % k)
        out.append(chosen)
    return out


W_C, H_C = float(os.environ.get("W_C", 5.90)), float(os.environ.get("H_C", 5.80))
fig, axes = plt.subplots(2, 2, figsize=(W_C / CM, H_C / CM))
fig.subplots_adjust(left=0.145, right=0.985, top=0.825, bottom=0.245, wspace=0.17, hspace=0.42)
YLO, YHI = -4.7, 5.5
AXIS_MARKS = []
# Where the thin-pool marker hangs.  "violin" = current (directly above the axis it qualifies).
# "label" (v17b) read as "both right-hand tiles are suspect"; "title" (v13/v16/v17a) read as a
# cohort-level caveat.  Both rejected variants stay reproducible so the two archived PDFs can be
# regenerated byte-for-byte from this file.
STAR_POS = os.environ.get("STAR_POS", "violin")
for gi, (gse, T, c1, cA, cB, r1, rA, rB, note) in enumerate(TILES):
    T = T[T.trait != "TRAIT"].dropna(subset=[c1, cA, cB]).reset_index(drop=True)
    ax = axes[gi // 2, gi % 2]
    open_spines(ax)                # statistical tiles: left + bottom only (author 2026-11-18)
    for k, (cc, rr) in enumerate([(c1, r1), (cA, rA), (cB, rB)]):
        v = robust_z(T[cc].values); vv = v[np.isfinite(v)]
        q1, q2, q3 = np.nanpercentile(v, [25, 50, 75])
        violin_clipped(ax, vv, k, 0.28, AXCOL[k], YLO, YHI)
        ax.add_patch(plt.Rectangle((k - 0.055, q1), 0.11, q3 - q1, facecolor="white",
                                   edgecolor=AXCOL[k], lw=0.5, alpha=0.92, zorder=4))
        ax.plot([k - 0.055, k + 0.055], [q2, q2], color=AXCOL[k], lw=1.0, zorder=5)
        # MS / SCZ markers.  When their robust z are close, the two "#rank" labels flank the same dot
        # and run into each other (house spec SS17.1: fitz merges flush words, the collide gate stays
        # silent).  Rule: keep the dots centred, keep MS right / SCZ left, and in a near-tie push the
        # MS label up and the SCZ label down by 0.15 cm so the two boxes cannot touch.
        pts = []
        for tr, col in [(MS, RED), (SCZ, BLU)]:
            row = T[T.trait == tr]
            if not len(row):
                continue
            rk = int(row[rr].iloc[0])
            zv = float(robust_z(T[cc].values)[T.trait.values == tr][0])
            pts.append((col, rk, zv))
        for col, rk, zv in pts:
            side = 1 if col == RED else -1
            dy = 6.2 if side > 0 else -6.2          # initial bearing; place_rank_labels re-searches
            ax.scatter([k], [zv], s=18, facecolor=col, edgecolor="white", lw=0.5, zorder=7)
            ax.annotate("#%d" % rk, (k, zv), xytext=(side * 72 / 2.54 * 0.30, dy),
                        textcoords="offset points", fontsize=S6, color=col,
                        ha="left" if side > 0 else "right", va="center", zorder=7, bbox=WB)
    _TILES_DONE = True
    ax.axhline(0, color=MUT, lw=0.5, ls=":", zorder=1)
    ax.set_xlim(-0.55, 2.55); ax.set_ylim(YLO, YHI)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels([])
    if note is not None and STAR_POS == "violin":
        AXIS_MARKS.append((ax, note))
    if gi // 2 == 1:
        for k, lab in enumerate(AXL):
            txt = lab + " *" if (note == k and STAR_POS == "label") else lab
            ax.text(rot_anchor(ax, k, txt, 0.25), YLO - 0.22, txt, ha="right", va="top",
                fontsize=S6, rotation=60, rotation_mode="anchor", color="black")
    ax.set_title((gse + " *") if (note is not None and STAR_POS == "title") else gse,
                 fontsize=S8, fontweight="normal", loc="left", pad=3.0)
    if gi % 2 == 0:
        ax.set_ylabel("robust z (\u0394)", fontsize=S8)
    else:
        ax.tick_params(axis="y", left=False, labelleft=False)
    ax.grid(True, axis="y", lw=0.3, color=GRID); ax.set_axisbelow(True)
place_rank_labels(fig, list(axes.ravel()))
_placed = place_axis_marks(fig, AXIS_MARKS)
print("axis marks: %d placed (STAR_POS=%s)" % (len(_placed), STAR_POS))
fig.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=RED, markeredgecolor="white",
                           markersize=3.6, label="MS"),
                    Line2D([0], [0], marker="o", color="none", markerfacecolor=BLU, markeredgecolor="white",
                           markersize=3.6, label="SCZ")],
           loc="upper right", bbox_to_anchor=(1.00, 1.005), ncol=2, fontsize=S6, frameon=False,
           handlelength=0.0, columnspacing=0.9, handletextpad=0.42)
savefig(fig, "fig7v17_C_threeaxes")

# %% PANEL D
T = TABD                 # vendored input tables (tables_v17/)

C = pd.concat([pd.read_csv(T + "/S13_two_contrasts.csv"),
               pd.read_csv(T + "/S23_new_targets_stats.csv")], ignore_index=True).drop_duplicates(["ds", "gene"])
RS = pd.read_csv(T + "/S26_rowset_C12.csv").rename(columns={"class": "klass"})
F2 = pd.read_csv(T + "/S29_flags_v2.csv").set_index("gene")
assert not (F2[["label_v2", "MOA_v2", "MS_v2"]] == 2).any().any(), "a half-filled flag is still pending"
LED = pd.read_csv(T + "/S26_target_ledger.csv").set_index("gene")
assert len(LED) == 83 and (LED.tier_audited != "X").sum() == 80, \
    "the caption says 18 rows survived audit out of 80 scored: 83 ledger rows minus 3 non-drug atoms"
n_scored = int((LED.tier_audited != "X").sum())
CELL = {(r.ds, r.gene): r for r in C.itertuples()}
DS4 = ["GSE118257", "GSE301908", "GSE231585", "GSE277430"]
DS3 = DS4[:3]
DET = {"A4": pd.read_csv(T + "/S15_det_GSE301908.csv", index_col=0).groupby(level=0).max()}
ND = {"A4": pd.read_csv(T + "/S15_n_GSE301908.csv", index_col=0)["n"]}
POP = list(DET["A4"].columns)
RATE = {d: (DET[d].div(ND[d].reindex(POP).values.astype(float), axis=1) * 100) for d in DET}
NONP = [p for p in POP if p not in ["OL", "OPC", "iOL", "Unclear", "QC"]]

ROWS_ = [(r.gene, "F" if r.floor else str(r.klass)[:1], r.tier_audited) for r in RS.itertuples()]
ORD = {"1": 0, "2": 1, "3": 2, "6": 3, "4": 4, "5": 5, "F": 6}
ROWS = sorted(ROWS_, key=lambda t: (ORD[t[1]], t[0]))
N = len(ROWS)

FW, FH = float(os.environ.get("FIG7D_FW", 16.2)), 7.90
FS, FSH, FST = 6.0, 6.0, 8.0   # house spec two tiers: 8.0 letter + 6.0 face
CW, GAPC, GAPB = .62, .04, .26
RWC, RGAP = .46, .025
X_TIER, LABR = .52, 2.30
X_C2 = 2.55
W_C2 = 4 * CW + 3 * GAPC
X_C1 = X_C2 + W_C2 + GAPB
W_C1 = 3 * CW + 2 * GAPC
X_GR = X_C1 + W_C1 + GAPB
EW = [.48, .48, .48]
FL_W = sum(EW) + .03 * 2
X_RATE = X_GR + FL_W + GAPB
RATE_W = len(POP) * RWC + (len(POP) - 1) * RGAP
X_DIV = X_RATE + RATE_W + .28
X_LEG = X_DIV + .18
LEG_W = FW - X_LEG - .25
TXT_W = LEG_W - .36
XS = [X_C2 + i * (CW + GAPC) for i in range(4)] + [X_C1 + i * (CW + GAPC) for i in range(3)]
YTOP, YBOT = 6.78, 1.72
RHP = min(.32, (YTOP - YBOT) / (N - .5))
assert X_RATE + RATE_W < X_DIV < X_LEG < FW

RCM = LinearSegmentedColormap.from_list("rp", ["#f4f2f8", "#d8cbe4", "#a585c4", "#6a3d9a", "#3a1a63"])
RNRM = Normalize(-2, 2)
RLOG = lambda v: np.clip(np.log10(max(v, .01)), -2, 2)
ORCM = LinearSegmentedColormap.from_list("or", ["#2166ac", "#67a9cf", "#f7f4e9", "#ef8a62", "#b2182b"])
NRM = Normalize(-3, 3)
pale = lambda col, f=0.26: tuple(1 - (1 - matplotlib.colors.to_rgb(col)[k]) * f for k in range(3))

fig = plt.figure(figsize=(FW / 2.54, FH / 2.54))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, FW); ax.set_ylim(0, FH); ax.set_aspect("equal"); ax.axis("off")
_rend = fig.canvas.get_renderer()


def tw(s, fs=FS, rot=0):
    t = ax.text(0, -99, s, fontsize=fs, rotation=rot)
    bb = t.get_window_extent(renderer=_rend)
    t.remove()
    return (bb.x1 - bb.x0) / fig.dpi * 2.54


# ---------------------------------------------------------------- titles
TIER_Y = 7.06
if os.environ.get("FIG7D_LETTER"):        # panel letters belong to journal assembly (house spec SS1); opt in with FIG7D_LETTER=1
    ax.text(.30, 7.55, "D", fontsize=FST, va="center", color=".1", weight="normal")
ax.text(X_C2 + W_C2 / 2, TIER_Y, "mature OL vs OPC", fontsize=FSH, ha="center", va="center", color=".1")
ax.text(X_C1 + W_C1 / 2, TIER_Y, "lineage vs non-lineage", fontsize=FSH, ha="center", va="center", color=".1")
ax.text(X_GR + FL_W / 2, TIER_Y, "drug\u2013target", fontsize=FSH, ha="center", va="center", color=".1")
ax.text(X_RATE + RATE_W / 2, TIER_Y, "% nuclei (log)", fontsize=FSH, ha="center", va="center", color=".1")
ax.text(X_TIER, TIER_Y, "tier", fontsize=FS, ha="left", va="center", color=".1")
ax.text(LABR, TIER_Y, "gene", fontsize=FS, ha="right", va="center", color=".1")
assert tw("drug\u2013target", FSH) <= FL_W + .5 and tw("% nuclei (log)", FSH) <= RATE_W
assert TIER_Y - .14 > YTOP - .02 + .06, "header line touches the first row"
assert tw("tier") + tw("gene") < LABR - X_TIER and .30 + tw("D", FST) < X_TIER

# ---------------------------------------------------------------- rows
for i, (g, cls, tier) in enumerate(ROWS):
    y = YTOP - RHP * (i + .5)
    ax.text(X_TIER, y, tier, fontsize=FS, ha="left", va="center",
            color=".1" if tier == "T1" else ".45", weight="bold" if tier == "T1" else "normal")
    ax.text(LABR, y, g, fontsize=FS, ha="right", va="center", color=".1",
            weight="bold" if g == "S1PR5" else "normal")
    for k, (kind, ds) in enumerate([("C2", d) for d in DS4] + [("C1", d) for d in DS3]):
        r = CELL.get((ds, g))
        x = XS[k]
        lor = None if r is None else (r.lor_OLvsOPC if kind == "C2" else r.lor_lineage)
        if r is None or (kind == "C1" and ds == "GSE277430") or (r is not None and not np.isfinite(lor)):
            ax.add_patch(plt.Rectangle((x, y - RHP / 2 + .02), CW, RHP - .04, facecolor="#ececec",
                                       edgecolor="0.72", lw=.3, hatch="///"))
            continue
        if cls == "F":
            ax.add_patch(plt.Rectangle((x, y - RHP / 2 + .02), CW, RHP - .04, facecolor="#dcdcdc",
                                       edgecolor="0.55", lw=.35))
            continue
        q = (r.q_OL_hi if r.lor_OLvsOPC > 0 else r.q_OL_lo) if kind == "C2" else \
            (r.q_lin_hi if r.lor_lineage > 0 else r.q_lin_lo)
        base = ORCM(NRM(np.clip(lor, -3, 3)))
        sig = np.isfinite(q) and q < 0.05
        ax.add_patch(plt.Rectangle((x, y - RHP / 2 + .02), CW, RHP - .04,
                                   facecolor=base if sig else pale(base), edgecolor="0.55" if sig else "0.8", lw=.35))
    for k, fl in enumerate(["label_v2", "MOA_v2", "MS_v2"]):
        st = int(F2.loc[g, fl])
        xx = X_GR + sum(EW[:k]) + .03 * k
        ax.add_patch(plt.Rectangle((xx + (EW[k] - .30) / 2, y - .055), .30, .11,
                                   facecolor="#08519c" if st else "white", edgecolor="0.55", lw=.35))
    if g in RATE["A4"].index:
        pct = RATE["A4"].loc[g].values.astype(float)
        jmax = int(np.nanargmax(pct))
        for j in range(len(POP)):
            ax.add_patch(plt.Rectangle((X_RATE + j * (RWC + RGAP), y - RHP / 2 + .02), RWC, RHP - .04,
                                       facecolor=RCM(RNRM(RLOG(pct[j]))),
                                       edgecolor=".1" if j == jmax else "0.8", lw=.85 if j == jmax else .25,
                                       zorder=6 if j == jmax else 2))
        ol, nn = pct[POP.index("OL")], np.nanmax([pct[POP.index(q2)] for q2 in NONP])
        if cls != "F" and np.isfinite(ol) and np.isfinite(nn) and ol > nn + 1.0 and ol > 1.5 * nn:
            ax.add_patch(plt.Rectangle((X_RATE - .10, y - RHP / 2 + .015), RATE_W + .20, RHP - .03,
                                       facecolor="none", edgecolor="#b2182b", lw=.9, zorder=7))
for i in range(1, N):
    if ROWS[i][1] != ROWS[i - 1][1]:
        ax.plot([.20, X_RATE + RATE_W + .05], [YTOP - RHP * i] * 2, color="0.75", lw=.4, zorder=6)

# ---------------------------------------------------------------- column labels below, rotated 45 deg
ylab = YBOT - RHP / 2 - .10
for k, ds in enumerate(DS4 + DS3):
    ax.text(XS[k] + CW / 2 + tw(ds) * .3536, ylab, ds, fontsize=FS, ha="right", va="top",
            rotation=45, rotation_mode="anchor", color=".15")
for k, w in enumerate(EW):
    _fl = ["label", "MOA", "MS"][k]
    ax.text(X_GR + sum(EW[:k]) + .03 * k + w / 2 + tw(_fl) * .3536, ylab, _fl, fontsize=FS,
            ha="right", va="top", rotation=45, rotation_mode="anchor", color=".15")
for j, pp in enumerate(POP):
    ax.text(X_RATE + j * (RWC + RGAP) + RWC / 2 + tw(pp) * .3536, ylab, pp, fontsize=FS,
            ha="right", va="top", rotation=45, rotation_mode="anchor", color=".15")
ax.text(X_RATE + RATE_W / 2, .78, "GSE301908", fontsize=FS, ha="center", va="center", color=".15")
assert ylab - tw("GSE118257") * .71 > .25, "rotated cohort labels run off the page"
assert .78 + .10 <= ylab - max(tw(p) for p in POP) * .71 - .04, "accession touches the population labels"
assert .78 - .10 > .20, "accession runs off the page"

# ---------------------------------------------------------------- legend column (no rule, no grey notes)


def leg(yy, s, fs=FS, x=None, col=".15"):
    assert tw(s, fs) <= TXT_W + .03, "legend text overruns: " + s
    ax.text(X_LEG + (0 if x is None else x), yy, s, fontsize=fs, va="center", color=col)


ly = TIER_Y
leg(ly, "symbols", FS, col=".1")
SYM = [("max", ["row max"]), ("frame", ["OL > all", "non-lineage"]), ("hatch", ["not in matrix"]),
       ("floor", ["below floor"]), ("pale", ["pale = n.s.", "(q \u2265 0.05)"])]
for k, (kind, lines) in enumerate(SYM):
    yy = ly - .34 - .32 * k
    if kind == "max":
        ax.add_patch(plt.Rectangle((X_LEG, yy - .075), .26, .15, facecolor="#a585c4", edgecolor=".1", lw=.85))
    elif kind == "frame":
        ax.add_patch(plt.Rectangle((X_LEG, yy - .075), .26, .15, facecolor="#d8cbe4", edgecolor="#b2182b", lw=.9))
    elif kind == "hatch":
        ax.add_patch(plt.Rectangle((X_LEG, yy - .075), .26, .15, facecolor="#ececec", edgecolor="0.72", lw=.35, hatch="///"))
    elif kind == "floor":
        ax.add_patch(plt.Rectangle((X_LEG, yy - .075), .26, .15, facecolor="#dcdcdc", edgecolor="0.55", lw=.35))
    else:
        ax.add_patch(plt.Rectangle((X_LEG, yy - .075), .26, .15, facecolor=pale("#ef8a62"), edgecolor="0.8", lw=.35))
    if len(lines) == 1:
        leg(yy, lines[0], x=.36)
    else:
        leg(yy + .085, lines[0], x=.36); leg(yy - .095, lines[1], x=.36)

y2 = ly - .34 - .32 * (len(SYM) - 1) - .62
leg(y2, "log2 OR", FS, col=".1")
bwa = LEG_W - .06
assert X_LEG + .02 + bwa <= FW - .22, "colour bar overruns the page"
cax = fig.add_axes([(X_LEG + .02) / FW, (y2 - .30) / FH, bwa / FW, .15 / FH]); cax.set_axis_off()
cax.imshow(np.linspace(-3, 3, 200).reshape(1, -1), aspect="auto", cmap=ORCM, vmin=-3, vmax=3)
def tickrow(spec, xc_fn, ytxt):
    """place colour-bar tick labels and assert >=0.12 cm of measured air between neighbours"""
    es = []
    for v, lab, ha in spec:
        x, w = xc_fn(v), tw(lab)
        ax.text(x, ytxt, lab, fontsize=FS, ha=ha, va="center", color=".15")
        es.append((x - w / 2, x + w / 2) if ha == "center"
                  else ((x, x + w) if ha == "left" else (x - w, x)))
    for a1, a2 in zip(es, es[1:]):
        assert a2[0] - a1[1] >= .12, "colour-bar ticks touch at y=%.2f (%.2f cm)" % (ytxt, a2[0] - a1[1])


tickrow([(-3, "\u22123", "left"), (0, "0", "center"), (3, "+3", "right")],
        lambda v: X_LEG + .02 + (v + 3) / 6 * bwa, y2 - .48)
y3 = y2 - 1.32
leg(y3, "% nuclei (log)", FS, col=".1")
cax2 = fig.add_axes([(X_LEG + .02) / FW, (y3 - .30) / FH, bwa / FW, .15 / FH]); cax2.set_axis_off()
cax2.imshow(np.linspace(-2, 2, 200).reshape(1, -1), aspect="auto", cmap=RCM, vmin=-2, vmax=2)
_tk = [(-2, "0.01", "left"), (0, "1", "center"), (2, "100", "right")]
tickrow(_tk, lambda v: X_LEG + .02 + (v + 2) / 4 * bwa, y3 - .48)
assert y3 - .48 > YBOT + .10, "legend column collides with the rotated labels"

savefig(fig, "fig7v17_D_targets", tight=False)
