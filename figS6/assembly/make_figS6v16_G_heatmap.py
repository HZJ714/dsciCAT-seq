# -*- coding: utf-8 -*-
"""Fig 6 · D-panel heatmap version (v21): one RNA and one ATAC heatmap, side by side.
 rows = the 74 traits (row order = ascending mean of the three chemistries' rank_rna_adj; only the expression channel orders, chromatin does not)
 each heatmap has 6 columns = 3 chemistries x [Δ1 column | rank column], with narrow gaps between chemistries
   column layout (8 slots): Δ0 rk0 | gap | Δ1 rk1 | gap | Δ2 rk2  - data column indices POS = (0,1)(3,4)(6,7)
 Δ1 column = within-column z score (coolwarm, clipped symmetrically at ±LIM) - Δ1 absolutes are never compared across chemistries, hence the within-column standardization
 rank column = the rank itself (single-hue sequential, rank 1 = darkest); ranks are the only thing allowed to be compared across chemistries
Only the three traits of figure C are labelled (MS / SCZ / HEIGHT); the other 71 rows carry no text.
This version follows the author's notes:
 - the grey "ladder" line across the top is deleted; content moves up (ρ/pass verdicts live only in the caption and the tables)
 - titles become regular Arial 8pt black: left "Gene expression", right "Gene score"
 - the three chemistry names move to the bottom of the heatmap, in full: dsciCAT-seq / 10xMultiome / HT-scCAT-seq
 - the small "Δ / rk" legend line is deleted (SUBRK=1 restores it) - the two columns' meanings are carried by the right colour-bar titles
 - every text in the figure is regular Arial 8pt black - no greys, no bold (in-cell ranks 6pt, column width only 0.52 cm)
Every number is read from tables/*.csv; nothing is recomputed in the figure (ρ/n are only formatted).
Environment: RKCMAP rank colour map / RKBIN rank binning levels (0=continuous); COMPACT=1 the main-figure compact row-pitch version; WD/WK/WG/RH column widths and row height (cm); LM/GAPC/RMC/TOP/BOT margins (cm);
 VERBOSE=1 adds n; MARKRANKS=1 adds #rna->#atac on the right; SUBRK=1 restores the Δ/rk small line; SUF output-name suffix; PNGOUT writes a PNG.
"""
import os, sys
# ---- portability (2026-11-18, GitHub cleanup): <work> was hard-coded here and would not run in the repo.
#       HERE = this script's directory; WORK = the analysis-tree root (DSCI_WORK_ROOT, inferred from the repo location by default).
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("DSCI_WORK_ROOT", os.path.abspath(os.path.join(HERE, "..", "..", "..")))
sys.path.insert(0, HERE)              # fig6_palette.py ships in the same directory as this script
if os.environ.get("DSCI_PYLIBS"):
    sys.path.insert(0, os.environ["DSCI_PYLIBS"])
sys.path.insert(0, os.path.join(WORK, "scripts", "MS_GWAS", "fig6_repro", "deps"))
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.rcParams["pdf.fonttype"] = 42
from matplotlib import font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

FDIR = os.environ.get("DSCI_FONTS", os.path.join(WORK, "fonts"))
for f in ("ARIAL.TTF", "ARIALBD.TTF"):
    if os.path.exists(os.path.join(FDIR, f)):
        fm.fontManager.addfont(os.path.join(FDIR, f))
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams.update({"axes.linewidth": 0.55, "font.size": 8.0})
CM = 2.54
_SHIPPED = os.path.join(os.path.dirname(HERE), "tables_v16")     # figS6/tables_v16
TABD = os.environ.get("TABD", _SHIPPED if os.path.exists(os.path.join(_SHIPPED, "fig6D_delta1_74traits_per_tech.csv"))
                       else os.path.join(WORK, "figures/routeA_out/D_crosschemistry/tables"))
OUT = os.environ.get("OUTDIR", os.environ.get("S6_OUT", os.path.join(WORK, "FY/oligo/fig6_rework_v2")))
COMPACT = os.environ.get("COMPACT", "0") == "1"
os.makedirs(OUT, exist_ok=True)

D = pd.read_csv(TABD + "/fig6D_delta1_74traits_per_tech.csv")
TECH = ["dsciCAT-seq (own)", "10x Multiome", "HT-scCAT-seq"]
NAME = ["dsciCAT\nseq", "10x\nMultiome", "HT-scCAT\nseq"]      # bottom full names (author-specified); one 8pt horizontal line would need 1.95 cm per group,
#  so D would grow to 15.5 cm and leave the main-figure row; hence two lines; at a 1.26 cm group pitch the widest line "Multiome" (1.16 cm) never overlaps
NAME1L = ["dsciCAT-seq", "10xMultiome", "HT-scCAT-seq"]
SHORT = ["dsciCAT-seq", "10x", "HT"]
from fig6_palette import TRAIT_COLOR   # shared with C: MS red / SCZ dark blue / HEIGHT grey
import matplotlib.patheffects as pe
# Submission hard gate: every visible glyph >= 5.5 pt (scripts/final_audit.py:206). The author asked to "shrink a bit",
# 5.0 pt would fit the column but break this gate, so the shrink stops at 5.5 pt (rank column 0.55 cm; a two-digit 5.5 pt number needs 0.22 cm, still comfortable).
FSR = float(os.environ.get("FSR", "5.5"))          # in-cell rank size (was 6.0; the author asked to shrink further; the 5.5 pt hard gate sets the floor)
# No path_effects: measured - the PDF backend then writes text as vector outlines, invisible to get_text()
# (the "every face number traces to a table row" acceptance gate would fail on the spot). Instead the same number is printed twice: a larger white layer underneath,
# then the colour layer: both are real text, white hugs only the glyph outline, and unlike a white box it does not hide the heatmap's thin cells.
FSW = float(os.environ.get("FSW", "6.8"))          # white underlay size = FSR + 1.3, ~0.65 pt visible all around
MARK = {"PASS_Multiple_sclerosis": "MS",
        "PASS_Schizophrenia_Pardinas2018": "SCZ",
        "UKB_460K.body_HEIGHTz": "HEIGHT"}
# ---- three placements for the rank digits (author 11-03 question: raw digits on the cell are illegible; a white backing is legible but hides the heatmap) ----
# The root cause is geometry: the delivered row height 0.058 cm = 1.64 pt, while the hard gate demands >= 5.5 pt type - the glyph is inevitably taller than the row,
# so "in the cell" and "not covering neighbouring rows" cannot both hold at equal row heights; only geometry or placement can change.
#   halo  = current: equal row heights, colour glyph + enlarged white underlay (presses on ~1 row above and below)
#   cell  = the 3 named rows grow to TALLH cm and the other 71 shrink proportionally (total heatmap height unchanged);
#           the named rows' rank cells stay blank and the digit prints inside its own row - no neighbouring row touched, no white underlay needed
#   none  = no digits: the heatmap untouched, all ranks go to the caption
NUMMODE = os.environ.get("NUMMODE", "none")   # author 11-03, seventh revision, decided: no ranks on the face; rankings go into the text (halo/cell kept as fallback)
TALLH = float(os.environ.get("TALLH", "0.22"))      # cell mode: physical height of a named row (cm, = 6.2 pt)
FSD = float(os.environ.get("FSD", "6.0"))           # cell mode: rank size (tall enough rows, larger than halo's 5.5)
# In cell mode rows are stretched; if the Δ cell stretched along, it would become a 4x-area block (area = visual weight, easily read as "bigger effect").
# KEEPD=1: a named row's Δ cell keeps the regular row height (centred), blank above and below - colour and area unchanged; only a reading strip is inserted into the row.
KEEPD = os.environ.get("KEEPD", "1") == "1"
FS = 8.0                                                        # size of every annotation in the figure (author: regular Arial 8pt black)

# ---- row order: only the expression channel (the reference channel) orders; chromatin does not participate ----
rk = D.pivot(index="trait", columns="tech", values="rank_rna_adj")[TECH]
order = rk.mean(axis=1).sort_values().index.tolist()
P = {t: D[D.tech == t].set_index("trait").loc[order] for t in TECH}
NT = len(order)

# ---- column layout ----
W_D = float(os.environ.get("WD", 0.55 if COMPACT else 0.92))   # Δ column width (cm)
W_K = float(os.environ.get("WK", 0.55 if COMPACT else 0.92))   # rank column width
W_G = float(os.environ.get("WG", 0.16 if COMPACT else 0.26))   # chemistry gap
RH = float(os.environ.get("RH", 0.058 if COMPACT else 0.112))  # row height (cm)
# Row height became "per-row physical height" so named rows can grow alone without changing the total heatmap height (the y axis then runs in cm, not rows)
MARKED = [i for i, tr in enumerate(order) if tr in MARK]
if NUMMODE == "cell" and MARKED:
    base = (NT * RH - len(MARKED) * TALLH) / (NT - len(MARKED))
    assert base > 0, "TALLH too large: the remaining rows would go negative"
    HROW = np.where(np.isin(np.arange(NT), MARKED), TALLH, base)
else:
    HROW = np.full(NT, RH)
YEDGE = np.concatenate([[0.0], np.cumsum(HROW)])     # row edges (cm)
YC = (YEDGE[:-1] + YEDGE[1:]) / 2.0                  # row centres (cm)
BASEH = float(np.min(HROW))                          # the regular row height (named Δ cells use this height under KEEPD)
xs = []
for i in range(3):
    xs += [W_D, W_K]
    if i < 2:
        xs.append(W_G)
NC = len(xs)                              # 8
xedge = np.concatenate([[0.0], np.cumsum(xs)])
XCEN = (xedge[:-1] + xedge[1:]) / 2.0
POS = [(3 * i, 3 * i + 1) for i in range(3)]        # (Δ column index, rank column index)
assert POS[2][1] == NC - 1, POS

# ---- layout: explicit cm placement. gridspec wspace is a column-width fraction; no margin you compute can fit an 8pt colour-bar title, so it inevitably spills into the heatmap ----
LM = float(os.environ.get("LM", 2.05 if COMPACT else 2.30))     # left trait names
GAPC = float(os.environ.get("GAPC", 0.60 if COMPACT else 0.80))  # between the two heatmaps
RMC = float(os.environ.get("RMC", 1.45 if COMPACT else 1.55))    # right colour bars
TOP = float(os.environ.get("TOP", 0.46))                          # top titles
BOT = float(os.environ.get("BOT", 0.68 if COMPACT else 0.72))      # bottom chemistry names (two lines of 8pt)
TW = 2 * xedge[-1] + LM + GAPC + RMC
TH = YEDGE[-1] + TOP + BOT
fig = plt.figure(figsize=(TW / CM, TH / CM), dpi=300)


def axpos(x0):
    return [x0 / TW, BOT / TH, xedge[-1] / TW, (TH - TOP - BOT) / TH]


def z(a):
    a = np.asarray(a, float)
    return (a - a.mean()) / a.std(ddof=0)


LIM, RMAX = 1.40, float(NT)
CD = mpl.colormaps["coolwarm"].copy(); CD.set_bad("none")
# ---- rank colour map (author 11-03: with the digits off the face, colours must convey “roughly how different”) ----
# The semantic stays: DARKER = BETTER rank (rank 1 darkest), so all built-in ramps use the `_r` variants.
# Why measured rather than “just pick something nice”: 74 rows in 4.29 cm - the colour difference between adjacent rows is all the information available.
# Sampling the ramp into 74 levels and computing adjacent-pair Lab ΔE (JND ~ 2.3): the current greyblue has median **1.09** and **0%** of pairs above JND -
# — “can't tell them apart” is a necessity, not a taste issue. Candidates: cividis_r 2.13 / viridis_r 3.13 / inferno_r 3.61.
# `RKBIN>0` additionally bins the 74 ranks (staircase): adjacent distinguishability comes from colour JUMPS instead of fine differences, at the cost of within-bin order.
RK_STOPS = {"greyblue": ["#F4F6F8", "#8FA0B8", "#1B2A45"],      # the previously delivered ramp (median adjacent ΔE 1.09)
            "ink": ["#FFFFFF", "#000000"],                       # pure grey, ΔE 1.39, no better than the current
            "cividis_r": "cividis_r", "viridis_r": "viridis_r", "inferno_r": "inferno_r",
            "magma_r": "magma_r"}
RKCMAP = os.environ.get("RKCMAP", "cividis_r")
RKBIN = int(os.environ.get("RKBIN", "0"))
_rk = RK_STOPS[RKCMAP]
CR = (mpl.colormaps[_rk].copy() if isinstance(_rk, str)
      else LinearSegmentedColormap.from_list("rk", _rk))
CR.set_bad("none")
if RKBIN:
    CR = CR.resampled(RKBIN)
    CRNORM = mpl.colors.BoundaryNorm(np.linspace(0, 1, RKBIN + 1), RKBIN)
else:
    CRNORM = mpl.colors.Normalize(0, 1)

for ch, ttl, ax in (("rna", "Gene expression", fig.add_axes(axpos(LM))),
                    ("atac", "Gene score", fig.add_axes(axpos(LM + xedge[-1] + GAPC)))):
    Cd = np.full((NT, NC), np.nan)
    Ck = np.full((NT, NC), np.nan)
    key = "d1_rna" if ch == "rna" else "d1_atac"
    rc = "rank_rna_adj" if ch == "rna" else "rank_atac_adj"
    for i, t in enumerate(TECH):
        Cd[:, POS[i][0]] = z(P[t][key].to_numpy())
        Ck[:, POS[i][1]] = 1.0 - (P[t][rc].to_numpy() - 1) / (RMAX - 1)
    if RKBIN:
        Ck[:, [POS[i][1] for i in range(3)]] = (np.floor(Ck[:, [POS[i][1] for i in range(3)]] * RKBIN)
                                                .clip(0, RKBIN - 1) + 0.5)   # bin centres, aligned with CRNORM
    if NUMMODE == "cell":
        # the named rows' rank cells stay blank: the digit is that cell's information; yielding to it loses nothing,
        # and without the white underlay no white block ever covers neighbouring rows (the Δ columns keep their colours intact)
        Ck[np.ix_(MARKED, [POS[i][1] for i in range(3)])] = np.nan
        if KEEPD:
            ZD = Cd.copy()
            Cd[np.ix_(MARKED, [POS[i][0] for i in range(3)])] = np.nan
    Xe, Ye = np.meshgrid(xedge, YEDGE)
    ax.pcolormesh(Xe, Ye, Cd, cmap=CD, vmin=-LIM, vmax=LIM, shading="flat", zorder=2)
    ax.pcolormesh(Xe, Ye, Ck, cmap=CR, norm=CRNORM, shading="flat", zorder=2)
    if NUMMODE == "cell" and KEEPD:
        # the named Δ cells are redrawn at the regular row height (colours still mapped by coolwarm over ±LIM, consistent with the other rows)
        for r in MARKED:
            for i in range(3):
                c = POS[i][0]
                ax.add_patch(Rectangle((xedge[c], YC[r] - BASEH / 2), xs[c], BASEH,
                                       facecolor=CD((np.clip(ZD[r, c], -LIM, LIM) + LIM) / (2 * LIM)),
                                       edgecolor="none", zorder=2.5))
    ax.set_xlim(0, xedge[-1]); ax.set_ylim(YEDGE[-1], 0.0)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    # each chemistry's (Δ, rank) pair gets a shared frame so the eye groups the two columns; chemistries separate by the narrow gap only
    for i in range(3):
        x0, x1 = xedge[POS[i][0]], xedge[POS[i][1] + 1]
        ax.add_patch(Rectangle((x0, YEDGE[0]), x1 - x0, YEDGE[-1] - YEDGE[0], fill=False, ec="#666666", lw=0.5,
                               clip_on=False, zorder=4))
    ax.add_patch(Rectangle((0, YEDGE[0]), xedge[-1], YEDGE[-1] - YEDGE[0], fill=False, ec="#333333", lw=0.7, clip_on=False))
    # offsets are converted to row units by physical height: when row heights change nothing crowds (compact was once caught with 2 overlaps)
    if os.environ.get("SUBRK") == "1":
        y_rk = -0.05
        for i in range(3):
            ax.text(XCEN[POS[i][0]], y_rk, "\u0394", ha="center", va="bottom", fontsize=FS, color="black")
            ax.text(XCEN[POS[i][1]], y_rk, "rk", ha="center", va="bottom", fontsize=FS, color="black")
        ax.set_title(ttl, fontsize=FS, color="black", pad=3 + 8 * FS / 72.0 * 2.83)
    else:
        ax.set_title(ttl, fontsize=FS, color="black", pad=3)
    # bottom chemistry names: a group pitch of only W_D+W_K+W_G ~ 1.2 cm - 8pt full names on one line inevitably overlap (measured 12-15 pt overlap),
    # so they are rotated 90°. Horizontal would need ~0.88 cm columns and grow D to 15 cm, over the main-figure row.
    y_tech = YEDGE[-1] + 0.09
    rot = os.environ.get("TECHROT", "0")
    for i in range(3):
        cx = (xedge[POS[i][0]] + xedge[POS[i][1] + 1]) / 2
        if rot == "90":
            ax.text(cx, y_tech, NAME1L[i], ha="center", va="top", rotation=90, fontsize=FS, color="black")
        else:
            ax.text(cx, y_tech, NAME[i], ha="center", va="top", fontsize=FS, color="black", linespacing=1.12)

# ---- only the three figure-C traits are labelled: left name + short dash; in-cell ranks use a small white box (row height 0.06 cm; a full-row frame would double up) ----
for tr, nm in MARK.items():
    if tr not in order:
        continue
    r = order.index(tr)
    for ch, ax in (("rna", fig.axes[0]), ("atac", fig.axes[1])):
        if NUMMODE == "none":
            break
        for i in range(3):
            v = int(P[TECH[i]]["rank_rna_adj" if ch == "rna" else "rank_atac_adj"].iloc[r])
            # author 11-03: rank digits take the trait colour, smaller, white box removed (it hid the heatmap's thin cells).
            # no backing at all, but a white outline via path effect: without it, #D51F26's "11"
            # would vanish on a dark-red cell - the outline costs 0.9 pt of glyph edge, hides no colour, and keeps the digit legible on any cell.
            _cx = (xedge[POS[i][1]] + xedge[POS[i][1] + 1]) / 2
            if NUMMODE == "cell":
                ax.text(_cx, YC[r], str(v), ha="center", va="center", fontsize=FSD,
                        color=TRAIT_COLOR[nm], zorder=7)
            else:
                ax.text(_cx, YC[r], str(v), ha="center", va="center", fontsize=FSW, color="white", zorder=6.5)
                ax.text(_cx, YC[r], str(v), ha="center", va="center", fontsize=FSR, color=TRAIT_COLOR[nm], zorder=7)
    # named traits / dashes use the same trait colours as figure C (fig6_palette is the single source; the acceptance script compares the printed ink on the PDF)
    tc = TRAIT_COLOR[nm]
    fig.axes[0].plot([-0.42, -0.16], [YC[r], YC[r]], color=tc, lw=1.0, clip_on=False, zorder=7)
    fig.axes[0].text(-0.50, YC[r], nm, ha="right", va="center", fontsize=FS, color=tc)
    if os.environ.get("MARKRANKS") == "1":
        fig.axes[1].text(xedge[-1] + 0.22, YC[r], "  ".join(
            "#%d\u2192#%d" % (P[t].rank_rna_adj.iloc[r], P[t].rank_atac_adj.iloc[r]) for t in TECH),
            ha="left", va="center", fontsize=FS, color="black")

# ---- right colour bars: two vertical bars stacked, ticks right-aligned, titles horizontal above (side-by-side 8pt titles would necessarily spill into the heatmap) ----
bx = (LM + 2 * xedge[-1] + GAPC + 0.42) / TW
bw = 0.13 / TW
for y0, h0, cmap, lab, ticks, fmt in (
        (0.575, 0.245, CD, "\u03941 (z)", [-LIM, 0, LIM], lambda v, p: "%+.1f" % v),
        (0.185, 0.245, CR, "Rank", [0, 0.5, 1], None)):
    cax = fig.add_axes([bx, y0, bw, h0])
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=(CRNORM if cmap is CR else
                                                mpl.colors.Normalize(-LIM, LIM)))
    sm.set_array([-LIM, LIM] if cmap is CD else [0, 1])
    mpl.colorbar.ColorbarBase(cax, mappable=sm, orientation="vertical", ticks=ticks)
    cax.yaxis.set_ticks_position("right"); cax.yaxis.set_label_position("right")
    if fmt:
        cax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(fmt))
    else:
        cax.set_yticklabels(["74", "37", "1"])
    cax.tick_params(labelsize=FS, length=1.6, width=0.5, pad=1.6)
    for s in cax.spines.values():
        s.set_linewidth(0.5)
    cax.text(0.5, 1.19, lab, transform=cax.transAxes, ha="center", va="bottom",
             fontsize=FS, color="black")

if os.environ.get("VERBOSE") == "1":
    LAD = pd.read_csv(TABD + "/fig6D_crosschem_ladder.csv").set_index("tech")
    for i, t in enumerate(TECH):
        fig.text(0.0075, 0.93 - 0.045 * i, "n %s: %s / %s" % (SHORT[i], f"{int(LAD.loc[t,'n_oligo']):,}",
                 f"{int(LAD.loc[t,'n_opc']):,}"), fontsize=FS, color="black", va="top", ha="left")

# ---- self-check: the boxed ranks and the Δ columns' z must trace digit-by-digit into the delivered tables ----
for tr, nm in MARK.items():
    if tr not in order:
        continue
    r = order.index(tr)
    chk = []
    for i, t in enumerate(TECH):
        dr, da = P[t].d1_rna.iloc[r], P[t].d1_atac.iloc[r]
        zr = (dr - P[t].d1_rna.mean()) / P[t].d1_rna.std(ddof=0)
        za = (da - P[t].d1_atac.mean()) / P[t].d1_atac.std(ddof=0)
        kr, ka = int(P[t].rank_rna_adj.iloc[r]), int(P[t].rank_atac_adj.iloc[r])
        assert abs(dr - float(D[(D.tech == t) & (D.trait == tr)].d1_rna.iloc[0])) < 1e-9, (nm, t, "d1_rna")
        assert abs(da - float(D[(D.tech == t) & (D.trait == tr)].d1_atac.iloc[0])) < 1e-9, (nm, t, "d1_atac")
        assert kr == int(D[(D.tech == t) & (D.trait == tr)].rank_rna_adj.iloc[0]), (nm, t, "rank_rna")
        assert ka == int(D[(D.tech == t) & (D.trait == tr)].rank_atac_adj.iloc[0]), (nm, t, "rank_atac")
        chk.append("%s #%d/#%d (z %+.2f/%+.2f)" % (SHORT[i], kr, ka, zr, za))
    print("[selfcheck] %-7s %s" % (nm, "  ".join(chk)))

suf = ("_compact" if COMPACT else "") + os.environ.get("SUF", "")
out = OUT + "/fig6D_heatmap_pair%s" % suf
fig.savefig(out + ".pdf", bbox_inches="tight", pad_inches=0.03)
if os.environ.get("PNGOUT"):
    fig.savefig(os.environ["PNGOUT"] + "/fig6D_heatmap_pair%s.png" % suf, dpi=300, bbox_inches="tight", pad_inches=0.03)
print("saved %s.pdf | %.2f x %.2f cm nominal | %d rows x %d cols (data %s) x 2 panels | NUMMODE=%s "
      "RKCMAP=%s RKBIN=%d heat=%.2f cm"
      % (out, TW, TH, NT, NC, POS, NUMMODE, RKCMAP, RKBIN, YEDGE[-1]))
