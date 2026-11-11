#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_figS7_v30_panels.py -- the Figure S7 panels of the FINAL (v29) delivery, one PDF each.

The sheet (assembled outside this file, see README "How the sheet is made") is

    row 1   A  74 x 8 rank ladder (as computed | trend removed)   B  tight-pool ladder
            C  human RNA vs ATAC, same donors                     D  mouse RNA vs ATAC, same nuclei
    band    rank colour bar + the three delta definitions          (sits on the ledger's matrix edge)
    row 2   E  74 x 21 axis ledger, full text width
    row 3   F  target x population detection map (GSE118257)      G  S1PR5 detection, per donor
    top     K  identity key band: the seven channels, by accession

inside a 17.00 x 14.20 cm two-column page.  This script draws panels only -- no letters, no composing,
no QA gates -- because those are journal assembly, not content (house spec §1, §9).

v29 vs the v27/v28 prototypes: the panel set, the data and every number are unchanged except
Panel C's MS annotation, which is corrected (see `same_donor_C`: the plotted axes give 2 -> 36, the
delivered face said 2 -> 21, and 21 belongs to a different ATAC series).  Panels C and D are also
newly DRAWN here: in the prototype chain they were carried over as the delivered v22 PDFs, which a
notebook cannot reproduce.  Everything else is the v28 geometry (ledger head margin 0.86 -> 0.24 cm,
detection legend stacked under its own heading).

Typography: the house two-tier scale, 8.0 pt (panel letters, which the composer draws) and 6.0 pt
(everything on the face), Arial, black except the trait/identity colours.  Panels are drawn at their
final width; nothing is rescaled after the fact (spec §9 forbids vector scaling).

Run from the figS7/ directory.  Environment:
  DSCI_WORK_ROOT   analysis tree (default: two levels above the repo) -- only needed for the Arial
                   TTF; every table this script reads is vendored under data/ or tables_v29/
  FIGS7_OUT        where the panel PDFs go (default figS7_out/panels_v30)
  FIGS7_DATA       vendored inputs (default data/)
  FIGS7_TABLES     small tables produced outside this repo (default tables_v29/)

Inputs and who made them (full provenance, including the pipeline parameters, is in README.md):
  data/FIG7_v1/tables/seven4_confound_ALL_cohorts.csv        upstream/confound_all_cohorts.py
  data/SUPPORT_v1/tables/fig7Sb_matched_window.csv           upstream/matched_window_table.py
  data/SUPPORT_v1/tables/fig7Sc_rank_matrix_4cohorts.csv     upstream/axis_rho_and_rank_matrix.py
  data/SUPPORT_v1/tables/S7_atac_crossaxis_human_full.csv    upstream/human_atac_threeaxis.py
  data/SUPPORT_v1/tables/S7_atac_crossaxis_mouse_full.csv    upstream/mouse_atac_threeaxis.py
  data/SUPPORT_v1/tables/S7_rna_crossaxis_mouse_full.csv     upstream/mouse_rna_threeaxis.py
  data/A6_GSE231584_atac_v1/tables/delta1_vs_peakdomain_allchannels.csv   upstream/delta1_vs_peakdomain.py
  data/SUPPORT_v1/tables/figS7c_crosschannel_pairs.csv       upstream/crosschannel_pairs_table.py
  data/SUPPORT_v1/tables/S7_mouse_same_nucleus_delta1.csv|.json  upstream/export_mouse_same_nucleus_delta1.py
  data/A3_GSE118257_v1/tables/A3_GSE118257_target_stats.csv  upstream/target_stats_GSE118257_v7.py
  tables_v29/S41_S1PR5_per_donor.csv, S41_summary_by_cohort.csv   handoff/donor_s1pr5_20261117/
"""
import json
import os
import sys

import numpy as np
import pandas as pd

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt                                    # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, ListedColormap   # noqa: E402
from scipy.stats import rankdata, spearmanr                        # noqa: E402
import fitz                                                        # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
OUT = os.environ.get("FIGS7_OUT", os.path.join(REPO, "figS7_out", "panels_v30"))
DATA = os.environ.get("FIGS7_DATA", os.path.join(REPO, "data"))
TABD = os.environ.get("FIGS7_TABLES", os.path.join(REPO, "tables_v29"))
os.makedirs(OUT, exist_ok=True)
_ROOT = os.environ.get("DSCI_WORK_ROOT", os.path.abspath(os.path.join(REPO, "..", "..")))

# ------------------------------------------------------------------ house style
FONT = os.environ.get("FIG7_FONT") or os.path.join(_ROOT, "fonts", "ARIAL.TTF")
mpl.font_manager.fontManager.addfont(FONT)
FAMILY = mpl.font_manager.FontProperties(fname=FONT).get_name()
mpl.rcParams.update({
    "font.family": FAMILY, "mathtext.fontset": "custom", "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 6.0, "font.weight": "normal", "axes.linewidth": 0.55, "lines.linewidth": 0.8,
    "xtick.labelsize": 6.0, "ytick.labelsize": 6.0, "legend.fontsize": 6.0,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2, "xtick.major.pad": 1.6, "ytick.major.pad": 1.6,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 600,
    "text.color": "0", "axes.labelcolor": "0", "xtick.color": "0", "ytick.color": "0",
})
LET, FS = 8.0, 6.0                    # the only two tiers on the face
CM, PT = 2.54, 72.0 / 2.54
GRN, RED, BLU = "#0B8457", "#D51F26", "#1F77B4"
AXCOL = {"d1": "#9467BD", "d2a": "#E69F00", "d2b": "#17BECF"}      # the three delta axes
cmap_rank = plt.get_cmap("viridis_r")                              # rank 1 (best) = yellow
MINUS = chr(0x2212)
D1 = chr(0x394) + "1"
COH = ["A2", "A3", "A4", "A7"]
GSE = {"A2": "GSE277430", "A3": "GSE118257", "A4": "GSE301908", "A7": "GSE231585"}
MS = "PASS_Multiple_sclerosis"
AD = "PASS_Alzheimers_Jansen2019"
SCZ = "PASS_Schizophrenia_Pardinas2018"
PL = "UKB_460K.blood_PLATELET_COUNT"
AVAIL, ML_SIDE, WIDE, GAP = 15.44, 0.80, 17.0, 0.22


def say(*a):
    print(" ".join(str(x) for x in a), flush=True)


def read(rel, **kw):
    p = os.path.join(DATA, rel)
    assert os.path.exists(p), "missing vendored input %s (see README, Data layout)" % p
    return pd.read_csv(p, **kw)


# ------------------------------------------------------------------ the seven channel colours
# Searched in ../prototype/palette_search7.py: maximise the smallest CIE76 distance between two
# channel colours over {normal, deuteranopia, protanopia}, subject to staying >= 18 from the six
# colours the figure already spends (three delta axes, three trait rule lines) and to a chroma >= 32
# so nothing comes out washed out.  The seven are then ORDERED so that neighbours in a strip row are
# as far apart as possible.  Order is ledger order: four human RNA cohorts, then the three other
# channels.  v25 scored these in the wrong colour space and reported distances about eight CIE points
# too generous; the check below recomputes them on every run and the render refuses to build.
COL7 = ["#E9A87C", "#3CBCEE", "#B1394C", "#24C6B1", "#3D519D", "#579333", "#EE7CD1"]
CHANNELS = ["GSE277430 RNA", "GSE118257 RNA", "GSE301908 RNA", "GSE231585 RNA",
            "GSE231584 ATAC", "GSE250589 RNA", "GSE250589 ATAC"]
GROUPCOL = dict(zip(CHANNELS, COL7))
SEMANTIC = [AXCOL["d1"], AXCOL["d2a"], AXCOL["d2b"], RED, GRN, BLU]
GC = [GROUPCOL[GSE[c] + " RNA"] for c in COH]


# the colour machinery is palette_search7.py's, copied rather than imported so that the render and
# the search cannot disagree about what a colour is.  Two things in it are load-bearing:
#   * the simulation matrices act on LINEAR rgb and return linear rgb.  An earlier version of the
#     check fed them through lin2enc on the way, which put gamma-encoded numbers into a linear-only
#     matrix and reported L* 81 for a colour that prints as dark slate - the v25 palette then looked
#     about eight CIE points safer than it is;
#   * everything is scored after a round trip through 8-bit sRGB, because that is what gets drawn.
XYZ2LIN = np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]])
H_DEUT = np.array([[1, 0, 0], [0.494207, 1.24827, 0], [-0.000962448, 0, 1]])
H_PROT = np.array([[1, 0, 0], [0, 1, 0], [-0.00862408, 0.0931816, 1]])


def hex2lin(h):
    h = h.lstrip("#")
    c = np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)]).astype(float) / 255.
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def lin2hex(lin):
    g = np.clip(lin, 0, 1)
    enc = np.where(g <= 0.0031308, 12.92 * g, 1.055 * g ** (1 / 2.4) - 0.055)
    return "#%02X%02X%02X" % tuple(np.clip(np.round(enc * 255), 0, 255).astype(int))


def lab1(lin):
    X, Y, Z = XYZ2LIN @ (np.clip(lin, 0, 1) * 100)
    f = lambda t: np.where(t > 0.008856, np.cbrt(t), 7.787 * t + 16 / 116)
    fx, fy, fz = f(X / 95.047), f(Y / 100.), f(Z / 108.883)
    return np.array([116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz)])


def variants(h):
    """the three colour-vision views of one hex, as CIE L*a*b*"""
    lin = hex2lin(lin2hex(hex2lin(h)))                  # 8-bit round trip, then score what is drawn
    return [lab1(lin), lab1(np.clip(H_DEUT @ lin, 0, 1)), lab1(np.clip(H_PROT @ lin, 0, 1))]


def check_palette():
    """adjacent = smallest distance between two strips that sit next to each other in a row (the
    comparison a reader actually has to make); pair = any two of the seven; semantic = any channel
    colour against one of the six colours the figure already spends.  Minimum over the three vision
    models.  The three floors below are the ones the search was run with; a later colour edit that
    breaks one stops the render instead of becoming a note in a markdown file."""
    V = [variants(h) for h in COL7]
    Sv = [variants(h) for h in SEMANTIC]
    adj = min(min(np.linalg.norm(V[i][k] - V[i + 1][k]) for k in range(3)) for i in range(6))
    pair = min(min(np.linalg.norm(V[i][k] - V[j][k]) for k in range(3))
               for i in range(7) for j in range(i + 1, 7))
    sem = min(min(np.linalg.norm(V[i][k] - Sv[j][k]) for k in range(3))
              for i in range(7) for j in range(6))
    say("palette: 7 channel colours - adjacent strips min CIE76 %.1f, any pair %.1f, vs the six "
        "semantic colours %.1f (worst of normal / deuteranopia / protanopia)" % (adj, pair, sem))
    assert sem >= 17.0, ("a channel colour is the same colour as an axis or rule line (%.1f)" % sem)
    assert pair >= 35.0, ("two channel colours are too close (%.1f)" % pair)
    assert adj >= 50.0, ("neighbouring strips are too close to tell apart (%.1f)" % adj)
    return adj, pair, sem


# ------------------------------------------------------------------ inputs, with the assertions the figure relies on
S4 = read("FIG7_v1/tables/seven4_confound_ALL_cohorts.csv")          # 74 traits x 4 cohorts
MW = read("SUPPORT_v1/tables/fig7Sb_matched_window.csv")             # +-15% peak-domain window
RM = read("SUPPORT_v1/tables/fig7Sc_rank_matrix_4cohorts.csv").set_index("trait")
SD = read("A6_GSE231584_atac_v1/tables/same_donor_RNA_vs_CHROMATIN_delta1_74traits.csv")
CAXF = read("SUPPORT_v1/tables/S7_atac_crossaxis_human_full.csv", index_col=0, header=[0, 1])
MX = read("SUPPORT_v1/tables/S7_atac_crossaxis_mouse_full.csv", index_col=0)
RX = read("SUPPORT_v1/tables/S7_rna_crossaxis_mouse_full.csv", index_col=0)
PK = read("A6_GSE231584_atac_v1/tables/delta1_vs_peakdomain_allchannels.csv")
CP = read("SUPPORT_v1/tables/figS7c_crosschannel_pairs.csv")
SNU = read("SUPPORT_v1/tables/S7_mouse_same_nucleus_delta1.csv")     # panel D: same nuclei, both channels
SNUJ = json.load(open(os.path.join(DATA, "SUPPORT_v1/tables/S7_mouse_same_nucleus_delta1.json")))
BUB = read("A3_GSE118257_v1/tables/A3_GSE118257_target_stats.csv")   # panel F
DN = pd.read_csv(os.path.join(TABD, "S41_S1PR5_per_donor.csv"))      # panel G, per donor
S41 = pd.read_csv(os.path.join(TABD, "S41_summary_by_cohort.csv"), index_col=0)

# A: the ladder rows are ordered by peak-domain size, exactly as the delivered sheet (37 = MS, 56 = PLATELET)
S = S4.sort_values(["n_peak_raw", "trait"], ascending=[False, True]).reset_index(drop=True)
assert len(S) == 74 and S.n_peak_raw.nunique() == 73
iMS, iPL = int(S.index[S.trait == MS][0]), int(S.index[S.trait == PL][0])
assert (iMS + 1, iPL + 1) == (37, 56), (iMS + 1, iPL + 1)
LO, HI = float(MW.window_lo.iloc[0]), float(MW.window_hi.iloc[0])
inw = ((S.n_peak_raw >= LO) & (S.n_peak_raw <= HI)).values
w0, w1 = int(np.where(inw)[0][0]), int(np.where(inw)[0][-1])
assert w1 - w0 + 1 == int(MW.n_in_window.iloc[0]) == 40 and inw[w0:w1 + 1].all()

RAW = np.array([S["rank_raw_%s" % c].values.astype(float) for c in COH]).T     # as computed
ADJ = np.array([S["rank_adj_%s" % c].values.astype(float) for c in COH]).T     # trend removed
ZRAD = np.column_stack([RAW, ADJ])
assert ZRAD.shape == (74, 8)
ms_raw = [int(v) for v in RAW[iMS]]
ms_adj = [int(v) for v in ADJ[iMS]]
assert ms_raw == [2, 4, 7, 2] and ms_adj == [1, 1, 3, 1], (ms_raw, ms_adj)
RHO_LN = {c: float(MW[MW.cohort == c].rho_lnp_vs_d1.iloc[0]) for c in COH}

# B: the tight pool = MS + its 15 nearest neighbours on ln(peak-domain size)
lnp = S.lnp.values
d = np.abs(lnp - lnp[iMS]); d[iMS] = np.inf
POOL = np.sort(np.append(np.sort(np.argsort(d)[:15]), iMS))
PR = {c: int(rankdata(-S[c].values[POOL].astype(float), method="min")[list(POOL).index(iMS)]) for c in COH}
assert [PR[c] for c in COH] == [1, 1, 1, 1] and len(POOL) == 16

# E: the ledger.  21 columns = four human RNA cohorts + human ATAC + mouse RNA + mouse ATAC, three
# delta axes each; rows are ordered by the mean of the four human cohorts' AS-COMPUTED delta-1 rank.
hATAC = {ax: CAXF[(chr(0x394), ax)].rank(ascending=False, method="min") for ax in ("d1", "d2a", "d2b")}
mATAC = {ax: MX[chr(0x394) + ax].rank(ascending=False, method="min") for ax in ("d1", "d2a", "d2b")}
mRNA = {ax: RX[chr(0x394) + ax].rank(ascending=False, method="min") for ax in ("d1", "d2a", "d2b")}
rms = {}
for ax in ("d1", "d2a", "d2b"):
    suf = "1" if ax == "d1" else "A" if ax == "d2a" else "B"
    rms[ax] = RM[[c for c in RM.columns if c.endswith("_" + suf)]].mean(axis=1)
JROWS = rms["d1"].sort_values().index.tolist()
GRP_J = [(GSE[c], [RM["c%d_1" % k], RM["c%d_A" % k], RM["c%d_B" % k]]) for k, c in enumerate(COH)]
GRP_J += [(chr(0x394) + "ATAC", [hATAC[a] for a in ("d1", "d2a", "d2b")]),
          (chr(0x394) + "ATAC", [mATAC[a] for a in ("d1", "d2a", "d2b")]),
          (chr(0x394) + "RNA", [mRNA[a] for a in ("d1", "d2a", "d2b")])]
ZJ = np.column_stack([g[1][si].reindex(JROWS).values.astype(float) for g in GRP_J for si in range(3)])
assert ZJ.shape == (74, 21)
# the order key is the mean of the RAW ranks (c*_1 == rank_raw_*), which is NOT the Fig 7B key
assert np.allclose(RM["c0_1"].astype(float).values, S4.set_index("trait")["rank_raw_A2"].astype(float)
                   .reindex(RM.index).values)
assert int(np.where(np.array(JROWS) == MS)[0][0]) + 1 == 2, "MS row position in the ledger moved"

# G: the donor table.  A population only enters a donor's comparison if that donor contributes
# MIN_POP nuclei of it; a donor is plotted if OL has >= MIN_POP nuclei and one comparator survives.
# the plotted order is margin ascending: the jitter is drawn from a fixed seed against THIS index,
# so the row order is part of the figure and must not be left to the table's own order
EV = DN[DN.evaluable].sort_values("margin_pp").reset_index(drop=True)
assert len(DN) == 28 and len(EV) == 27 == int(S41.loc["ALL", "donors_evaluable"])
assert int(EV.rule.sum()) == int(S41.loc["ALL", "rule_met"]) == 22, int(EV.rule.sum())
assert [int(EV[EV.cohort == g].rule.sum()) for g in ("GSE301908", "GSE118257", "GSE231585")] == [8, 7, 7]


# ------------------------------------------------------------------ canvas mechanics
fig = None
FR = {}
INK = [0.0, 0.0, 0.0, 0.0]


def new_page(w, h):
    global fig
    fig = plt.figure(figsize=(w / CM, h / CM))
    return fig


def ax_cm(x0, y0, w, h):
    """x0, y0 = top-left corner of the axes, in cm on this panel's own canvas."""
    W, H = fig.get_size_inches() * CM
    return fig.add_axes([x0 / W, (H - y0 - h) / H, w / W, h / H])


def txt(x, y, s, **kw):
    W, H = fig.get_size_inches() * CM
    kw.setdefault("fontsize", FS)
    kw.setdefault("color", "0")
    return fig.text(x / W, (H - y) / H, s, **kw)


def rect(x, y, w, h, fc, ec="none"):
    """a figure-space rectangle, top-left corner in cm"""
    W, H = fig.get_size_inches() * CM
    fig.add_artist(plt.Rectangle((x / W, 1 - (y + h) / H), w / W, h / H,
                                 fc=fc, ec=ec, zorder=6, lw=0.3 if ec != "none" else 0))


def line(x0, y0, x1, y1, lw=0.5, color="0"):
    W, H = fig.get_size_inches() * CM
    fig.add_artist(plt.Line2D([x0 / W, x1 / W], [1 - y0 / H, 1 - y1 / H],
                              lw=lw, color=color, zorder=3))


def tw(s, fs=FS):
    """measured width of a string in cm on the current canvas - never estimated from the character
    count: every fit-to-width assert in this file (and the rotated-head projection rule) needs the
    real ink width, and 6 pt Arial with digits and capitals differs from any rule of thumb by 10%."""
    t = txt(0, 0, s, fontsize=fs)
    fig.canvas.draw()
    w = t.get_window_extent(renderer=fig.canvas.get_renderer()).width * CM / fig.dpi
    t.remove()
    return w


def swatch(x, y, col, label, sw=0.20, txtofs=None, tail=0.30):
    """a colour strip and its label, returning the x after it.  Two geometries are in use: the
    narrow one (0.20 strip, text 0.29 along) for the delta key, where three entries share a row with
    the rank bar, and the wide one (0.26 strip, text 0.34 along) for the identity key band, where the
    strip has to read as the same object as a ledger column strip."""
    txtofs = sw + 0.09 if txtofs is None else txtofs
    rect(x, y - 0.065, sw, 0.13, col, ec="0.55")
    txt(x + txtofs, y, label, ha="left", va="center")
    return x + txtofs + tw(label) + tail


def head45(x, y, s, ang=45):
    """angled column head, right-anchored at the top (spec §14.2).  The horizontal projection of a
    head is len(s) * cos(theta), and it must fit the column pitch: 45 deg clears every pitch used
    here; only the withdrawn panel F needed 60 deg for "GSE250589 ATAC" (1.73 cm -> 1.22 cm against a
    1.10 cm pitch)."""
    txt(x - 0.02, y, s, ha="right", va="top", rotation=ang, rotation_mode="anchor")
    return tw(s)


def ink_rect(path):
    """the union of text, drawing and image ink of a rendered PDF, in pt"""
    global INK
    d = fitz.open(path); p = d[0]
    PW, PH = p.rect.width, p.rect.height
    X0, Y0, X1, Y1 = [], [], [], []
    for b in p.get_text("dict")["blocks"]:
        for ln in b.get("lines", []):
            for sp in ln.get("spans", []):
                r = sp["bbox"]; X0.append(r[0]); Y0.append(r[1]); X1.append(r[2]); Y1.append(r[3])
    for dr in p.get_drawings():
        r = dr["rect"]
        if r.is_valid and not r.is_empty and r.get_area() < 0.6 * PW * PH:
            X0.append(r.x0); Y0.append(r.y0); X1.append(r.x1); Y1.append(r.y1)
    for im in p.get_image_info():
        bb = fitz.Rect(im["bbox"])
        if bb.is_valid and not bb.is_empty:
            X0.append(bb.x0); Y0.append(bb.y0); X1.append(bb.x1); Y1.append(bb.y1)
    d.close()
    assert X0, "no ink"
    pad = 1.0
    INK = (min(X0), min(Y0), max(X1), max(Y1))
    return fitz.Rect(max(0.0, min(X0) - pad), max(0.0, min(Y0) - pad),
                     min(PW, max(X1) + pad), min(PH, max(Y1) + pad))


def export(key, fname):
    """save, refuse to export anything whose ink touches the canvas edge (that is how a clipped label
    shows up), then crop the page to the ink.  FR[key] carries the cropped size in cm and the ink
    offset, which is what the composer needs to align a panel's matrix edge."""
    tmp = os.path.join(OUT, "_tmp_%s.pdf" % key)
    fig.savefig(tmp)
    rect = ink_rect(tmp)
    d = fitz.open(tmp); pw, ph = d[0].rect.width, d[0].rect.height; d.close()
    ix0, iy0, ix1, iy1 = INK
    assert min(ix0, iy0) > 0.2 and ix1 < pw - 0.2 and iy1 < ph - 0.2, (
        "%s: ink reaches the canvas edge, something is clipped (ink %.1f-%.1f x %.1f-%.1f pt on a "
        "%.1f x %.1f pt canvas)" % (key, ix0, ix1, iy0, iy1, pw, ph))
    s = fitz.open(tmp)
    o = fitz.open()
    pg = o.new_page(width=rect.width, height=rect.height)
    pg.show_pdf_page(pg.rect, s, 0, clip=rect)
    dst = os.path.join(OUT, fname)
    o.save(dst, deflate=True)
    o.close(); s.close(); os.remove(tmp); plt.close(fig)
    FR[key] = dict(file=fname, w=rect.width / PT, h=rect.height / PT,
                   ink_x0=INK[0] / PT, ink_y0=INK[1] / PT)
    say("panel %-3s -> %-30s %5.2f x %5.2f cm" % (key, fname, rect.width / PT, rect.height / PT))
    return dst


# %% PANEL K
# ================================================================ K  identity key band (legend artwork)
def key_band():
    """K: every colour the figure spends, in one row at the top, grouped by assay channel.

    The author's two complaints about earlier keys were "one row only" and "mouse datasets must
    follow the legend too".  Both are satisfied here: seven entries, one row, and each entry is an
    accession - no "human RNA" / "mouse ATAC" wording anywhere on the sheet.  The channel is stated
    once per block instead of once per entry, which is what buys the room: GSE250589 legitimately
    appears in both blocks, and the headings are what tell its two columns apart.
    """
    W, H = AVAIL, 0.46
    new_page(W, H)
    x, y = 0.05, 0.24
    for head_, names in (("RNA", [CHANNELS[0], CHANNELS[1], CHANNELS[2], CHANNELS[3], CHANNELS[5]]),
                         ("ATAC", [CHANNELS[4], CHANNELS[6]])):
        txt(x, y, head_, ha="left", va="center")
        x += tw(head_) + 0.26
        for nm in names:
            x = swatch(x, y, GROUPCOL[nm], nm.replace(" RNA", "").replace(" ATAC", ""), sw=0.26,
                       txtofs=0.34)
        x += 0.34
    x -= 0.34 + 0.30
    assert x <= W + 0.01, ("one-row key needs %.2f cm of %.2f" % (x, W))
    say("key band: one row, %d entries, %.2f cm of the %.2f cm text column" % (len(CHANNELS), x, W))
    export("K", "S7v30_K_keyband.pdf")


# %% PANEL A
# ================================================================ A  74 x 8 rank ladder
def ladder_A():
    """A: the 74 traits x four human cohorts, twice - the rank as scDRS produced it, and the rank
    after the peak-domain-size trend is removed.  Rows are ordered by peak-domain size (n_peak_raw
    descending), the delivered sheet's key, NOT by rank; MS is row 37, PLATELET row 56.  The eight
    columns are the same four cohorts in both halves, so the strips repeat the cohort colours."""
    W, gy, gh = 4.55, 0.86, 2.97
    gx, gw = 1.19, 3.20
    H = gy + gh + 0.40
    new_page(W, H)
    pitch = gw / 8.0
    a = fig.add_axes([gx / W, 1 - (gy + gh) / H, gw / W, gh / H])
    a.imshow(ZRAD, aspect="auto", cmap=cmap_rank, vmin=1, vmax=74)
    a.set_xlim(-0.5, 7.5); a.set_ylim(73.5, -0.5); a.set_xticks([]); a.set_yticks([])
    a.tick_params(length=0)
    a.axhline(iMS, color=RED, lw=0.9, zorder=7)
    a.axhline(iPL, color="0", lw=0.9, zorder=7)
    a.axvline(3.5, color="0", lw=0.8, zorder=8)
    halfw = 4 * pitch
    for half, lab in zip((0, 4), ("as computed", "trend removed")):
        wcm = tw(lab)
        assert wcm <= halfw - 0.08, ("half label %r is %.2f cm, half is %.2f" % (lab, wcm, halfw))
        xc = gx + (half + 2) * pitch
        txt(xc, 0.28, lab, ha="center", va="center")
        x0, x1 = gx + half * pitch + 0.06, gx + (half + 4) * pitch - 0.06
        line(x0, 0.46, x1, 0.46)                       # the brace over each half
        line(x0, 0.46, x0, 0.53); line(x1, 0.46, x1, 0.53)
    for k in range(8):
        rect(gx + k * pitch + 0.015, 0.66, pitch - 0.03, 0.13, GC[k % 4], ec="0.55")
    txt(gx - 0.07, gy + gh * (iMS + 0.5) / 74.0, "MS", ha="right", va="center", color=RED)
    txt(gx - 0.07, gy + gh * (iPL + 0.5) / 74.0, "PLATELET", ha="right", va="center")
    yl = gy + gh + 0.17
    txt(gx - 0.07, yl, "MS rank", ha="right", va="center")
    for k in range(4):
        txt(gx + (k + 0.5) * pitch, yl, str(ms_raw[k]), ha="center", va="center",
            color=RED if ms_raw[k] <= 2 else "0")
        txt(gx + (4 + k + 0.5) * pitch, yl, str(ms_adj[k]), ha="center", va="center", color=RED)
    export("A", "S7v30_A_ladder.pdf")


# %% PANEL B
# ================================================================ B  tight-pool ladder
def tight_pool_B():
    """B: the same four cohorts, but ranked inside a tight pool - MS plus the 15 traits whose
    peak-domain size is nearest to MS's (|ln n_peak| nearest neighbours), 16 traits in all.  Inside
    that pool MS is rank 1 in every cohort, which is the size-matched version of the claim; the
    +-15% window and its rho live in Table S7, not on the face."""
    POOLRANK = np.full((74, 4), np.nan)
    for k, c in enumerate(COH):
        POOLRANK[POOL, k] = rankdata(-S[c].values[POOL].astype(float), method="min")
    W, gy, gh = 3.30, 0.86, 2.97
    bx, bw = 0.04, 3.20
    H = gy + gh + 0.16
    new_page(W, H)
    a = fig.add_axes([bx / W, 1 - (gy + gh) / H, bw / W, gh / H])
    cmap_pool = ListedColormap(cmap_rank(np.linspace(0, 1, 17))[:16]); cmap_pool.set_bad("0.88")
    a.imshow(np.ma.masked_invalid(POOLRANK), aspect="auto", cmap=cmap_pool, vmin=1, vmax=16)
    a.set_xlim(-0.5, 3.5); a.set_ylim(73.5, -0.5); a.set_xticks([]); a.set_yticks([])
    a.tick_params(length=0); a.set_facecolor("0.88")
    a.axhline(iMS, color=RED, lw=0.9, zorder=7)
    a.axhline(iPL, color="0", lw=0.9, zorder=7)
    pitch = bw / 4.0
    for k in range(4):
        rect(bx + k * pitch + 0.02, 0.66, pitch - 0.04, 0.13, GC[k], ec="0.55")
    export("B", "S7v30_B_tightpool.pdf")


def _scatter_pair(a, x, y, mspt, col, lbl, lines, xlim, ylim):
    """the shared body of C and D: 74 trait points, dashed zero cross, MS highlighted, the two
    statistics lines in the upper-left corner.

    The view is not chosen here: xlim / ylim are the delivered panels' own limits, recovered by
    regressing the 74 dot centres against the tick marks of the signed PDFs, and verify_v29_ink.py
    re-checks every dot to <0.05 pt.  Autoscaling would give a slightly tighter frame and would move
    the statistics text onto the cloud; D in particular keeps 1.1 cm of headroom above its points
    because that is where its two annotation lines live.  The 0.25 cm of canvas the callers leave
    above the axes is not decoration either - the top y tick label is centred on its tick and would
    otherwise hang off the page, which export() rightly refuses to export."""
    a.scatter(x, y, s=6, color=col, lw=0, alpha=0.85, zorder=3)
    DASH = (0, (2.0, 1.6))
    a.axhline(0, color="0", lw=0.6, zorder=2, ls=DASH)
    a.axvline(0, color="0", lw=0.6, zorder=2, ls=DASH)
    a.scatter([mspt[0]], [mspt[1]], s=26, facecolor=RED, edgecolor="black", lw=0.5, zorder=5)
    # spec 14.8: a highlighted point must say who it is, next to itself
    a.annotate(lbl, xy=mspt, xytext=(0, -9), textcoords="offset points", ha="center", va="center",
               fontsize=FS, color=RED, zorder=8,
               bbox=dict(boxstyle="square,pad=0.10", fc="white", ec="none", alpha=0.92))
    a.set_xlim(*xlim); a.set_ylim(*ylim)
    a.text(0.035, 0.965, "\n".join(lines), transform=a.transAxes, fontsize=FS, ha="left", va="top",
           linespacing=1.35, zorder=6)
    a.tick_params(labelsize=FS, length=1.6, colors="0")
    a.spines[["top", "right"]].set_visible(False)


# %% PANEL C
# ================================================================ C  human RNA vs ATAC, same donors
def same_donor_C():
    """C: does the RNA delta-1 order predict the ATAC delta-1 order in the same people?

    x = GSE231585 RNA delta-1 (the delivered A7 sweep, all samples), y = GSE231584 ATAC delta-1
    (size-matched libraries).  The two assays were built separately from the same donors, which is
    why the pair is labelled "separately profiled" in figS7c_crosschannel_pairs.csv.  The rho on the
    face is the Spearman of exactly these two columns (0.722, and 0.137 after residualising both on
    ln n_peak, P = 0.24) - it is read from the pairs table, never typed.

    v29 change, and it is a correction: the delivered face said "MS 2 -> 21".  2 is MS's rank in the
    plotted x and 21 is its rank in a DIFFERENT series (the donor-paired CHR column of
    same_donor_RNA_vs_CHROMATIN_delta1_74traits.csv, which is not what is plotted here).  On the
    plotted axes MS is rank 2 in RNA and rank 36 in ATAC, so the arrow now says 2 -> 36.  Both
    numbers are computed from the plotted columns and asserted.
    """
    HUM = S4.merge(PK[["trait", "median_groupdiff"]], on="trait")
    assert len(HUM) == 74
    i = int(HUM.index[HUM.trait == MS][0])
    x, y = HUM["A7"].values, HUM["median_groupdiff"].values
    rank_x, rank_y = int((x > x[i]).sum()) + 1, int((y > y[i]).sum()) + 1
    assert (rank_x, rank_y) == (2, 36), (rank_x, rank_y)
    u = CP[CP.pair.astype(str).str.contains("GSE231585 " + chr(0x2194) + " GSE231584", regex=False)].iloc[0]
    assert abs(float(u.rho_raw) - spearmanr(x, y).statistic) < 1e-4, "the plotted cloud is not the rho"
    W, Hh = 3.45, 4.00
    aw, ah, pad_l, pad_t = 2.198, 2.850, 1.147, 0.25
    new_page(W, Hh)
    a = ax_cm(pad_l, pad_t, aw, ah)
    # the dot colour is the delivered one (BLU here, GRN in D): the cloud is a cross-CHANNEL pair,
    # so it does not belong to a cohort and must not borrow a cohort colour from the key band
    _scatter_pair(a, x, y, (x[i], y[i]), BLU, ("MS %d " + chr(0x2192) + " %d") % (rank_x, rank_y),
                  [(chr(0x3C1) + " %.3f " + chr(0x2192) + " %.3f") % (u.rho_raw, u.rho_adj),
                   "P %.2f" % u.P_adj],
                  xlim=(-6.689, 3.874), ylim=(-0.5335, 0.6999))
    a.set_xlabel("RNA " + D1 + " (GSE231585)", fontsize=FS)
    a.set_ylabel("ATAC " + D1 + " (GSE231584)", fontsize=FS)
    a.set_xticks([-5, 0]); a.set_yticks([-0.50, -0.25, 0.00, 0.25, 0.50])
    export("C", "S7v30_C_human_samedonor.pdf")


# %% PANEL D
# ================================================================ D  mouse EAE, same nucleus
def same_nucleus_D():
    """D: the mouse EAE multiome, where RNA and ATAC are measured in the SAME nucleus, so the
    comparison needs no donor pairing at all.  x = RNA delta-1, y = ATAC delta-1, 74 traits, 14,928
    nuclei in 25 libraries (the numbers come from the JSON that the export script wrote).  rho 0.248
    with P 0.033: much weaker than the human same-donor pair, which is the point - the two channels
    of one nucleus disagree more than two channels of one donor do.
    """
    assert len(SNU) == 74 and SNUJ["n_nuclei"] == 14928 and SNUJ["n_libs"] == 25
    i = int(SNU.index[SNU.trait == MS][0])
    x, y = SNU["delta1_rna"].values, SNU["delta1_atac"].values
    rho = spearmanr(x, y)
    assert abs(rho.statistic - SNUJ["spearman_rho"]) < 1e-6 and abs(rho.pvalue - SNUJ["spearman_p"]) < 1e-9
    rank_x = int((x > x[i]).sum()) + 1
    rank_y = int((y > y[i]).sum()) + 1
    assert (rank_x, rank_y) == (SNUJ["ms_rank_rna"], SNUJ["ms_rank_atac"]) == (11, 36)
    W, Hh = 3.00, 4.00
    aw, ah, pad_l, pad_t = 2.198, 2.850, 0.677, 0.25
    new_page(W, Hh)
    a = ax_cm(pad_l, pad_t, aw, ah)
    _scatter_pair(a, x, y, (x[i], y[i]), GRN, ("MS %d " + chr(0x2192) + " %d") % (rank_x, rank_y),
                  [(chr(0x3C1) + " %.3f, P %.2g") % (rho.statistic, rho.pvalue), "same nucleus"],
                  xlim=(-4.396, 3.779), ylim=(-1.800, 1.5046))
    a.set_xlabel("RNA " + D1 + " (GSE250589)", fontsize=FS)
    a.set_ylabel("ATAC " + D1, fontsize=FS, labelpad=1.5)
    a.set_yticks([-1.0, 0.0, 1.0])
    export("D", "S7v30_D_mouse_samenucleus.pdf")




# %% PANEL E
# ================================================================ E  74 x 21 axis ledger
def ledger_E(width=AVAIL):
    """E: the whole rank picture in one table - 74 traits x 21 columns (four human RNA cohorts,
    human ATAC, mouse RNA, mouse ATAC, each on the three delta axes), every cell a rank out of 74.

    Rows are ordered by the mean of the four human cohorts' AS-COMPUTED delta-1 rank (MS lands on
    row 2).  That is deliberately not the Fig 7B key: Fig 7B orders its rows by the mean
    size-ADJUSTED rank, where MS is row 1.  This panel exists to show the un-adjusted basis, so it
    is ordered by it, and the ladder A - which is ordered by peak-domain size - is a third key.  The
    legend says all three; the face cannot.

    The head margin is 0.24 cm, not the 0.86 cm the earlier renderer reserved: 0.62 cm of that was
    for rotated cohort heads that the key band replaced, and leaving it in put 1.43 cm of dead air
    between this panel and the band that decodes it (spec §18.2).
    """
    h_head, h_mat, pad_l, pad_r = 0.24, 2.55, 0.62, 0.10
    Hh = h_head + h_mat + 0.40
    gw = width - pad_l - pad_r
    cpitch = gw / 21.0
    new_page(width, Hh)
    a = fig.add_axes([pad_l / width, 1 - (h_head + h_mat) / Hh, gw / width, h_mat / Hh])
    a.imshow(np.ma.masked_invalid(ZJ), aspect="auto", cmap=cmap_rank, vmin=1, vmax=74)
    a.set_xlim(-0.5, 20.5); a.set_ylim(73.5, -0.5); a.set_xticks([]); a.set_yticks([])
    a.tick_params(length=0); a.set_facecolor("0.88")
    for gi in range(1, 7):
        a.axvline(3 * gi - 0.5, color="0", lw=0.5, zorder=6)          # one line per channel group
    a.axhline(JROWS.index(MS), color=RED, lw=0.9, zorder=7)
    for nm, col in ((AD, GRN), (SCZ, BLU)):
        if nm in JROWS:
            a.axhline(JROWS.index(nm), color=col, lw=0.7, zorder=7)
    for gi, col in enumerate(COL7):                                   # channel strips above the matrix
        rect(pad_l + 3 * gi * cpitch + 0.02, 0.04, 3 * cpitch - 0.04, 0.13, col)
    txt(pad_l - 0.06, h_head + h_mat * (JROWS.index(MS) + 0.5) / 74.0, "MS",
        ha="right", va="center", color=RED)
    ys = h_head + h_mat + 0.06                                        # delta strips below the matrix
    for jj in range(7):
        for si, axk in enumerate(("d1", "d2a", "d2b")):
            rect(pad_l + (3 * jj + si) * cpitch + 0.015, ys, cpitch - 0.03, 0.12, AXCOL[axk])
    export("E", "S7v30_E_ledger.pdf")
    assert abs(0.04 + 0.13 + 0.07 - h_head) < 1e-9, "strip-to-matrix pitch moved"
    say("E  matrix begins %.3f cm inside the panel; band above it will align to that edge"
        % matrix_x0("E"))


def matrix_x0(key):
    """Where a panel's matrix actually begins, in its own cropped canvas coordinates.

    Needed because export() crops to ink, and the ink that sets the left edge is a *text* bbox: the
    ledger's right-aligned "MS" row label.  "The matrix starts pad_l inside the panel" is therefore
    not a statement you can derive from the source constants; you have to measure the produced PDF
    (spec §18.4)."""
    d = fitz.open(os.path.join(OUT, FR[key]["file"]))
    ib = [fitz.Rect(i["bbox"]) for i in d[0].get_image_info()]
    d.close()
    assert ib, key
    return max(ib, key=lambda r: r.width * r.height).x0 / PT


# %% PANEL RB
# ================================================================ RB  the sheet's only rank bar
def rank_bar():
    """RB: the viridis_r scale the three rank tables share, 1 (best) to 74."""
    W, H = 4.44, 0.59
    new_page(W, H)
    lab = "rank (1 = highest of the 74)"
    txt(0.05, 0.20, lab, ha="left", va="center")
    xb = 0.05 + tw(lab) + 0.14
    bw = W - xb - 0.40
    assert bw >= 1.00, ("no room for the bar, %.2f cm" % bw)
    a = ax_cm(xb, 0.135, bw, 0.13)
    a.imshow(np.linspace(1, 74, 256)[None, :], cmap=cmap_rank, vmin=1, vmax=74, aspect="auto")
    a.set_xticks([]); a.set_yticks([]); a.tick_params(length=0)
    for sp in a.spines.values():
        sp.set_linewidth(0.4)
    txt(xb, 0.36, "1", ha="left", va="top")
    txt(xb + bw, 0.36, "74", ha="right", va="top")
    export("RB", "S7v30_RB_rankbar.pdf")


# %% PANEL L2
# ================================================================ L2  the three delta definitions
def axis_key(width):
    """L2: what delta-1 / delta-2A / delta-2B subtract, in one row.

    Spelled from the pipeline's own definitions (A7_GSE231585_delta2_full74.py:13, and the A2/A3/A4
    siblings), not from the loosest possible paraphrase: delta-2A subtracts MG + LYM, i.e. the immune
    pool, not microglia alone, and delta-2B subtracts every non-neuronal nucleus - the manuscript
    abbreviation table still calls it "neuron minus glia", which is a naming problem, not a
    computation problem (see NEEDS_ABBREVIATIONS.md).

    v30 (author 2026-11-18): the key now uses the SAME three strings as the Fig 7C column labels,
    OL/MOL - OPC, AST - (MG+LYM), NEU - non-NEU.  One figure family, one vocabulary.  This costs two
    new rows in the manuscript abbreviation table (NEU, LYM) - they were never missing from the
    pipeline, only from the table; NEEDS_ABBREVIATIONS.md §1 lists exactly what to add."""
    items = [("d1", D1 + " = OL/MOL " + MINUS + " OPC"),
             ("d2a", D1.replace("1", "2A") + " = AST " + MINUS + " (MG+LYM)"),
             ("d2b", D1.replace("1", "2B") + " = NEU " + MINUS + " non-NEU")]
    W, H = width, 0.46
    new_page(W, H)
    x = 0.05
    for k, lab in items:
        x = swatch(x, 0.24, AXCOL[k], lab) + 0.16
    need = x - 0.16 - 0.30
    assert need <= W + 0.01, ("delta key needs %.2f cm, band is %.2f" % (need, W))
    say("delta key: one row, %.2f cm of %.2f" % (need, W))
    export("L2", "S7v30_L2_axiskey.pdf")


# %% PANEL F
# ================================================================ F  target x population detection
def detection_map():
    """F: which cells actually carry the five Fig 7D targets, in GSE118257 (17,799 nuclei, 9 donors
    from 20 libraries).  Area = % of nuclei in that population detecting the gene; fill = mean
    log2(CP10k+1), depth-corrected expression.  The population n under each column is the denominator
    the percentage is taken over, so it is on the face; the diagnosis split is not (Table S7).

    v28 moved the whole legend block up under its own "% detecting" heading: it had been sitting
    2.1 cm below it, in the strip once reserved for the rotated population labels (spec §18.3).  The
    1.62 cm of label room at the bottom is measured - the longest rotated "name  n" label is 1.47 cm.
    """
    COLS = ["OPCs", "COPs", "ImOlGs", "Oligo1", "Oligo2", "Oligo3", "Oligo4", "Oligo5", "Oligo6",
            "Astro", "Micro/Mac", "Immune", "Endo", "Pericyte/SMC", "Neuron"]
    CORE = ["S1PR5", "S1PR1", "NFE2L2", "ITGA4", "CD52"]
    CMAP_EXP = LinearSegmentedColormap.from_list("w2r", ["#FFFFFF", "#FFD9B3", "#FF7F0E", "#D51F26"])
    px, pad_l, py = 0.36, 0.95, 0.34
    gw, leg_w = px * len(COLS), 1.85
    iy, head_y = 0.44, 0.22
    ph = py * len(CORE)
    W_, Hh = pad_l + gw + 0.35 + leg_w, iy + ph + 1.62
    new_page(W_, Hh)
    a = ax_cm(pad_l, iy, gw, ph)
    a.set_xlim(0, gw); a.set_ylim(0, ph); a.axis("off")
    NI = {}
    for ri, gs in enumerate(CORE):
        for ci, k in enumerate(COLS):
            row = BUB[(BUB.group == k) & (BUB.gene == gs)]
            assert len(row) == 1, (gs, k)
            p_, e = float(row.pct.iloc[0]), float(row.mean_exp.iloc[0])
            x, y = px * (ci + 0.5), py * (len(CORE) - ri - 0.5)
            if p_ <= 0:
                a.plot(x, y, marker="o", ms=1.4, mfc="none", mec="0.45", mew=0.4, zorder=3)
            else:
                a.scatter(x, y, s=4.0 + 52.0 * p_ / 100.0, c=[min(max(e, 2.2), 4.4)],
                          cmap=CMAP_EXP, vmin=2.2, vmax=4.4, edgecolors="0.25",
                          linewidths=0.25, zorder=5)
    for ri, gs in enumerate(CORE):
        txt(pad_l - 0.08, iy + py * (ri + 0.5), gs, ha="right", va="center")
    # v30: the populations that carry an axis term are printed with the axis symbol (AST / MG / LYM /
    # NEU), so the column heads of F match the delta key and Fig 7C.  The oligodendrocyte subtypes and
    # the vascular labels keep the annotation table's own names - they are a vocabulary of source
    # labels, not part of any contrast, and no abbreviation for them exists in the manuscript table.
    SHORT = {"Micro/Mac": "MG", "Immune": "LYM", "Astro": "AST", "Neuron": "NEU",
             "Pericyte/SMC": "Peri"}          # full names in the legend
    ylab = iy + ph + 0.10
    for ci, k in enumerate(COLS):
        NI[k] = int(BUB[(BUB.group == k) & (BUB.gene == "S1PR5")].cells.iloc[0])
        txt(pad_l + px * (ci + 0.5), ylab, "%s  %s" % (SHORT.get(k, k), "{:,}".format(NI[k])),
            ha="center", va="top", rotation=90)
    # the accession is a column-group head over the matrix, not a title beside the panel letter
    txt(pad_l + gw / 2.0, head_y, "GSE118257", ha="center", va="center")
    LX = pad_l + gw + 0.35
    txt(LX, head_y, "% detecting", ha="left", va="center")
    gy_c = 0.66
    for i, pv in enumerate([1, 10, 50]):
        a.scatter(LX + i * 0.62 + 0.22 - pad_l, iy + ph - gy_c, s=4.0 + 52.0 * pv / 100.0, c="0.55",
                  edgecolors="0.25", linewidths=0.25, zorder=5, clip_on=False)
        txt(LX + i * 0.62 + 0.22, gy_c + 0.20, "%d%%" % pv, ha="center", va="top")
    txt(LX, 1.22, "fill = mean", ha="left", va="center")
    txt(LX, 1.42, "log2(CP10k+1)", ha="left", va="center")
    aCB = ax_cm(LX, 1.60, 1.20, 0.13)
    aCB.imshow(np.linspace(2.2, 4.4, 256)[None, :], cmap=CMAP_EXP, vmin=2.2, vmax=4.4, aspect="auto")
    aCB.set_xticks([]); aCB.set_yticks([]); aCB.tick_params(length=0)
    for sp in aCB.spines.values():
        sp.set_linewidth(0.4)
    txt(LX, 1.77, "2.2", ha="left", va="top")
    txt(LX + 1.20, 1.77, "4.4", ha="right", va="top")
    assert sum(NI.values()) == 17799, sum(NI.values())
    export("F", "S7v30_F_detection.pdf")
    return NI


# %% PANEL G
# ================================================================ G  S1PR5 detection, per donor
def donor_s1pr5(width=6.40):
    """G: the Fig 7D detection rule re-tested one donor at a time, 27 evaluable donors in three
    cohorts.  Each dot is a donor; y is that donor's OWN within-donor contrast - the percentage of
    OL nuclei detecting S1PR5 minus the percentage in that donor's strongest OTHER population - so
    the axis is a margin, not a group difference.

    rule met (filled) = margin >= 1.0 pp AND (OL >= 1.5x the comparator OR the comparator is 0);
    open circles fail it, and all five failures are thin-comparator donors (31-410 nuclei in the
    population that beat or nearly beat OL), which is why the rule has a floor at all.  GSE277430 is
    not here: its released object is an OL-enriched subset with no non-lineage comparator, and its
    scores were built on normalised counts.  Diagnosis is deliberately not encoded (16 MS vs 11
    Ctrl, pooled Mann-Whitney P = 0.336, and the two cohorts trend in opposite directions): the
    table and the legend carry it.

    Dot area = the comparator population's nucleus count; the thick bar is the cohort median; the
    red dashed line is the +1.0 pp arm of the rule and the tinted band above it is where that arm
    holds.  Colour = cohort, matching the key band and every other panel.
    """
    H = 4.12
    new_page(width, H)
    pad_l, pad_r, gy, gb = 1.06, 0.22, 0.62, 1.30
    aw = width - pad_l - pad_r
    a = fig.add_axes([pad_l / width, gb / H, aw / width, (H - gy - gb) / H])
    a.axhspan(1.0, 16.5, color="#FDF3E4", lw=0, zorder=0)
    a.axhline(0, color="0.45", ls=":", lw=0.7, zorder=2)
    a.axhline(1.0, color=RED, ls=(0, (2.0, 1.6)), lw=0.7, zorder=3)
    a.grid(True, axis="y", lw=0.32, color="0.88", zorder=1)
    a.set_axisbelow(True)
    a.set_xlim(-0.75, 2.75); a.set_ylim(-1.9, 16.5)
    a.set_yticks([0, 4, 8, 12])
    a.tick_params(axis="y", labelsize=FS, length=1.6)
    a.set_xticks(range(3)); a.set_xticklabels([""] * 3)
    a.set_ylabel("S1PR5: OL " + MINUS + " strongest\nother population (pp)", fontsize=FS, linespacing=1.3)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
    rng = np.random.default_rng(11)                       # the jitter is fixed, so the figure is stable
    sz = 2.2 + 7.6 * np.sqrt(EV.n_binding.values) / np.sqrt(float(EV.n_binding.max()))
    jit = rng.uniform(-0.17, 0.17, len(EV))
    groups = ["GSE301908", "GSE118257", "GSE231585"]
    for gi, g in enumerate(groups):
        sub = EV[EV.cohort == g]
        assert len(sub) == int(S41.loc[g, "donors_evaluable"]), (g, len(sub))
        col = GROUPCOL[g + " RNA"]
        for _, r in sub.iterrows():
            k = int(EV.index[EV.donor == r.donor][0])
            a.scatter(gi + jit[k], r.margin_pp, s=float(sz[k]), c=col if r.rule else "white",
                      edgecolors="0.1" if not r.rule else col, linewidths=0.7,
                      zorder=6 if r.rule else 7)
        a.plot([gi - 0.21, gi + 0.21], [float(np.median(sub.margin_pp))] * 2, color=col, lw=1.6, zorder=8)
        xc = pad_l + aw * (gi + 0.5) / 3.0
        txt(xc, gy + 0.10, "%d/%d" % (int(sub.rule.sum()), len(sub)), ha="center", va="center",
            color="0.15")
        head45(xc, H - gb + 0.10, g)
    txt(pad_l + 0.03, gy + 0.42, "+1.0 pp", ha="left", va="center", color=RED)
    export("G", "S7v30_G_donor.pdf")


# ------------------------------------------------------------------ run everything
if __name__ == "__main__":
    check_palette()
    key_band()
    ladder_A()
    tight_pool_B()
    same_donor_C()
    same_nucleus_D()
    ledger_E(AVAIL)
    rank_bar()
    axis_key(AVAIL - matrix_x0("E") - FR["RB"]["w"] - 2 * GAP + 0.30)
    detection_map()
    donor_s1pr5(6.40)
    json.dump({k: {a: (round(b, 4) if isinstance(b, float) else b) for a, b in v.items()}
               for k, v in FR.items()}, open(os.path.join(OUT, "panels_v30_frames.json"), "w"), indent=1)
    say("panels:", " ".join(sorted(FR)))
