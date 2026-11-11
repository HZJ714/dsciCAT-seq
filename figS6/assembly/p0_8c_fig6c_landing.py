#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p0_8c_fig6c_landing.py - Fig6C v31: the "landing" heatmap of 38 *cis*-window risk genes x 6 OL subtypes (main-figure version)

v30 -> v31 compliance fixes (against skills/figure-house-spec/SKILL.md)
  1. §8.3/§2  all five face explanation lines deleted: reading-essential marks go into the swatch legend (§12); the rest (Δ definition, bootstrap,
     ARI, permutation P, detection-rate counts) all go into the caption (FIG6C_v31_CAPTION.md).
  2. §14.1    the 0.42 cm gap between the RNA/ATAC blocks was a direct FAIL, so the columns now touch end-to-end, separated by a thin line only;
     the freed width goes to the columns (0.55 -> 0.70 cm/column).
  3. §2        all face text is pure BLACK (the #333/#666 greys were violations); enriched subtypes no longer use gold text,
     replaced by a GOLD identity bar under the column (§14.1 permits a column-identity bar; it is a mark, not a text colour).
  4. §2        titles lose bold (face defaults to regular black text, no bold, no italics).
  5. the colour scale moved from a right-hand vertical bar to a horizontal gradient bar INSIDE the legend, eliminating overlap with the graphic; the end ticks centre exactly on the bar ends x (§14.11).

Display = the tested quantity, and the tested quantity is checkable by eye
  landing column L[g,s] = mean_log2(g,s) − mean_log2(g,COP/NFOL) (the COP/NFOL column is identically 0, narrowed as the explicit baseline)
  Δ column Δ[g]   = (1006·L[MOL3/6] + 2574·L[MOL6])/3580 − 0 = the weighted mean of the two enriched columns
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap

CM = 1 / 2.54
HERE = os.path.dirname(os.path.abspath(__file__))
# Default = three levels above the repo GitHub/figS6/assembly (the analysis-tree root); override with DSCI_WORK_ROOT.
WORK = os.environ.get("DSCI_WORK_ROOT", os.path.abspath(os.path.join(HERE, "..", "..", "..")))
BASE = WORK + "/figures/routeA_out/P0_8C"
# In the repo the notebook points P08C_TABLES/S6_OUT at tables_v16 and fig6_out; the defaults keep the analysis-tree paths.
TAB = os.environ.get("P08C_TABLES", BASE + "/tables")
PAN = os.environ.get("S6_OUT", BASE + "/panels"); os.makedirs(PAN, exist_ok=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arial_style import use_arial, use_two_tier
import fitz as _fz          # column labels need measured text widths; fitz is only imported for real in the legend section, so alias it here


def tw(t, fs=6.0):
    return _fz.get_text_length(t, fontname="helv", fontsize=fs) * 2.54 / 72.0
use_arial(); use_two_tier()

SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
ENA, ENB = "MOL3/6", "MOL6"
# MOL3/6=1006 is the full RNA AnnData used by the C1 pre-check; the 991 on the Fig 6A axis is the subset with scDRS scores (15 cells fewer)
NCELL = dict(zip(SUB, [203, 490, 1016, 688, 1006, 2574]))
LO, HI = -1.2, 1.2
DETECT_MIN = 10.0
# z display mode: same convention as p0_8c_fig6c_top1000.py (see its comments); both blocks share one colour bar / legend.
ZMODE = os.environ.get("PANEL_ZMODE", "0") == "1"
DRAW_LEG = os.environ.get("PANEL_LEGEND", "1") == "1"
ZCLIP = float(os.environ.get("PANEL_ZCLIP", 2.0))
CONTRAST = os.environ.get("PANEL_CONTRAST", "pool")   # mol6 = the main-figure convention: Δ = MOL6 − COP/NFOL
JC = 5 if CONTRAST == "mol6" else None                 # None = cell-count weighting of the two enriched columns (old convention)
HAS_D = os.environ.get("PANEL_DCOLUMN", "1") == "1"    # the rightmost Δ column; the author decided to drop it once the single-column convention was set
SHOW_NAMES = os.environ.get("PANEL_NAMES", "1") == "1" # D (the null panel) stays unnamed
NSUB = 6
GOLD = "#b8860b"          # channel/subtype identity colour, marks only, never text
BLACK = "#000000"

# ---------------------------------------------------------------- data (single source: the cell-level cached npz)
S0 = pd.read_csv(TAB + "/p08c_c1_rna_atac_by_subtype.csv")[["gene", "mouse", "top_subtype"]]
zz = np.load(TAB + "/p08c_c2_celllevel38.npz", allow_pickle=True)
XG, AG, key, cg = zz["XG"], zz["AG"], zz["key"], [str(g) for g in zz["genes"]]
rows = {s: np.flatnonzero(key == s) for s in SUB}
col = np.array([cg.index(m) for m in S0.mouse])


def summarize(M):
    Mn = np.column_stack([M[rows[s]][:, col].mean(0) for s in SUB]).astype(float)
    Pc = np.column_stack([(M[rows[s]][:, col] > 0).mean(0) * 100 for s in SUB]).astype(float)
    return Mn, Pc


Mn_r, pct_r = summarize(XG)
Mn_a, pct_a = summarize(AG)
Lr, La = Mn_r - Mn_r[:, :1], Mn_a - Mn_a[:, :1]
wA, wB = NCELL[ENA], NCELL[ENB]
if JC is None:
    d_r = (wA * Mn_r[:, 4] + wB * Mn_r[:, 5]) / (wA + wB) - Mn_r[:, 0]
    d_a = (wA * Mn_a[:, 4] + wB * Mn_a[:, 5]) / (wA + wB) - Mn_a[:, 0]
    assert np.allclose(d_r, (wA * Lr[:, 4] + wB * Lr[:, 5]) / (wA + wB), atol=1e-9)  # Δ is a linear combination of the landing columns
else:            # main-figure convention: only compared with MOL6, so cell weights are no longer needed
    d_r, d_a = Lr[:, JC].copy(), La[:, JC].copy()

rng = np.random.default_rng(20261116); br, ba = [], []
for _ in range(500):
    idx = {s: rng.choice(rows[s], len(rows[s]), replace=True) for s in SUB}
    for M, acc in ((XG, br), (AG, ba)):
        p = np.column_stack([M[idx[s]][:, col].mean(0) for s in SUB])
        acc.append((p[:, JC] - p[:, 0]) if JC is not None else
                   (wA * p[:, 4] + wB * p[:, 5]) / (wA + wB) - p[:, 0])
SE_R, SE_A = np.array(br).std(0), np.array(ba).std(0)
res_r, res_a = np.abs(d_r) > 2 * SE_R, np.abs(d_a) > 2 * SE_A

comp = d_r / d_r.std() + d_a / d_a.std()
order = np.argsort(-comp)
S = S0.iloc[order].reset_index(drop=True)
Lr, La, d_r, d_a = Lr[order], La[order], d_r[order], d_a[order]
# [BUG FIX] comp must be reordered together with order. Previously only S/L/d/res/pct were reordered, so inside pick_labels
# `key=lambda x: -abs(comp[x])` looked up the UNreordered array by row x - i.e. ranked priority by another gene's effect size.
# Consequence: pool, row order and name-to-row mapping were ALL correct (no face mislabel), but "who deserves the name" was wrong -
# 4 of the 10 names were not the ones the rule should have picked. Neither the three gates nor visual review can catch this.
comp = comp[order]
res_r, res_a, pct_r, pct_a = res_r[order], res_a[order], pct_r[order], pct_a[order]
print("sign-resolvable RNA %d/38 ATAC %d/38 | <10%%-detected landing cells RNA %d ATAC %d (of 228)"
      % (res_r.sum(), res_a.sum(), (pct_r < DETECT_MIN).sum(), (pct_a < DETECT_MIN).sum()))

# ---------------------------------------------------------------- layout
# Panel size is env-overridable (the FigS6 three-row version needs 8.38 wide); unset = the delivered main-figure size as-is.
PW = float(os.environ.get("PANEL_W", 10.72)); PH = float(os.environ.get("PANEL_H", 5.50))
NR = 38
NC = NSUB + (1 if HAS_D else 0)
LBLW = 1.14
x_lbl, x_body = 0.30, 1.54
TITLE_ON = os.environ.get("PANEL_TITLE", "1") == "1"   # as in D: title goes to the caption, the saved height to the body
YA, YG = PH - 0.27, PH - 0.55
BODY_T = (PH - 0.62) if TITLE_ON else PH - 0.34
NARROW = PW < 10.0          # FigS6 narrow version: subtype labels hang vertically (45° would cut diagonally into the gene-name column)
BODY_B = float(os.environ.get("PANEL_BODY_B", 2.05 if PW >= 10.0 else 2.60))
CW0 = float(os.environ.get("PANEL_CW0",
            float(os.environ.get("PANEL_CW", 0.70)) if (ZMODE and not HAS_D) else 0.26))
CW = float(os.environ.get("PANEL_CW", 0.70))
XOFF = [0.0] + [CW0 + i * CW for i in range(NC - 1)]
BW = XOFF[-1] + CW
x_blk = [x_body, x_body + BW]                          # the two blocks touch end-to-end, no gap (§14.1)
RH = (BODY_T - BODY_B) / NR
# With the title externalized the body moves up, so the channel names "RNA"/"ATAC" must rise above the body - otherwise they sit on
# the first row's cells (verify_text_on_data would report DOT-UNDER-TEXT 'RNA').
if not TITLE_ON:
    # With va="baseline" the box BOTTOM edge sits ~0.045 cm below the baseline (the PDF emits boxes by the font's ascender/descender,
    # taller than matplotlib's extent, see SKILL §17.2); +0.03 would be pressed by the first row; top edge = BODY_T+0.25,
    # so the body needs at least 0.34 cm headroom to the page edge, else show_pdf_page trims at the CropBox (§17.6).
    YG = BODY_T + 0.06
LEG_Y = 0.52
assert x_blk[1] + BW <= PW - 0.20, x_blk[1] + BW


def pick_labels(gap_cm=0.26):
    """Naming rule for rows (stated in the caption, never hand-picked): max detection >= 30% and RNA sign-resolvable (|Δ|>2SE),
    then greedy naming by descending |composite Δ|, named rows at least `need` rows apart (so 6pt labels never collide)."""
    need = max(1, int(np.ceil(gap_cm / RH)))
    cand = np.flatnonzero((pct_r.max(1) >= 30) & res_r)
    keep = []
    for i in sorted(cand, key=lambda x: -abs(comp[x])):
        if all(abs(i - j) >= need for j in keep):
            keep.append(i)
    return sorted(keep)


PRESPEC = [g.strip() for g in os.environ.get("PANEL_PRESPEC", "").split(",") if g.strip()]
if PRESPEC:
    """The prespecified list wins (author 2026-11-17). The largest-|Δ| genes are picked BY THE TESTED QUANTITY ITSELF,
    using them as evidence is circular; instead, only the genes that panel S6A resolved independently by locus are named -
    the names' justification lives OFF the face (GWAS lead SNP ±500 kb), orthogonal to the panel statistics.
    The spacing rule yields to the prespecified list: colliding labels are pushed up and pointed back to their row with a short leader."""
    have = set(S.gene)
    miss = [g for g in PRESPEC if g not in have]
    assert not miss, ("prespecified gene not among the panel rows", miss)
    lab_rows = sorted(int(np.flatnonzero(S.gene.values == g)[0]) for g in PRESPEC)
else:
    lab_rows = pick_labels()
print("naming %d:" % len(lab_rows), ", ".join(S.gene.iloc[lab_rows]))

fig = plt.figure(figsize=(PW / 2.54, PH / 2.54))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, PW); ax.set_ylim(0, PH); ax.axis("off")
CMAP = LinearSegmentedColormap.from_list("by", ["#2166ac", "#f7f7f7", "#b2182b"])
if TITLE_ON:
    ax.text(x_body, YA, "Risk-set landing across OL subtypes", fontsize=8.0, ha="left",
            va="baseline", color=BLACK)


def cell_y(i):
    return BODY_T - (i + 1) * RH


def rowz_disp(Ld):
    """Same function as the top1000 panel: the denominator adds a shrinkage term s0; Δ applies the same cell weights to the display.
    Row order still follows the delivered "composite Δ" (the raw quantity) - the caption states that rule; it does not change with the display."""
    # A gene undetected everywhere is NaN for the whole row (greyed below by pct<DETECT_MIN), so centering/spread/shrink
    # must all use the nan* versions: one NaN would drag np.median of s0 to NaN and whiten the whole heatmap.
    mu = np.nanmean(Ld, 1, keepdims=True); sd = np.nanstd(Ld, 1, ddof=0, keepdims=True)
    s0 = float(np.nanmedian(sd))
    assert np.isfinite(s0) and s0 > 0, ("shrinkage s0 not finite", s0)
    Z = (Ld - mu) / np.sqrt(sd ** 2 + s0 ** 2)
    fz = np.array([NCELL[s] for s in SUB], float); fz = fz / fz.sum()
    return Z, (Z[:, JC] if JC is not None else fz[4] * Z[:, 4] + fz[5] * Z[:, 5]), s0


if ZMODE:
    Lr, d_r, _s0r = rowz_disp(Lr)
    La, d_a, _s0a = rowz_disp(La)
    LO, HI = -ZCLIP, ZCLIP
    print("[zmode] RNA s0 = %.4g | ATAC s0 = %.4g" % (_s0r, _s0a), flush=True)

_placed = []
LBL_H = 0.22   # the real 6pt gene-name glyph height (ascender+descender): the push-up collision step
for b, (Mat, x0, dlt, res, pct, nm) in enumerate(((Lr, x_blk[0], d_r, res_r, pct_r, "RNA"),
                                                  (La, x_blk[1], d_a, res_a, pct_a, "ATAC"))):
    for i in range(NR):
        y, h = cell_y(i), RH - 0.004
        for j in range(NSUB):
            v, wx, ww = Mat[i, j], x0 + XOFF[j], (CW0 if j == 0 else CW)
            if pct[i, j] < DETECT_MIN:
                ax.add_patch(Rectangle((wx, y), ww, h, fc="#d9d9d9", ec="none"))
            else:
                t = max(min(v, HI), LO)
                ax.add_patch(Rectangle((wx, y), ww, h, fc=CMAP((t - LO) / (HI - LO)), ec="none"))
                if abs(v) > HI:   # the reference column (j=0) is marked too: under z it is no longer identically 0 and all off-scale cells sit there (the old j>0 was a leftover that swallowed every mark)
                    ax.plot(wx + (ww - 0.07 if v > 0 else 0.07), y + h / 2, marker=">" if v > 0 else "<",
                            ms=1.5, color="#ffffff", mew=0.4, zorder=6)
        if not HAS_D:
            continue
        wx = x0 + XOFF[NC - 1]
        if res[i]:
            t = max(min(dlt[i], HI), LO)
            ax.add_patch(Rectangle((wx, y), CW, h, fc=CMAP((t - LO) / (HI - LO)), ec="none"))
            if abs(dlt[i]) > HI:
                ax.plot(wx + (CW - 0.07 if dlt[i] > 0 else 0.07), y + h / 2,
                        marker=">" if dlt[i] > 0 else "<", ms=1.5, color="#ffffff", mew=0.4, zorder=6)
        else:
            ax.add_patch(Rectangle((wx, y), CW, h, fc="#d9d9d9", ec="none"))
    _n = int(((np.abs(Mat) > HI) & (pct >= DETECT_MIN)).sum())      # this panel's display is called Mat, not L
    print("  %s block off-scale triangles %d (of which reference column COP/NFOL %d) | |z|>%.1f and not greyed"
          % (nm, _n, int(((np.abs(Mat[:, 0]) > HI) & (pct[:, 0] >= DETECT_MIN)).sum()), HI), flush=True)
    ax.add_patch(Rectangle((x0, BODY_B), BW, BODY_T - BODY_B, fc="none", ec="#333333", lw=0.4))
    if HAS_D:
        ax.add_patch(Rectangle((x0 + XOFF[NC - 1], BODY_B), CW, BODY_T - BODY_B, fc="none", ec=BLACK, lw=0.7))
    ax.text(x0 + BW / 2, YG, nm, fontsize=6.0, ha="center", va="baseline", color=BLACK)
    for j, s in enumerate(SUB + (["\u0394"] if HAS_D else [])):
        ww = CW0 if j == 0 else CW
        if NARROW:
            # The MEASURED semantics of rotation=90 + ha=right + va=center: the anchor lands at the text's vertical centre,
            # horizontal right edge. To hang it top-aligned under the body, drop the anchor by L/2 per label's measured length.
            ax.text(x0 + XOFF[j] + ww / 2 + 0.115, BODY_B - 0.12 - tw(s) / 2.0, s,
                    fontsize=6.0, rotation=90, ha="right", va="center", color=BLACK)
        else:
            ax.text(x0 + XOFF[j] + ww / 2, BODY_B - 0.07, s, fontsize=6.0, rotation=45,
                    ha="right", va="top", color=BLACK)
    # the gold identity bar under the enriched subtype (a mark, not a text colour)
    for k in ((JC,) if JC is not None else (4, 5)):
        ax.add_patch(Rectangle((x0 + XOFF[k] + 0.06, BODY_B - 0.055), CW - 0.12, 0.035, fc=GOLD, ec="none"))
    if not HAS_D:   # with no Δ column, "which is the numerator and which the denominator" rests entirely on these two identity bars
        ax.add_patch(Rectangle((x0 + XOFF[0] + 0.06, BODY_B - 0.055), CW0 - 0.12, 0.035, fc=BLACK, ec="none"))

for i in sorted(lab_rows, reverse=True) if SHOW_NAMES else []:   # bottom-up
    yc = cell_y(i) + (RH - 0.004) / 2
    y = max(yc, _placed[-1] + LBL_H) if _placed else yc
    _placed.append(y)
    ax.text(x_lbl + LBLW - 0.04, y, S.gene[i], fontsize=6.0, ha="right", va="center", color=BLACK)
    if abs(y - yc) > 0.015:                                  # leader: name's right end -> its own row's left edge
        ax.plot([x_lbl + LBLW - 0.01, x_body - 0.05], [y, yc], lw=0.3, color="#333333", zorder=5)

# ---------------------------------------------------------------- swatch legend (§12: marks and their text share one coordinate system)
import fitz
def tw(t, fs=6.0):
    return fitz.get_text_length(t, fontname="helv", fontsize=fs) * 2.54 / 72.0


def legend(x0, x1, y, items, gap=0.42):
    """items = [(kind, text)]. With kind='grad' the text sits ABOVE the gradient bar (not consuming horizontal width),
    the end-tick text centres exactly on the bar-end x (§14.11). Too wide fails hard; never solved by shrinking type."""
    MKW = {"grad": 1.45, "sw": 0.26, "car": 0.30, "bar": 0.42}
    wid = [MKW[k] + (0.0 if k == "grad" else 0.06 + tw(t)) for k, t in items]
    tot = sum(wid) + gap * (len(items) - 1)
    if x0 + tot > x1:
        raise SystemExit("legend too wide %.2f > %.2f - please shorten the wording" % (x0 + tot, x1))
    x = x0 + (x1 - x0 - tot) / 2.0
    for (kind, t), w in zip(items, wid):
        if kind == "grad":
            ax.imshow(np.linspace(HI, LO, 64)[None, :], extent=(x, x + MKW[kind], y - 0.075, y + 0.075),
                      cmap=CMAP, aspect="auto", zorder=3)
            ax.add_patch(Rectangle((x, y - 0.075), MKW[kind], 0.15, fc="none", ec="#333333", lw=0.3))
            for px, lab in ((x, "%+.1f" % HI), (x + MKW[kind] / 2, "0"),
                            (x + MKW[kind], "\u2212%.1f" % abs(LO))):
                ax.text(px, y - 0.115, lab, fontsize=6.0, ha="center", va="top", color=BLACK)
            ax.text(x + MKW[kind] / 2, y + 0.115, t, fontsize=6.0, ha="center", va="bottom", color=BLACK)
        elif kind == "sw":
            ax.add_patch(Rectangle((x, y - 0.075), MKW[kind], 0.15, fc="#d9d9d9", ec="#8c8c8c", lw=0.3))
            ax.text(x + MKW[kind] + 0.06, y, t, fontsize=6.0, ha="left", va="center", color=BLACK)
        elif kind == "car":
            ax.plot(x + 0.07, y, marker="<", ms=2.2, color="#2166ac", mew=0.3)
            ax.plot(x + 0.23, y, marker=">", ms=2.2, color="#b2182b", mew=0.3)
            ax.text(x + MKW[kind] + 0.06, y, t, fontsize=6.0, ha="left", va="center", color=BLACK)
        else:
            ax.add_patch(Rectangle((x, y - 0.055), MKW[kind], 0.075, fc=GOLD, ec="none"))
            ax.text(x + MKW[kind] + 0.06, y, t, fontsize=6.0, ha="left", va="center", color=BLACK)
        x += w + gap
    return x - gap


if DRAW_LEG:
  _it = [("grad", "z" if ZMODE else "log2 vs COP/NFOL"), ("sw", "not detected"),
         ("car", "off scale"), ("bar", "enriched subtypes")]
  if PW < 10.0:      # narrow version (FigS6 8.38 cm): the legend wraps to two rows; wording kept, type not shrunk (spec §3.3 breathes first)
      # when the two gradient bars sit side by side, the end ticks "−1.2" and "+0.35" read as one string (gates silent, eyes misread) -> widen the gap
      xe = legend(x_body, x_blk[1] + BW, LEG_Y, _it[2:], gap=0.90)
      # each gradient bar's end ticks centre on its bar end: if the gap is smaller than the sum of the half-widths of the two labels it reads "−1.2+0.35" (gates silent, eyes misread)
      xe = legend(x_body, x_blk[1] + BW, LEG_Y + 0.44, _it[:2], gap=1.10)
  else:
      xe = legend(x_body, x_blk[1] + BW, LEG_Y, _it)
  print("legend occupies %.2f-%.2f (available %.2f-%.2f)" % (x_body, xe, x_body, x_blk[1] + BW))
else:
    xe = x_body
    print("legend handed to the shared band (PANEL_LEGEND=0)", flush=True)

out = PAN + "/fig6v31_C_landing%s.pdf" % os.environ.get("PANEL_TAG", "")
fig.savefig(out, dpi=600)
C = 2.54 / 72
d = fitz.open(out); p = d[0]
W, H = p.rect.width * C, p.rect.height * C
ws = p.get_text("words")
worst = max((w[2] * C for w in ws), default=0); bot = max((w[3] * C for w in ws), default=0)
sp = [s for b in p.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]]
grey = sorted({"%06x" % s["color"] for s in sp if s["color"] not in (0, 0xffffff)})
print("size %.2f x %.2f | text right %.2f (limit %.2f) | bottom %.2f | non-black text %s" % (W, H, worst, W - 0.20, bot, grey))
d.close()
