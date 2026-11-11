# -*- coding: utf-8 -*-
"""Fig6 final: A | B | D in one row (17.00 x ~5.9 cm), two-tier type, zero vector scaling
=================================================================
Per FIG6_S6_LAYOUT_FINAL_20261117.md §3:
  - the gene-level C is pulled from the main figure (no defensible positive) -> only A (cell-level score) / B (UMAP localisation) / D (the 74-trait dual-omics ladder) remain;
  - hard page-width gate 17.00 cm: A 5.80 + B 5.20 + D 5.50 + 2 x 0.24 gaps = 16.98;
  - all three panels are RE-RENDERED from scratch (§3.3 forbids vector scaling - type would drift);
  - two-tier type: 8.0 (letters / axis titles / row heads) and 6.0 (ticks / corner notes / subtype names / statistics), no 5.5/7.0 third tier;
  - B's UMAP is RE-PROJECTED onto the new frame (set_aspect locks the shape; the scatter is not stretched); B's face title is deleted in the final,
    (house spec: no title next to the panel letter; the content goes into the caption);
  - A keeps its four-line statistics table (final content), with 5.5/5.6 pt raised to 6.0 pt.

Usage: DSCI_WORK_ROOT=<work> python3 fig6v36_abd_final.py
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import blended_transform_factory as btf
import anndata as ad
from scipy.stats import mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("DSCI_WORK_ROOT", os.path.abspath(os.path.join(HERE, "..", "..", "..")))
sys.path.insert(0, HERE)
sys.path.insert(0, WORK + "/GitHub/figS6/assembly")
from arial_style import use_arial, FONTS_PDF
use_arial({"axes.linewidth": 0.6, "savefig.dpi": 600})
import fitz
from arial_style import patch_fitz_arial
patch_fitz_arial(fitz)

CM = 2.54
PT = 72 / 2.54
OUT = os.environ.get("FIG6_OUT", os.path.join(WORK, "GitHub/fig6/fig6_out"))
DST = os.environ.get("FIG6_PANELS", os.path.join(WORK, "figures/routeA_out/P0_8C/panels"))
FY = os.environ.get("DSCI_FY_ROOT", os.path.join(WORK, "FY"))
# Per-nucleus table cache: A/B only need the four columns cell_type + UMAP_1/2 + MS norm_score.
# Read the cache when present (the entry point for the repo notebook running alone); otherwise recompute from h5ad + .score.gz and write back.
CACHE = os.environ.get("FIG6_CACHE", os.path.join(OUT, "caches"))
os.makedirs(CACHE, exist_ok=True)
H5_R = os.path.join(FY, "oligo", "3k_dsci_oligo_subset.h5ad")
MS_SCORE = os.path.join(FY, "oligo", "scdrs_ms_out", "PASS_Multiple_sclerosis.score.gz")
GS_GROUP = os.path.join(FY, "oligo", "scdrs_ms_out", "PASS_Multiple_sclerosis.scdrs_group.cell_type")
SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
STALLION = {"COP/NFOL": "#D51F26", "MFOL": "#272E6A", "MOL2": "#208A42",
            "MOL3": "#89288F", "MOL3/6": "#F47D2B", "MOL6": "#FEE500"}
FS_LAB, FS_TITLE = 6.0, 8.0                     # the two-tier policy (§8)
WB = dict(fc="white", ec="none", alpha=0.9, pad=0.2)


def cliffs(x, y):
    x, y = np.asarray(x), np.asarray(y)
    return float((x[:, None] > y[None, :]).mean() - (x[:, None] < y[None, :]).mean())


def bh(p):
    p = np.asarray(p, float); q = np.full(len(p), np.nan)
    idx = np.where(~np.isnan(p))[0]; m = len(idx)
    if m == 0:
        return q
    o = idx[np.argsort(p[idx])]
    qs = p[o] * m / np.arange(1, m + 1)
    q[o] = np.minimum(np.minimum.accumulate(qs[::-1])[::-1], 1.0)
    return q


def fmtp(p):
    return "%.3f" % p if p >= 0.001 else "%.1fe%d" % (p / 10 ** int(np.floor(np.log10(p))),
                                                      int(np.floor(np.log10(p))))


def dark(c, thr=0.55):
    rgb = np.array(mplc.to_rgb(c)); lum = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
    return tuple(rgb * (thr / lum)) if lum > thr else tuple(rgb)


# ------------------------------------------------------------------ data (the same chain as the final)
_cc = os.path.join(CACHE, "fig6_cellscore.csv")
if os.path.exists(_cc):
    df = pd.read_csv(_cc)
    print("[cache] per-nucleus score table %s (%d nuclei)" % (os.path.basename(_cc), len(df)))
else:
    obs = ad.read_h5ad(H5_R, backed="r").obs
    ms = pd.read_csv(MS_SCORE, sep="\t", index_col=0)["norm_score"]
    df = obs[["cell_type", "UMAP_1", "UMAP_2"]].join(ms.to_frame("norm_score")).dropna()
    df["cell_type"] = df.cell_type.astype(str)
    df.to_csv(_cc, index=False)
grp = pd.read_csv(os.environ.get("FIG6_GROUP", GS_GROUP), sep=None, engine="python", index_col="group")
ref = df.loc[df.cell_type == "COP/NFOL", "norm_score"].values
rows = []
for ct in SUB:
    v = df.loc[df.cell_type == ct, "norm_score"].values
    rows.append(dict(cell_type=ct, n=len(v), mean=v.mean(), median=float(np.median(v)),
                     sem=v.std(ddof=1) / np.sqrt(len(v)),
                     mwu_p_vs_COPNFOL=np.nan if ct == "COP/NFOL"
                     else float(mannwhitneyu(v, ref, alternative="two-sided").pvalue),
                     group_P=float(grp.loc[ct, "assoc_mcp"]),
                     cliffs_delta=0.0 if ct == "COP/NFOL" else cliffs(v, ref)))
A = pd.DataFrame(rows)
A["group_q_BH5"] = bh(A.group_P.where(A.cell_type != "COP/NFOL").values)
OLD = pd.read_csv(os.path.join(OUT, "fig6A_stats.csv"))
assert np.allclose(A["n"].values, OLD["n"].values) and np.allclose(A["median"].values, OLD["median"].values) \
    and np.allclose(A["cliffs_delta"].values, OLD["cliffs_delta"].values, atol=1e-6), "statistics must equal the delivered fig6A_stats.csv"
print("nuclei %d | consistent with the delivered fig6A_stats.csv" % len(df))
data = [df.loc[df.cell_type == ct, "norm_score"].values for ct in SUB]


# ------------------------------------------------------------------ A: box plot + the four-line statistics table
def render_A(w, h, path):
    fig, ax = plt.subplots(figsize=(w / CM, h / CM))
    fig.subplots_adjust(left=0.175, right=0.985, top=0.975, bottom=0.255)
    bp = ax.boxplot(data, showfliers=False, whis=[5, 95], patch_artist=True, widths=0.55,
                    medianprops=dict(color="black", lw=1.0), whiskerprops=dict(lw=0.7),
                    capprops=dict(lw=0.7))
    for i, (patch, ct) in enumerate(zip(bp["boxes"], SUB)):
        patch.set_facecolor(STALLION[ct]); patch.set_alpha(0.55); patch.set_linewidth(0.6)
        patch.set_edgecolor(dark(STALLION[ct], 0.45))
        xj = np.random.RandomState(7).normal(i + 1, 0.07, len(data[i]))
        ax.scatter(xj, data[i], s=1.5, color=dark(STALLION[ct]), alpha=0.35, lw=0, zorder=3)
    ax.axhline(0, color="0.55", ls="--", lw=0.7)
    # The statistics table must sit ABOVE ALL jittered points: the delivered v21 started at the whisker end (wmax), letting points beyond the 95th percentile poke into
    # the table area (verify_text_on_data reported 18 DOT-UNDER-TEXT on the delivered page); here it starts at the data's true upper bound.
    wmax = max(l.get_ydata().max() for l in bp["whiskers"])
    dmax = max(float(np.max(d)) for d in data)
    top = max(wmax, dmax)
    # The small table keeps only ONE line, MWU P (author 2026-11-17): the four-line table pushed ylim up an extra 1.88 score units
    # (~17% of the axis height) - the very reason A/D grew to 5.9 cm and the aspect-locked B could not fill the row.
    # group P / group q / Cliff delta all move off the face into the caption (the numbers still live in fig6A_stats.csv; the caption is ledger-generated).
    R1 = top + 0.16
    for i, ct in enumerate(SUB):
        x = i + 1
        if ct == "COP/NFOL":
            # Row labels RIGHT-align at 1.32 (0.23 cm right of the first column's centre): centred, "group q" (0.75 cm) would end
            # 0.10 cm from MFOL's number, reading as "group q0.452"; right-aligned, the row labels and the number column sit 0.42 cm apart.
            for y, lab in ((R1, "MWU P"),):
                ax.text(1.32, y, lab, ha="right", va="bottom", fontsize=FS_LAB, color="0.15",
                        bbox=WB, zorder=8)
            continue
        r = A.loc[A.cell_type == ct].iloc[0]
        for y, txt, col in ((R1, fmtp(r.mwu_p_vs_COPNFOL), "black"),):
            ax.text(x, y, txt, ha="center", va="bottom", fontsize=FS_LAB, color=col, bbox=WB, zorder=8)
        v = data[i]
        ax.errorbar(x, v.mean(), yerr=v.std(ddof=1) / np.sqrt(len(v)), fmt="o", ms=2.4, mfc="white",
                    mec="black", mew=0.6, capsize=1.5, capthick=0.6, lw=0.6, ecolor="black", zorder=6)
    ax.set_xticks(range(1, len(SUB) + 1))
    # The x axis carries subtype names only (author 2026-11-17): n lives in the text/caption already, repeating it on the face is redundant;
    # the labels shrink from "COP/NFOL  n=203" (1.70 cm) to "COP/NFOL" (0.93 cm), returning the saved bottom margin to the data area.
    ax.set_xticklabels(list(SUB), fontsize=FS_LAB, rotation=90, ha="center")
    for lbl, ct in zip(ax.get_xticklabels(), SUB):
        lbl.set_color(dark(STALLION[ct]))
    # 0.62 score units (~0.3 cm) of breathing above the table row: at 0.34, "MWU P"'s left edge sat 0.17 cm from the panel letter,
    # no gate fires but it reads as cursive. The row height comes from B's measurement; the extra height here merely compresses the data area 3%, no type shrink.
    ax.set_ylim(min(-2.9, df.norm_score.min() - 0.1), R1 + 0.62)
    ax.set_yticks([t for t in [-2, -1, 0, 1, 2, 3] if t < R1 + 0.62])
    ax.set_ylabel("scDRS norm_score\n(MS GWAS)", linespacing=1.35, fontsize=FS_TITLE)
    ax.tick_params(labelsize=FS_LAB, length=2.0, width=0.5, pad=1.2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.02); plt.close(fig)


# ------------------------------------------------------------------ B: UMAP (aspect-locked, unstretched)
BLUERED = LinearSegmentedColormap.from_list("bluered", ["#9ECAE1", "#FFFFFF", "#E60000"])


def render_B(w, h, path):
    """Aspect-locked UMAP + the colour bar HORIZONTAL below the map.
    Why: the UMAP's 1.63 length/width ratio is hard geometry - in a three-column row, B's column width cannot stretch to the 5.9 cm row height,
    and a vertical colour bar splits that height difference into two bands of blank (ugly). Horizontal, the bar's ink bottom lands at the column floor,
    the blank merges into one bottom band, and the hexbin stays isotropic, never stretched by the frame."""
    fig = plt.figure(figsize=(w / CM, h / CM))
    LP, RP, TP = 0.84, 0.05, 0.05          # left: UMAP_2 label+ticks; right; top (cm)
    XLAB, BARGAP, BARH, CTICK = 0.40, 0.16, 0.17, 0.34   # x-label band / bar gap / bar height / bar ticks
    x0, x1 = np.percentile(df["UMAP_1"], [0.2, 99.8])
    y0, y1 = np.percentile(df["UMAP_2"], [0.2, 99.8])
    ASP = (x1 - x0) / (y1 - y0)                   # the isotropy ratio given by the data itself
    bw = w - LP - RP
    bh = bw / ASP                                 # frame height derived backwards so nothing stretches
    ax = fig.add_axes([LP / w, (h - TP - bh) / h, bw / w, bh / h])
    hb = ax.hexbin(df["UMAP_1"], df["UMAP_2"], C=df["norm_score"], reduce_C_function=np.mean,
                   gridsize=25, mincnt=5, cmap=BLUERED, linewidths=0, edgecolors="none")
    hb.set_rasterized(True)
    vals = hb.get_array()
    v = max(np.nanpercentile(vals, 99), abs(np.nanpercentile(vals, 1)))
    hb.set_clim(-v, v)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_aspect("equal")                        # re-projection: isotropic scatter, never stretched by the frame
    xl, yl = ax.get_xlim(), ax.get_ylim()
    for nm, du, dv, ha, col in [("COP/NFOL", 0.06, 0.10, "left", "#084594"),
                                 ("MOL6", -0.04, 0.06, "right", "#B30000"),
                                 ("MOL2", -0.06, -0.08, "right", "#B30000")]:
        s = df[df.cell_type == nm]
        px, py = s.UMAP_1.median(), s.UMAP_2.median()
        u, vv = (px - xl[0]) / np.diff(xl), (py - yl[0]) / np.diff(yl)
        ax.annotate(nm, xy=(px, py), xytext=(u + du, vv + dv), xycoords="data",
                    textcoords="axes fraction", ha=ha, va="center", fontsize=FS_LAB,
                    color=col, bbox=WB,
                    arrowprops=dict(arrowstyle="->", color=col, lw=0.8, shrinkA=1, shrinkB=1.5))
    ax.set_xlabel("UMAP_1", fontsize=FS_TITLE, labelpad=1.5)
    ax.set_ylabel("UMAP_2", fontsize=FS_TITLE, labelpad=1.5)
    # horizontal colour bar: same width as the UMAP frame, tucked under the x labels
    byc = (h - TP - bh - XLAB - BARGAP - BARH) / h
    axc = fig.add_axes([LP / w, byc, bw / w, BARH / h])
    cb = fig.colorbar(hb, cax=axc, orientation="horizontal")
    cb.outline.set_linewidth(0.5)
    cb.set_ticks([round(-v, 2), 0, round(v, 2)])
    cb.ax.tick_params(labelsize=FS_LAB, length=2, width=0.5, pad=1.2, bottom=True,
                      top=False, left=False, right=False)
    cb.set_label("norm_score", fontsize=FS_LAB, labelpad=2.0)   # horizontal bar: the label sits at the left end
    fig.savefig(path, bbox_inches="tight", pad_inches=0.02); plt.close(fig)


# ------------------------------------------------------------------ C/D: the 74-trait dual-omics ladder
LT = pd.read_csv(OUT + "/fig6C_delta1_own74.csv")
ST = pd.read_csv(OUT + "/fig6C_stats.csv").iloc[0]
RHO = float(ST.rho_delta1); P_ADD = (1 + 5) / 10001
assert abs(RHO - 0.819) < 1e-3 and abs(P_ADD - 6.0e-4) < 1e-5
MARK = {"PASS_Multiple_sclerosis": "MS", "PASS_Schizophrenia_Pardinas2018": "SCZ",
        "UKB_460K.body_HEIGHTz": "HEIGHT"}
RED, BLU, MUT, GRID = "#D51F26", "#272E6A", "0.42", "0.85"


def render_C(w, h, path):
    fig, axc = plt.subplots(figsize=(w / CM, h / CM))
    fig.subplots_adjust(left=0.175, right=0.985, top=0.965, bottom=0.145)
    axc.grid(True, lw=0.32, color=GRID, zorder=0); axc.set_axisbelow(True)
    x, y = LT.d1_atac.values, LT.d1_rna.values
    lo, hi = min(x.min(), y.min()) * 1.12 - 0.06, max(x.max(), y.max()) * 1.12 + 0.06
    axc.scatter(x, y, s=8, facecolor="0.75", edgecolor="none", zorder=2)
    axc.plot([lo, hi], [lo, hi], ls="--", lw=0.5, color=MUT, zorder=1)
    axc.axhline(0, color=GRID, lw=0.5, zorder=1); axc.axvline(0, color=GRID, lw=0.5, zorder=1)
    axc.set_xlim(lo, hi); axc.set_ylim(lo, hi)
    # The ladder does not lock 1:1 (consistent with the delivered v23): both axes share units, but this is a "sides of a reference line" plot,
    # points are circular symbols, immune to stretching; the shape that must stay locked is B's UMAP embedding geometry.
    # the dashed y=x is not 45° on the page; that goes into the caption rather than faking coincidence with an aspect lock.
    # the x title wraps to two lines: one "Gene score (MOL6 − COP/NFOL)" line measures 4.4 cm, wider than the axis body (4.1 cm),
    # and the axis centre sits at 58% of the page width, so its right end would cross the CropBox by 0.13 cm and show_pdf_page would clip it (word-level gates cannot see it,
    # as they only walk the on-page words). Matching the ylabel's two-line wrap is the stable choice.
    axc.set_xlabel("Gene score\n(MOL6 \u2212 COP/NFOL)", fontsize=FS_TITLE, labelpad=2.0,
                   linespacing=1.35)
    axc.set_ylabel("Gene expression\n(MOL6 \u2212 COP/NFOL)", fontsize=FS_TITLE, linespacing=1.35,
                   labelpad=1.5)
    axc.tick_params(labelsize=FS_LAB, length=2.0, width=0.5, pad=1.2)
    for sp in axc.spines.values():
        sp.set_linewidth(0.5)
    fig.canvas.draw(); rnd = fig.canvas.get_renderer()
    ar = axc.get_window_extent(); px_ = axc.transData.transform(np.c_[x, y])
    pts = [(LT.trait.iloc[i], px_[i, 0], px_[i, 1]) for i in range(len(px_))]

    # box probing uses one resident Text + get_window_extent(renderer); the earlier approach drew fig.canvas per candidate -
    # on a 0.01-step grid that was 24k full-figure redraws (>8 min); this is O(1) text measurement instead.
    _probe = axc.text(0, 0, "", transform=axc.transAxes, fontsize=FS_LAB, linespacing=1.25,
                      va="top", ha="left")
    fig.canvas.draw()

    def box_of(txt, fx, fy):
        _probe.set_text(txt); _probe.set_position((fx, fy))
        return mpl.transforms.Bbox(_probe.get_window_extent(rnd).get_points())
    cor = [axc.text(0.975, 0.105, "P = %.0e" % P_ADD, transform=axc.transAxes, fontsize=FS_LAB,
                    color="0.15", ha="right", va="baseline"),
           axc.text(0.975, 0.185, "\u03c1 = %.3f" % RHO, transform=axc.transAxes, fontsize=FS_LAB,
                    color="0.15", ha="right", va="baseline")]
    reserved = [t.get_window_extent(rnd) for t in cor]
    # the single-row frame is smaller, so the leader budget and search step widen accordingly (§3.3: crowding is solved by breathing room, never type shrink)
    PAD, MAX_LEAD = 2.5, 150.0
    LEADRET = 7.0   # leader endpoint retreats along the outward normal/diagonal (display px; at dpi=100, 7 px = 5 pt = 0.18 cm)

    def seg_hits_rect(a, b, bb, inset=0.0):
        """Liang-Barsky: does segment ab enter the interior of rectangle bb (inset pulls the box in, so a corner-touching line only touching the boundary does not count)."""
        x0, y0 = a[0], a[1]; dx, dy = b[0] - a[0], b[1] - a[1]
        lx, ly, rx, ry = bb.x0 + inset, bb.y0 + inset, bb.x1 - inset, bb.y1 - inset
        u1, u2 = 0.0, 1.0
        for pk, qk in ((-dx, x0 - lx), (dx, rx - x0), (-dy, y0 - ly), (dy, ry - y0)):
            if pk == 0:
                if qk < 0:
                    return False
            else:
                u = qk / pk
                if pk < 0:
                    u1 = max(u1, u)
                else:
                    u2 = min(u2, u)
                if u1 > u2:
                    return False
        return True
    # the y=x reference line (display coordinates): label boxes must not sit on it. SCZ's points are on the line to begin with,
    # so a search that only checks "not on points / not on other text" would place the label on the line; the gate flags LINE-ACROSS-TEXT at 7%
    # (it crosses the dashed line, not a leader - for two rounds I misread it as the leader).
    REF0, REF1 = axc.transData.transform((lo, lo)), axc.transData.transform((hi, hi))
    grid = [(fx, fy) for fy in np.arange(0.98, 0.05, -0.01) for fx in np.arange(0.02, 0.88, 0.01)]
    # each anchor keeps its own v23 anchor position (transAxes). The earlier scoring pulled all three labels toward MS's (0.43, 0.9),
    # making SCZ's leader span half the figure and flinging HEIGHT's to the top-right - not "nearest", a wrong search target.
    ANCHOR0 = {"Multiple_sclerosis": (0.430, 0.960), "HEIGHTz": (0.660, 0.720),
               "Schizophrenia": (0.600, 0.570)}
    for key, col in [("Multiple_sclerosis", RED), ("Schizophrenia", BLU), ("HEIGHTz", MUT)]:
        r = LT[LT.trait.str.contains(key, regex=False)].iloc[0]
        nm = MARK[r.trait]
        axc.scatter([r.d1_atac], [r.d1_rna], s=20, facecolor=col, edgecolor="white", lw=0.4, zorder=4)
        ptx, pty = axc.transData.transform((r.d1_atac, r.d1_rna))
        txt = "%s\nRNA #%d  ATAC #%d" % (nm, r.rank_rna_raw, r.rank_atac_raw)
        best = None
        for fx, fy in grid:
            bx = box_of(txt, fx, fy)
            if (bx.x0 < ar.x0 + 0.5 or bx.x1 > ar.x1 - 0.5 or bx.y0 < ar.y0 + 0.5 or bx.y1 > ar.y1 - 0.5):
                continue
            if any(tr != r.trait and bx.x0 - PAD <= ppx <= bx.x1 + PAD
                   and bx.y0 - PAD <= ppy <= bx.y1 + PAD for tr, ppx, ppy in pts):
                continue
            if any(bx.overlaps(rb) for rb in reserved):
                continue
            if seg_hits_rect(REF0, REF1, bx, inset=1.0):
                continue
            # try ALL FOUR corners, take the feasible shortest. Trying only the "nearest corner" fails for cases like HEIGHT where point and box are nearly the same height:
            # the nearest corner is top-right, and after the retreat the leader still grazes parallel to the box's top edge, so every inset test rejects it.
            ccx, ccy = (bx.x0 + bx.x1) / 2.0, (bx.y0 + bx.y1) / 2.0
            # attachment candidates = 4 corners + the perpendicular projection onto 4 edges. Corners only, SCZ's "point just above-left of the box"
            # geometry draws a line grazing parallel to the box's top edge (gate measured 7%, the leader itself crossing);
            # an edge projection gives a short stub perpendicular to the edge, inherently non-grazing.
            px_ = min(max(ptx, bx.x0), bx.x1); py_ = min(max(pty, bx.y0), bx.y1)
            att = [((bx.x0, bx.y1), (bx.x0 - ccx, bx.y1 - ccy)), ((bx.x1, bx.y1), (bx.x1 - ccx, bx.y1 - ccy)),
                   ((bx.x0, bx.y0), (bx.x0 - ccx, bx.y0 - ccy)), ((bx.x1, bx.y0), (bx.x1 - ccx, bx.y0 - ccy)),
                   ((px_, bx.y1), (0.0, 1.0)), ((px_, bx.y0), (0.0, -1.0)),
                   ((bx.x1, py_), (1.0, 0.0)), ((bx.x0, py_), (-1.0, 0.0))]
            cnd = None
            for (ax_, ay_), (ox, oy) in att:
                ld = np.hypot(ax_ - ptx, ay_ - pty)
                if ld > MAX_LEAD:
                    continue
                _n0 = max(np.hypot(ox, oy), 1e-6)
                ex, ey = ax_ + LEADRET * ox / _n0, ay_ + LEADRET * oy / _n0
                # the self-check box must be OUTSET by 3.5 px: matplotlib's text extent measures its own ascent,
                # while the PDF word bbox follows the embedded font's ascender/descender, ~0.04 cm taller, measured;
                # add the 0.5 pt outline and round caps - without the outset, "geometry clean, gate dirty" (SCZ lost three rounds).
                if seg_hits_rect((ptx, pty), (ex, ey), bx, inset=-3.5):
                    continue
                if any(seg_hits_rect((ptx, pty), (ex, ey), rb, inset=1.0) for rb in reserved):
                    continue
                if cnd is None or ld < cnd[0]:
                    cnd = (ld, ex, ey)
            if cnd is None:
                continue
            lead = cnd[0]
            ax0, ay0 = ANCHOR0[key]
            score = lead + 0.6 * np.hypot((fx - ax0) * ar.width, (fy - ay0) * ar.height)
            if best is None or score < best[0]:
                best = (score, fx, fy, bx, cnd)
        assert best is not None, "anchor %s does not fit the new frame (in-frame / off-points / leader <=%.0fpx all impossible)" % (nm, MAX_LEAD)
        _, fx, fy, bx, cnd = best
        lead, ex, ey = cnd
        dx, dy = axc.transData.inverted().transform((ex, ey))
        axc.plot([r.d1_atac, dx], [r.d1_rna, dy], lw=0.5, color=col, zorder=5,
                 solid_capstyle="round")
        reserved.append(bx)
        axc.text(fx, fy, txt, transform=axc.transAxes, fontsize=FS_LAB, color=col,
                 linespacing=1.25, va="top", ha="left", zorder=6)
        reserved.append(bx)
    _probe.set_visible(False)      # the probe only measures; left visible it ghosts at the last candidate position
    fig.savefig(path, bbox_inches="tight", pad_inches=0.02); plt.close(fig)


# ------------------------------------------------------------------ iterate to the target frames + compose
def size_of(p):
    d = fitz.open(p); r = d[0].rect; d.close(); return r.width * CM / 72, r.height * CM / 72


def fit(fn, wtarget, htarget, name, free_h=False):
    """free_h=True: the height is set by the in-figure aspect lock (B's UMAP data ratio 1.63); only the width converges."""
    path = os.path.join(DST, name); w, h = wtarget, htarget
    for it in range(5):
        fn(w, h, path)
        aw, ah = size_of(path)
        if abs(aw - wtarget) <= 0.04 and (free_h or abs(ah - htarget) <= 0.04):
            break
        w *= wtarget / aw
        if not free_h:
            h *= htarget / ah
    aw, ah = size_of(path)
    print("  %-34s %.2f x %.2f cm (target %.2f x %.2f)" % (name, aw, ah, wtarget, htarget))
    return path, aw, ah


GAP = 0.24
# Row height no longer comes from A/D's "content stacks" (A's small table is down to one line); it is DERIVED FROM B's aspect lock:
#   B ink height = column width / 1.63 + 0.90 (0.90 = x labels 0.40 + bar gap 0.16 + colour bar 0.17 + ticks 0.34, measured)
#   to fill the row height B needs column width = 1.63*(H-0.90); the remaining width goes to A and C, the three exactly spanning 17.00.
# Targets keep 0.12 of slack: fit's convergence tolerance is ±0.04 per panel, and three panels drifting the same way would break the 17.00 hard gate (the previous version failed at 17.04).
WA, WD, WBx = 5.80, 4.85, 5.75
assert WA + WBx + WD + 2 * GAP <= 17.00 - 0.10, WA + WBx + WD + 2 * GAP
os.makedirs(DST, exist_ok=True)
# Render B first to get its natural height (aspect-locked UMAP + horizontal bar, height set by geometry), then hand the row height to A/C to align to.
# The analytic H = WB/1.63 + 0.90 measured 0.13 short (the bar's tick row height varies with type size), so no guessing: measure directly.
pb = fit(render_B, WBx, 5.00, "fig6v36_B_umap.pdf", free_h=True)
H = pb[2]
pa = fit(render_A, WA, H, "fig6v36_A_score.pdf")
pc = fit(render_C, WD, H, "fig6v36_C_ladder.pdf")
if os.environ.get("FIG6_PANELS_ONLY", "0") != "1":
        # Journal-sheet composition (the row, the letters, the page-width gate) is the delivered
        # driver's job, not the notebook's: the author's reader contract is "inputs -> the panels",
        # with the composed page kept as the driver's product. The notebook exports
        # FIG6_PANELS_ONLY=1 and stops here with the three panel PDFs at the delivered sizes.
    items = [pa, pb, pc]
    print("Fig6 final A %.2f + B %.2f + C %.2f + 2 x %.2f = %.2f cm (target cap 17.00) | row height measured from B = %.2f"
          % (WA, WBx, WD, GAP, WA + WBx + WD + 2 * GAP, H))
    HP = max(i[2] for i in items)
    # Letter strip: A's ylabel is a two-line vertical text; after shortening the page its top edge hits the letter (collide_gate measured
    # 'norm_score' x 'A' 0.15 x 0.15 cm). Move all three panels down by LSTRIP so the letters own a top strip,
    # and the three letters share one baseline - steadier than nudging labels panel by panel; the only cost is +0.34 cm of page height.
    LSTRIP = 0.34
    doc = fitz.open()
    # Page width lands exactly on 17.00: house-spec gates two-column at <=17.00, but the typesetter asks "so the total is 17 cm",
    # so the slack is spread into the gaps and the page comes out exactly 17.00 - NO panel is enlarged (§3.3 zero vector scaling stands),
    # the three panels simply sit 0.035 cm further apart. After spreading, the gap itself is still sanity-checked (too wide = scattered, too narrow = glued).
    GAP = (17.00 - sum(i[1] for i in items)) / (len(items) - 1)
    assert 0.20 <= GAP <= 0.40, "spread gap %.3f cm is unreasonable (three panels total %.2f)" % (
        GAP, sum(i[1] for i in items))
    tw = sum(i[1] for i in items) + (len(items) - 1) * GAP
    print("gap auto-spread to %.3f cm => page width %.2f" % (GAP, tw))
    page = doc.new_page(width=tw * PT, height=(HP + LSTRIP) * PT)
    x = 0.0
    for k, (p, pw, ph) in zip("ABC", items):
        # top-aligned: the three columns' letters and ink tops share a line; B's aspect-locked height difference is absorbed into a blank band at the column bottoms,
        # which reads as "intentional whitespace" more than splitting the difference top and bottom (spec §3.3: add breathing room).
        ytop = LSTRIP + (0.0 if k == "B" else (HP - ph) / 2.0)
        src = fitz.open(p)
        page.show_pdf_page(fitz.Rect(x * PT, ytop * PT, (x + pw) * PT, (ytop + ph) * PT), src, 0)
        src.close()
        page.insert_text((x * PT + 1.1, (ytop - LSTRIP) * PT + 9.5), k, fontname="Arial-Regular",
                         fontfile=FONTS_PDF["normal"], fontsize=8, color=(0, 0, 0))
        x += pw + GAP
    OUTPDF = os.path.join(DST, "Figure_6_v36_ABC_row.pdf")
    doc.save(OUTPDF)
    d = fitz.open(OUTPDF); p = d[0]
    C = 2.54 / 72
    print("composed %s: %.2f x %.2f cm | page-width gate 17.00 -> %s"
          % (os.path.basename(OUTPDF), p.rect.width * C, p.rect.height * C,
             "PASS" if p.rect.width * C <= 17.005 else "FAIL"))
    ws = p.get_text("words")
    print("text rightmost %.2f / %.2f, bottom %.2f / %.2f"
          % (max(w[2] * C for w in ws), p.rect.width * C, max(w[3] * C for w in ws), p.rect.height * C))
    szs = sorted({round(s["size"], 2) for bl in p.get_text("dict")["blocks"] for ln in bl.get("lines", [])
                  for s in ln["spans"] if s["text"].strip()})
    print("face type sizes:", szs)
    p.get_pixmap(dpi=300).save(OUTPDF.replace(".pdf", ".png"))
else:
    print("FIG6_PANELS_ONLY=1: panels written to %s (row composition stays in the driver)" % DST)
