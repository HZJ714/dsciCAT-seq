# -*- coding: utf-8 -*-
"""FigS6 v16 panel re-render (three rows, six cells, 8.38 cm wide each)
====================================================
Per FIG6_S6_LAYOUT_FINAL_20261117.md §4.1:
  row 1  A four MS loci tracks  | B 331-platform-gene window permutation
  row 2  C top-1000 evidence chain (223 rows) | D 38-gene landing  <- produced by p0_8c_fig6c_{top1000,landing}.py via env
  row 3  E RNA per-cell norm_score, three chemistries | F ATAC per-cell norm_score, three chemistries
This script only renders A/B/E/F (C/D reuse the main figure's data chain, so a second implementation cannot drift the numbers).

Discipline (figure-house-spec):
  * Two font sizes only: 8.0 (panel letters / axis titles / row headers) and 6.0 (ticks, corner notes, box names, statistics).
    v10's 5.5/5.7/5.8/6.2/6.8 sizes were all raised; crowding is never solved by shrinking type (spec §3.3: breathe first, wrap second, trim last).
  * All face text is pure black (v10's 0.2/0.35/0.45 grey and red/blue emphasis text was blackened; emphasis colour stays on the marks themselves).
  * Panel letters are added by the compose layer; no title sentence inside a cell.
Usage: DSCI_WORK_ROOT=<work> python3 figS6v16_panels.py
"""
import os, sys, glob
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("DSCI_WORK_ROOT", os.path.abspath(os.path.join(HERE, "..", "..", "..")))
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, WORK + "/scripts")
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "upstream"))   # repo layout: ../upstream
from support_common import R, say
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import fitz
from arial_style import use_arial, patch_fitz_arial
use_arial({"axes.linewidth": 0.5, "savefig.dpi": 600})
patch_fitz_arial(fitz)

# ---- portability contract (the repo notebook and the local pipeline share one code path; only env vars differ) ----
#   S6_OUT   panel PDF output dir       default = analysis tree .../P0_8C/panels (historical behaviour kept)
#   S6_V8    delivered tables_v8       default = .../SUPPORT_v2/tables_v8
#   S6_CACHE small re-input caches (gene coords / UMAP table / per-nucleus scores) - the run also works without them,
#            it just has to recompute from GTF / h5ad / 6 GB .score.gz; results are written back to the cache.
#   S6_SCMS  scDRS output tree (locus hit table + GENCODE GTF)
PAN = os.environ.get("S6_OUT", R + "/P0_8C/panels"); os.makedirs(PAN, exist_ok=True)
SUP2 = R + "/SUPPORT_v2"; V8 = os.environ.get("S6_V8", SUP2 + "/tables_v8")
CACHE = os.environ.get("S6_CACHE", os.path.join(PAN, "..", "caches")); os.makedirs(CACHE, exist_ok=True)
SCMS = os.environ.get("S6_SCMS", WORK + "/FY/oligo/scdrs_ms_out")
R_OUT = R
CM = 2.54
FS_LAB, FS_TITLE = 6.0, 8.0                      # the two-size policy
RED, OLIGO, OPC = "#D51F26", "#1F77B4", "#FF7F0E"
ARMS = ["dsciCAT-seq (own)", "10x Multiome", "HT-scCAT-seq"]
ARM_LAB = {"dsciCAT-seq (own)": "dsciCAT-seq", "10x Multiome": "10x Multiome",
           "HT-scCAT-seq": "HT-scCAT-seq"}
W_CELL = 8.38
BLUERED = mcolors.LinearSegmentedColormap.from_list("bluered", ["#9ECAE1", "#FFFFFF", "#E60000"])


def tw(t, fs=FS_LAB):
    """Arial(helv) text width in cm: wrap/trim decisions are always measured, never estimated."""
    return fitz.get_text_length(t, fontname="helv", fontsize=fs) * CM / 72.0


def savefig(fig, name):
    out = os.path.join(PAN, name + ".pdf")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    d = fitz.open(out); p = d[0]; C = CM / 72
    w, h = p.rect.width * C, p.rect.height * C
    ws = p.get_text("words")
    right = max((x[2] * C for x in ws), default=0); bot = max((x[3] * C for x in ws), default=0)
    # fitz 1.24 span["color"] is an sRGB int (not a tuple): only pure black/white allowed on the face
    cols = sorted({s["color"] for bl in p.get_text("dict")["blocks"]
                   for ln in bl.get("lines", []) for s in ln["spans"] if s["text"].strip()})
    nonbw = ["#%06x" % c for c in cols if c not in (0x000000, 0xffffff)]
    szs = sorted({round(s["size"], 2) for bl in p.get_text("dict")["blocks"]
                  for ln in bl.get("lines", []) for s in ln["spans"] if s["text"].strip()})
    d.close()
    print("  %-30s %.2f x %.2f cm | ink right %.2f/%.2f bottom %.2f/%.2f | sizes %s | non-bw %s"
          % (name + ".pdf", w, h, right, w, bot, h, szs, nonbw))
    assert right <= w - 0.01 and bot <= h - 0.01, name
    return out


# ================================================================ A: four MS loci tracks
def panel_A():
    import ms_locus_overlap as mlo
    HIT = pd.read_csv(SCMS + "/ms_locus_overlap_table.csv")
    gc_ = os.path.join(CACHE, "s6a_gene_coords.csv")
    if os.path.exists(gc_):                       # the cache stores only (chrom, tss, name)
        _g = pd.read_csv(gc_)
        GENES = list(zip(_g.chrom.astype(str), _g.tss.astype(int), _g.name.astype(str)))
        say("[A] gene coordinates from cache %s (%d entries)" % (os.path.basename(gc_), len(GENES)))
    else:
        GENES = mlo.load_genes(SCMS + "/gencode.v45.basic.annotation.gtf.gz")
        pd.DataFrame({"chrom": [g[0] for g in GENES], "tss": [g[1] for g in GENES],
                      "name": [g[2] for g in GENES]}).to_csv(gc_, index=False)
    SEL4 = ["PHLDB1", "PHGDH", "BCAS1", "PLLP"]
    fig = plt.figure(figsize=(11.36 / CM, 3.40 / CM))   # 3.40 canvas -> 2.75 ink after tight, matching B's 2.73 height
    # Keeps the original S6 row-1 layout: A full 11.36 + B 5.40 (a histogram never needs half a page), both cells 3.00 high.
    # Top margin 1.15: two title lines + two in-cell corner-note lines; the second title line (chr11) would touch the first corner line (GWAS) by 0.04 cm, measured.
    gs = fig.add_gridspec(1, 4, left=0.055, right=0.995, top=1 - 1.15 / 3.40, bottom=0.42 / 3.40,
                          wspace=0.16)
    for k, gname in enumerate(SEL4):
        h = HIT[HIT.gene == gname].iloc[0]
        ax = fig.add_subplot(gs[0, k])
        ch = str(h["chr"]); lo, hi = int(h.snp_pos) - 500_000, int(h.snp_pos) + 500_000
        for c2, t, nm in GENES:
            if c2 == ch and lo <= t <= hi and nm.upper() != gname:
                ax.plot([t, t], [0.42, 0.58], color="grey", lw=0.8)
        ax.plot([h.gene_tss, h.gene_tss], [0.42, 0.58], color=RED, lw=2.0)
        ax.plot([h.snp_pos, h.snp_pos], [0.70, 0.86], color="black", lw=1.2)
        ax.plot([h.gene_tss, h.gene_tss], [0.40, 0.28], color=RED, lw=0.5)
        ha = "left" if h.gene_tss > hi - 0.18 * (hi - lo) else (
            "right" if h.gene_tss < lo + 0.18 * (hi - lo) else "center")
        ax.text(h.gene_tss, 0.17, gname, fontsize=FS_LAB, ha=ha, va="center", color="black")
        # Two-line titles: a single "PHLDB1  chr11" line @8pt measures 1.85 cm > facet width 1.79 cm and would collide
        ax.set_title("%s\nchr%s" % (gname, ch), fontsize=FS_TITLE, pad=1.2, linespacing=1.25)
        ax.text(0.02, 0.93, "GWAS lead-SNP\nP = %.0e" % h.pvalue, transform=ax.transAxes,
                ha="left", va="top", fontsize=FS_LAB, linespacing=1.3, color="black")
        ax.text(0.02, 0.66, "dist %+.0f kb" % (h.distance_bp / 1000.0), transform=ax.transAxes,
                ha="left", va="top", fontsize=FS_LAB, color="black")
        ax.set_xlim(lo, hi); ax.set_ylim(0, 1.0)
        ax.set_xticks([lo, hi])
        ax.set_xticklabels(["%.1f" % (lo / 1e6), "%.1f" % (hi / 1e6)], fontsize=FS_LAB)
        for lb, al in zip(ax.get_xticklabels(), ["left", "right"]):
            lb.set_ha(al)
        ax.set_yticks([])
        for sp in ["top", "right", "left"]:
            ax.spines[sp].set_visible(False)
        for sp in ax.spines.values():
            sp.set_linewidth(0.5)
        ax.tick_params(width=0.5, length=2, pad=0.8)
    return savefig(fig, "figS6v16_A_loci")


# ================================================================ B: 331 platform genes, window permutation
CM = 2.54   # shared-kernel hygiene: the D/E panel cells overwrite CM with 1/2.54 (cm->inch inverse); re-asserted per section, a no-op when the whole file runs
def panel_B():
    WPERM = pd.read_csv(V8 + "/figS6v8_window_permutation.csv").set_index("stat")
    WNULL = pd.read_csv(V8 + "/figS6v8_window_permutation_null.csv")
    fig, ax = plt.subplots(figsize=(5.40 / CM, 3.00 / CM))
    fig.subplots_adjust(left=0.175, right=0.985, top=0.90, bottom=0.235)
    obs_pct = float(WPERM.loc["observed (%)", "value"])
    ax.hist(WNULL.null_pct.values, bins=np.arange(6, 17, 0.5), color=(0.69, 0.769, 0.871),
            edgecolor=(0.28, 0.40, 0.56), linewidth=0.4, zorder=3)
    ax.axvline(obs_pct, color=RED, lw=1.1, zorder=4)
    ax.set_xlim(5.4, 16.6); ax.set_xticks([6, 8, 10, 12, 14, 16])
    # A single 8pt line measures 4.9 cm > the 4.4 cm available and would run past the page right (the self-check caught 5.49>5.38) -> two-line titles, no font shrink
    ax.set_xlabel("% of the 331 platform genes\ninside window", fontsize=FS_TITLE,
                  labelpad=2.0, linespacing=1.3)
    ax.set_ylabel("permutations", fontsize=FS_TITLE, labelpad=2.0)
    ax.tick_params(labelsize=FS_LAB, length=2, width=0.5, pad=1.0)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    pperm = float(WPERM.loc["permutation P (one-sided, observed vs null)", "value"])
    # 2026-11-18 (author, 2nd pass): the top-LEFT placement still sat on the histogram - the left
    # bins run taller than the naive estimate, the two long lines landed on 32 % ink density and
    # the text was unreadable. The empty region is the TOP-RIGHT corner (measured 0 ink at
    # x >= 15.0 cm page, y 0.4-1.6 cm: the right bins top out at ~0.13 of ymax). The red observed
    # line is at page x 14.34, so only 1.89 cm of width is available right of it; the widest of the
    # original two lines is 2.10 cm -> the note is re-typeset as FOUR short lines (widest 1.047 cm,
    # left edge 15.26, 0.88 cm clear of the line), no white pad needed. Measured, not eyeballed.
    ax.text(0.985, 0.965,
            "obs 38/331\n= %.1f%%\nn = 10,000\nP = %.3f" % (obs_pct, pperm),
            transform=ax.transAxes, fontsize=FS_LAB, va="top", ha="right", linespacing=1.35,
            color="black", zorder=6)
    return savefig(fig, "figS6v16_B_perm")


# ================================================================ C: bin-size UMAP (the Fig 6B supplement)
CM = 2.54   # shared-kernel hygiene: the D/E panel cells overwrite CM with 1/2.54 (cm->inch inverse); re-asserted per section, a no-op when the whole file runs
# Original S6 row-2 layout: RNA in four cells (single cells / gridsize 15 / 25 / 40) + one colour bar on the right, full width 17.00.
# Only two compliance fixes vs the original: (1) set_aspect("equal") locks the shape (the original's 3.55x3.20 flattened the 1.63
# UMAP); (2) all type 8.0/6.0, bold and grey text removed, annotation colour kept on marks only.
def panel_C():
    import anndata as ad
    CW, CH = 17.00, 3.00
    fig = plt.figure(figsize=(CW / CM, CH / CM))
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)      # exact cm placement, no tight bbox
    cc_ = os.path.join(CACHE, "s6c_umap_score.csv")
    if os.path.exists(cc_):                           # the four columns are this panel's entire input
        dU = pd.read_csv(cc_)
        print("  [C] UMAP+score from cache %s" % os.path.basename(cc_))
    else:
        obs = ad.read_h5ad(WORK + "/FY/oligo/3k_dsci_oligo_subset.h5ad", backed="r").obs
        CT = pd.read_csv(R + "/O1_own_multiome_paired/tables/cell_level_table.csv", low_memory=False)
        dU = obs[["UMAP_1", "UMAP_2", "cell_type"]].join(
            CT.set_index("cell_id")[["PASS_Multiple_sclerosis_rna"]], how="inner")
        dU = dU.astype({"UMAP_1": float, "UMAP_2": float, "PASS_Multiple_sclerosis_rna": float})
        dU["cell_type"] = dU["cell_type"].astype(str)     # cell_type is categorical; astype(float) on the whole table would blow up
        dU.to_csv(cc_, index=False)
    say = lambda m: print("  " + m)
    say("[C] UMAP reconciled on %d nuclei" % len(dU))
    xr = max(np.nanpercentile(dU.PASS_Multiple_sclerosis_rna, 99),
             abs(np.nanpercentile(dU.PASS_Multiple_sclerosis_rna, 1)))
    x0, x1 = np.percentile(dU.UMAP_1, [0.2, 99.8]); y0, y1 = np.percentile(dU.UMAP_2, [0.2, 99.8])
    ASP = (x1 - x0) / (y1 - y0)
    TW, TY = 3.55, 0.46
    TH = TW / ASP                                   # locked aspect: cell height from the data's own length/width ratio
    TX = [0.62, 4.52, 8.42, 12.32]
    TIT = ["single cells", "gridsize 15", "gridsize 25 (Fig 6B)", "gridsize 40"]
    def tile(i):
        a = fig.add_axes([TX[i] / CW, 1 - (TY + TH) / CH, TW / CW, TH / CH])
        a.set_title(TIT[i], fontsize=FS_TITLE, pad=1.8, loc="left")
        a.set_xticks([]); a.set_yticks([]); a.set_xlim(x0, x1); a.set_ylim(y0, y1)
        a.set_xlabel("UMAP_1", fontsize=FS_LAB, labelpad=1.2)
        a.set_aspect("equal")
        for sp in a.spines.values(): sp.set_linewidth(0.5)
        return a
    ax1 = tile(0)
    ax1.scatter(dU.UMAP_1, dU.UMAP_2, s=0.35, c=dU.PASS_Multiple_sclerosis_rna, cmap=BLUERED,
                vmin=-xr, vmax=xr, linewidths=0, rasterized=True)
    HB = []
    for i, gs_ in ((1, 15), (2, 25), (3, 40)):
        a = tile(i)
        hb = a.hexbin(dU.UMAP_1, dU.UMAP_2, C=dU.PASS_Multiple_sclerosis_rna,
                      reduce_C_function=np.mean, gridsize=gs_, mincnt=1, cmap=BLUERED,
                      linewidths=0, edgecolors="none", zorder=2)
        hb.set_clim(-xr, xr); hb.set_rasterized(True); HB.append(hb)
    sp1 = dU.UMAP_1.max() - dU.UMAP_1.min(); sp2 = dU.UMAP_2.max() - dU.UMAP_2.min()
    mol = dU[dU.cell_type == "MOL6"][["UMAP_1", "UMAP_2"]].median()
    cop = dU[dU.cell_type == "COP/NFOL"][["UMAP_1", "UMAP_2"]].median()
    # The original ha="center" let the arrow pierce the middle of the text box (text-on-data gate, LINE-ACROSS-TEXT 22%).
    # Now anchored at the box's lower-left corner + white pad: the leader ends at the box corner, shrinkB backs off another 2 pt. Offsets compensate half a cell per label, by measured width.
    WB = dict(facecolor="white", edgecolor="none", alpha=0.85, pad=0.6)
    for a in (ax1, HB[0].axes):
        a.annotate("MOL6", xy=(mol.UMAP_1, mol.UMAP_2),
                   xytext=(mol.UMAP_1 + 0.02 * sp1, mol.UMAP_2 + 0.26 * sp2), fontsize=FS_LAB,
                   color=RED, ha="left", va="bottom", bbox=WB,
                   arrowprops=dict(arrowstyle="->", color=RED, lw=0.7, shrinkA=1, shrinkB=2))
        a.annotate("COP/NFOL", xy=(cop.UMAP_1, cop.UMAP_2),
                   xytext=(cop.UMAP_1 + 0.44 * sp1, cop.UMAP_2 + 0.30 * sp2), fontsize=FS_LAB,
                   color="#084594", ha="center", va="bottom", bbox=WB,
                   arrowprops=dict(arrowstyle="->", color="#084594", lw=0.7, shrinkA=1, shrinkB=2))
    cbax = fig.add_axes([16.05 / CW, 1 - (TY + TH) / CH, 0.30 / CW, TH / CH])
    cbR = matplotlib.colorbar.ColorbarBase(cbax, cmap=BLUERED,
                                           norm=matplotlib.colors.Normalize(-xr, xr))
    cbR.outline.set_linewidth(0.4)
    cbR.ax.tick_params(labelsize=FS_LAB, length=1.8, width=0.4, pad=1.0)
    cbR.set_label("norm_score", fontsize=FS_LAB, labelpad=1.5)
    out = os.path.join(PAN, "figS6v16_C_umapbins.pdf")
    fig.savefig(out); plt.close(fig)
    d = fitz.open(out); p = d[0]; C = CM / 72
    ws = p.get_text("words")
    szs = sorted({round(s["size"], 2) for bl in p.get_text("dict")["blocks"]
                  for ln in bl.get("lines", []) for s in ln["spans"] if s["text"].strip()})
    cols = sorted({s["color"] for bl in p.get_text("dict")["blocks"]
                   for ln in bl.get("lines", []) for s in ln["spans"] if s["text"].strip()})
    print("  %-30s %.2f x %.2f cm | ink right %.2f/%.2f | sizes %s | non-bw %s"
          % ("figS6v16_C_umapbins.pdf", p.rect.width * C, p.rect.height * C,
             max(x[2] * C for x in ws), p.rect.width * C, szs,
             ["#%06x" % c for c in cols if c not in (0x000000, 0xffffff)]))
    d.close()
    return out


# ================================================================ E/F: per-cell norm_score box plots across three chemistries
CM = 2.54   # shared-kernel hygiene: the D/E panel cells overwrite CM with 1/2.54 (cm->inch inverse); re-asserted per section, a no-op when the whole file runs
CELLDIR = os.environ.get("S6_CELLS", R_OUT + "/D_crosschemistry/tables_log1p_20261104")
CELLF = {"dsciCAT-seq (own)": CELLDIR + "/cells_dsciCAT_seq.csv",
         "10x Multiome": CELLDIR + "/cells_10x.csv",
         "HT-scCAT-seq": CELLDIR + "/cells_HT_scCAT_seq_HT.csv"}
DIRF = {"dsciCAT-seq (own)": ("dsciCAT_seq_rna", "dsciCAT_seq_atac"),
        "10x Multiome": ("10x_rna", "10x_atac"), "HT-scCAT-seq": ("HT_scCAT_seq_HT_rna", "HT_scCAT_seq_HT_atac")}
_C = {}


def ms_scores(tech, chan):
    """Per-nucleus MS scDRS norm_score + the Fig 2 label. The panel only uses these two columns, so the cache is exactly those two;
    when missing, recompute from .score.gz (part of the ~6 GB per-tech per-channel tree) and write back."""
    key = (tech, chan)
    if key in _C:
        return _C[key]
    slug = DIRF[tech][0 if chan == "rna" else 1]
    fc_ = os.path.join(CACHE, "s6f_%s_%s.csv" % (slug, chan))
    if os.path.exists(fc_):
        df = pd.read_csv(fc_).set_index("cell")
        say("[F] %-16s %-4s per-nucleus scores from cache %s (%d nuclei)" % (tech, chan, os.path.basename(fc_), len(df)))
        _C[key] = df
        return df
    rn, an = DIRF[tech]
    base = R_OUT + "/D_crosschemistry/scdrs/" + rn if chan == "rna" else \
        R_OUT + "/D_crosschemistry/scdrs_log1p_20261104/" + an
    fs = sorted(glob.glob(base + "/*/PASS_Multiple_sclerosis.full_score.gz"))
    if not fs:
        raise SystemExit("[F] no .score.gz under %s and no cache %s: export DSCI_WORK_ROOT, "
                         "or run upstream/S6v16_export_caches.py first to build the caches" % (base, fc_))
    df = pd.concat([pd.read_csv(f, sep="\t", index_col=0, low_memory=False) for f in fs])
    df.index.name = "cell"
    lb = pd.read_csv(CELLF[tech]).set_index("cell")["label"]
    df = df[df.index.isin(lb.index)].copy()
    df["label"] = lb.loc[df.index].values
    df[["norm_score", "label"]].to_csv(fc_)
    df = df[["norm_score", "label"]]
    _C[key] = df
    return df


MWU = pd.read_csv(V8 + "/figS6v10_G_opc_oligo_mwu.csv", dtype={"p_fig": str, "p_label": str})
GPT = {"dsciCAT-seq (own)": dict(s=0.9, alpha=0.20, lighten=0.55),
       "10x Multiome": dict(s=1.6, alpha=0.30, lighten=0.0),
       "HT-scCAT-seq": dict(s=1.6, alpha=0.30, lighten=0.0)}


def _dark(c, thr=0.55):
    rgb = np.array(mcolors.to_rgb(c)); lum = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
    return tuple(rgb * (thr / lum)) if lum > thr else tuple(rgb)


def _lighten(c, f):
    rgb = np.array(mcolors.to_rgb(c)); return tuple(rgb * f + (1 - f))


def panel_boxes(chan, title_side):
    fig = plt.figure(figsize=(W_CELL / CM, 4.20 / CM))
    gs = fig.add_gridspec(1, 3, left=0.135, right=0.985, top=1 - 0.62 / 4.60, bottom=0.155,
                          wspace=0.20)
    rng = np.random.default_rng(20261117)
    for ci, tech in enumerate(ARMS):
        ax = fig.add_subplot(gs[0, ci])
        df = ms_scores(tech, chan)
        mo = df[df.label == "Oligo"].norm_score.values; op = df[df.label == "OPC"].norm_score.values
        st = GPT[tech]
        for pos, arr, col in [(0, op, OPC), (0.6, mo, OLIGO)]:
            ax.scatter(pos + rng.normal(0, 0.07, len(arr)), arr, s=st["s"],
                       color=_lighten(col, 1 - st["lighten"]) if st["lighten"] else col,
                       alpha=st["alpha"], linewidths=0, rasterized=True, zorder=1)
        bp = ax.boxplot([op, mo], positions=[0, 0.6], widths=0.5, patch_artist=True,
                        showfliers=False, whis=[2, 98], medianprops=dict(lw=1.3),
                        boxprops=dict(lw=0.8), whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7))
        for ln in bp["caps"]:
            c = float(np.mean(ln.get_xdata())); y = float(ln.get_ydata()[0])
            ln.set_data([c - 0.13, c + 0.13], [y, y])
        for b, col in zip(bp["boxes"], [OPC, OLIGO]):
            b.set_facecolor(col); b.set_alpha(0.30); b.set_edgecolor(_dark(col)); b.set_zorder(2)
        for ln, col in zip(bp["medians"], [OPC, OLIGO]):
            ln.set_color(_dark(col)); ln.set_zorder(3)
        for key in ("whiskers", "caps"):
            for ln in bp[key]:
                ln.set_color(OLIGO if ln.get_xdata()[0] > 0.2 else OPC); ln.set_zorder(2)
        r = MWU[(MWU.tech == ARM_LAB[tech]) & (MWU.channel == chan)].iloc[0]
        # Ledger consistency: the face label must be exactly p_fig's "P = " form. If the table is ever rewritten without p_label following,
        # this explodes immediately - no silent "P on the face disagrees with Table S6" drift.
        assert str(r.p_label) == "P = " + str(r.p_fig), (tech, chan, r.p_label, r.p_fig)
        ax.text(0.97, 0.05, str(r.p_label), fontsize=FS_LAB, ha="right", va="bottom",
                color="black", transform=ax.transAxes)
        ax.set_ylim(-1.62, 1.4); ax.set_xlim(-0.5, 1.1)
        ax.axhline(0, color="0.80", lw=0.5, zorder=1)
        ax.set_xticks([0, 0.6]); ax.set_xticklabels(["OPC", "Oligo"], fontsize=FS_LAB)
        ax.set_yticks([-1, 0, 1])
        if ci:
            ax.set_yticklabels([])
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
        for sp in ["left", "bottom"]:
            ax.spines[sp].set_linewidth(0.5)
        ax.tick_params(width=0.4, length=1.8, pad=0.8, labelsize=FS_LAB)
        ax.set_title(ARM_LAB[tech], fontsize=FS_TITLE, pad=2.0)
        if ci == 0:
            ax.set_ylabel("per-cell norm_score, %s" % ("RNA" if chan == "rna" else "ATAC"),
                          fontsize=FS_TITLE, labelpad=1.5)
    return savefig(fig, "figS6v16_%s_%sboxes" % (title_side, "RNA" if chan == "rna" else "ATAC"))


# ================================================================ F: both channels merged into one cell
def _boxes_cell(ax, tech, chan, rng):
    df = ms_scores(tech, chan)
    mo = df[df.label == "Oligo"].norm_score.values; op = df[df.label == "OPC"].norm_score.values
    st = GPT[tech]
    for pos, arr, col in [(0, op, OPC), (0.6, mo, OLIGO)]:
        ax.scatter(pos + rng.normal(0, 0.07, len(arr)), arr, s=st["s"],
                   color=_lighten(col, 1 - st["lighten"]) if st["lighten"] else col,
                   alpha=st["alpha"], linewidths=0, rasterized=True, zorder=1)
    bp = ax.boxplot([op, mo], positions=[0, 0.6], widths=0.5, patch_artist=True,
                    showfliers=False, whis=[2, 98], medianprops=dict(lw=1.3),
                    boxprops=dict(lw=0.8), whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7))
    for ln in bp["caps"]:
        c = float(np.mean(ln.get_xdata())); y = float(ln.get_ydata()[0])
        ln.set_data([c - 0.13, c + 0.13], [y, y])
    for b, col in zip(bp["boxes"], [OPC, OLIGO]):
        b.set_facecolor(col); b.set_alpha(0.30); b.set_edgecolor(_dark(col)); b.set_zorder(2)
    for ln, col in zip(bp["medians"], [OPC, OLIGO]):
        ln.set_color(_dark(col)); ln.set_zorder(3)
    for key in ("whiskers", "caps"):
        for ln in bp[key]:
            ln.set_color(OLIGO if ln.get_xdata()[0] > 0.2 else OPC); ln.set_zorder(2)
    r = MWU[(MWU.tech == ARM_LAB[tech]) & (MWU.channel == chan)].iloc[0]
    assert str(r.p_label) == "P = " + str(r.p_fig), (tech, chan, r.p_label, r.p_fig)
    # The P label sits lower-right, but Oligo's lower whisker runs right through the digits (verify_text_on_data threshold 6%,
    # which does not fire; by eye it reads as "1.8|e-18"). The white pad hides the whisker - a whisker is a continuous quantity, no misreading.
    ax.text(0.97, 0.05, str(r.p_label), fontsize=FS_LAB, ha="right", va="bottom",
            color="black", transform=ax.transAxes, zorder=7,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.8))
    ax.set_ylim(-1.62, 1.4); ax.set_xlim(-0.5, 1.1)
    ax.axhline(0, color="0.80", lw=0.5, zorder=1)
    ax.set_xticks([0, 0.6]); ax.set_xticklabels(["OPC", "Oligo"], fontsize=FS_LAB)
    ax.set_yticks([-1, 0, 1])
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    for sp in ["left", "bottom"]:
        ax.spines[sp].set_linewidth(0.5)
    ax.tick_params(width=0.4, length=1.8, pad=0.8, labelsize=FS_LAB)


def panel_F_combined():
    """F: RNA and ATAC merged into ONE cell (author 2026-11-17: the narrative first reads "MS risk nucleus by nucleus",
    the two channels are two sides of one thing; splitting into F and G would read one thing as two).
    2 rows = channels, 3 columns = chemistries; column names only on the top row, channel names in a two-line ylabel."""
    H0 = 6.00
    fig = plt.figure(figsize=(W_CELL / CM, H0 / CM))
    gs = fig.add_gridspec(2, 3, left=0.148, right=0.985, top=1 - 0.50 / H0, bottom=0.070,
                          hspace=0.44, wspace=0.20)
    rng = np.random.default_rng(20261117)
    for ri, chan in enumerate(["rna", "atac"]):
        for ci, tech in enumerate(ARMS):
            ax = fig.add_subplot(gs[ri, ci])
            _boxes_cell(ax, tech, chan, rng)
            if ri == 0:
                ax.set_title(ARM_LAB[tech], fontsize=FS_TITLE, pad=2.0)
            if ci == 0:
                ax.set_ylabel("per-cell norm_score\n%s" % ("RNA" if chan == "rna" else "ATAC"),
                              fontsize=FS_TITLE, labelpad=1.5, linespacing=1.30)
            else:
                ax.set_yticklabels([])
    return savefig(fig, "figS6v16_F_cellboxes")


if __name__ == "__main__":
    # 2026-11-18 fix: the old __main__ still called panel_boxes("rna","F")/panel_boxes("atac","G"),
    # the abandoned split-into-two-cells plan (round 5, rejected by the author: one thing cannot read as two). The delivered F uses
    # panel_F_combined() (2 channels x 3 chemistries, one cell). Running this script now produces the delivered four panels.
    print("S6 panels re-rendered (two-size type / pure-black face text)")
    panel_A(); panel_B(); panel_C(); panel_F_combined()
