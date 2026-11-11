#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build Figure_6.ipynb — panel-by-panel reproduction of the journal Figure 6
(MS scDRS selectivity along the oligodendrocyte lineage), organised exactly
like Figure_4.ipynb. Run:  python3 build_fig6_notebook.py
"""
import os
import nbformat as nbf

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "fig6")


def md(s):
    return nbf.v4.new_markdown_cell(s)


def code(s):
    return nbf.v4.new_code_cell(s)


cells = []

# =============================================================================
# Section 0 — environment, fonts, versions
# =============================================================================
cells.append(md("""# Figure 6 — MS scDRS selectivity along the oligodendrocyte lineage

This notebook reproduces the four panels of **Figure 6** from the upstream
analysis objects, panel by panel:

| Panel | File | Size | Content |
|---|---|---|---|
| A | `fig6a_boxplot.pdf` | 8.30 × 6.40 cm | scDRS MS norm-score boxplot by subtype; MWU P vs the lowest-median subtype (COP/NFOL); whiskers 5–95 % |
| B | `fig6b_umap.pdf` | 8.30 × 5.60 cm | UMAP hexbin of binned mean score (gridsize 25, ≥5 cells/bin, ±0.74) |
| C | `fig6c_manhattan.pdf` | 10.40 × 4.60 cm | MS genome-wide Manhattan (178 lead SNPs, PMID 31604244); hit genes ±500 kb labelled |
| D | `fig6d_74trait_strip.pdf` | 6.32 × 4.60 cm | MOL6 selectivity across 74 GWAS gene sets (row z; MS rank 1/74, empirical P) |

Each section ends with its own panel PDF. The UMAP embedding is taken as a
*display axis* from the h5ad (upstream Seurat coordinates); the content of
the figure is the scDRS score × annotated `cell_type`."""))
cells.append(md("""---

## 0. Environment, paths and versions

Set the paths once; everything below reads only these objects. Python ≥3.8
with the packages printed in Section 0.2; no R is needed for Figure 6."""))
cells.append(code('''import os, sys, gzip, bisect, random, re, glob
import numpy as np, pandas as pd, scipy, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import font_manager as fm
import scanpy as sc
# ---- Arial-first font registration + project rcParams (use_arial equivalent) ----
FONT_DIR = os.environ.get("FONT_DIR", os.path.join(os.environ.get("DSCI_WORK_ROOT", "."), "fonts"))
for _f in sorted(glob.glob(os.path.join(FONT_DIR, "*.ttf"))) + \\
         sorted(glob.glob(os.path.join(FONT_DIR, "*.TTF"))) + \\
         sorted(glob.glob(os.path.join(FONT_DIR, "*.otf"))):
    try:
        fm.fontManager.addfont(_f)
    except Exception:
        pass
_avail = {f.name for f in fm.fontManager.ttflist}
_FAMILY = next((n for n in ("Arial", "Liberation Sans") if n in _avail), "DejaVu Sans")
matplotlib.rcParams.update({"font.family": _FAMILY, "mathtext.fontset": "custom",
                             "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
                             "font.size": 7.0, "axes.titlesize": 7.0, "axes.labelsize": 7.0,
                             "xtick.labelsize": 6.0, "ytick.labelsize": 6.0,
                             "legend.fontsize": 6.0, "axes.titleweight": "normal",
                             "axes.linewidth": 0.6, "xtick.major.width": 0.5,
                             "ytick.major.width": 0.5, "xtick.minor.width": 0.4,
                             "ytick.minor.width": 0.4, "xtick.major.size": 2.0,
                             "ytick.major.size": 2.0, "xtick.minor.size": 1.0,
                             "ytick.minor.size": 1.0, "xtick.major.pad": 1.6,
                             "ytick.major.pad": 1.6, "lines.linewidth": 0.8,
                             "patch.linewidth": 0.5, "figure.dpi": 120, "savefig.dpi": 600})
print("font family:", _FAMILY)
# ---- paths ----
DATA = os.environ.get("S6_DATA", os.path.join(os.environ.get("DSCI_WORK_ROOT", "."), "FY", "oligo", "scdrs_ms_out"))
H5AD = os.environ.get("S6_H5AD", os.path.join(os.environ.get("DSCI_WORK_ROOT", "."), "FY", "oligo", "3k_dsci_oligo_subset.h5ad"))
MISC = os.environ.get("S6_MISC", os.path.join(os.environ.get("DSCI_WORK_ROOT", "."), "FY", "oligo", "scdrs_traits_draft", "full74"))
ARCH = os.environ.get("S6_ARCH", os.path.join(os.environ.get("DSCI_WORK_ROOT", "."), "FY", "oligo", "scdrs_crossplatform_20261020"))
OUT  = os.environ.get("FIG6_OUT", "fig6_out")
os.makedirs(OUT, exist_ok=True)
CM = 2.54   # cm → inch
# ---- layout helpers ----
def open_spines(ax, top=True, right=True):
    if top:
        ax.spines["top"].set_visible(False)
    if right:
        ax.spines["right"].set_visible(False)
def frame(ax, lw=0.6):
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_linewidth(lw)
def align_xlabel(ax, fig, y_fig=0.040):
    from matplotlib.transforms import blended_transform_factory as btf
    ax.xaxis.set_label_coords(0.5, y_fig, transform=btf(ax.transAxes, fig.transFigure))
# ---- GWAS Catalog / GENCODE parsers (shared by C and S6A/S6E) ----
def load_ms_snps(assoc_zip, pmid="31604244", p_thresh=5e-8):
    """MS genome-wide-significant lead SNPs from the GWAS Catalog (PMID 31604244)."""
    rows = []
    with zipfile.ZipFile(assoc_zip) as z:
        name = z.namelist()[0]
        with z.open(name) as raw:
            f = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
            header = next(f).rstrip("\\n").split("\\t")
            idx = {h: i for i, h in enumerate(header)}
            for line in f:
                cols = line.rstrip("\\n").split("\\t")
                if len(cols) <= idx["CHR_POS"]:
                    continue
                if cols[idx["PUBMEDID"]] != pmid:
                    continue
                try:
                    p = float(cols[idx["P-VALUE"]])
                except ValueError:
                    continue
                if p > p_thresh:
                    continue
                if not cols[idx["CHR_POS"]].strip() or not cols[idx["CHR_ID"]].strip():
                    continue
                rows.append({"snp": cols[idx["SNP_ID_CURRENT"]],
                             "chr": str(cols[idx["CHR_ID"]]).replace("chr", ""),
                             "pos": int(cols[idx["CHR_POS"]]), "p": p,
                             "mapped_gene": cols[idx["MAPPED_GENE"]]})
    df = pd.DataFrame(rows).drop_duplicates(subset=["snp", "chr", "pos"])
    return df
def load_genes(gtf_gz):
    """GENCODE gene features -> list of (chr, tss, name)."""
    genes = []
    with gzip.open(gtf_gz, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\\n").split("\\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            m = re.search(r'gene_name "([^"]+)"', cols[8])
            if not m:
                continue
            chrom = cols[0].replace("chr", ""); start, end = int(cols[3]), int(cols[4])
            genes.append((chrom, start if cols[6] == "+" else end, m.group(1)))
    return genes'''))

cells.append(md("""### 0.1 Input files

| Variable | Upstream object | Used by |
|---|---|---|
| `DATA` | `scdrs_ms_out/` — `PASS_Multiple_sclerosis.score.gz`, `MS_scDRS_stats.csv`, GWAS Catalog zip, GENCODE gtf, `ms_locus_overlap_table.csv` | A–C |
| `H5AD` | `3k_dsci_oligo_subset.h5ad` (5,977 cells, 6 subtypes, UMAP) | A, B |
| `MISC` | `scdrs_traits_draft/full74/` — 74 trait `norm_score` files | D |

### 0.2 Software versions"""))
cells.append(code('''for m in (np, pd, scipy, matplotlib, sc):
    print(f"  {m.__name__:<12} {m.__version__}")
print("  Python", sys.version.split()[0])'''))

# =============================================================================
# Section 1 — shared data
# =============================================================================
cells.append(md("""---

## 1. Shared loading (one cell)

Loads the annotated h5ad + MS scDRS scores (→ `df` with `cell_type`,
`UMAP_1/2`, `norm_score`), the pre-computed MWU P values (→ `mw_p`), the 178
MS lead SNPs (→ `ms_snps`), the ±500 kb overlap table (→ `hit`), and GENCODE
TSS positions (→ `genes`)."""))
cells.append(code('''import zipfile, io
adata = sc.read_h5ad(H5AD)
score = pd.read_csv(os.path.join(DATA, "PASS_Multiple_sclerosis.score.gz"),
                    sep="\\t", index_col=0)
if "UMAP_1" in adata.obs.columns:
    umap = adata.obs[["UMAP_1", "UMAP_2"]]
else:
    umap = pd.DataFrame(adata.obsm["X_umap"][:, :2], index=adata.obs_names,
                        columns=["UMAP_1", "UMAP_2"])
df = adata.obs[["cell_type"]].join(umap).join(score, how="inner")
stats = pd.read_csv(os.path.join(DATA, "MS_scDRS_stats.csv"))
mw_p = dict(zip(stats["cell_type"], stats["MWU_p_vs_lowest_median"]))
order = stats["cell_type"].tolist()
pal = dict(zip(order, plt.cm.tab10.colors))
def fmtp(p):
    if p >= 0.001:
        return "%.3f" % p
    e = int(np.floor(np.log10(p)))
    return "%.1fe%d" % (p / 10 ** e, e)
ms_snps = load_ms_snps(os.path.join(DATA, "gwas-catalog-associations-full.zip"))
hit = pd.read_csv(os.path.join(DATA, "ms_locus_overlap_table.csv"))
ms_snps = ms_snps.sort_values(["chr", "pos"])
genes = load_genes(os.path.join(DATA, "gencode.v45.basic.annotation.gtf.gz"))
print(f"cells {len(df)}  | MS lead SNPs {len(ms_snps)}  | overlap pairs {len(hit)}  | genes {len(genes)}")'''))

# =============================================================================
# Section 2 — Panel A (boxplot)
# =============================================================================
cells.append(md("""---

## 2. Panel A — scDRS boxplot by OL-lineage subtype (8.30 × 6.40 cm)

**Data chain**: `norm_score` of `PASS_Multiple_sclerosis.score.gz` joined to
the h5ad `cell_type`. **Core parameters**: whiskers = 5–95 percentiles,
`showfliers=False`; MWU P of each subtype vs the lowest-median subtype
(COP/NFOL), taken from `MS_scDRS_stats.csv` (computed upstream on the full
cell pool; the notebook only *reads* it to label the panel)."""))
cells.append(code('''matplotlib.rcParams.update({"xtick.labelsize": 6.5, "ytick.labelsize": 6.5})
fig, ax = plt.subplots(figsize=(8.3 / CM, 6.4 / CM))
fig.subplots_adjust(left=0.150, right=0.985, top=0.975, bottom=0.255)
data = [df.loc[df["cell_type"] == ct, "norm_score"].values for ct in order]
bp = ax.boxplot(data, showfliers=False, whis=[5, 95], patch_artist=True,
                widths=0.55, medianprops=dict(color="black", lw=1.0),
                whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7))
for i, (patch, ct) in enumerate(zip(bp["boxes"], order)):
    patch.set_facecolor(pal[ct]); patch.set_alpha(0.55); patch.set_linewidth(0.6)
    x = np.random.RandomState(7).normal(i + 1, 0.07, len(data[i]))
    ax.scatter(x, data[i], s=1.5, color=pal[ct], alpha=0.35, lw=0, zorder=3)
ax.axhline(0, color="grey", ls="--", lw=0.7)
ax.set_xticks(range(1, len(order) + 1))
ax.set_xticklabels(order, rotation=40, ha="right", fontsize=7)
for lbl, ct in zip(ax.get_xticklabels(), order):
    lbl.set_color([float(c) for c in pal[ct][:3]]); lbl.set_fontweight("bold")
ax.set_ylabel("scDRS norm_score\\n(MS GWAS)", linespacing=1.15)
wmax = max(l.get_ydata().max() for l in bp["whiskers"])
y1, y2 = wmax + 0.22, wmax + 1.00
for i, ct in enumerate(order):
    p = mw_p.get(ct, np.nan)
    if pd.isna(p):
        continue
    ax.text(i + 1, y1 if i % 2 == 0 else y2, fmtp(p), ha="center",
            va="bottom", fontsize=6, color="crimson" if p < 0.05 else "grey",
            fontweight="bold")
ax.set_ylim(-2.6, y2 + 0.62)
ax.set_yticks([t for t in [-2, -1, 0, 1, 2, 3] if t < y2 + 0.62])
open_spines(ax)
fig.savefig(os.path.join(OUT, "fig6a_boxplot.pdf"))
plt.close(fig)
print("A done")'''))

# =============================================================================
# Section 3 — Panel B (UMAP hexbin)
# =============================================================================
cells.append(md("""---

## 3. Panel B — UMAP hexbin of binned mean score (8.30 × 5.60 cm)

**Core parameters**: hexbin `gridsize=25`, `mincnt=5` (drop bins with <5
cells), color scale ±`vB` = ±max |binned mean| (= 0.74): low = light blue
`#9ECAE1`, 0 = white, high = red. The COP/NFOL precursor tail maps to the
negative end, the mature MOL6 cloud to the positive end."""))
cells.append(code('''BLUERED = LinearSegmentedColormap.from_list(
    "bluered", ["#9ECAE1", "#FFFFFF", "#E60000"])
fig = plt.figure(figsize=(8.3 / CM, 5.6 / CM))
gs = fig.add_gridspec(1, 2, width_ratios=[22, 0.75], wspace=0.045,
                      left=0.085, right=0.87, top=0.875, bottom=0.155)
ax = fig.add_subplot(gs[0]); axc = fig.add_subplot(gs[1])
hb = ax.hexbin(df["UMAP_1"], df["UMAP_2"], C=df["norm_score"],
               reduce_C_function=np.mean, gridsize=25, mincnt=5,
               cmap=BLUERED, linewidths=0, edgecolors="none")
hb.set_rasterized(True)
vals = hb.get_array()
vB = max(np.nanpercentile(vals, 99), abs(np.nanpercentile(vals, 1)))
hb.set_clim(-vB, vB)
ax.set_xlabel("UMAP_1"); ax.set_ylabel("UMAP_2")
align_xlabel(ax, fig)
ax.set_xticks([]); ax.set_yticks([])
ax.set_title("MS scDRS score on UMAP", fontsize=7, pad=3, loc="left")
cop = df[df["cell_type"] == "COP/NFOL"]
mol6 = df[df["cell_type"] == "MOL6"]
ax.annotate("COP/NFOL",
            xy=(cop["UMAP_1"].median(), cop["UMAP_2"].median()),
            xytext=(cop["UMAP_1"].median() + 1.6, cop["UMAP_2"].median() + 2.6),
            fontsize=6.0, fontweight="bold", color="#084594",
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.4),
            arrowprops=dict(arrowstyle="->", color="#084594", lw=0.8,
                            shrinkA=1, shrinkB=1.5))
ax.annotate("MOL6",
            xy=(mol6["UMAP_1"].median(), mol6["UMAP_2"].median()),
            xytext=(mol6["UMAP_1"].median() + 2.4, mol6["UMAP_2"].median() + 2.2),
            fontsize=6.0, fontweight="bold", color="#B30000",
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.4),
            arrowprops=dict(arrowstyle="->", color="#B30000", lw=0.8,
                            shrinkA=1, shrinkB=1.5))
cb = fig.colorbar(hb, cax=axc)
cb.set_label("norm_score", fontsize=6, labelpad=1.5)
cb.ax.tick_params(labelsize=6.0, length=2, width=0.5, pad=1.5)
cb.set_ticks([round(-vB, 2), 0, round(vB, 2)])
cb.outline.set_linewidth(0.5)
frame(ax)
fig.savefig(os.path.join(OUT, "fig6b_umap.pdf"))
plt.close(fig)
print(f"B done (vB={vB:.2f}, bins={len(vals)})")'''))

# =============================================================================
# Section 4 — Panel C (Manhattan)
# =============================================================================
cells.append(md("""---

## 4. Panel C — compact MS Manhattan (10.40 × 4.60 cm)

**Data chain**: 178 genome-wide-significant MS lead SNPs (GWAS Catalog,
PMID 31604244), plotted per chromosome with a 1.08× inter-chromosome gap;
y = `−log10 P`, capped at 80. **Core parameters**: hit window ±500 kb (labels
only for genes in the `ms_locus_overlap_table.csv` overlap table); 16 labels
placed greedily with leader lines; genome-wide threshold 5e-8."""))
cells.append(code('''matplotlib.rcParams.update({"xtick.labelsize": 6.0, "ytick.labelsize": 6.5})
# ---- single-font true sub/superscript (Arial lacks U+2080/U+207B; mathtext
# ---- would fall back to DejaVu), identical to the project's arial_style -------
def sup_text(ax, x, y, base, sup, fontsize=None, transform=None, **kw):
    fs = fontsize or plt.rcParams["font.size"]
    tr = transform or ax.transData
    t1 = ax.text(x, y, base, fontsize=fs, va="baseline", transform=tr, **kw)
    fig = ax.figure
    fig.canvas.draw()
    bb = t1.get_window_extent(renderer=fig.canvas.get_renderer())
    xp, yp = tr.inverted().transform((bb.x1, bb.y0 + 0.74 * (bb.y1 - bb.y0)))
    t2 = ax.text(xp, yp, sup, fontsize=fs * 0.70, va="baseline", transform=tr, **kw)
    return t1, t2
def ten_pow_sup(ax, x, y, coef="5", exp="8", prefix="p = ", fontsize=7, **kw):
    return sup_text(ax, x, y, "%s%s \u00d7 10" % (prefix, coef), "\u2212%s" % exp,
                    fontsize=fontsize, **kw)
def vsub_label(ax, y, base, sub, tail="", fontsize=7, pad_pt=6.0, drop_pt=None,
              dx_pt=None, sub_scale=0.70, **kw):
    """Vertical-axis title with a true subscript: base + lowered small sub + tail."""
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bb = ax.get_window_extent(renderer=r)
    fs = fontsize
    gap = r.points_to_pixels(0.10 * fs)
    drop = r.points_to_pixels(0.30 * fs) if drop_pt is None else r.points_to_pixels(drop_pt)
    if dx_pt is not None:
        x = r.points_to_pixels(dx_pt) / bb.width
    else:
        x = (r.points_to_pixels(pad_pt) - bb.x0) / bb.width
    def adv(txt, size):
        t = ax.text(0, 0, txt, fontsize=size, rotation=90, ha="center", va="center",
                    transform=ax.transAxes, **kw)
        fig.canvas.draw()
        h = t.get_window_extent(renderer=r).height
        t.remove()
        return h
    a1, a2 = adv(base, fs), adv(sub, fs * sub_scale)
    a3 = adv(tail, fs) if tail else 0.0
    total = a1 + gap + a2 + (gap + a3 if tail else 0.0)
    yc = y - total / (2.0 * bb.height) + a1 / (2.0 * bb.height)
    out = [ax.text(x, yc, base, fontsize=fs, rotation=90, ha="center", va="center",
                   transform=ax.transAxes, **kw)]
    yc += (a1 + gap) / bb.height
    out.append(ax.text(x + drop / bb.width, yc, sub, fontsize=fs * sub_scale, rotation=90,
                       ha="center", va="center", transform=ax.transAxes, **kw))
    if tail:
        yc += (a2 + gap) / bb.height
        out.append(ax.text(x, yc, tail, fontsize=fs, rotation=90, ha="center", va="center",
                           transform=ax.transAxes, **kw))
    fig.canvas.draw()
    return out
W, H = 10.4, 4.6  # cm
chrs = [c for c in sorted(ms_snps["chr"].unique(), key=lambda x: (len(x), x))]
offsets, cum = {}, 0
for c in chrs:
    offsets[c] = cum
    cum += ms_snps.loc[ms_snps["chr"] == c, "pos"].max() * 1.08
ms_snps["x"] = ms_snps.apply(lambda r: offsets[r["chr"]] + r["pos"], axis=1)
LABEL = ["PHGDH", "PHLDB1", "PRR5L", "ZBTB38", "RRAS2", "SBF1", "TUBB4A",
         "CNP", "SEZ6L2", "BCAS1", "PLLP", "OLIG1", "APOE", "ATP1A1",
         "ADAMTS4", "CAMSAP2"]
PRIO = {g: 0 for g in ["PHGDH", "PHLDB1", "PRR5L", "ZBTB38"]}
PRIO.update({g: 1 for g in ["RRAS2", "SBF1"]})
PRIO.update({g: 2 for g in LABEL if g not in PRIO})
lab_map = {}
for _, h in hit.iterrows():
    if h["gene"] in LABEL:
        key = (str(h["chr"]), int(h["snp_pos"]))
        cur = lab_map.get(key)
        if cur is None or PRIO[h["gene"]] < PRIO[cur]:
            lab_map[key] = h["gene"]
YMAX = 80.0
fig, ax = plt.subplots(figsize=(W / CM, H / CM))
fig.subplots_adjust(left=0.085, right=0.995, top=0.925, bottom=0.13)
colors = plt.cm.tab20(np.linspace(0, 1, len(chrs)))
for i, c in enumerate(chrs):
    sub = ms_snps[ms_snps["chr"] == c]
    ax.scatter(sub["x"], -np.log10(sub["p"]), s=4, c=[colors[i % 20]], lw=0)
ax.set_ylim(0, YMAX); ax.set_xlim(0, cum)
def _ov(a, b, padx=1.2, pady=0.8):
    return not (a.x1 + padx < b.x0 or b.x1 + padx < a.x0 or
                a.y1 + pady < b.y0 or b.y1 + pady < a.y0)
items = []
for (ch, pos), gene in lab_map.items():
    r = ms_snps[(ms_snps["chr"] == ch) & (ms_snps["pos"] == pos)]
    if len(r):
        items.append((r["x"].iloc[0], -np.log10(r["p"].iloc[0]), gene))
items.sort()
fig.canvas.draw()
rend = fig.canvas.get_renderer()
inv = ax.transData.inverted()
ax_bb = ax.get_window_extent(rend)
line_px = rend.points_to_pixels(5.5 * 1.45)
line_data = abs(inv.transform((0, line_px))[1] - inv.transform((0, 0))[1])
px_per_data = ax_bb.width / cum
placed = []
for x, y, gene in items:
    ax.scatter([x], [y], s=11, c="red", lw=0, zorder=5)
    t = ax.text(x, y, gene, fontsize=6.0, color="red", ha="center",
                va="bottom", zorder=7, clip_on=False)
    k = 0
    while k < 14:
        t.set_position((x, y + 0.6 + k * line_data * 1.15))
        bb = t.get_window_extent(rend)
        if bb.x0 < ax_bb.x0 + 1:
            t.set_x((ax_bb.x0 + 1 - bb.x0) / px_per_data + x)
        elif bb.x1 > ax_bb.x1 - 1:
            t.set_x(x - (bb.x1 - ax_bb.x1 + 1) / px_per_data)
        bb = t.get_window_extent(rend)
        if all(not _ov(bb, p) for p in placed):
            break
        k += 1
    placed.append(t.get_window_extent(rend))
    y_lab = t.get_position()[1]
    if y_lab - y > line_data * 1.3:
        ax.plot([x, x], [y + 0.8, y_lab - 0.2], color="indianred", lw=0.4,
                zorder=4, clip_on=False)
thr = -np.log10(5e-8)
ax.axhline(thr, color="grey", lw=0.7, ls="--", zorder=2)
ten_pow_sup(ax, 0.30, 0.93, coef="5", exp="8", prefix="grey dashed: ",
            fontsize=8.0, transform=ax.transAxes, ha="left", color="#444444")
xt, xl = [], []
for c in chrs:
    sub = ms_snps[ms_snps["chr"] == c]
    xt.append(offsets[c] + sub["pos"].max() / 2); xl.append(c)
ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=5.5, rotation=90)
ax.tick_params(axis="x", length=2, pad=1)
ax.set_yticks([0, 10, 20, 30, 40, 50, 60, 70])
vsub_label(ax, 0.5, "\u2212log", "10", "(P)", fontsize=7.9)
ax.set_xlabel("")
for sp in ["top", "right"]:
    ax.spines[sp].set_visible(False)
fig.savefig(os.path.join(OUT, "fig6c_manhattan.pdf"))
plt.close(fig)
print(f"C done ({W} x {H} cm)")'''))

# =============================================================================
# Section 5 — Panel D (74-trait strip)
# =============================================================================
cells.append(md("""---

## 5. Panel D — MS selectivity across 74 GWAS gene sets (6.32 × 4.60 cm)

**Data chain**: 74 trait scDRS score files (`MISC`, one `norm_score` per cell
each) → 6-subtype mean per trait → row-wise z. **Core parameters**: the
pre-registered statistic is the **MOL6 z**; MS ranks 1/74 (z = +0.863),
one-sided empirical P = (1 + #{z ≥ z_MS})/74 = 0.0135. The strip shows all 74
traits (MS red; UKB traits grey, PASS disease traits steel), with three
reference traits labelled by phenotype."""))
cells.append(code('''matplotlib.rcParams.update({"font.size": 6.0, "axes.linewidth": 0.5})  # analyze_full74 values
FULL = MISC
H_FIG6D = 5.30                  # calibrated so the tight bbox lands on 6.32 x 4.60 cm
MS = "PASS_Multiple_sclerosis"
ORDER = ["MOL6", "MOL2", "MOL3/6", "MFOL", "MOL3", "COP/NFOL"]
ct = pd.read_csv(os.path.join(DATA, "cell_to_type.csv")).set_index("cell_id")["cell_type"]
files = {os.path.basename(f).replace(".score.gz", ""): f
         for f in glob.glob(f"{FULL}/*/*.score.gz") + glob.glob(f"{FULL}/*.score.gz")
         if "full_score" not in f and "/scores/" not in f}
print(f"found {len(files)} trait score files")
assert len(files) >= 2, "too few trait score files"
rows = {}
for t, f in sorted(files.items()):
    s = pd.read_csv(f, sep="\\t", index_col=0)["norm_score"]
    lab = ct.reindex(s.index)
    rows[t] = {c: s[lab == c].mean() for c in ORDER}
P = pd.DataFrame(rows).T[ORDER]
assert MS in P.index, "MS not scored yet"
Z = P.sub(P.mean(axis=1), axis=0).div(P.std(axis=1), axis=0)
T = Z["MOL6"].rename("MOL6_z").sort_values(ascending=False)
v = Z["MOL6"]
rank = int((v > v[MS]).sum()) + 1
p_emp = float((v >= v[MS]).sum()) / len(v)
print(f"MOL6_z(MS) = {v[MS]:+.3f} -> rank {rank}/{len(v)}, empirical P = {p_emp:.4f}")
short = lambda t: t.replace("PASS_", "").replace("UKB_460K.", "")
MS_C, GRAY = "#c0392b", "0.55"
fig, ax = plt.subplots(figsize=(6.8 / 2.54, H_FIG6D / 2.54), dpi=300)
x = np.arange(len(T)); y = T.values
cols = [MS_C if t == MS else (GRAY if t.startswith("UKB_") else "#7f8fa6") for t in T.index]
sz = np.where(np.asarray(T.index) == MS, 30, 9)
lo, hi = float(np.nanmin(y)), float(np.nanmax(y))
ax.set_ylim(lo - 0.55 * (hi - lo), hi + 0.42 * (hi - lo))
ax.scatter(x, y, c=cols, s=sz, zorder=3, edgecolor="none")
ax.axhline(0, color="0.8", lw=0.5, zorder=1)
i = list(T.index).index(MS)
ax.annotate(f"MS: rank {rank}/{len(T)}\\n{rank}/{len(T)} nominal (descriptive)", (i, y[i]),
            xytext=(i + 2.2, hi + 0.30 * (hi - lo)), fontsize=5.5, color=MS_C,
            ha="left", va="top", linespacing=1.25,
            arrowprops=dict(arrowstyle="-", color=MS_C, lw=0.5, shrinkA=2, shrinkB=3))
refs = {"PASS_Alzheimers_Jansen2019": 0.0, "UKB_460K.body_HEIGHTz": -0.30,
        "PASS_Schizophrenia_Pardinas2018": 0.0}
REFNAME = {"PASS_Alzheimers_Jansen2019": "Alzheimer disease",
           "UKB_460K.body_HEIGHTz": "Adult height (z-scored GWAS)",
           "PASS_Schizophrenia_Pardinas2018": "Schizophrenia"}
for t, dy in refs.items():
    if t in T.index:
        j = list(T.index).index(t)
        ax.annotate(REFNAME[t], (j, y[j]), xytext=(j - 0.6, y[j] - 0.18 + dy),
                    fontsize=5.5, color="0.4", ha="left", va="top", rotation=32)
ax.set_xticks([]); ax.set_xlim(-2.5, len(T) + 3)
ax.set_yticks([-2, -1, 0, 1])
ax.set_ylabel("MOL6 selectivity\\n(row z of mean norm_score)", fontsize=5.8)
ax.set_title(f"{len(T)} GWAS gene sets, same cells", fontsize=6.4, pad=4)
ax.tick_params(labelsize=5.5, width=0.4, length=1.8)
for sp in ["top", "right"]:
    ax.spines[sp].set_visible(False)
fig.savefig(f"{OUT}/fig6d_74trait_strip.pdf", bbox_inches="tight", pad_inches=0.03)
plt.close(fig)
import fitz
d = fitz.open(f"{OUT}/fig6d_74trait_strip.pdf")[0]
print(f"fig6d_74trait_strip.pdf  {d.rect.width/28.3465:.2f} x {d.rect.height/28.3465:.2f} cm  rank={rank}/{len(T)} P={p_emp:.4f}")'''))

nb = nbf.v4.new_notebook(metadata={
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.8"}})
nb["cells"] = cells
out = os.path.join(OUT, "Figure_6.ipynb")
nbf.write(nb, out)
print("wrote", out, len(cells), "cells")
