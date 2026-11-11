#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build Figure_S4.ipynb — self-contained reproduction of the journal
Supplementary Figure 4 (Fig S4), organised panel-by-panel exactly like
Figure_4.ipynb. Run:  python3 build_figS4_notebook.py
"""
import os
import nbformat as nbf

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figS4")


def md(s):
    return nbf.v4.new_markdown_cell(s)


def code(s):
    return nbf.v4.new_code_cell(s)


cells = []

# =============================================================================
# Section 0 — setup, fonts, helpers
# =============================================================================
cells.append(md("""# Figure S4 — time axes, peak-to-gene links and gene-set structure behind Fig 4

This notebook reproduces the four panels of Supplementary Figure 4 **from the
original data objects** (the only read-back step parses an intermediate PDF
that this notebook itself produces — the MultiVelo model-1 heat map raster).

* **A** — MultiVelo latent time on the Fig 4A embedding + cell-type composition
  of the latent-time axis.
* **B** — the 15,607 peak-to-gene links (r, distance, peaks per gene).
* **C** — gene-set funnel (1,809 → 372 / 114, overlap 32, not nested), the
  model-1 / model-2 heat maps, and two direction examples (Nfasc, Ndrg1).
* **D** — three-dimensional s–u–c phase portraits for eight genes.

Each section ends with its own panel PDF."""))
cells.append(code('''import os, subprocess, sys, glob
import numpy as np, pandas as pd, scipy, matplotlib
import anndata as ad
from matplotlib import font_manager as fm
# Arial-first font registration (self-contained; falls back to Liberation Sans
# whose metrics match Arial). Note the file names are ARIAL.TTF (uppercase .TTF).
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
# rcParams identical to the project's fig_common.py (import-time use_arial values)
matplotlib.rcParams.update({"font.size": 7.0,
                             "axes.titlesize": 7.0, "axes.labelsize": 7.0,
                             "xtick.labelsize": 6.0, "ytick.labelsize": 6.0,
                             "legend.fontsize": 6.0, "axes.titleweight": "normal",
                             "figure.dpi": 120, "savefig.dpi": 600,
                             "mathtext.fontset": "custom",
                             "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                             "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
                             "xtick.major.size": 2.0, "ytick.major.size": 2.0,
                             "xtick.minor.size": 1.0, "ytick.minor.size": 1.0,
                             "xtick.major.pad": 1.6, "ytick.major.pad": 1.6,
                             "lines.linewidth": 0.8, "patch.linewidth": 0.5,
                             "axes.linewidth": 0.6, "pdf.fonttype": 42})
print("font family:", _FAMILY)
# R interpreter (only used for the blueYellow ramp, same as Figure 4)
RSCRIPT = os.environ.get("RSCRIPT", "Rscript")
def rrun(code_str, script="panel_s4.R", env=None):
    header = 'DATA <- Sys.getenv("DATA"); OUT <- Sys.getenv("OUT")\\n'
    with open(script, "w") as f:
        f.write(header + code_str)
    e = dict(os.environ)
    if env:
        e.update(env)
    subprocess.run([RSCRIPT, script], check=True, env=e)
DATA = os.environ.get("DATA", os.path.join(os.environ.get("DSCI_WORK_ROOT", "."), "FY", "oligo"))
OUT  = os.environ.get("FIG4_OUT", "fig4_out")
os.makedirs(OUT, exist_ok=True)
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
plt.rcParams.update({"font.family": _FAMILY, "font.size": 7.0,
                     "axes.titlesize": 7.0, "axes.labelsize": 7.0,
                     "xtick.labelsize": 6.0, "ytick.labelsize": 6.0,
                     "legend.fontsize": 6.0, "axes.titleweight": "normal",
                     "figure.dpi": 120, "savefig.dpi": 600,
                     "mathtext.fontset": "custom",
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                     "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
                     "xtick.major.size": 2.0, "ytick.major.size": 2.0,
                     "xtick.minor.size": 1.0, "ytick.minor.size": 1.0,
                     "xtick.major.pad": 1.6, "ytick.major.pad": 1.6,
                     "lines.linewidth": 0.8, "patch.linewidth": 0.5,
                     "axes.linewidth": 0.6})
CM = 1 / 2.54
# ---- upstream data (all drawn from the original analysis objects) ----
D = ad.read_h5ad(f"{DATA}/multivelo/results/multivelo_dynamical.h5ad")
OBS = D.obs[["cell_type", "latent_time"]].copy()
OBS["umap1"] = D.obsm["X_umap"][:, 0]
OBS["umap2"] = D.obsm["X_umap"][:, 1]
CT    = list(D.obs.cell_type.cat.categories)
CT_COLS = [D.uns["cell_type_colors"][i] for i in range(len(CT))]
vt = D.var[[c for c in D.var.columns if c.startswith("fit_")] + ["velo_s_genes"]].copy()
vt.insert(0, "gene", D.var_names)
dec = pd.read_csv(f"{DATA}/multivelo/results/decoupling.tsv", sep="\\t", index_col=0)
dir_df = pd.read_csv(f"{DATA}/multivelo/results/decoupling_direction.tsv", sep="\\t")
CAT = dict(zip(dir_df.gene, dir_df.cat))
SHORT = {"priming (chromatin-leading)": "priming",
         "lagging shutdown (chromatin-open)": "lagging",
         "coupled": "coupled", "decoupled-other": "decoupled"}
def _cat(g):
    return SHORT.get(CAT.get(g), "not high-conf.")
TRAJ8 = ["Anln", "Mbp", "Mal", "Car2", "Tnr", "Plp1", "Ptgds", "Tcf7l2"]
def _rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))
SC = np.array([_rgb(h) for h in CT_COLS])[np.searchsorted(np.array(CT),
                                                         D.obs.cell_type.values)]
# ---- blueYellow ramp (same palette as Fig 4A) exported by R ----
rrun(r"""
suppressPackageStartupMessages({library(scMEGA); library(ArchR)})
pal <- paletteContinuous(set="blueYellow", n=100)
write.table(pal, file.path(OUT, "blueYellow.csv"),
            sep=",", col.names=FALSE, row.names=FALSE)
cat("ramp rows:", length(pal), "\\n")
""", "ramp.R", env={"DATA": DATA, "OUT": OUT})
pal = pd.read_csv(f"{OUT}/blueYellow.csv", header=None).values[:, 0]
RAMP = np.array([(int(h[1:3], 16) / 255, int(h[3:5], 16) / 255,
                  int(h[5:7], 16) / 255) for h in pal])
RAMPCM = ListedColormap(RAMP)'''))

# -----------------------------------------------------------------------------
# Versions (Python + R) printed up front for reproducibility
# -----------------------------------------------------------------------------
cells.append(md("""### 0.1 Input files

| Object | Produced by | Used in |
|---|---|---|
| `multivelo_dynamical.h5ad` | MultiVelo `recover_dynamics_chrom` (5977 × 1809) | A–D |
| `multivelo/results/decoupling.tsv` + `decoupling_direction.tsv` | LRT refit | C |
| `scmega_p2g.csv` + `scmega_p2g_S1markers.csv` | `SelectGenes()` (15,607 links) | B |
| `multivelo/plots/heatmap__model1__1798genes.pdf` | MultiVelo model-1 heat map (wordless raster) | C2 |
| blueYellow ramp | `paletteContinuous("blueYellow")` (R) | A |

### 0.2 Software versions

Prints every relevant package version (Python for panels A–D, R for the ramp)."""))
cells.append(code('''print("==== Python packages ====")
for m in (np, pd, scipy, matplotlib, ad):
    print(f"  {m.__name__:<12} {m.__version__}")
print("==== R (blueYellow ramp) ====")
rrun(r"""
ip <- installed.packages()
for (p in c("scMEGA","ArchR","Seurat","Signac")) {
  if (p %in% rownames(ip)) cat(sprintf("  %-22s %s\\n", p, ip[p, "Version"]))
}
cat("  R ", R.version.string, "\\n", sep="")
""", "versions.R", env={"DATA": DATA, "OUT": OUT})'''))

# =============================================================================
# Section 1 — Panel A (latent-time UMAP + composition)
# =============================================================================
cells.append(md("""---

## 1. Panel A — latent time + cell-type composition (17.00 × 4.20 cm)

**Data chain** (all in Python):

1. `multivelo_dynamical.h5ad` — `obs[latent_time]` is the learned dynamical
   time (0–1) of each cell, on the same RNA embedding as Fig 4A
   (`obsm[X_umap]`, identical to `rna_umap.tsv`).
2. Left: the 5,977 cells are drawn on that embedding, coloured by latent time
   with the same blueYellow ramp as Fig 4A, so the two time axes can be
   compared directly.
3. Right: latent time is binned into 20 equal-width bins; the stacked bars give
   the six-subtype composition of each bin — this is exactly the horizontal
   axis of the Fig S4C heat maps.

> Latent time is only weakly correlated with the scMEGA pseudotime of Fig 4A
> (cell-level Spearman ρ = 0.20) and does not order the mature subtypes
> (MOL3/6 0.747, MOL2 0.665, MOL6 0.659) — the two axes are not interchangeable."""))
cells.append(code('''from matplotlib.patches import Rectangle
fig = plt.figure(figsize=(17.00 * CM, 4.20 * CM))
rect = [0.055, 0.115, 0.330, 0.700]
ax = fig.add_axes(rect)
ax.scatter(OBS.umap1, OBS.umap2, s=0.75, linewidths=0, c=OBS.latent_time,
           cmap=RAMPCM, vmin=0, vmax=1, rasterized=True)
ax.set_aspect("equal"); ax.axis("off")
ax.set_title("MultiVelo latent time on the Fig 4A embedding", fontsize=6.4, pad=3)
cb = fig.add_axes([rect[0] + rect[2] + 0.012, rect[1] + 0.06, 0.0135, rect[3] - 0.12])
cb.imshow(np.linspace(0, 1, 256).reshape(-1, 1), aspect="auto", origin="lower",
          cmap=RAMPCM, extent=[0, 1, 0, 1])
cb.yaxis.set_ticks_position("right"); cb.set_xticks([])
cb.set_yticks([0, 0.5, 1]); cb.set_yticklabels(["0", "0.5", "1"])
cb.set_ylim(0, 1)
cb.tick_params(labelsize=5.5, width=0.4, length=1.3, pad=1.0)
cb.yaxis.set_label_position("right")
cb.set_ylabel("Latent time", fontsize=5.8, labelpad=1.5)
for s_ in cb.spines.values():
    s_.set_linewidth(0.4)
ax = fig.add_axes([0.505, 0.150, 0.460, 0.560])
lt = OBS.latent_time.values
nb = 20
ltb = np.clip((lt * nb).astype(int), 0, nb - 1)
ctv = OBS.cell_type.values
comp = np.zeros((len(CT), nb))
for j in range(nb):
    m = ltb == j
    for k, t in enumerate(CT):
        comp[k, j] = (ctv[m] == t).mean() if m.sum() else 0.0
bottom = np.zeros(nb)
for k, t in enumerate(CT):
    ax.bar(range(nb), comp[k], bottom=bottom, color=CT_COLS[k], width=0.99,
           linewidth=0)
    bottom += comp[k]
ax.set_xticks(range(0, nb, 2))
ax.set_xticklabels([str(i * 5) for i in range(0, nb, 2)], fontsize=5.5)
ax.set_xlabel("latent-time bin (heatmap x axis, %)", fontsize=5.8, labelpad=2)
ax.set_ylabel("cell fraction", fontsize=5.8)
ax.set_ylim(0, 1.0); ax.set_xlim(-0.6, nb - 0.4)
ax.tick_params(labelsize=5.6, width=0.4, length=1.4)
ax.set_title("What the latent-time axis is made of", fontsize=6.4, pad=3)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
axl = fig.add_axes([0.505, 0.755, 0.460, 0.115]); axl.axis("off")
axl.legend([Rectangle((0, 0), 1, 1, fc=CT_COLS[k]) for k in range(len(CT))],
           list(CT), loc="center", ncol=6, frameon=False, handletextpad=0.3,
           columnspacing=0.85, fontsize=5.5, handlelength=0.85,
           handleheight=0.85, borderaxespad=0.0)
fig.savefig(f"{OUT}/figS4v3a_umap.pdf")
print("A done")'''))

# =============================================================================
# Section 2 — Panel B (peak-to-gene links)
# =============================================================================
cells.append(md("""---

## 2. Panel B — peak-to-gene links (17.00 × 3.40 cm)

**Data chain**: `scmega_p2g.csv` is the *kept* peak-to-gene link set produced
by scMEGA (`SelectGenes()`: correlation r ≥ 0.40, distance ≤ 250 kb,
FDR ≤ 1e-4 — 15,607 links, 2,486 genes, 13,010 peaks). Three histograms:
link correlation (floor 0.40, median ≈ 0.56), peak–TSS distance (cap 250 kb,
median ≈ 107 kb), and peaks per gene (median 5, max 27), with the 331 S1
marker genes overlaid in orange. The axes are truncated at the link-selection
cut-offs, so no enrichment statement is made."""))
cells.append(code('''from matplotlib.gridspec import GridSpec
p = pd.read_csv(f"{DATA}/scmega_p2g.csv")
s1 = pd.read_csv(f"{DATA}/scmega_p2g_S1markers.csv")["gene"].unique()
fig = plt.figure(figsize=(17.00 * CM, 3.40 * CM))
gs = GridSpec(1, 3, wspace=0.34, left=0.055, right=0.985, top=0.805, bottom=0.185)
ax = fig.add_subplot(gs[0, 0])
ax.hist(p["Correlation"], bins=np.arange(0.40, 1.001, 0.04), color="#4C72B0",
        edgecolor="white", lw=0.45)
ax.axvline(p["Correlation"].min(), color="crimson", lw=0.9)
ax.set_xlabel("peak-to-gene r", fontsize=5.8, labelpad=1)
ax.set_ylabel("links", fontsize=5.8)
ax.set_title("r floor %.2f\\n(median %.2f)" % (p["Correlation"].min(),
                                              p["Correlation"].median()),
             fontsize=6.2, pad=3)
ax.tick_params(labelsize=5.8, width=0.4, length=1.4)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax = fig.add_subplot(gs[0, 1])
ax.hist(p["distance"] / 1000.0, bins=np.arange(0, 260, 10), color="#4C72B0",
        edgecolor="white", lw=0.45)
ax.axvline(p["distance"].max() / 1000.0, color="crimson", lw=0.9)
ax.set_xlabel("distance to TSS (kb)", fontsize=5.8, labelpad=1)
ax.set_title("cap 250 kb\\n(median %.0f kb)" % (p["distance"].median() / 1000),
             fontsize=6.2, pad=3)
ax.tick_params(labelsize=5.8, width=0.4, length=1.4)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax = fig.add_subplot(gs[0, 2])
cnt = p.groupby("gene")["peak"].nunique()
b = np.arange(0, 28, 2)
ax.hist(cnt.values, bins=b, color="#4C72B0", edgecolor="white", lw=0.45,
        label="all linked genes")
ax.hist(cnt[cnt.index.isin(s1)].values, bins=b, color="#DD8452",
        edgecolor="white", lw=0.45, label="S1 markers")
ax.set_xlabel("peaks linked per gene", fontsize=5.8, labelpad=1)
ax.set_ylabel("genes", fontsize=5.8)
ax.set_title("%d genes, median %.0f\\n(max %d)" % (cnt.size, cnt.median(), cnt.max()),
             fontsize=6.2, pad=3)
ax.legend(fontsize=5.6, frameon=False, loc="upper right")
ax.tick_params(labelsize=5.8, width=0.4, length=1.4)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.savefig(f"{OUT}/figS4v3b_p2g.pdf")
print("B done")'''))

# =============================================================================
# Section 3 — Panel C (funnel + model heatmaps + direction examples)
# =============================================================================
cells.append(md("""---

## 3. Panel C — gene-set structure and programmes (17.00 × 7.55 + blocks)

**C1 — funnel** (17.00 × 2.30 cm): the 1,809 fitted genes give two *partially
overlapping* subsets — 372 LRT-decoupled genes (FDR < 0.05, Fig 4C volcano)
and 114 high-confidence velocity genes (82 coupled · 31 priming · 1 lagging
shutdown) — which intersect in only 32 genes; the funnel therefore **forks**
(372 and 114 both come directly from 1,809), it is not a chain.

**C2 — model heat maps** (17.00 × 7.55 cm): model-1 (1,798 genes, 99.4%) uses
the raster produced by MultiVelo (`plots/heatmap__model1__1798genes.pdf` — a
wordless, row-ordered image; the four k-means programmes are C0 n=576 early,
C1 n=491 MFOL onset, C2 n=433 mature sustained, C3 n=298 mid-transient) next
to the model-2 heat map (11 genes, 0.6%, including Plp1 and Tcf7l2) redrawn
here from the h5ad `Ms` layer, z-scored per gene and binned over latent time
(600 bins) with the same cell-type top bar.

**C3 — direction examples** (17.00 × 4.20 cm): the single lagging-shutdown
gene Nfasc (model 2, LRT 700, FDR 2e-153; t_s50 0.4 < t_c50 9.8) and the
coupled gene Ndrg1 (model 1, LRT −6,755, FDR 1), with t_c50 / t_s50 dashed."""))
cells.append(code('''from matplotlib.patches import Rectangle as Rect
# ---- C1 funnel ----
dd = pd.read_csv(f"{DATA}/multivelo/results/decoupling.tsv", sep="\\t")
dg = pd.read_csv(f"{DATA}/multivelo/results/decoupling_direction.tsv", sep="\\t")
n_fit, n_dec, n_hi = len(dd), int((dd.FDR < 0.05).sum()), len(dg)
both = len(set(dd.loc[dd.FDR < 0.05, "gene"]) & set(dg.gene))
cnt = dg.cat.value_counts()
W, H = 17.00, 2.30
fig = plt.figure(figsize=(W * CM, H * CM))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0]); ax.axis("off")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
HP = H * CM * 72.0
dy = lambda pt: pt / HP
def box(cx, cy, w, big, lines, col, hpt=40):
    h = dy(hpt)
    ax.add_patch(Rect((cx - w / 2, cy - h / 2), w, h, fc=col, ec="0.35",
                      lw=0.5, zorder=1))
    ax.text(cx - w / 2 + dy(6), cy + dy(11), big, fontsize=8.0, weight="bold",
            ha="left", va="center", color="black", zorder=2)
    y = cy - dy(1)
    for txt, fs, wt, cc in lines:
        ax.text(cx - w / 2 + dy(6), y, txt, fontsize=fs, ha="left",
                va="center", color=cc, weight=wt, zorder=2)
        y -= dy(7.2)
box(0.105, 0.50, 0.190, f"{n_fit:,}",
    [("fitted by MultiVelo", 5.5, "normal", "0.35")], "#DCE6F1")
yfork = 0.93
ax.plot([0.200, 0.2225], [0.50, 0.50], color="0.45", lw=0.7, zorder=0)
ax.plot([0.2225, 0.2225], [0.50, yfork], color="0.45", lw=0.7, zorder=0)
ax.plot([0.2225, 0.7575], [yfork, yfork], color="0.45", lw=0.7, zorder=0)
ax.plot([0.7575, 0.7575], [yfork, 0.815], color="0.45", lw=0.7, zorder=0)
ax.annotate("", xy=(0.243, 0.50), xytext=(0.2225, 0.50),
            arrowprops=dict(arrowstyle="-|>", color="0.45", lw=0.7), zorder=0)
ax.annotate("", xy=(0.7575, 0.788), xytext=(0.7575, 0.815),
            arrowprops=dict(arrowstyle="-|>", color="0.45", lw=0.7), zorder=0)
box(0.350, 0.50, 0.205, f"{n_dec}",
    [("decoupled", 5.6, "bold", "black"),
     ("FDR < 0.05", 5.5, "normal", "0.35"),
     ("Fig 4C volcano", 5.5, "normal", "0.35")], "#BFD3E6")
box(0.7575, 0.50, 0.415, f"{n_hi}",
    [("high-confidence velocity genes", 5.6, "bold", "black"),
     ("%d coupled \\u00b7 %d priming \\u00b7 %d lagging"
      % (cnt.get("coupled", 0), cnt.get("priming (chromatin-leading)", 0),
         cnt.get("lagging shutdown (chromatin-open)", 0)),
      5.5, "normal", "black"),
     ("Fig 4C right panel, and the two curves below", 5.5, "normal", "0.35")],
    "#9DB9D5")
ax.text(0.5035, 0.585, f"overlap {both}", fontsize=5.6, weight="bold",
        ha="center", va="center", color="black")
ax.text(0.5035, 0.385, "not nested", fontsize=5.5, ha="center", va="center",
        color="0.3")
fig.savefig(f"{OUT}/figS4v3c_funnel.pdf")
print("C1 funnel done")'''))

# ---- C2 model heatmaps (fitz composition) ----
cells.append(md("""**C2 model heat maps** — the model-1 block is the raster
produced by MultiVelo (`plots/heatmap__model1__1798genes.pdf`, wordless image,
rows in k-means programme order, columns = 5,977 cells by latent time); the
model-2 block is redrawn here from the `Ms` layer of the h5ad, z-scored per
gene, binned over latent time, with a cell-type top bar (the same latent-time
order as the model-1 map)."""))
cells.append(code('''import fitz, os
PT = 72 / 2.54
# ---- embed Arial for fitz text (helv/hebo -> Arial-Regular/Arial-Bold) ----
_AR = os.path.join(os.environ.get("FONT_DIR", os.path.join(os.environ.get("DSCI_WORK_ROOT", "."), "fonts")), "ARIAL.TTF")
_AB = os.path.join(os.environ.get("FONT_DIR", os.path.join(os.environ.get("DSCI_WORK_ROOT", "."), "fonts")), "ARIALBD.TTF")
_ins = fitz.Page.insert_text
def _resolve(fontname, fontfile):
    if fontfile is None and fontname in ("helv", "hebo", "Helvetica"):
        bold = fontname == "hebo"
        return ("Arial-Bold" if bold else "Arial-Regular"), (_AB if bold else _AR)
    return fontname, fontfile
def _insert(self, point, text, fontname="helv", fontfile=None, **kw):
    fontname, fontfile = _resolve(fontname, fontfile)
    return _ins(self, point, text, fontname=fontname, fontfile=fontfile, **kw)
fitz.Page.insert_text = _insert
_fcache = {}
def _glen(text, fontname="helv", fontsize=11, **kw):
    bold = fontname == "hebo"
    if bold not in _fcache:
        _fcache[bold] = fitz.Font(fontfile=(_AB if bold else _AR))
    return _fcache[bold].text_length(text, fontsize=fontsize)
fitz.get_text_length = _glen
# the model-2 block is a *pre-existing* project artefact drawn with the
# fig_common rcParams (tick pad 1.2, size 1.6); restore them for this block.
plt.rcParams.update({"font.size": 6.0, "axes.titlesize": 6.2,
                     "axes.labelsize": 6.0, "xtick.labelsize": 5.5,
                     "ytick.labelsize": 5.5, "legend.fontsize": 5.5,
                     "xtick.major.size": 1.6, "ytick.major.size": 1.6,
                     "xtick.major.pad": 1.2, "ytick.major.pad": 1.2,
                     "axes.linewidth": 0.5})
# ---- model-2 heat map (redrawn from h5ad Ms layer) ----
Ms2 = D[:, D.var["fit_model"].values == 2].layers["Ms"]
Ms2 = np.asarray(Ms2.todense() if hasattr(Ms2, "todense") else Ms2)  # 5977 x 11
genes2 = list(D.var_names[D.var["fit_model"].values == 2])
order = np.argsort(D.obs["latent_time"].values)
z = (Ms2 - Ms2.mean(0, keepdims=True)) / (Ms2.std(0, keepdims=True) + 1e-9)
nb = 600
bins = np.array_split(order, nb)
M2 = np.column_stack([z[b].mean(0) for b in bins])  # 11 x 600
strip = np.asarray([np.bincount([CT.index(x) for x in D.obs.cell_type.values[b]],
                                minlength=len(CT)).argmax() for b in bins])
fig2 = plt.figure(figsize=(8.35 * CM, 5.60 * CM))
gs2 = fig2.add_gridspec(2, 1, height_ratios=[0.06, 1.0], hspace=0.06,
                        left=0.135, right=0.905, top=0.965, bottom=0.115)
axs = fig2.add_subplot(gs2[0])
axs.imshow(strip.reshape(1, -1), aspect="auto", cmap=ListedColormap(CT_COLS),
           vmin=-0.5, vmax=len(CT) - 0.5)
axs.set_xticks([]); axs.set_yticks([])
ax = fig2.add_subplot(gs2[1])
ax.imshow(M2, aspect="auto", cmap="viridis", vmin=-2, vmax=2, interpolation="nearest")
ax.set_yticks(range(len(genes2)))
ax.set_yticklabels(genes2, fontsize=5.6, style="italic")
ax.set_xticks([0, nb / 2, nb - 1])
ax.set_xticklabels(["COP/NFOL", "latent time", "MOL6"], fontsize=5.5)
ax.set_xlim(-0.5, nb - 0.5)
cb = fig2.add_axes([0.935, 0.32, 0.014, 0.36])
cb.imshow(np.linspace(-2, 2, 64).reshape(-1, 1), aspect="auto", cmap="viridis",
          origin="lower", extent=[0, 1, -2, 2])
cb.set_ylim(-2, 2); cb.set_xlim(0, 1); cb.set_xticks([])
cb.set_yticks([-2, 0, 2]); cb.set_yticklabels(["-2", "0", "+2"], fontsize=5.5)
cb.yaxis.set_ticks_position("right"); cb.yaxis.set_label_position("right")
cb.tick_params(width=0.4, length=1.2, pad=1.2)
for s_ in cb.spines.values():
    s_.set_linewidth(0.4)
fig2.savefig(f"{OUT}/figS4a2_model2_heatmap.pdf")
plt.close(fig2)
# restore the BASE rcParams used by the other S4 panels
plt.rcParams.update({"font.size": 7.0, "axes.titlesize": 7.0,
                     "axes.labelsize": 7.0, "xtick.labelsize": 6.0,
                     "ytick.labelsize": 6.0, "legend.fontsize": 6.0,
                     "xtick.major.size": 2.0, "ytick.major.size": 2.0,
                     "xtick.major.pad": 1.6, "ytick.major.pad": 1.6,
                     "axes.linewidth": 0.6})
# ---- compose the C2 block (model-1 raster + model-2 redraw + captions) ----
W2, H2 = 17.0, 7.55
src = fitz.open(f"{DATA}/multivelo/plots/heatmap__model1__1798genes.pdf")
sp = src[0]
sw, sh = float(sp.rect.width), float(sp.rect.height)
hm_w = 8.35
hm_h = sh / sw * hm_w
x0, y0 = 0.02, 1.15
x2 = W2 - x0 - 8.35
doc = fitz.open()
page = doc.new_page(width=W2 * PT, height=H2 * PT)
page.insert_text((46, 8.6), "Dynamical models: gene programmes along the lineage",
                 fontname="hebo", fontsize=7.0)
page.insert_text((x0 * PT, 18.5), "model 1 \\u00b7 1798 genes (99.4%)",
                 fontname="helv", fontsize=5.6)
page.insert_text((x2 * PT, 18.5), "model 2 \\u00b7 11 genes (0.6%)",
                 fontname="helv", fontsize=5.6)
lx = x0 * PT
for i, t in enumerate(CT):
    r = fitz.Rect(lx, 24.0, lx + 5.0, 29.0)
    page.draw_rect(r, color=None, fill=_rgb(CT_COLS[i]), width=0)
    page.insert_text((lx + 6.5, 28.4), t, fontname="helv", fontsize=5.5)
    lx += 6.5 + fitz.get_text_length(t, fontname="helv", fontsize=5.5) + 7
page.show_pdf_page(fitz.Rect(x0 * PT, y0 * PT, (x0 + hm_w) * PT, (y0 + hm_h) * PT),
                   src, 0)
src.close()
m2 = fitz.open(f"{OUT}/figS4a2_model2_heatmap.pdf")
p2 = m2[0]
w2p, h2p = float(p2.rect.width) / PT, float(p2.rect.height) / PT
x2p = W2 - x0 - w2p
page.show_pdf_page(fitz.Rect(x2p * PT, y0 * PT, (x2p + w2p) * PT, (y0 + h2p) * PT),
                   m2, 0)
m2.close()
page.insert_text((x0 * PT, (y0 + max(hm_h, h2p) + 0.32) * PT),
                 "Rows: smoothed spliced RNA (Ms), z-scored per gene (1798 model-1 "
                 "rows, in k-means programme order).  Columns: 5,977 cells "
                 "ordered by latent time (COP/NFOL \\u2192 MOL6); top bar = cell type.",
                 fontname="helv", fontsize=5.5, color=(0.35, 0.35, 0.35))
page.insert_text((x0 * PT, (y0 + max(hm_h, h2p) + 0.62) * PT),
                 "model 1 = chromatin-led shutdown (RNA falls while chromatin "
                 "closes);  model 2 = transcription-led, chromatin stays open "
                 "(includes Plp1, Tcf7l2).",
                 fontname="helv", fontsize=5.5, color=(0.35, 0.35, 0.35))
doc.save(f"{OUT}/figS4a_model_heatmaps.pdf", garbage=4, deflate=True)
doc.close()
print("C2 model heatmaps done  (%.2f x %.2f cm; model1 block %.2f x %.2f cm)"
      % (W2, H2, hm_w, hm_h))'''))

# ---- C3 direction examples ----
cells.append(md("""**C3 direction examples** — the c(t) / s(t) curves are
rebuilt from each gene's fitted parameters (`fit_alpha_c`, `fit_alpha`,
`fit_beta`, `fit_gamma`, `fit_scale_cc`, `fit_t_sw1..3`) using MultiVelo's own
`predict_exp`, exactly as the dynamical model defines them (400 points over
0–20)."""))
cells.append(code('''from multivelo.dynamical_chrom_func import predict_exp
def gene_traj(gene, T=20.0, nt=400):
    v = vt.set_index("gene").loc[gene]
    a_c, a, b, g_, scc = (float(v.fit_alpha_c), float(v.fit_alpha),
                          float(v.fit_beta), float(v.fit_gamma),
                          float(v.fit_scale_cc))
    sw = [float(v.fit_t_sw1), float(v.fit_t_sw2), float(v.fit_t_sw3)]
    model = int(v.fit_model)
    bounds = [0.0] + [min(max(x, 0.0), T) for x in sw] + [T]
    if model == 0:
        segs = [(a, True), (a, False), (a, False), (0.0, False)]
    elif model == 1:
        segs = [(a, True), (a, True), (a, False), (0.0, False)]
    else:
        segs = [(a, True), (a, True), (0.0, True), (0.0, False)]
    ts, cs, ss = [], [], []
    c0 = u0 = s0 = 0.0
    for (t0, t1), (al, op) in zip(zip(bounds[:-1], bounds[1:]), segs):
        tau = np.linspace(0, max(t1 - t0, 1e-6), max(int(nt * (t1 - t0) / T), 2))
        pred = predict_exp(tau, c0, u0, s0, a_c, al, b, g_, scale_cc=scc,
                           pred_r=True, chrom_open=op)
        c0, u0, s0 = pred[-1]
        ts.append(t0 + tau); cs.append(pred[:, 0]); ss.append(pred[:, 2])
    t = np.concatenate(ts); c = np.concatenate(cs); s = np.concatenate(ss)
    keep = np.r_[True, np.diff(t) > 1e-9]
    t, c, s = t[keep], c[keep], s[keep]
    i_c = np.argmax(c >= 0.5 * c.max()); i_s = np.argmax(s >= 0.5 * s.max())
    tc50 = t[i_c] if c[i_c] >= 0.5 * c.max() else t[-1]
    ts50 = t[i_s] if s[i_s] >= 0.5 * s.max() else t[-1]
    return t, c, s, tc50, ts50, model, sw
def _fmt(x, pat="{:.0f}"):
    try:
        return pat.format(x)
    except Exception:
        return "-"
fig = plt.figure(figsize=(17.00 * CM, 4.20 * CM))
gs = fig.add_gridspec(1, 2, wspace=0.24, left=0.062, right=0.985,
                      top=0.815, bottom=0.155)
for k, g in enumerate(["Nfasc", "Ndrg1"]):
    ax = fig.add_subplot(gs[0, k])
    t, c, s, tc50, ts50, model, sw = gene_traj(g)
    ax.plot(t, c / c.max(), color="steelblue", lw=1.1,
            label="c (chromatin, norm)")
    ax.plot(t, s / s.max(), color="crimson", lw=1.1,
            label="s (spliced RNA, norm)")
    for x in sw:
        if 0 <= x <= 20:
            ax.axvline(x, ls=":", c="0.55", lw=0.6)
    ax.axvline(tc50, ls="--", c="steelblue", lw=0.7)
    ax.axvline(ts50, ls="--", c="crimson", lw=0.7)
    row = dec.loc[g]
    ax.set_title("%s \\u00b7 %s\\nmodel %d \\u00b7 LRT %s \\u00b7 FDR %s\\n"
                 "tc50 %.1f \\u00b7 ts50 %.1f"
                 % (g, _cat(g), model, _fmt(row.LRT), _fmt(row.FDR, "{:.0g}"),
                    tc50, ts50),
                 fontsize=5.5, pad=2, linespacing=1.3)
    ax.set_xlabel("trajectory time", fontsize=5.8, labelpad=1)
    ax.set_ylabel("normalized value", fontsize=5.8)
    ax.set_xlim(0, 20); ax.set_ylim(0, 1.12)
    ax.set_xticks([0, 5, 10, 15, 20]); ax.set_yticks([0, 0.5, 1.0])
    ax.tick_params(labelsize=5.6, width=0.4, length=1.4)
    if k == 0:
        ax.legend(fontsize=5.5, frameon=False, loc="upper right")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.savefig(f"{OUT}/figS4v3c_dec.pdf")
print("C3 dec done")'''))

# =============================================================================
# Section 4 — Panel D (3D phase portraits)
# =============================================================================

cells.append(md("""---

## 4. Panel D — 3D s–u–c phase portraits (17.00 × 8.40 cm)

**Data chain**: for each of the eight genes (Anln, Mbp, Mal, Car2, Tnr, Plp1,
Ptgds, Tcf7l2), the h5ad layers `Ms`, `Mu` and `ATAC` give the per-cell
s / u / c coordinates (coloured by subtype); `varm["fit_anchor_{s,u,c}"]`
over `[fit_anchor_min_idx, fit_anchor_max_idx]` gives the fitted trajectory;
`varm["fit_anchor_{s,u,c}_sw"]` gives the switch-time markers (magenta).
Headers carry model / LRT / FDR / t_c50 / t_s50 (the latter two from the
`predict_exp` reconstruction of C3). The view is identical in all eight
panels."""))
cells.append(code('''from matplotlib.ticker import MaxNLocator
def dense(x):
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x).ravel()
def _3d(ax, s, u, c, cols, traj=None, swpts=None):
    ax.scatter(s, u, c, c=cols, s=1.0, linewidths=0, alpha=0.6, rasterized=True,
               depthshade=False)
    if traj is not None:
        ax.plot(traj[0], traj[1], traj[2], color="black", lw=0.8, alpha=0.9,
                rasterized=True)
    if swpts is not None:
        for k, mk in enumerate(["o", "X", "^"]):
            ax.plot([swpts[0][k]], [swpts[1][k]], [swpts[2][k]], mk, ms=2.6,
                    mfc="m", mec="m", mew=0.5, rasterized=True, zorder=5)
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.set_major_locator(MaxNLocator(2, prune="lower"))
        a.set_major_formatter(plt.FuncFormatter(lambda v, _: "%.4g" % v))
        a.pane.set_alpha(0.0)
        a.grid(False)
    ax.view_init(elev=16, azim=-68)
    ax.tick_params(labelsize=5.5, width=0.4, length=1.2, pad=0.5)
    ax.set_box_aspect((1, 1, 1))
    ax.dist = 7.8
LBL = {"s": (0.155, -0.085, "center", "top"),
       "u": (0.865, -0.115, "center", "top"),
       "c": (0.945, 0.480, "left", "center")}
def AXLBL(pos, fig):
    for k, (dx, dyy, ha, va) in LBL.items():
        fig.text(pos.x0 + dx * pos.width, pos.y0 + dyy * pos.height, k,
                 fontsize=5.5, ha=ha, va=va, color="0.15", weight="bold")
def _figw(fig, t):
    fig.canvas.draw()
    return t.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.bbox.width
from matplotlib.lines import Line2D
fig = plt.figure(figsize=(17.00 * CM, 8.40 * CM))
gs = fig.add_gridspec(2, 4, hspace=0.92, wspace=0.06, left=0.015, right=0.99,
                      top=0.795, bottom=0.150)
for i, g in enumerate(TRAJ8):
    sub = D[:, g]
    c = dense(sub.layers["ATAC"]); u = dense(sub.layers["Mu"])
    ss = dense(sub.layers["Ms"])
    out = D.uns["velo_s_params"].get("outlier", 99)
    keep = (c <= np.percentile(c, out)) & (u <= np.percentile(u, out)) & \\
           (ss <= np.percentile(ss, out))
    msk = keep & np.isfinite(c) & np.isfinite(u) & np.isfinite(ss)
    mi = int(sub.var["fit_anchor_min_idx"][0])
    ma = int(sub.var["fit_anchor_max_idx"][0])
    traj = tuple(dense(sub.varm[f"fit_anchor_{k}"])[mi:ma + 1] for k in ("s", "u", "c"))
    swp = tuple(dense(sub.varm[f"fit_anchor_{k}_sw"]) for k in ("s", "u", "c"))
    ax = fig.add_subplot(gs[i // 4, i % 4], projection="3d")
    _3d(ax, ss[msk], u[msk], c[msk], SC[msk], traj=traj, swpts=swp)
    pos = gs[i // 4, i % 4].get_position(fig)
    AXLBL(pos, fig)
    *_, tc50, ts50, model, sw = gene_traj(g)
    row = dec.loc[g]
    x0 = pos.x0 + 0.030 * pos.width
    y1 = pos.y1 + 0.032
    t1 = fig.text(x0, y1, g, fontsize=6.4, style="italic", weight="bold",
                  va="bottom", ha="left")
    fig.text(x0 + _figw(fig, t1) + 0.006, y1, "|  " + _cat(g), fontsize=5.5,
             ha="left", va="bottom", color="0.3")
    y2 = pos.y1 + 0.012
    fig.text(x0, y2, "model %d \\u00b7 LRT %s \\u00b7 FDR %s"
             % (model, _fmt(row.LRT), _fmt(row.FDR, "{:.0g}")),
             fontsize=5.5, ha="left", va="bottom", color="0.35")
    fig.text(x0, pos.y1 - 0.008, "tc50/ts50 = %.1f / %.1f" % (tc50, ts50),
             fontsize=5.5, ha="left", va="bottom", color="0.35")
fig.text(0.015, 0.010, "s = spliced RNA, u = unspliced, c = chromatin "
                       "accessibility (raw model units; view identical in all "
                       "eight). Black line = fitted trajectory; magenta = "
                       "switch times t_sw1/2/3.",
         fontsize=5.5, color="0.35", ha="left", va="bottom")
axl = fig.add_axes([0.02, 0.905, 0.9, 0.05]); axl.axis("off")
handles = [Line2D([0], [0], marker="o", ls="", ms=3.0, mfc=CT_COLS[k],
                  mec="none", label=CT[k]) for k in range(len(CT))]
handles.append(Line2D([0], [0], marker="o", ls="", ms=3.0, mfc="m", mec="none",
                      label="switch times t_sw1/2/3"))
axl.legend(handles=handles, loc="center left", ncol=len(CT) + 1, frameon=False,
           handletextpad=0.25, columnspacing=1.0, fontsize=5.5)
fig.savefig(f"{OUT}/figS4v3d_3d.pdf")
print("D 3d done")'''))

nb = nbf.v4.new_notebook(metadata={
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.8"}})
nb["cells"] = cells
out = os.path.join(OUT, "Figure_S4.ipynb")
nbf.write(nb, out)
print("wrote", out, len(cells), "cells")
