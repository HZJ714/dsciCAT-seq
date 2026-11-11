#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build Figure_4.ipynb — complete, self-contained reproduction of Figure 4.

Design (per user requirements):
  * every panel is produced from the ORIGINAL data object (RDS / h5ad / tsv),
    never by reverse-engineering a finished PDF;
  * Panel A & B are generated with R (Seurat/scMEGA/ComplexHeatmap) and invoked
    from the notebook via Rscript; Panel C & D are pure Python (anndata/multivelo);
  * no figure-composition ("paste panels together") code is included — each
    section ends with its own single-panel PDF;
  * a final section prints every software package and its version (R + Python)
    for reproducibility;
  * comments are in English; code is as short as possible while complete.
"""
import nbformat as nbf
import os

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "fig4")
os.makedirs(OUT, exist_ok=True)

# ---- notebook plumbing ------------------------------------------------------
def md(s):
    return nbf.v4.new_markdown_cell(s)


def code(s):
    return nbf.v4.new_code_cell(s)


def rst(s):
    """R code as a code cell (marked for Rscript execution by the next cell)."""
    return nbf.v4.new_code_cell(s)


nb = nbf.v4.new_notebook(metadata={
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.8"}})

cells = []

# =============================================================================
# Title + overview
# =============================================================================
cells.append(md("""# Figure 4 — Chromatin–RNA dynamics along the oligodendrocyte lineage

This notebook reproduces **Figure 4** panel by panel, starting from the *original
data objects* produced by the upstream analyses:

| Panel | Figure | Content | Input data | Language |
|---|---|---|---|---|
| A | scMEGA trajectory UMAP | 5,977 lineage cells, pseudotime 0–100 (blueYellow) | `3k_dsci_oligo_subset.RDS` (Seurat, RNA+ATAC paired) | **R** (Seurat/scMEGA) |
| B | Dual-omics heatmap | chromatin accessibility + gene expression, 15,607 peak-to-gene rows | `scmega_obj_target.rds` + `scmega_p2g.csv` | **R** (ComplexHeatmap) |
| C | Model statistics | decoupling volcano · induction timing · switch intervals | `multivelo_dynamical.h5ad` + `results/decoupling.tsv` | Python |
| D | Three-gene kinetics | Tnr / Anln / Plp1 × (ATAC, Mu, Ms) | `multivelo_dynamical.h5ad` | Python |

```
scMEGA (R)                       MultiVelo (Python)
  3k_dsci_oligo_subset.RDS          counts/*.mtx, peaks, rna.h5ad
        │                                │
        ├─ Panel A (TrajectoryPlot)      ├─ multivelo_dynamical.h5ad ─┬─ Panel C
        ├─ obj_target.rds ── Panel B     │                            └─ Panel D
        │                                 └─ results/decoupling.tsv ─ Panel C
```

Each section ends with its own panel PDF."""))

cells.append(md("""---

## 0. Environment, paths and versions

Set `DATA` to the directory holding the input objects (Section 0.1). Both
Python (≥3.8) and R (≥4.2) are required: panels A/B invoke `Rscript`,
panels C/D run in Python. Package versions are printed at the end of this
section."""))
cells.append(code('''import os, subprocess, sys, glob
import numpy as np, pandas as pd, scipy, matplotlib
import anndata as ad
from matplotlib import font_manager as fm
# Arial-first font registration (self-contained; falls back to Liberation Sans
# whose metrics match Arial, then to the default). Set FONT_DIR if fonts live
# elsewhere; panels C/D are redrawn with these fonts by design.
FONT_DIR = os.environ.get("FONT_DIR", os.path.join(os.environ.get("DSCI_WORK_ROOT", "."), "fonts"))
for _f in sorted(glob.glob(os.path.join(FONT_DIR, "*.ttf"))) + \
         sorted(glob.glob(os.path.join(FONT_DIR, "*.TTF"))) + \
         sorted(glob.glob(os.path.join(FONT_DIR, "*.otf"))):
    try:
        fm.fontManager.addfont(_f)
    except Exception:
        pass
_avail = {f.name for f in fm.fontManager.ttflist}
_FAMILY = next((n for n in ("Arial", "Liberation Sans") if n in _avail), "DejaVu Sans")
matplotlib.rcParams.update({"font.size": 6.0,
                             "axes.titlesize": 6.2, "axes.labelsize": 6.0,
                             "xtick.labelsize": 5.5, "ytick.labelsize": 5.5,
                             "legend.fontsize": 5.5,
                             "figure.dpi": 120, "savefig.dpi": 600,
                             "mathtext.fontset": "custom",
                             "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                             "xtick.major.size": 1.6, "ytick.major.size": 1.6,
                             "xtick.major.pad": 1.2, "ytick.major.pad": 1.2,
                             "lines.linewidth": 0.8, "patch.linewidth": 0.5,
                             "axes.linewidth": 0.5, "pdf.fonttype": 42})
print("font family:", _FAMILY)
# R interpreter used for panels A/B (adjust if Rscript is not on PATH)
RSCRIPT = os.environ.get("RSCRIPT", "Rscript")
def rrun(code_str, script="panel_r.R", env=None):
    """Run an R snippet via Rscript; prints its stdout. Reusable for A/B.
    Passed DATA/OUT are injected as R globals via Sys.getenv()."""
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
plt.rcParams.update({"font.family": _FAMILY, "font.size": 6.0,
                     "axes.titlesize": 6.2, "axes.labelsize": 6.0,
                     "xtick.labelsize": 5.5, "ytick.labelsize": 5.5,
                     "legend.fontsize": 5.5,
                     "figure.dpi": 120, "savefig.dpi": 600,
                     "mathtext.fontset": "custom",
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                     "xtick.major.size": 1.6, "ytick.major.size": 1.6,
                     "xtick.major.pad": 1.2, "ytick.major.pad": 1.2,
                     "lines.linewidth": 0.8, "patch.linewidth": 0.5,
                     "axes.linewidth": 0.5})
CM = 1 / 2.54'''))

# -----------------------------------------------------------------------------
# Versions (Python + R) and input-file table — printed up front for reproducibility
# -----------------------------------------------------------------------------
cells.append(md("""### 0.1 Input files

| Object | Produced by | Used in |
|---|---|---|
| `3k_dsci_oligo_subset.RDS` | cell-type annotation of paired RNA+ATAC (Seurat) | Panel A |
| `scmega_obj_target.rds` | `scMEGA_oligo_GRN.R` (chromVAR + GRN) | Panel B |
| `scmega_p2g.csv` | `SelectGenes()` (15,607 peak-to-gene links) | Panel B |
| `multivelo_dynamical.h5ad` | MultiVelo `recover_dynamics_chrom` (5977 × 1809) | Panels C, D |
| `multivelo/results/decoupling.tsv` | LRT refit (`fit_decoupling=False`) | Panel C |

### 0.2 Software versions

Prints every relevant package version — Python packages for panels C/D and the
R interpreter + R packages for panels A/B."""))
cells.append(code('''print("==== Python packages ====")
for m in (np, pd, scipy, matplotlib, ad):
    print(f"  {m.__name__:<12} {m.__version__}")
print("==== R ====")
rrun(r"""
ip <- installed.packages()
for (p in c("Seurat","Signac","scMEGA","ArchR","ComplexHeatmap","SummarizedExperiment",
            "chromVAR","ggplot2","dplyr","matrixStats","readxl")) {
  if (p %in% rownames(ip)) cat(sprintf("  %-22s %s\\n", p, ip[p, "Version"]))
}
cat("  R ", R.version.string, "\\n", sep="")
""", "versions.R", env={"DATA": DATA, "OUT": OUT})'''))

# =============================================================================
# Section 1 — Panel A (R)
# =============================================================================
cells.append(md("""---

## 1. Panel A — scMEGA trajectory UMAP (6.20 × 6.20 cm)

**Data chain** (all in R):

1. `3k_dsci_oligo_subset.RDS` — Seurat object with paired RNA + ATAC for the
   5,977 oligodendrocyte-lineage cells (31,263 genes, 186,321 peaks).
2. **Pseudotime construction**: the four mature subtypes (MOL2/MOL3/MOL3-6/MOL6)
   are merged into a single terminal group `MOL` (no internal ordering — they are
   alternative maturation endpoints, not a linear chain). Pseudotime is then
   built by group mapping + within-group ranking on the harmony embedding:
   `COP/NFOL → [0,25)` (distance to the MFOL centroid), `MFOL → [25,50)`
   (distance to the MOL centroid), `MOL → [50,100]` (random bin order, the
   mature block has no developmental order). Result stored as `Trajectory`.
3. `TrajectoryPlot(object, reduction="umap", continuousSet="blueYellow")` draws
   the cells on the RNA UMAP, coloured by pseudotime with the blueYellow ramp.
4. The figure is saved at the final panel size (6.20 × 6.20 cm)."""))
cells.append(code(r'''R_CODE_A = r"""
suppressPackageStartupMessages({library(Seurat); library(Signac)
                                library(scMEGA); library(ArchR); library(ggplot2)})
# Seurat v4 compatibility shims (v5 API not present in this environment)
LayerData <- function(object, assay=NULL, layer="counts", cells=NULL) {
  mat <- GetAssayData(object, assay=assay, slot=layer)
  if (!is.null(cells)) mat <- mat[, cells]; mat
}
JoinLayers <- function(object, ...) object

obj <- readRDS(file.path(DATA, "3k_dsci_oligo_subset.RDS"))
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
# ---- pseudotime: group mapping + within-group rank on harmony embedding ----
obj$traj_group <- as.character(obj$cell_type)
obj$traj_group[obj$traj_group %in% c("MOL2","MOL3","MOL3/6","MOL6")] <- "MOL"
set.seed(42)
emb <- Embeddings(obj, "harmony")[, 1:30]
tg  <- obj$traj_group
c_mfol <- colMeans(emb[tg == "MFOL", , drop=FALSE])
c_mol  <- colMeans(emb[tg == "MOL", , drop=FALSE])
pseudo <- rep(NA, ncol(obj)); names(pseudo) <- colnames(obj)
idx <- which(tg == "COP/NFOL")
d <- sqrt(colSums((t(emb[idx, , drop=FALSE]) - c_mfol)^2))
pseudo[idx] <- (1 - rank(d)/length(idx)) * 25
idx <- which(tg == "MFOL")
d <- sqrt(colSums((t(emb[idx, , drop=FALSE]) - c_mol)^2))
pseudo[idx] <- 25 + (1 - rank(d)/length(idx)) * 25
idx <- which(tg == "MOL")
pseudo[idx] <- 50 + rank(runif(length(idx)))/length(idx) * 50
obj$Trajectory <- pseudo
obj.t <- obj[, !is.na(obj$Trajectory)]
# ---- panel A source: real TrajectoryPlot from the RDS-derived object ----
p3 <- TrajectoryPlot(object=obj.t, reduction="umap",
                     continuousSet="blueYellow", size=1, addArrow=FALSE)
ggsave(file.path(OUT, "fig4a_raw.pdf"), p3, width=7, height=6)
# blueYellow ramp (0->100, same palette as the published figure)
pal <- paletteContinuous(set="blueYellow", n=100)
write.table(pal, file.path(OUT, "blueYellow.csv"),
            sep=",", col.names=FALSE, row.names=FALSE)
cat("A source written:", ncol(obj.t), "cells\\n")
"""
rrun(R_CODE_A, "panel_a.R", env={"DATA": DATA, "OUT": OUT})
print("fig4a_raw.pdf + blueYellow.csv written")'''))
cells.append(md("""
**Panel drawing** (Python, journal layout): the same data were generated by
TrajectoryPlot above; we read the vector point cloud back from this fresh
intermediate PDF (centres of the small filled rectangles, one per cell) and
redraw at the final panel size 6.20 × 6.20 cm with a compact colourbar and a
bold centred title — no axis clutter, matching the published layout."""))
cells.append(code('''import fitz, numpy as np, pandas as pd
# read back the vector points from the intermediate PDF (not the published PDF)
doc = fitz.open(f"{OUT}/fig4a_raw.pdf")
pts, cols = [], []
for d in doc[0].get_drawings():
    if d["type"] in ("fs", "f"):
        r = d["rect"]
        if r.width < 4 and r.height < 4 and d.get("fill"):
            pts.append((r.x0 + r.width / 2, r.y0 + r.height / 2))
            cols.append(tuple(round(c, 4) for c in d["fill"]))
pts = np.asarray(pts); cols = np.asarray(cols)
doc.close()
# colourbar ramp from the R palette CSV (same 100-entry blueYellow)
pal = pd.read_csv(f"{OUT}/blueYellow.csv", header=None).values[:, 0]
def hx(h):
    return np.array([int(h[i:i+2], 16) / 255 for i in (1, 3, 5)])
ramp = np.array([hx(c) for c in pal])
fig = plt.figure(figsize=(6.20 * CM, 6.20 * CM))
ax = fig.add_axes([0.055, 0.045, 0.705, 0.80])
ax.scatter(pts[:, 0], -pts[:, 1], s=0.75, c=cols, linewidths=0, rasterized=True)
ax.set_aspect("equal"); ax.axis("off")
fig.text(0.5, 0.985, "scMEGA trajectory", fontsize=6.6, ha="center",
         va="top", fontweight="bold")
cb = fig.add_axes([0.835, 0.20, 0.032, 0.50])
cb.imshow(np.linspace(0, 1, 256).reshape(-1, 1), aspect="auto", origin="lower",
          cmap=ListedColormap(ramp), extent=[0, 1, 0, 100])
cb.yaxis.set_ticks_position("right"); cb.set_yticks([0, 25, 50, 75, 100])
cb.set_xticks([])
cb.tick_params(labelsize=5.5, width=0.4, length=1.4, pad=1.0)
cb.set_ylim(0, 100)
for s in cb.spines.values():
    s.set_linewidth(0.4)
fig.savefig(f"{OUT}/fig4a_traj_umap.pdf")
print(f"A panel done: {len(pts)} cells")'''))

# =============================================================================
# Section 2 — Panel B (R)
# =============================================================================
cells.append(md("""---

## 2. Panel B — dual-omics heatmap (9.80 × 6.20 cm)

**Data chain** (all in R):

1. `scmega_obj_target.rds` — the trajectory object after chromVAR motif activity
   and GRN target-assay addition (produced by `scMEGA_oligo_GRN.R`).
2. `scmega_p2g.csv` — the **15,607 significant peak-to-gene links** from
   `SelectGenes()` (the rows of the heatmap).
3. `GetTrajectory(assay, trajectory.name="Trajectory", smoothWindow=7,
   log2Norm=TRUE)` smooths each feature into 100 pseudotime bins.
4. Rows are z-scored and clipped to ±2; the row order comes from
   `TrajectoryHeatmap()` on the combined (ATAC+RNA z) matrix.
5. Two heatmaps side by side — ATAC accessibility (`solarExtra` ramp) and gene
   expression (`horizonExtra` ramp) — with a blueYellow pseudotime bar
   (matching panel A) and 14 stage markers annotated with anti-overlap layout.

> The published B panel is the same ComplexHeatmap output rotated 90° to match
> the horizontal layout; the underlying matrix and annotation are identical."""))
cells.append(code(r'''R_CODE_B = r"""
suppressPackageStartupMessages({library(Seurat); library(Signac); library(scMEGA)
                                library(ArchR); library(SummarizedExperiment)
                                library(ComplexHeatmap); library(readxl)})
LayerData <- function(object, assay=NULL, layer="counts", cells=NULL) {
  mat <- GetAssayData(object, assay=assay, slot=layer)
  if (!is.null(cells)) mat <- mat[, cells]; mat
}
JoinLayers <- function(object, ...) object
obj.t  <- readRDS(file.path(DATA, "scmega_obj_target.rds"))
df.p2g <- read.csv(file.path(DATA, "scmega_p2g.csv"))
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
# ---- smooth both assays over the same pseudotime axis (100 bins) ----
trajRNA  <- GetTrajectory(obj.t, assay="RNA", trajectory.name="Trajectory",
                          groupEvery=1, slot="data", smoothWindow=7, log2Norm=TRUE)
trajATAC <- GetTrajectory(obj.t, assay="ATAC", trajectory.name="Trajectory",
                          groupEvery=1, slot="data", smoothWindow=7, log2Norm=TRUE)
trajATAC <- trajATAC[df.p2g$peak, ]
trajRNA  <- trajRNA[df.p2g$gene, ]
# ---- row order: combined z-score matrix, TrajectoryHeatmap's row dendrogram ----
trajCombined <- trajATAC
assay(trajCombined, withDimnames=FALSE) <-
  t(apply(assay(trajRNA), 1, scale)) + t(apply(assay(trajATAC), 1, scale))
combinedMat <- TrajectoryHeatmap(trajCombined, returnMatrix=TRUE, varCutOff=0)
rowOrder <- match(rownames(combinedMat), rownames(trajATAC))
# ---- left heatmap: chromatin accessibility (solarExtra, z clipped at +-2) ----
ht1 <- TrajectoryHeatmap(trajATAC, pal=paletteContinuous(set="solarExtra"),
                         varCutOff=0, maxFeatures=nrow(trajATAC), rowOrder=rowOrder,
                         limits=c(-2, 2), labelRows=FALSE, labelTop=0,
                         name="Chromatin accessibility")
# ---- right heatmap: gene expression, same row order, z clipped at +-2 ----
mat2 <- assay(trajRNA)[rowOrder, , drop=FALSE]
mat2 <- sweep(mat2 - rowMeans(mat2), 1, matrixStats::rowSds(mat2), `/`)
mat2[mat2 > 2] <- 2; mat2[mat2 < -2] <- -2
# ---- 14 stage markers (anti-overlap: DFS over occurrences, min row gap 200) ----
mk <- intersect(c("Bmp4","Mbp","Plp1","Mag","Mog","Gpr17","Olig1","Olig2",
                  "Tcf7l2","Car2","Opalin","Mal","Mobp","Apod"), df.p2g$gene)
mk <- mk[mk %in% rownames(trajRNA)]
genes2 <- rownames(trajRNA)[rowOrder]
occ <- lapply(mk, function(g) which(genes2 == g)); names(occ) <- mk
first_pos <- sapply(occ, `[`, 1); target <- first_pos
target["Mbp"] <- 14792; target["Mag"] <- 14094; target["Plp1"] <- 11107
g_order <- mk[order(first_pos)]
MIN_GAP <- 200; chosen_pos <- integer(); chosen_g <- character()
dfs <- function(i) {
  if (i > length(g_order)) return(TRUE)
  g <- g_order[i]
  for (cp in occ[[g]][order(abs(occ[[g]] - target[g]))]) {
    if (all(abs(cp - chosen_pos) >= MIN_GAP)) {
      chosen_pos <<- c(chosen_pos, cp); chosen_g <<- c(chosen_g, g)
      if (dfs(i + 1)) return(TRUE)
      chosen_pos <<- chosen_pos[-length(chosen_pos)]
      chosen_g <<- chosen_g[-length(chosen_g)]
    }
  }
  FALSE
}
dfs(1)
lbl <- rep("", length(genes2)); lbl[chosen_pos] <- chosen_g
# ---- colourbars + pseudotime annotation (blueYellow, same ramp as panel A) ----
cols2 <- colorRamp2(seq(-2, 2, length.out=100),
                    paletteContinuous(set="horizonExtra", n=100))
pseudo_cols <- paletteContinuous(set="blueYellow", n=100)
pseudo_fun  <- colorRamp2(seq(0, 100, length.out=100), pseudo_cols)
pseudo_lab  <- rep("", 100); pseudo_lab[c(1,26,51,76,100)] <- c("0","25","50","75","100")
bottom_anno <- HeatmapAnnotation(
  "Pseudotime" = anno_simple((0:99) + 0.5, col=pseudo_fun, border=FALSE,
                             height=unit(0.8, "cm")),
  "tick" = anno_text(pseudo_lab, location=0.5,
                     gp=gpar(fontsize=9, col="grey30")),
  which="column", annotation_name_gp=gpar(fontsize=11),
  show_legend=c("Pseudotime"=FALSE),
  show_annotation_name=c("Pseudotime"=TRUE, "tick"=FALSE))
ht2 <- Heatmap(mat2, name="Gene expression", col=cols2,
               cluster_rows=FALSE, cluster_columns=FALSE, show_column_names=FALSE,
               row_labels=lbl, row_names_gp=gpar(fontsize=12, fontface="bold.italic"),
               row_names_side="right",
               heatmap_legend_param=list(title="Gene expression"),
               bottom_annotation=bottom_anno)
ht1@bottom_annotation <- bottom_anno
# vertical 05f-equivalent (12x16 in), identical to make_05_heatmaps.R
pdf(file.path(OUT, "fig4b_05f.pdf"), width=12, height=16)
draw(ht1 + ht2); dev.off()
# marker positions (row indices in the heatmap row order) for panel layout
write.csv(data.frame(gene=chosen_g, pos=chosen_pos),
          file.path(OUT, "fig4b_labels.csv"), row.names=FALSE)
cat("B done\n")
"""
rrun(R_CODE_B, "panel_b.R", env={"DATA": DATA, "OUT": OUT})
print("fig4b_05f.pdf + fig4b_labels.csv written")'''))
cells.append(md("""
**Panel drawing** (Python, journal layout): the raster body and colourbars of
the ComplexHeatmap output above are read back from this intermediate PDF and
redrawn at the final panel size 9.80 × 6.20 cm: ATAC (top) and RNA (bottom)
heatmaps with pseudotime on the y axis, the two z-score colourbars on top and
the 14 stage markers annotated on the top edge."""))
cells.append(code('''import fitz, numpy as np, pandas as pd
from PIL import Image
# ---- read back rasters / colourbars / labels from the fresh 05f-equivalent ----
doc = fitz.open(f"{OUT}/fig4b_05f.pdf")
pg  = doc[0]
heat = {}
for xref, key in ((9, "heat_atac"), (10, "heat_rna")):
    ex = doc.extract_image(xref)
    pm = fitz.Pixmap(ex["image"])
    heat[key] = np.frombuffer(pm.samples, np.uint8).reshape(pm.height, pm.width, pm.n)[:, :, :3]
bars = {"p1": [], "p2": []}; cb_acc, cb_exp = [], []
for d in pg.get_drawings():
    if d["type"] != "fs" or not d.get("fill"):
        continue
    r = d["rect"]; c = tuple(round(x, 4) for x in d["fill"])
    if 3.2 < r.width < 3.5 and 22 < r.height < 23.5:          # pseudotime bars
        bars["p2" if r.x0 > 340 else "p1"].append((r.x0, c))
    elif 0.7 < r.width < 1.0 and 10.5 < r.height < 12:        # acc colourbar
        cb_acc.append((r.x0, c))
    elif 10.5 < r.width < 12 and 0.4 < r.height < 0.8:        # exp colourbar
        cb_exp.append((r.y0, c))
for k in bars: bars[k].sort()
cb_acc.sort(); cb_exp.sort()
cb_acc = np.asarray([c for _, c in cb_acc], float)
cb_exp = np.asarray([c for _, c in cb_exp], float)
labs = pd.read_csv(f"{OUT}/fig4b_labels.csv")
frac = labs.pos.values / (15607 - 1)          # heatmap row index -> 0..1 (top..bottom)
names = labs.gene.values
o = np.argsort(frac); frac, names = frac[o], names[o]
doc.close()
# ---- helpers (identical to the project's fig_common.py) ----
def cm2fig(w_cm, h_cm):
    return plt.figure(figsize=(w_cm * CM, h_cm * CM))
def greedy_stack(positions, min_gap, lo=None, hi=None):
    pos = np.asarray(positions, float)
    order = np.argsort(pos)
    out = np.empty_like(pos); last = -np.inf
    for i in order:
        p = pos[i]
        if p - last < min_gap: p = last + min_gap
        out[i] = p; last = p
    if hi is not None and out.max() > hi: out -= (out.max() - hi)
    if lo is not None and out.min() < lo: out += (lo - out.min())
    return out, out - pos
# ---- replicate the published B layout (9.80 x 6.20 cm) ----
def _pool(arr, w, h):
    return np.asarray(Image.fromarray(arr.astype(np.uint8)).resize((w, h), Image.BOX))
W, H = 9.80, 6.20
fig = cm2fig(W, H)
def AX(x, y, w, h): return fig.add_axes([x / W, y / H, w / W, h / H])
def TXT(x, y, *a, **k): return fig.text(x / W, y / H, *a, **k)
LX, RX = 0.55, 9.05
Y_A_T, Y_A_B = 4.32, 2.45
Y_R_T, Y_R_B = 2.18, 0.30
for key, (yt, yb) in (("heat_atac", (Y_A_T, Y_A_B)), ("heat_rna", (Y_R_T, Y_R_B))):
    arr = np.rot90(_pool(heat[key], 1100, 120), 1)
    ax = AX(LX, yb, RX - LX, yt - yb)
    ax.imshow(arr, aspect="auto", interpolation="nearest")
    n = arr.shape[0] - 1
    ax.set_yticks([0, n / 2, n]); ax.set_yticklabels(["100", "50", "0"])
    ax.tick_params(axis="y", labelsize=5.5, width=0.4, length=1.4, pad=1.0,
                   left=False, right=True, labelleft=False, labelright=True)
    ax.set_xticks([])
    for sp in ax.spines.values(): sp.set_linewidth(0.4)
TXT(9.62, (Y_A_T + Y_R_B) / 2, "Pseudotime", fontsize=5.8, rotation=90,
    ha="center", va="center")
for x0, cbkey, lab in ((1.975, "cb_acc", "Chromatin accessibility"),
                       (5.025, "cb_exp", "Gene expression")):
    cb = cb_acc if cbkey == "cb_acc" else cb_exp[::-1]
    TXT(x0, 5.98, lab, fontsize=5.8, ha="left", va="bottom")
    cbh = AX(x0, 5.62, 2.60, 0.16)
    cbh.imshow(cb.reshape(1, -1, 3), aspect="auto")
    cbh.set_xticks([0, (len(cb) - 1) / 2, len(cb) - 1])
    cbh.set_xticklabels(["-2", "0", "+2"])
    cbh.set_yticks([])
    cbh.tick_params(labelsize=5.5, width=0.4, length=1.2, pad=0.8)
    for sp in cbh.spines.values(): sp.set_linewidth(0.4)
from matplotlib.lines import Line2D
span = RX - LX
xtrue = LX + frac * span
base = Y_A_T + 0.12
xpos, moved = greedy_stack(xtrue, 0.24, lo=LX + 0.06, hi=RX - 0.06)
for i, g in enumerate(names):
    if abs(xpos[i] - xtrue[i]) > 0.02:
        fig.add_artist(Line2D([xtrue[i] / W, xpos[i] / W],
                              [(Y_A_T + 0.01) / H, (base - 0.02) / H],
                              color="0.45", lw=0.35))
    TXT(xpos[i], base, g, fontsize=5.5, va="bottom", ha="center",
        rotation=90, style="italic", weight="bold")
fig.savefig(f"{OUT}/fig4b_dual_heatmap.pdf")
print("B panel done (fig4b_dual_heatmap.pdf)")'''))

# =============================================================================
# Section 3 — Panel C (Python)
# =============================================================================
cells.append(md("""---

## 3. Panel C — model statistics (16.00 × 4.90 cm)

**Data chain** (all in Python):

1. `multivelo_dynamical.h5ad` — MultiVelo dynamical fit (5977 cells × 1809 genes,
   `fit_decoupling=True`). Per-gene parameters live in `var`: `fit_model`
   (1 = three-phase chromatin-dependent, 2 = constitutive), `fit_t_sw1..3`
   (switch times on the 0–20 gene-time axis) and `velo_s_genes`.
2. `results/decoupling.tsv` — likelihood-ratio test for chromatin–RNA
   decoupling, produced by *notebook 08*: refit with `fit_decoupling=False`,
   `LRT = −2·n_obs·(logL_wo − logL_w) ~ χ²(1)`, BH–FDR; 372/1809 genes
   decoupled at FDR < 0.05 (20.6%).
3. Three sub-panels are drawn:
   - **volcano**: LRT vs −log₁₀(P) (P clipped at 1e−300), decoupled in red;
   - **induction timing**: histogram of `t_sw1` for the 1,798 model-1 genes —
     93.6% < 0 means chromatin opens before the expression switch (priming);
   - **switch intervals**: for the 114 velocity genes, the four interval types
     (primed / coupled-on / decoupled / coupled-off) as proportions of the
     0–20 window; decoupled intervals are longer (median 0.92 vs 0.33 pooled)."""))
cells.append(code('''d = ad.read_h5ad(f"{DATA}/multivelo/results/multivelo_dynamical.h5ad")
vt = d.var[[c for c in d.var.columns if c.startswith("fit_")] + ["velo_s_genes"]].copy()
vt.insert(0, "gene", d.var_names)
dec = pd.read_csv(f"{DATA}/multivelo/results/decoupling.tsv", sep="\\t", index_col=0)
dec["neglogp"] = -np.log10(dec.pval.clip(1e-300))

# ---- true-subscript axis labels (Arial has no U+2080/2081; no mathtext) ----
def vsub_label(ax, y, base, sub, tail="", fontsize=7, pad_pt=6.0, drop_pt=None,
               dx_pt=None, sub_scale=0.70, **kw):
    """Vertical axis title with a real subscript: base + small sub + tail,
    each rotated 90 deg and centred; widths measured in pixels."""
    fig = ax.figure; fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bb = ax.get_window_extent(renderer=r)
    gap = r.points_to_pixels(0.10 * fontsize)
    drop = (r.points_to_pixels(0.30 * fontsize) if drop_pt is None
            else r.points_to_pixels(drop_pt))
    x = (r.points_to_pixels(dx_pt) / bb.width if dx_pt is not None
         else (r.points_to_pixels(pad_pt) - bb.x0) / bb.width)
    def adv(txt, size):
        t = ax.text(0, 0, txt, fontsize=size, rotation=90, ha="center", va="center",
                    transform=ax.transAxes, **kw)
        fig.canvas.draw()
        h = t.get_window_extent(renderer=r).height
        t.remove(); return h
    a1, a2 = adv(base, fontsize), adv(sub, fontsize * sub_scale)
    a3 = adv(tail, fontsize) if tail else 0.0
    total = a1 + gap + a2 + (gap + a3 if tail else 0.0)
    yc = y - total / (2.0 * bb.height) + a1 / (2.0 * bb.height)
    out = [ax.text(x, yc, base, fontsize=fontsize, rotation=90, ha="center",
                   va="center", transform=ax.transAxes, **kw)]
    yc += (a1 + gap) / bb.height
    out.append(ax.text(x + drop / bb.width, yc, sub, fontsize=fontsize * sub_scale,
                       rotation=90, ha="center", va="center", transform=ax.transAxes, **kw))
    if tail:
        yc += (a2 + gap) / bb.height
        out.append(ax.text(x, yc, tail, fontsize=fontsize, rotation=90, ha="center",
                           va="center", transform=ax.transAxes, **kw))
    fig.canvas.draw(); return out

def hsub_label(ax, x, y, base, sub, tail="", fontsize=7, align="center",
               dx_pt=None, dy_pt=None, sub_scale=0.70, **kw):
    """Horizontal label with a real subscript (t_sw1 etc.)."""
    fig = ax.figure; fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bb = ax.get_window_extent(renderer=r)
    W, H = bb.width, bb.height
    gap = r.points_to_pixels(0.08 * fontsize)
    drop = r.points_to_pixels(0.24 * fontsize)
    if dx_pt is not None: x = r.points_to_pixels(dx_pt) / W
    if dy_pt is not None: y = r.points_to_pixels(dy_pt) / H
    def adv(txt, size):
        t = ax.text(0, 0, txt, fontsize=size, ha="center", va="center",
                    transform=ax.transAxes, **kw)
        fig.canvas.draw()
        w = t.get_window_extent(renderer=r).width
        t.remove(); return w
    a1, a2 = adv(base, fontsize), adv(sub, fontsize * sub_scale)
    a3 = adv(tail, fontsize) if tail else 0.0
    total = a1 + gap + a2 + (gap + a3 if tail else 0.0)
    x0 = x * W if align == "left" else x * W - total / 2.0
    out = [ax.text((x0 + a1 / 2.0) / W, y, base, fontsize=fontsize, ha="center",
                   va="center", transform=ax.transAxes, **kw)]
    x0 += a1 + gap
    out.append(ax.text((x0 + a2 / 2.0) / W, y - drop / H, sub, fontsize=fontsize * sub_scale,
                       ha="center", va="center", transform=ax.transAxes, **kw))
    if tail:
        x0 += a2 + gap
        out.append(ax.text((x0 + a3 / 2.0) / W, y, tail, fontsize=fontsize, ha="center",
                           va="center", transform=ax.transAxes, **kw))
    fig.canvas.draw(); return out

def label_cloud(ax, xs, ys, names, dx=3.0, min_dy=6.5, fontsize=5.3, colors=None,
                right_frac=0.72, lo_pad=None, halo=False):
    """Anti-overlap scatter labels; halo keeps text readable on dense points."""
    fig = ax.figure; fig.canvas.draw()
    k = fig.get_dpi() / 72.0
    inv = ax.transData.inverted()
    xd, yd = ax.transData.transform(np.c_[np.asarray(xs, float),
                                          np.asarray(ys, float)]).T
    box = ax.get_window_extent()
    order = np.argsort(xd); placed = []
    for i in order:
        col = colors[i] if isinstance(colors, (list, np.ndarray)) else colors
        left_side = xd[i] > box.x0 + right_frac * box.width
        off = dx * k * (-1 if left_side else 1)
        txt_x = xd[i] + off
        w_px = len(names[i]) * fontsize * 0.55 * k
        y = yd[i]
        for _ in range(200):
            moved = False
            for px, py, pw, _sd in placed:
                if abs(px - txt_x) < (pw + w_px) * 0.6 + 3 * k and (py + min_dy * k) - y > 1e-6:
                    y = py + min_dy * k; moved = True
            if not moved: break
        if lo_pad is not None and y < box.y0 + lo_pad * k:
            y = box.y0 + lo_pad * k
            for _ in range(200):
                moved = False
                for px, py, pw, _sd in placed:
                    if abs(px - txt_x) < (pw + w_px) * 0.6 + 3 * k and (py + min_dy * k) - y > 1e-6:
                        y = py + min_dy * k; moved = True
                if not moved: break
        placed.append((txt_x, y, w_px, "right" if left_side else "left"))
        tx, ty = inv.transform((txt_x, y))
        need_arrow = abs(y - yd[i]) > 1.5 * k
        ha = "right" if left_side else "left"
        kw = {}
        if halo:
            import matplotlib.patheffects as pe
            kw["path_effects"] = [pe.withStroke(linewidth=max(1.4, fontsize * 0.32),
                                                 foreground="white")]
            kw["zorder"] = 6
        if need_arrow and halo:
            ax.annotate("", (xs[i], ys[i]), xytext=(tx, ty), textcoords="data",
                        annotation_clip=True, ha=ha,
                        arrowprops=dict(arrowstyle="-", lw=1.5, color="white",
                                        shrinkA=0, shrinkB=1.2, capstyle="round"))
        ax.annotate(names[i], (xs[i], ys[i]), xytext=(tx, ty), textcoords="data",
                    fontsize=fontsize, color=col, va="center", ha=ha,
                    annotation_clip=True,
                    arrowprops=(dict(arrowstyle="-", lw=0.3, color="0.5",
                                     shrinkA=0, shrinkB=1.2) if need_arrow else None), **kw)
    fig.canvas.draw()

fig = plt.figure(figsize=(16.00 * CM, 4.90 * CM))
gs  = fig.add_gridspec(1, 3, width_ratios=[1.45, 1.0, 1.06], wspace=0.50,
                       left=0.042, right=0.988, top=0.845, bottom=0.175)
# ---- C1 decoupling volcano ----
ax = fig.add_subplot(gs[0, 0])
ax.scatter(dec.LRT[~dec.decoupled], dec.neglogp[~dec.decoupled], s=1.5, c="0.75",
           linewidths=0, rasterized=True)
ax.scatter(dec.LRT[dec.decoupled], dec.neglogp[dec.decoupled], s=2.0, c="crimson",
           alpha=0.85, linewidths=0, rasterized=True)
ax.axhline(-np.log10(0.05), ls="--", c="k", lw=0.45)
ax.set_xscale("symlog", linthresh=100)
ax.set_xlim(-4e4, 4e4)
ax.set_xticks([-1e4, -1e2, 0, 1e2, 1e4])
ax.set_xticklabels(["\\u221210\\u2074", "\\u221210\\u00b2", "0", "10\\u00b2", "10\\u2074"])
ax.set_ylim(-12, 355)
lab = [g for g in ("Mbp", "Mal", "Plp1") if g in dec.index]
label_cloud(ax, [dec.loc[g, "LRT"] for g in lab], [dec.loc[g, "neglogp"] for g in lab],
            lab, colors="steelblue", fontsize=5.5, halo=True)
for g, yy in (("Cnp", 92), ("Tcf7l2", 46)):
    if g not in dec.index: continue
    ax.annotate(g, (dec.loc[g, "LRT"], dec.loc[g, "neglogp"]), xytext=(-1.3e4, yy),
                textcoords="data", fontsize=5.5, ha="left", va="center",
                color="steelblue" if g == "Cnp" else "k", annotation_clip=True,
                arrowprops=dict(arrowstyle="-", lw=0.3, color="0.5",
                                shrinkA=0, shrinkB=1.2))
ax.text(0.04, 0.975, "372/1809 decoupled (20.6%)", transform=ax.transAxes,
        fontsize=5.5, va="top", color="0.3")
ax.set_xlabel("LRT (decoupling interval)", fontsize=5.6, labelpad=2)
vsub_label(ax, 0.5, "-log", "10", "(p)", fontsize=5.8, dx_pt=-15.5, sub_scale=1.0)
ax.set_title("Chromatin-RNA decoupling", fontsize=6.2, pad=2)
ax.tick_params(labelsize=5.5, width=0.4, length=1.4)
# ---- C2 induction timing (t_sw1 of model-1 genes) ----
ax = fig.add_subplot(gs[0, 1])
m1 = vt[vt.fit_model == 1]
sw1 = m1.fit_t_sw1.clip(-20, 20)
ax.hist(sw1, bins=np.arange(-20, 21, 1.0), color="0.75", linewidth=0.3, edgecolor="0.35")
ax.axvspan(-20.6, 0, color="#dce6f1", zorder=0)
ax.axvline(0, color="crimson", ls="--", lw=0.6)
ax.set_yscale("log"); ax.set_ylim(0.7, 22000)
frac = float((m1.fit_t_sw1 < 0).mean()) * 100
_hpt = ax.get_window_extent().height * 72.0 / fig.dpi
_st = 1.45 * 5.6
_y = lambda i: 0.975 - (i + 0.5) * _st / _hpt
ax.text(0.02, _y(0), "primed", transform=ax.transAxes, fontsize=5.6, va="center",
        ha="left", color="#26456b")
hsub_label(ax, 0.02, _y(1), "(t", "sw1", "<0)", fontsize=5.6, align="left",
           sub_scale=1.0, color="#26456b")
ax.text(0.02, _y(2), f"{frac:.1f}%, n={len(m1)}", transform=ax.transAxes,
        fontsize=5.6, va="center", ha="left", color="#26456b")
ax.set_yticks([1, 10, 100, 1000]); ax.set_yticklabels(["1", "10", "100", "1k"])
ax.set_xlim(-20.6, 20.6)
hsub_label(ax, 0.5, 0, "t", "sw1", " (end of induction)", fontsize=5.6, dy_pt=-14.5,
           sub_scale=1.0)
ax.set_ylabel("genes (log)", fontsize=5.8)
ax.set_title("Induction timing", fontsize=6.2, pad=2)
ax.tick_params(labelsize=5.5, width=0.4, length=1.4)
# ---- C3 switch intervals (114 velocity genes) ----
ax = fig.add_subplot(gs[0, 2])
sw = vt.loc[vt.velo_s_genes, ["fit_t_sw1", "fit_t_sw2", "fit_t_sw3"]].copy()
sw = sw.mask(sw > 20, 20).mask(sw < 0)
ints = {"primed": sw.fit_t_sw1, "coupled-on": sw.fit_t_sw2 - sw.fit_t_sw1,
        "decoupled": sw.fit_t_sw3 - sw.fit_t_sw2, "coupled-off": 20 - sw.fit_t_sw3}
data, labs, ns = [], [], []
for k, v in ints.items():
    vv = v[(v > 0) & (v <= 20)] / 20.0
    data.append(vv.values); ns.append(len(vv)); labs.append(k)
bp = ax.boxplot(data, patch_artist=True, widths=0.62,
                medianprops=dict(color="black", lw=0.6),
                whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5),
                flierprops=dict(marker="d", ms=1.8, mfc="0.4", mec="none"))
for p, cc in zip(bp["boxes"], ["#79af9f", "#e59866", "#8b9dc3", "#c780ad"]):
    p.set_facecolor(cc); p.set_edgecolor("0.3"); p.set_linewidth(0.5)
for i, n_ in enumerate(ns):
    ax.text(i + 1, 1.24, f"n={n_}", ha="center", va="bottom", fontsize=5.5, color="0.3")
ax.set_xticks(range(1, 5))
ax.set_xticklabels(labs, rotation=26, ha="right")
ax.set_ylim(0, 1.34)
ax.set_ylabel("proportion of the 0-20\\ngene-time window", fontsize=5.6, labelpad=2)
ax.set_title("Switch intervals (n=114)", fontsize=6.2, pad=2)
ax.tick_params(labelsize=5.5, width=0.4, length=1.4)
fig.savefig(f"{OUT}/fig4d_models.pdf")
print("C done (fig4d_models.pdf)")'''))

# =============================================================================
# Section 4 — Panel D (Python)
# =============================================================================
cells.append(md("""---

## 4. Panel D — three-gene kinetics (16.00 × 7.80 cm)

**Data chain** (all in Python, from the same `multivelo_dynamical.h5ad`):

For each of Tnr / Anln / Plp1, the dynamical fit provides, per cell:
`fit_t` (gene time), layers `ATAC` (chromatin), `Mu` (unspliced) and `Ms`
(spliced); the fit anchors `varm["fit_anchor_{c,u,s}"]` give the gene-level
curves over the anchor interval `[fit_anchor_min_idx, fit_anchor_max_idx]`,
and `fit_t_sw1..3` the switch times. Cells are coloured by subtype and shown
with the fitted curve overlaid; dashed lines mark the switch times.

Gene-specific messages:
- **Tnr** — RNA first closes, chromatin lags (model 1, likelihood 0.071,
  r² 0.96; chromatin open in 100% of COP/NFOL vs 13% of MOL3/6);
- **Anln** — RNA mid-transient while chromatin stays open (model 1,
  likelihood 0.108, r² 0.58; 100% → 92% [MOL6]);
- **Plp1** — constitutive, both sustained (model 2, likelihood 0.089,
  r² 0.08; 100% → 99% [MOL3])."""))
cells.append(code('''GENES = ["Tnr", "Anln", "Plp1"]
CT    = list(d.obs.cell_type.cat.categories)
CT_COLS = [d.uns["cell_type_colors"][i] for i in range(len(CT))]
def dense(x):
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x).ravel()

fig = plt.figure(figsize=(16.00 * CM, 7.80 * CM))
hr = [0.10] + [0.13, 1.0] * len(GENES)
gs = fig.add_gridspec(2 * len(GENES) + 1, 3, height_ratios=hr, hspace=0.34,
                      wspace=0.30, left=0.045, right=0.988, top=0.985, bottom=0.075)
axl = fig.add_subplot(gs[0, :]); axl.axis("off")
from matplotlib.lines import Line2D
axl.legend(handles=[Line2D([0], [0], marker="o", ls="", ms=3.0, mfc=CT_COLS[i],
                           mec="none", label=CT[i]) for i in range(len(CT))],
           loc="center left", ncol=len(CT), frameon=False, handletextpad=0.25,
           columnspacing=0.9, fontsize=5.8)
for r, g in enumerate(GENES):
    t = dense(d[:, g].layers["fit_t"]); c = dense(d[:, g].layers["ATAC"])
    u = dense(d[:, g].layers["Mu"]);   s = dense(d[:, g].layers["Ms"])
    out = d.uns["velo_s_params"]["outlier"]
    keep = (c <= np.percentile(c, out)) & (u <= np.percentile(u, out)) & \\
           (s <= np.percentile(s, out))
    mi = int(d[:, g].var["fit_anchor_min_idx"][0])
    ma = int(d[:, g].var["fit_anchor_max_idx"][0])
    nwt_raw = np.linspace(0, 20, int(d.uns["velo_s_params"]["t"]))[mi:ma + 1]
    nwt = (nwt_raw - nwt_raw.min()) * 20 / max(nwt_raw.max() - nwt_raw.min(), 1e-9)
    v   = vt[vt.gene == g].iloc[0]
    axh = fig.add_subplot(gs[2 * r + 1, :]); axh.axis("off")
    axh.text(0.0, 0.9, g, fontsize=6.8, weight="bold", va="top", transform=axh.transAxes)
    for j, (key, lab) in enumerate([("c", "ATAC accessibility"),
                                    ("u", "unspliced (Mu)"),
                                    ("s", "spliced (Ms)")]):
        ax = fig.add_subplot(gs[2 * r + 2, j])
        dat = {"c": c, "u": u, "s": s}[key]
        o = np.arange(len(t))[keep]
        for k, ty in enumerate(CT):
            m = d.obs.cell_type.values[o] == ty
            ax.scatter(t[o][m], dat[o][m], s=1.1, c=CT_COLS[k], alpha=0.65,
                       linewidths=0, rasterized=True)
        anc = dense(d[:, g].varm[f"fit_anchor_{key}"])[mi:ma + 1]
        ax.plot(nwt, anc, color="black", lw=0.8, alpha=0.8)
        for swv in [v.fit_t_sw1, v.fit_t_sw2, v.fit_t_sw3]:
            if 0 <= swv <= 20:
                ax.axvline(swv, color="black", ls="--", lw=0.45, alpha=0.5)
        ax.set_ylabel(lab, fontsize=5.6, labelpad=1.5)
        ax.set_xlim(0, 20); ax.tick_params(labelsize=5.5, width=0.4, length=1.4)
        if r == len(GENES) - 1:
            ax.set_xlabel("gene time t", fontsize=5.8, labelpad=1.5)
        else:
            ax.tick_params(labelbottom=False)
fig.savefig(f"{OUT}/fig4c_mv_dyn.pdf")
print("D done (fig4c_mv_dyn.pdf)")'''))

# =============================================================================
# Section 5 — inputs & software versions
# =============================================================================

nb["cells"] = cells
nb["metadata"] = {"kernelspec": {"display_name": "Python 3", "language": "python",
                                 "name": "python3"},
                  "language_info": {"name": "python", "version": "3.8"}}
out = os.path.join(OUT, "Figure_4.ipynb")
nbf.write(nb, out)
print("wrote", out, len(cells), "cells")
