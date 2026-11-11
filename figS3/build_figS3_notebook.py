#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build Figure_S3.ipynb — reproduction of Figure S3 (single-omic views of the same
three libraries, cross-library cluster matching, and how much of the oligoden-
drocyte lineage COP/NFOL needs to be seen at all).

Design (the rules the author set for Figure 1 / 2 / 3 / 5):
  * reproduction first: every panel states the input files it consumes and the
    shortest complete code that turns them into the panel;
  * styling is deliberately minimal — no journal geometry, no panel letters, no
    compositing (those live in figures/scripts/compose_figS3*.py and
    figures/scripts/s3_cop_stability.py and are not part of the science);
  * §1 reduces the three per-library objects to the two things every panel needs
    (a cell-type call + average expression over the genes the three libraries
    share), and reads the ATAC embedding from Figure 3's cache when it exists;
  * §4 is the expensive measurement of the figure: the same delivered clustering
    parameters applied to lineage subsets down-sampled to eight depths, 20 seeds
    each. Its two known pitfalls are written into the code as comments, because
    both changed the published numbers;
  * comments in English, code as short as possible while complete.

Run:  python build_figS3_notebook.py
Out:  <repo>/figS3/Figure_S3.ipynb
"""
import os
import nbformat as nbf

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
nb = nbf.v4.new_notebook(metadata={
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.9"}})
C = []
md = lambda s: C.append(nbf.v4.new_markdown_cell(s))
code = lambda s: C.append(nbf.v4.new_code_cell(s.strip("\n")))

# ===========================================================================
md("""# Figure S3 — Same cells, single-omic views: what matches across libraries, and how deep you have to sequence to see COP/NFOL

| Panel | What it shows | Input | Code |
|---|---|---|---|
| A | RNA UMAP of each library, coloured by its own cell-type call | 3 per-library RDS (§1) | R/Seurat (§2) |
| B | ATAC UMAP of the same cells | joint multiomic RDS (or Figure 3's cache) | R/Seurat (§2) |
| C–E | cluster–cluster correlation of average expression: 50k vs HT, 100k vs HT, 50k vs 100k | `S3_avgexpr_*.csv` (§1) | Python (§3) |
| F | can a COP/NFOL cluster be found at all, as a function of how many lineage cells enter the clustering (8 conditions × 20 seeds) | the 4 lineage subsets (§4.1) | R/Seurat + Python (§4.2, §5.1) |
| G | the seven-marker identity of the COP/NFOL / OPC / MOL clusters found in each library | same run | Python (§5.2) |

```
        3k_dsci_ht_10x_gs_subset.RDS  (Figure 2's joint object, 70,643 nuclei)
              │  one subset() per library - written out and checked in Figure 3 §1.0:
              │    EL14_subset.RDS     batch == 'Cortex_EL14'       ->  4,029 nuclei
              │    3k_FY13_subset.RDS  orig.ident == 'dsciCAT-seq-50k'  -> 17,093
              │    3k_FY14_subset.RDS  orig.ident == 'dsciCAT-seq-100k' -> 30,350
              │    HT_subset_harmony.RDS orig.ident == 'HT-scCAT-seq-HT' -> 14,762
              │  (+ each one's own library-level reductions, graphs, RNA_snn_res.1, cell_type)
              ▼
EL14_subset.RDS   3k_FY13_subset.RDS   3k_FY14_subset.RDS        (Figure 3's fig3_light.RDS)
   ct labels,        ct labels,           ct labels,                umap.atac of all 51,472
   RNA umap,         RNA umap,              RNA umap                    nuclei, barcodes 100%
   RNA data          RNA data               RNA data                        │
        └───────────────────────┬───────────────────────────────────────────┘
                        §1  S3_light.RDS + S3_avgexpr_{EL14,FY13,FY14}.csv
                        (cell types over the 5,597 genes the three libraries share)
                                │
              ┌─────────────────┴──────────────────┐
              §2  panels A/B (RNA + ATAC UMAPs)     §3  panels C/D/E (cluster-cluster r)

HT_subset_harmony.RDS + the 3 subsets ── §4.1 lineage cache (Oligo + OPC only, RNA counts)
        └─ §4.2  same clustering as the delivered pipeline, 8 conditions × 20 seeds
                 -> S3_u13_clusters.csv, S3_marker_table.csv, S3_strict4029.csv
        └─ §5.1  panel F        └─ §5.2  panel G
```""")

# ===========================================================================
md("""---

## 0. Environment, inputs, outputs

| Variable | Default | Contents |
|---|---|---|
| `FIGS3_DATA` | `<work>/FY` | the four per-library Seurat objects of §0.1 — all four are `subset()`s of Figure 2's joint object, the same selection Figure 3 §1.0 runs and verifies |
| `FIGS3_JOINT` | `FIGS3_DATA/3k_dsci_ht_10x_gs_subset.RDS` | the joint object they are cut from (only read when Figure 3's cache is absent, for `umap.atac`) |
| `FIGS3_LIGHT` | `../fig3/fig3_out/fig3_light.RDS` | Figure 3's derived cache; if present it supplies the ATAC embedding and §1 skips the 13 GB object |
| `FIGS3_OUT` | `figS3_out` | created here; every panel PDF/PNG/CSV lands here |
| `FIGS3_FORCE` | unset | `1` re-runs §1 and §4.1 even when their caches exist |

Python ≥3.8 with `numpy pandas scipy matplotlib`. §1, §2 and §4 need `Rscript`
with `Seurat`, `ggplot2`, `dplyr`, `tidyr`, `ggrepel`, `pals`, `grid`.""")

code("""
import os, sys, subprocess, time
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as _fm

DATA = os.environ.get('FIGS3_DATA', os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')))
LIGHT = os.environ.get('FIGS3_LIGHT', os.path.join('..', 'fig3', 'fig3_out', 'fig3_light.RDS'))
def show(*paths):        # report the saved file only -- nothing is displayed inline
    for q in paths:
        print('%s  (%.0f KB)' % (q, os.path.getsize(q) / 1024.0))

OUT  = os.environ.get('FIGS3_OUT', 'figS3_out')
RSCRIPT = os.environ.get('RSCRIPT', 'Rscript')
os.makedirs(OUT, exist_ok=True)

_FONT_DIR = os.environ.get('FONT_DIR', os.environ.get('DSCI_FONTS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'fonts')))
if os.path.isdir(_FONT_DIR):
    for _f in sorted(os.listdir(_FONT_DIR)):
        if _f.lower().endswith(('.ttf', '.otf')):
            try: _fm.fontManager.addfont(os.path.join(_FONT_DIR, _f))
            except Exception: pass
_avail = {f.name for f in _fm.fontManager.ttflist}
plt.rcParams.update({'font.family': next((n for n in ('Arial', 'Liberation Sans')
                                          if n in _avail), 'DejaVu Sans'),
                     'pdf.fonttype': 42, 'axes.linewidth': 0.5,
                     'xtick.labelsize': 6, 'ytick.labelsize': 6,
                     'axes.titlesize': 8, 'axes.labelsize': 7, 'legend.fontsize': 6,
                     'savefig.dpi': 300})

def write_r(name, src):
    p = os.path.join(OUT, name)
    with open(p, 'w') as fh: fh.write(src)
    return p

def run_r(path, tag=''):
    t0 = time.time()
    rc = subprocess.run([RSCRIPT, path], env=dict(os.environ, FIGS3_DATA=DATA,
                                                  FIGS3_LIGHT=LIGHT, FIGS3_OUT=OUT))
    assert rc.returncode == 0, '%s failed (exit %d)' % (path, rc.returncode)
    print('%s done in %.1f min' % (tag or os.path.basename(path), (time.time() - t0) / 60))

# the three libraries, in figure order, and the cell-type order of panels C-E
LIB = dict(EL14='HT-scCAT-seq (n = 4,029)', FY13='dsciCAT-seq-50k (n = 17,093)',
           FY14='dsciCAT-seq-100k (n = 30,350)')
GLIA = ['COP/NFOL', 'OPC', 'Oligo', 'AST', 'MG', 'EC', 'Ependyma']
COP = 'COP/NFOL'
print('DATA =', DATA, '| LIGHT =', LIGHT, '| OUT =', OUT, '| font =', plt.rcParams['font.family'])""")

md("""### 0.1 Input files

| Used by | Path | What it is |
|---|---|---|
| §1 | `FIGS3_DATA/HT/EL14_subset.RDS` (0.6 GB) | HT-scCAT-seq, single library, 4,029 nuclei = `subset(joint, batch == 'Cortex_EL14')` + its own library-level analysis (`cell_type`, `RNA_snn_res.1`, `umap`) |
| §1 | `FIGS3_DATA/FY13/3k_FY13_subset.RDS` (1.7 GB) | dsciCAT-seq-50k, 17,093 nuclei = `subset(joint, orig.ident == 'dsciCAT-seq-50k')`, same layout |
| §1 | `FIGS3_DATA/FY14/3k_FY14_subset.RDS` (2.6 GB) | dsciCAT-seq-100k, 30,350 nuclei = `subset(joint, orig.ident == 'dsciCAT-seq-100k')`, same layout |
| §4.1 | `FIGS3_DATA/HT/HT_subset_harmony.RDS` (1.8 GB) | the four HT-scCAT-seq libraries pooled, 14,762 nuclei = `subset(joint, orig.ident == 'HT-scCAT-seq-HT')`, pooled-library-level analysis |
| §1 | `FIGS3_DATA/3k_dsci_ht_10x_gs_subset.RDS` (13 GB) | the joint multiomic object those four are cut from — read only for `umap.atac`, and only if Figure 3's cache is absent |

**None of the four is separate data.** They are the four single-library /
pooled-library selections of the joint object that Figure 2 produces, and
Figure 3 §1.0 is where the selection is written out and checked: same barcodes
in the same order, same RNA count matrix, same QC metadata, cell for cell (the
two dsciCAT-seq objects keep 33 ATAC peaks more than the joint object, which is
the only place their metadata differs). What each object adds on top of its
subset is the analysis run at that scope — its own PCA/UMAP/WNN reductions,
neighbour graphs, `RNA_snn_res.1` and, for the three single-library ones, the
hand-written `cell_type` table. That is what panels A–E measure and what §4.1
re-clusters, so the files are read rather than rebuilt; Figure 3 writes the pure
subsets too when it runs with `FIG3_SAVE_DERIVED=1`.

Every object is the **analysis object of the published figure**, not a re-run of
it: the same cells, the same reductions, the same `cell_type` column. The
COP/NFOL re-annotation applied to FY13/FY14 in §1 is the one Figure 2/3 use
(`RNA_snn_res.1` cluster 29 / 31).""")

md("""### 0.2 Provenance — read this before comparing with the printed figure

* **The four input objects are selections of one joint object** (§0.1): the cell
  sets, count matrices and QC metadata are Figure 2's, verified cell for cell in
  Figure 3 §1.0; only the library-level analysis is added per object.
* **A/B are the published embeddings** — RNA space is each library's own `umap`
  reduction, ATAC space is the joint object's `umap.atac`, both coloured by that
  library's cell-type call. Nothing is re-clustered.
* **C–E are new measurements on §1's output**, not on a saved matrix: the cell
  types of two libraries are matched by the Pearson correlation of their average
  expression vectors over the 5,597 high-variable genes the three libraries share
  (union of the three `VariableFeatures` lists ∩ common genes — the same
  definition the delivered panel used). Colour is the z-score **along the query
  row**, so the brightest cell of a row is "which type of the other library this
  type looks most like". The support numbers (`300/300`) are a gene bootstrap.
* **F and G are the delivered analysis, re-run.** §4.2 applies exactly the
  clustering the delivered objects went through (NormalizeData → 3,000 HVG vst →
  PCA 30 → 10 neighbours → SLM `algorithm = 3`, resolutions 0.5/1/1.5/2, primary
  `RNA_snn_res.1`) to lineage subsets at eight depths, 20 down-sampling seeds
  each, and calls a cluster COP-like with a rule derived from the positive
  controls (FY13 cluster 29, FY14 cluster 31) and then applied blind:
  `Bmp4 ≥ 0.25 & Gpr17 ≥ 0.5 & %Gpr17 ≥ 30 & C1ql1 ≤ 0.2 & Ptprz1 ≤ 1.5 & n ≥ 10`.
* **Two bugs that changed published numbers are kept visible** (§4.2 comments):
  cluster means must be divided column-wise (`sweep(..., 2, n, "/")`), and a
  down-sample must be computed into a variable before `subset()`, because Seurat
  reuses the promise and 20 "different" seeds silently drew the same 20 times.
* **The delivered F/G layout is artwork** (swim-lanes, tick marks drawn as
  strokes, hatch for partial hits). §5 draws the same numbers as a plain
  per-condition strip: what is plotted is exactly `S3_u13_clusters.csv`.""")

md("""### 0.3 Input check""")

code("""
need = [('§1 EL14',    os.path.join(DATA, 'HT/EL14_subset.RDS')),
        ('§1 FY13',    os.path.join(DATA, 'FY13/3k_FY13_subset.RDS')),
        ('§1 FY14',    os.path.join(DATA, 'FY14/3k_FY14_subset.RDS')),
        ('§4.1 HTpool',os.path.join(DATA, 'HT/HT_subset_harmony.RDS'))]
for k, p in need:
    print('%-12s %-9s %s' % (k, 'ok' if os.path.exists(p) else 'MISSING', p))
assert all(os.path.exists(p) for _, p in need), 'input missing - set FIGS3_DATA'
print('ATAC embedding will come from Figure 3 cache: %s (exists: %s)'
      % (LIGHT, os.path.exists(LIGHT)))""")

md("""**Expected** — four objects present; the Figure 3 cache is optional (it only
saves the 13 GB read in §1).""")

md("""### 0.4 Software""")

code("""
import scipy
print('python   ', sys.version.split()[0])
for m in (numpy := np, pandas := pd, scipy, matplotlib):
    print('%-9s' % m.__name__.split('.')[0], m.__version__)
subprocess.run([RSCRIPT, '--version'])""")

# ===========================================================================
md("""---

## 1. The step every panel depends on — three libraries, one annotation, average expression over the genes they share

**Input**: the three per-library objects of §0.1 (+ Figure 3's cache if present) —
the same cells Figure 3 §1.0 subsets out of the joint object, carrying the
library-level analysis on top. **Output**: `figS3_out/S3_light.RDS`, `S3_avgexpr_{EL14,FY13,FY14}.csv`,
`S3_counts_{EL14,FY13,FY14}.csv`, `S3_summary.txt`.

The cell-type call is the figure's: `cell_type`, except that the COP/NFOL cluster
is taken by name from `RNA_snn_res.1` (29 in FY13, 31 in FY14; EL14 has none).
Average expression is taken over **all** genes and then cut to the shared
high-variable set — union of the three libraries' `VariableFeatures`, intersected
with the genes all three carry. That cut is what makes the cross-library
correlations of §3 comparable: every type vector has the same 5,597 dimensions.""")

code("""
R1 = r'''
suppressPackageStartupMessages({library(Seurat); library(dplyr)})
DATA  <- Sys.getenv('FIGS3_DATA', os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')))
LIGHT <- Sys.getenv('FIGS3_LIGHT', '')
OUT   <- Sys.getenv('FIGS3_OUT', 'figS3_out')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sub <- c(EL14 = 'HT/EL14_subset.RDS', FY13 = 'FY13/3k_FY13_subset.RDS',
         FY14 = 'FY14/3k_FY14_subset.RDS')
copclust <- c(EL14 = NA_integer_, FY13 = 29L, FY14 = 31L)   # the COP/NFOL cluster of each library
libs <- list(); log <- character()
for (nm in names(sub)) {
  t0 <- proc.time()[3]
  o <- readRDS(file.path(DATA, sub[[nm]])); DefaultAssay(o) <- 'RNA'
  ct <- as.character(o$cell_type)
  if (!is.na(copclust[[nm]])) ct[o$RNA_snn_res.1 == copclust[[nm]]] <- 'COP/NFOL'
  d <- GetAssayData(o, assay = 'RNA', slot = 'data')
  cts <- sort(unique(ct))
  avg <- t(sapply(cts, function(k) rowMeans(as.matrix(d[, ct == k, drop = FALSE]))))
  libs[[nm]] <- list(n = length(ct), ct = ct, cells = colnames(o),
                     rna = as.data.frame(Embeddings(o, 'umap')), avg = avg,
                     genes = rownames(o), hvg = VariableFeatures(o),
                     counts = table(ct)[cts])
  log <- c(log, sprintf('%s cells=%d types=%d genes=%d hvg=%d loaded %.1f min', nm, ncol(o),
                        length(cts), length(libs[[nm]]$genes), length(libs[[nm]]$hvg),
                        (proc.time()[3] - t0) / 60))
  cat(tail(log, 1), '\\n'); rm(o, d); gc(verbose = FALSE)
}
common <- Reduce(intersect, lapply(libs, function(s) s$genes))
hv <- intersect(unique(Reduce(union, lapply(libs, function(s) s$hvg))), common)
log <- c(log, sprintf('common genes %d | shared HVG used %d', length(common), length(hv)))
cat(tail(log, 1), '\\n')
for (nm in names(libs)) {
  s <- libs[[nm]]
  s$avg <- s$avg[, hv, drop = FALSE]
  write.csv(round(s$avg, 5), file.path(OUT, paste0('S3_avgexpr_', nm, '.csv')))
  write.csv(as.data.frame(s$counts), file.path(OUT, paste0('S3_counts_', nm, '.csv')))
  s$genes <- s$hvg <- NULL; libs[[nm]] <- s
}
# ATAC embedding: Figure 3's cache already holds it (same cells, matched by barcode)
if (!is.na(LIGHT) && nzchar(LIGHT) && file.exists(LIGHT)) {
  S3 <- readRDS(LIGHT)
  for (nm in names(libs)) {
    idx <- match(libs[[nm]]$cells, S3$libs[[nm]]$cells)
    libs[[nm]]$atac <- S3$libs[[nm]]$atac[idx, , drop = FALSE]
    log <- c(log, sprintf('%s ATAC from Figure 3 cache: %d / %d cells', nm, sum(!is.na(idx)),
                          libs[[nm]]$n))
  }
} else {
  ht <- readRDS(file.path(DATA, '3k_dsci_ht_10x_gs_subset.RDS'))
  E <- Embeddings(ht, 'umap.atac')
  for (nm in names(libs)) {
    idx <- match(libs[[nm]]$cells, colnames(ht))
    libs[[nm]]$atac <- as.data.frame(E[idx, , drop = FALSE])
    log <- c(log, sprintf('%s ATAC from joint object: %d / %d cells', nm, sum(!is.na(idx)),
                          libs[[nm]]$n))
  }
  rm(ht, E); gc(verbose = FALSE)
}
cat(paste(log, collapse = '\\n'), '\\n')
types <- sort(unique(unlist(lapply(libs, function(s) s$ct))))
cols <- unique(c(as.vector(pals::glasbey(32)), as.vector(pals::alphabet(26))))[1:length(types)]
names(cols) <- types
saveRDS(list(libs = libs, cols = cols, types = types, hv = hv), file.path(OUT, 'S3_light.RDS'))
writeLines(log, file.path(OUT, 'S3_summary.txt'))
cat('wrote', file.path(OUT, 'S3_light.RDS'), '\\n')
'''
cache = os.path.join(OUT, 'S3_light.RDS')
if os.path.exists(cache) and os.environ.get('FIGS3_FORCE') != '1':
    print('cache exists, skipped §1:', cache, '(%.1f MB) - delete it or set FIGS3_FORCE=1'
          % (os.path.getsize(cache) / 1e6))
    print(open(os.path.join(OUT, 'S3_summary.txt')).read())
else:
    run_r(write_r('figS3_s1_object.R', R1), '§1')""")

md("""**Expected** — EL14 4,029 cells / 22 types, FY13 17,093 / 32, FY14 30,350 / 33;
common genes ≈ 31,263 and **5,597 shared HVG**; ATAC embeddings for all 51,472
nuclei (100 % matched either from Figure 3's cache or from the joint object).
Measured on this machine, 2026-10-22: the three objects load in ≈ 2.5 min total,
the whole step ≈ 3 min with Figure 3's cache (≈ 8 min without it).""")

# ===========================================================================
md("""---

## 2. Panels A and B — the same cells in RNA space and in ATAC space

**Input**: `S3_light.RDS`. **Output**: `figS3_out/figS3ab_umaps.pdf` (+ `.png`),
a 2 × 3 grid: top row = each library's own RNA UMAP, bottom row = the joint
object's ATAC UMAP for the same cells, both coloured by that library's cell-type
call with the shared palette. Nothing is re-clustered; the panels differ only in
which embedding the published coordinates came from.""")

code("""
R2 = r'''
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(ggrepel); library(grid)})
OUT <- Sys.getenv('FIGS3_OUT', 'figS3_out')
S   <- readRDS(file.path(OUT, 'S3_light.RDS')); libs <- S$libs; cols <- S$cols
pt  <- c(EL14 = 0.07, FY13 = 0.055, FY14 = 0.05)
ttl <- c(EL14 = 'HT-scCAT-seq (n = 4,029)', FY13 = 'dsciCAT-seq-50k (n = 17,093)',
         FY14 = 'dsciCAT-seq-100k (n = 30,350)')
norm2 <- function(coords, span = 2.4) {
  x <- coords[, 1]; y <- coords[, 2]; s <- max(diff(range(x)), diff(range(y)))
  data.frame(x = (x - mean(range(x))) / s * span, y = (y - mean(range(y))) / s * span)
}
panel <- function(df, t) {
  lab <- df %>% group_by(cell_type) %>% summarise(x = median(x), y = median(y), .groups = 'drop')
  ggplot(df, aes(x, y)) +
    geom_point(aes(color = cell_type), size = df$pt[1], stroke = 0, alpha = 0.95) +
    scale_color_manual(values = cols, guide = 'none') +
    geom_text_repel(data = lab, aes(x, y, label = cell_type, color = cell_type),
                    size = 5 / 2.845276, box.padding = unit(0.25, 'lines'),
                    point.padding = unit(0.2, 'lines'), min.segment.length = 0.1,
                    segment.size = 0.3, segment.alpha = 0.7, max.overlaps = Inf,
                    force = 4, force_pull = 1.2, bg.color = 'white', bg.r = 0.15, seed = 42) +
    coord_fixed() + xlim(-1.7, 1.7) + ylim(-1.7, 1.7) + theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0)) +
    annotate('text', x = -1.68, y = 1.66, label = t, size = 5 / 2.845276, hjust = 0)
}
ps <- list()
for (emb in c('rna', 'atac')) for (nm in names(libs)) {
  df <- norm2(libs[[nm]][[emb]]); df$cell_type <- libs[[nm]]$ct; df$pt <- pt[[nm]]
  ps[[length(ps) + 1]] <- panel(df, paste(ttl[[nm]], toupper(emb)))
}
draw <- function(g, x, y, w, h) { pushViewport(viewport(x = x, y = y, width = w, height = h,
  just = c('left', 'bottom'), default.units = 'cm')); grid.draw(g); popViewport() }
for (dev in c('pdf', 'png')) {
  if (dev == 'pdf') pdf(file.path(OUT, 'figS3ab_umaps.pdf'), width = 15 / 2.54, height = 10 / 2.54)
  else png(file.path(OUT, 'figS3ab_umaps.png'), width = 15, height = 10, units = 'cm', res = 300)
  grid.newpage()
  for (i in seq_along(ps)) draw(ggplotGrob(ps[[i]]), ((i - 1) %% 3) * 5,
                                (1 - (i - 1) %/% 3) * 5, 5, 5)
  dev.off()
}
cat('wrote figS3ab_umaps.pdf + .png (2 rows x 3 columns, 5 x 5 cm each)\\n')
'''
run_r(write_r('figS3_s2_umaps.R', R2), '§2')""")

code("""show(os.path.join(OUT, 'figS3ab_umaps.png'))""")

md("""**Expected** — six 5 × 5 cm panels; the type sets differ per library exactly as
in Figure 3 (22 / 32 / 33 types, COP/NFOL present in the two dsciCAT-seq panels
and absent from HT-scCAT-seq), and the ATAC row keeps the same cells in a
different layout.""")

# ===========================================================================
md("""---

## 3. Panels C–E — which type of one library does a type look like in another

**Input**: `S3_avgexpr_{EL14,FY13,FY14}.csv` (§1). **Output**:
`figS3_out/figS3cde_xcorr.pdf` (+ `.png`) and `figS3_xcorr_stats.csv`.

For a pair of libraries, every cell type of one is compared with every cell type
of the other by the Pearson correlation of its 5,597-gene average-expression
vector. The colour is the z-score **along the query row**, so the brightest cell
of a row is that type's best match; the query's own row is boxed in red and the
`COP/NFOL` row is the one the figure is about. The support number under each
panel is a gene bootstrap: 300 resamples of the 5,597 genes, how often the best
match of `COP/NFOL` stays the same.""")

code("""
A = {k: pd.read_csv(os.path.join(OUT, 'S3_avgexpr_%s.csv' % k), index_col=0)
     for k in ('EL14', 'FY13', 'FY14')}
# if the matrix the delivered panels were drawn from is reachable, check against it
REF = os.environ.get('S3_REF', os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'scripts/cache'))
for k, a in A.items():
    p = os.path.join(REF, 'avgexpr_%s.csv' % k)
    if os.path.exists(p):
        b = pd.read_csv(p, index_col=0).reindex(columns=a.columns, index=a.index)
        print('%s vs delivered cache %s: %d x %d | max |diff| = %.1e'
              % (k, 'same shape' if b.notna().all().all() else 'shape differs',
                 a.shape[0], a.shape[1], np.nanmax(np.abs(a.values - b.values))))

def order(names):                                     # glia first, then alphabetical
    return [x for x in GLIA if x in names] + sorted(n for n in names if n not in GLIA)
def xcorr(a, b):
    X, Y = A[a].values, A[b].values
    R = np.corrcoef(np.vstack([X, Y]))[:len(X), len(X):]
    M = pd.DataFrame(R, index=list(A[a].index), columns=list(A[b].index))
    return M.loc[order(M.index), order(M.columns)]     # reorder rows AND columns together

PAIRS = [('FY13', 'EL14', 'C', 'dsciCAT-seq-50k  vs  HT-scCAT-seq'),
         ('FY14', 'EL14', 'D', 'dsciCAT-seq-100k  vs  HT-scCAT-seq'),
         ('FY13', 'FY14', 'E', 'dsciCAT-seq-50k  vs  dsciCAT-seq-100k')]
rng = np.random.default_rng(1)
NB = 300
stats = []
fig, axes = plt.subplots(1, 3, figsize=(17.4, 5.6))
for ax, (q, t, letter, ttl) in zip(axes, PAIRS):
    M = xcorr(q, t)
    Z = M.sub(M.mean(1), axis=0).div(M.std(1, ddof=0), axis=0).clip(-2.6, 2.6)
    ax.imshow(Z.values, cmap='bwr_r' if False else matplotlib.colors.LinearSegmentedColormap.from_list(
        'bz', ['#2166ac', '#8fb4dc', '#ffffff', '#f0a58a', '#b2182b']),
        vmin=-2.6, vmax=2.6, aspect='auto', interpolation='nearest')
    ax.set_xticks(np.arange(M.shape[1])); ax.set_yticks(np.arange(M.shape[0]))
    ax.set_xticklabels(M.columns, rotation=90, fontsize=4.2, ha='center', va='bottom')
    ax.set_yticklabels(M.index, fontsize=4.2, ha='right', va='center')
    ax.xaxis.set_ticks_position('top'); ax.tick_params(length=1.5, pad=1)
    for s in ax.spines.values(): s.set_visible(False)
    i = list(M.index).index(COP)
    ax.add_patch(plt.Rectangle((-0.5, i - 0.5), M.shape[1], 1, fill=False,
                               ec='#8c1015', lw=1.0))            # the COP/NFOL row
    best = M.loc[COP].idxmax(); r = M.loc[COP, best]
    # gene bootstrap of the best match
    Xi, Yi = A[q].loc[M.index].values, A[t].loc[M.columns].values
    cnt = {}
    for _ in range(NB):
        g = rng.integers(0, Xi.shape[1], Xi.shape[1])
        Rb = np.corrcoef(np.vstack([Xi[:, g], Yi[:, g]]))[:len(Xi), len(Xi):]
        nm = M.columns[int(np.argmax(Rb[i]))]; cnt[nm] = cnt.get(nm, 0) + 1
    top, ntop = max(cnt.items(), key=lambda kv: kv[1])
    ax.set_xlabel('%s: %s → %s  r = %.2f   ·   %s/%d bootstraps' %
                  (COP, COP, best, r, ntop, NB), fontsize=5.5, labelpad=4)
    ax.set_title('%s  %s' % (letter, ttl), fontsize=6.5)
    stats.append(dict(panel=letter, query_library=q, target_library=t, query=COP, best=best,
                      r=round(r, 4), bootstrap=top, bootstrap_n=ntop, bootstrap_total=NB,
                      opc_opc=round(M.loc['OPC', 'OPC'], 4), oligo_oligo=round(M.loc['Oligo', 'Oligo'], 4),
                      n_query_types=M.shape[0], n_target_types=M.shape[1]))
    M.to_csv(os.path.join(OUT, 'figS3_xcorr_%s_%s.csv' % (q, t)))
S = pd.DataFrame(stats)
S.to_csv(os.path.join(OUT, 'figS3_xcorr_stats.csv'), index=False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'figS3cde_xcorr.pdf')); fig.savefig(
    os.path.join(OUT, 'figS3cde_xcorr.png'), dpi=300); plt.close(fig)
print('wrote figS3cde_xcorr.pdf | shared HVG per vector:', A['EL14'].shape[1])
print(S.to_string(index=False))""")

code("""show(os.path.join(OUT, 'figS3cde_xcorr.png'))""")

md("""**Expected** — 5,597 genes per vector; HT has 22 types (no COP/NFOL), 50k 32,
100k 33. When `S3_REF` points at the cache the delivered panels were drawn from,
§1's matrices match it to `max |diff| = 0` (only the gene column order differs,
which is irrelevant to a correlation). The measured result, identical to the
delivered annotation:

| Panel | pair | COP/NFOL best match | r | OPC → OPC | bootstrap |
|---|---|---|---|---|---|
| C | 50k vs HT | **OPC** | 0.81 | 0.95 | 300/300 |
| D | 100k vs HT | **OPC** | 0.80 | 0.96 | 300/300 |
| E | 50k vs 100k | **COP/NFOL** | 0.99 | 1.00 | 300/300 |

That is the figure's argument in one line: a dsciCAT-seq COP/NFOL cluster matches
an **HT OPC** cluster almost as well as it matches the other dsciCAT library, so
the state is not a technology artefact — but the HT library never *splits* it.""")

# ===========================================================================
md("""---

## 4. Panel F's measurement — how much lineage has to enter the clustering before COP/NFOL can be found

### 4.1 The lineage subsets

**Input**: the four objects of §0.1 plus `HT_subset_harmony.RDS` (HT pooled).
**Output**: `figS3_out/S3_lineage_{FY13,FY14,HTp,EL14}.rds` — the **RNA assay** of
the cells whose `cell_type` is `Oligo`, `OPC` or `COP/NFOL`, i.e. the lineage
cells that the delivered sub-clustering was run on, with its counts, its
log-normalised data *and* the `scale.data` the whole-brain object carried.

Caching the assay rather than a counts matrix is what makes §4.2 the delivered
computation: Seurat's `RunPCA` reads `scale.data`, and when the subset's own
variable genes are not all in it, it silently uses the intersection. That
intersection is small — 425 of 3,000 for the 144 EL14 lineage cells — and it is
the space the published clusters were found in, so re-scaling here would move
boundaries and change the panel. §4.2 prints the number for every condition.

This is the only place in the notebook that reads the HT-pooled object, and the
only expensive one. Measured here: the lineage entering the sub-clustering is
2,944 of 17,093 nuclei (17.2 %) in FY13, 4,557 of 30,350 (15.0 %) in FY14, 571 of
14,762 (3.9 %) in the pooled HT object and 144 of 4,029 (3.6 %) in EL14 — those
four numbers are the x-axis of panel F.""")

code("""
R41 = r'''
suppressPackageStartupMessages({library(Seurat)})
DATA <- Sys.getenv('FIGS3_DATA', os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY'))); OUT <- Sys.getenv('FIGS3_OUT', 'figS3_out')
SRC <- c(FY13 = 'FY13/3k_FY13_subset.RDS', FY14 = 'FY14/3k_FY14_subset.RDS',
         HTp  = 'HT/HT_subset_harmony.RDS', EL14 = 'HT/EL14_subset.RDS')
for (k in names(SRC)) {
  f <- file.path(OUT, paste0('S3_lineage_', k, '.rds'))
  if (file.exists(f) && Sys.getenv('FIGS3_FORCE') != '1') { cat(k, 'cached\\n'); next }
  o <- readRDS(file.path(DATA, SRC[[k]])); nall <- ncol(o)
  keep <- as.character(o$cell_type) %in% c('Oligo', 'OPC', 'COP_NFOL', 'COP/NFOL')
  o <- subset(o, cells = colnames(o)[keep]); DefaultAssay(o) <- 'RNA'
  saveRDS(o[['RNA']], f)
  cat(sprintf('%-5s whole brain %6d | lineage %5d (%.2f%%)\\n', k, nall, sum(keep), 100 * sum(keep) / nall))
  rm(o); gc(verbose = FALSE)
}
'''
run_r(write_r('figS3_s41_lineage.R', R41), '§4.1')""")

md("""### 4.2 One clustering, eight depths, twenty seeds

**Input**: the four lineage files. **Output**: `figS3_out/S3_u13_clusters.csv`
(every cluster of every condition/seed/resolution, with its marker means and the
COP-like call), `S3_marker_table.csv` (the four full conditions, seven markers,
mean + % per cluster — panel G), `S3_strict4029.csv` (the strict
whole-brain→4,029-nuclei draws), `S3_summary.txt`.

The clustering is the delivered one, and the primary resolution is
`RNA_snn_res.1` — the column the published COP/NFOL labels were actually read
from. Two implementation notes are part of the result and are therefore in the
code, not in a footnote: cluster means must be divided **column-wise**, and a
down-sample must be materialised before `subset()`, because Seurat reuses the
promise and 20 "different" seeds silently drew the same 20 times.

One detail worth knowing before reading the hit counts: `ScaleData` is **not**
run here. The cached assay carries the whole-brain `scale.data`, `RunPCA` uses
the genes the subset's variable list shares with it, and that is the space in
which the published clusters were found — re-scaling on 144 cells would move
boundaries and change the panel.""")

code("""
R42 = r'''
suppressPackageStartupMessages({library(Seurat)})
DATA <- Sys.getenv('FIGS3_DATA', os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY'))); OUT <- Sys.getenv('FIGS3_OUT', 'figS3_out')
MK   <- c('Gpr17', 'Bmp4', 'C1ql1', 'Tmem144', 'Ptprz1', 'Plp1', 'Mbp')
RES  <- c(0.5, 1.0, 1.5, 2.0); PRIMARY <- 'RNA_snn_res.1'
NSEED <- as.integer(Sys.getenv('NSEED', '20'))
loadk <- function(k) CreateSeuratObject(counts = readRDS(file.path(OUT, paste0('S3_lineage_', k, '.rds'))),
                                        project = k, assay = 'RNA')
recly <- function(o) { set.seed(42)                       # same seed for every condition
  o <- NormalizeData(o, verbose = FALSE)
  o <- suppressWarnings(FindVariableFeatures(o, selection.method = 'vst', nfeatures = 3000, verbose = FALSE))
  # no ScaleData here on purpose: the cached assay brings the whole-brain scale.data,
  # RunPCA then uses the genes the two HVG lists share - exactly what the delivered
  # run did. The count is printed so the reader can see how small that space is.
  cat(sprintf('   %s cells: %d vst genes, %d of them scaled\n', ncol(o),
              length(VariableFeatures(o)),
              sum(VariableFeatures(o) %in% rownames(GetAssayData(o, assay = 'RNA',
                                                                slot = 'scale.data')))))
  o <- RunPCA(o, npcs = 30, verbose = FALSE)
  o <- FindNeighbors(o, reduction = 'pca', dims = 1:10, graph.name = 'RNA_snn', verbose = FALSE)
  FindClusters(o, algorithm = 3, resolution = RES, verbose = FALSE) }
tab <- function(o, cond, seed = NA) {
  gs <- intersect(MK, rownames(o))
  ex <- as.matrix(GetAssayData(o, assay = 'RNA', slot = 'data')[gs, , drop = FALSE]); storage.mode(ex) <- 'double'
  out <- NULL
  for (cn in paste0('RNA_snn_res.', RES)) {
    cl <- as.character(o@meta.data[[cn]]); lv <- sort(unique(cl))
    W <- sapply(lv, function(k) as.numeric(cl == k)); n <- colSums(W)
    # column-wise division: (ex %*% W) / n would recycle element-wise and is WRONG
    m <- sweep(ex %*% W, 2, n, '/'); p <- 100 * sweep((ex > 0) %*% W, 2, n, '/')
    d <- data.frame(cond = cond, seed = ifelse(is.na(seed), -1, seed), res = cn, cl = lv,
                    n = as.integer(n), nlin = sum(n), stringsAsFactors = FALSE)
    for (k in MK) d[[paste0('mean_', k)]] <- as.numeric(m[k, ])
    for (k in MK) d[[paste0('pct_', k)]]  <- as.numeric(p[k, ])
    d$cop <- with(d, mean_Bmp4 >= 0.25 & mean_Gpr17 >= 0.5 & pct_Gpr17 >= 30 &
                          mean_C1ql1 <= 0.2 & mean_Ptprz1 <= 1.5 & n >= 10)
    d$celllevel_copfrac <- 100 * mean((ex['Gpr17', ] > 0) & (ex['Bmp4', ] > 0) &
                                      (ex['C1ql1', ] == 0) & (ex['Ptprz1', ] < 1))
    out <- rbind(out, d) }
  stopifnot(all(out$pct_Gpr17 <= 100 + 1e-6, na.rm = TRUE))
  out }
subN <- function(o, n, seed) { set.seed(seed)
  # materialise the sample: subset(o, cells = sample(...)) lets Seurat reuse the
  # promise, and every seed then draws the SAME cells (this is what turned 20
  # independent draws into one and changed the published x/20 numbers)
  sel <- sample(colnames(o), min(n, ncol(o))); subset(o, cells = sel) }
ALL <- NULL; t0 <- proc.time()[3]
lin <- lapply(c('FY13', 'FY14', 'HTp', 'EL14'), loadk); names(lin) <- c('FY13', 'FY14', 'HTp', 'EL14')
for (k in names(lin)) { d <- tab(recly(lin[[k]]), paste0(k, '_full_lineage')); ALL <- rbind(ALL, d)
  cat(sprintf('%-22s lineage %5d | COP cluster at %s: %s\\n', paste0(k, '_full_lineage'), d$nlin[1],
              PRIMARY, ifelse(any(d$cop & d$res == PRIMARY), 'YES', 'no'))) }
N4029 <- round(ncol(lin$FY13) * 4029 / 17093)   # lineage cells left when 50k is thinned to 4,029 nuclei
cat(sprintf('(dsciCAT whole brain at 4,029 nuclei = lineage %d cells)\\n', N4029))
SPECS <- list(c('FY13_wholebrain_to_4029n', 'FY13', N4029),
              c('FY13_lineage_to_HTpooledN', 'FY13', ncol(lin$HTp)),
              c('FY13_lineage_to_EL14N', 'FY13', ncol(lin$EL14)),
              c('HTpooled_lineage_to_EL14N', 'HTp', ncol(lin$EL14)))
for (sd in seq_len(NSEED)) {
  for (sp in SPECS) {
    d <- tab(recly(subN(lin[[sp[2]]], as.integer(sp[3]), sd)), sp[1], sd); ALL <- rbind(ALL, d) }
  cat('.', ''); flush(stdout()) }
cat('\\n')
write.csv(ALL, file.path(OUT, 'S3_u13_clusters.csv'), row.names = FALSE)
# panel G's table is the four full conditions at the primary resolution
mt <- ALL[grepl('_full_lineage', ALL$cond) & ALL$res == PRIMARY, ]
mt$group <- sub('_lineage', '', mt$cond)
write.csv(mt, file.path(OUT, 'S3_marker_table.csv'), row.names = FALSE)
# the strict depth: thin the WHOLE BRAIN to 4,029 nuclei first, then take the lineage
fy <- readRDS(file.path(DATA, 'FY13/3k_FY13_subset.RDS')); DefaultAssay(fy) <- 'RNA'
S <- NULL
for (sd in seq_len(NSEED)) { set.seed(sd)
  sel <- sample(colnames(fy), min(4029, ncol(fy))); w <- subset(fy, cells = sel)
  keep <- as.character(w$cell_type) %in% c('Oligo', 'OPC')
  w2 <- subset(w, cells = colnames(w)[keep])
  d <- tab(recly(w2), 'FY13_strict_4029nuclei', sd); S <- rbind(S, d)
  cat(sprintf('strict seed=%2d lineage %d cells | COP cluster %s\\n', sd, sum(d$n),
              any(d$cop & d$res == PRIMARY))) }
write.csv(S, file.path(OUT, 'S3_strict4029.csv'), row.names = FALSE)
cat(sprintf('%.0f s | %d condition rows, %d marker rows, %d strict rows\\n',
            (proc.time()[3] - t0), nrow(ALL), nrow(mt), nrow(S)))
'''
res = os.path.join(OUT, 'S3_u13_clusters.csv')
if os.path.exists(res) and os.environ.get('FIGS3_FORCE') != '1':
    print('results exist, skipped §4.2:', res, '- delete it or set FIGS3_FORCE=1 to re-run')
else:
    run_r(write_r('figS3_s42_u13.R', R42), '§4.2')""")

md("""**Expected** — lineage cells entering the clustering: FY13 2,944 · FY14 4,557 ·
HT pooled 571 · EL14 144, and the down-sampling ladder 144 / 571 / 694. At
`RNA_snn_res.1` with the blind rule: the three full conditions and HT-pooled give
a COP-like cluster **every** time (1/1), the thinned dsciCAT conditions only
sometimes, and 144 lineage cells **never** (0/20). Measured on this machine,
2026-10-22: §4.1 ≈ 4 min (the four objects), §4.2 ≈ 4 min for 84 re-clusterings +
20 strict draws.""")

# ===========================================================================
md("""---

## 5. Panels F and G — the two pictures drawn from §4.2

### 5.1 Panel F — resolvability as a function of depth

**Input**: `S3_u13_clusters.csv`, `S3_strict4029.csv`. **Output**:
`figS3_out/figS3f_resolvability.pdf` (+ `.png`) and `figS3_f_conditions.csv`.""")

code("""
cl = pd.read_csv(os.path.join(OUT, 'S3_u13_clusters.csv'))
st = pd.read_csv(os.path.join(OUT, 'S3_strict4029.csv'))
PRIM = 'RNA_snn_res.1'
def summarise(d, cond, kind, plat, note):
    d = d[(d.cond == cond) & (d.res == PRIM)]
    seeds = sorted(d.seed.unique()); nh = [s for s in seeds if d[d.seed == s].cop.any()]
    sizes = [int(d[(d.seed == s) & d.cop].n.max()) for s in nh]
    return dict(panel='F', condition=cond, kind=kind, platform=plat, label=note,
                n_seed=len(seeds), n_hit=len(nh),
                lineage_cells=int(round(d.nlin.median())) if len(d) else 0,
                cells_lo=min(sizes) if sizes else 0, cells_hi=max(sizes) if sizes else 0,
                bmp4=float(d[d.cop].mean_Bmp4.median()) if sizes else np.nan,
                gpr17=float(d[d.cop].mean_Gpr17.median()) if sizes else np.nan)
F = [summarise(cl, 'FY13_lineage_to_EL14N',      'part', 'dsciCAT-seq', 'lineage → 144'),
     summarise(cl, 'FY13_lineage_to_HTpooledN',  'part', 'dsciCAT-seq', 'lineage → 571'),
     summarise(st, 'FY13_strict_4029nuclei',     'part', 'dsciCAT-seq', 'whole brain → 4,029'),
     summarise(cl, 'FY13_wholebrain_to_4029n',   'part', 'dsciCAT-seq', 'proxy: lineage → 694'),
     summarise(cl, 'FY13_full_lineage',          'full', 'dsciCAT-seq', '50k full · 17,093 nuclei'),
     summarise(cl, 'FY14_full_lineage',          'full', 'dsciCAT-seq', '100k full · 30,350 nuclei'),
     summarise(cl, 'EL14_full_lineage',          'none', 'HT-scCAT-seq', 'single lib · 4,029 nuclei'),
     summarise(cl, 'HTp_full_lineage',           'full', 'HT-scCAT-seq', '4 lib pooled · 14,762 nuclei'),
     summarise(cl, 'HTpooled_lineage_to_EL14N',  'part', 'HT-scCAT-seq', 'pooled lineage → 144')]
F = pd.DataFrame(F)
F.to_csv(os.path.join(OUT, 'figS3_f_conditions.csv'), index=False)
print(F[['label', 'lineage_cells', 'n_hit', 'n_seed', 'cells_lo', 'cells_hi',
         'bmp4', 'gpr17']].to_string(index=False))""")

code("""
# one row per condition; x = lineage cells entering (log); the mark is the outcome
fig, ax = plt.subplots(figsize=(7.2, 3.4))
COL = {'dsciCAT-seq': '#238a4d', 'HT-scCAT-seq': '#2f4b7c'}
for k, r in enumerate(F.itertuples()):
    y = len(F) - 1 - k
    col = COL[r.platform]
    hit = r.n_hit / r.n_seed
    if hit == 0:                                             # never separated
        ax.scatter([r.lineage_cells], [y], marker='x', color=col, s=95, linewidths=1.3, zorder=4)
    else:                                                    # hatched = sometimes, filled = always
        ax.scatter([r.lineage_cells], [y], s=95, marker='o', facecolors=(col if hit == 1 else 'white'),
                   edgecolors=col, linewidths=1.0, hatch=('' if hit == 1 else '///'), zorder=4)
    note = ('%d/%d draws' % (r.n_hit, r.n_seed)) if r.n_seed > 1 else ('found' if hit else 'not found')
    if r.cells_hi:
        note = note + ' · ' + (str(r.cells_lo) if r.cells_lo == r.cells_hi
                               else '%d-%d' % (r.cells_lo, r.cells_hi)) + ' cells'
    ax.text(r.lineage_cells * 1.4, y, note, va='center', fontsize=5.5, color=col)
ax.set_yticks(np.arange(len(F))[::-1])
ax.set_yticklabels(['%s · %s' % (r.platform, r.label) for r in F.itertuples()], fontsize=5.5)
ax.set_xscale('log'); ax.set_xlim(110, 40000); ax.set_ylim(-0.7, len(F) - 0.3)
ax.set_xticks([144, 571, 2944, 4557])          # 694 is the nominal proxy depth; the
ax.set_xticklabels(['144', '571', '2,944', '4,557'], fontsize=5.5)   # strict draws land at 706 (median)
ax.set_xlabel('oligodendrocyte-lineage cells entering the delivered sub-clustering (Oligo + OPC)',
              fontsize=5.8)
ax.tick_params(length=2, labelsize=5.5)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
ax.axvline(500, color='#cccccc', lw=0.6, zorder=1)
ax.grid(axis='x', lw=0.3, color='#eeeeee', zorder=0)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'figS3f_resolvability.pdf')); fig.savefig(
    os.path.join(OUT, 'figS3f_resolvability.png'), dpi=300); plt.close(fig)
show(os.path.join(OUT, 'figS3f_resolvability.png'))""")

md("""**Expected** — 144 lineage cells: 0/20 draws find a COP-like cluster (both
platforms); 571: ~6/20; 694: ~7–9/20; every *full* condition: 1/1. The point of
the panel is the left-to-right ordering — the same chemistry succeeds or fails
depending on how many lineage cells enter the clustering, and the pooled HT
library (571 cells) succeeds where the thinned dsciCAT library at the same depth
only sometimes does.""")

md("""### 5.2 Panel G — the marker identity of the clusters that were found

**Input**: `S3_marker_table.csv`. **Output**: `figS3_out/figS3g_marker_columns.pdf`
(+ `.png`).

For each of the four full conditions the best COP-like cluster is taken (the rule
itself, or, if it never fires, the cluster with the highest `Bmp4 × Gpr17`), then
the OPC cluster (highest `C1ql1` among the rest) and the MOL cluster (highest
`Plp1` among the rest, n ≥ 20). Dot size = % expressing, colour = the marker's
z-score across the twelve columns — exactly the two channels of the delivered
panel.""")

code("""
MK7 = ['Gpr17', 'Bmp4', 'C1ql1', 'Tmem144', 'Ptprz1', 'Plp1', 'Mbp']
mt = pd.read_csv(os.path.join(OUT, 'S3_marker_table.csv'))
GROUPS = [('dsciCAT-seq-50k', 'FY13_full'), ('dsciCAT-seq-100k', 'FY14_full'),
          ('HT-scCAT-seq, 4 lib pooled', 'HTp_full'), ('HT-scCAT-seq, single lib', 'EL14_full')]
def cols_for(g):
    t = mt[mt.group == g]
    cop = t[t.cop].sort_values('n', ascending=False).iloc[0] if t.cop.any() else \\
          t.loc[(t.mean_Bmp4 * t.mean_Gpr17).idxmax()]
    opc = t[t.cl != cop.cl].sort_values('mean_C1ql1', ascending=False).iloc[0]
    mol = t[(t.cl != cop.cl) & (t.cl != opc.cl) & (t.n >= 20)].sort_values('mean_Plp1',
                                                                           ascending=False).iloc[0]
    return [dict(kind='COP/NFOL', r=cop, ok=bool(cop.cop)), dict(kind='OPC', r=opc, ok=None),
            dict(kind='MOL', r=mol, ok=None)]
G = [dict(c, grp=nm) for nm, g in GROUPS for c in cols_for(g)]
mat  = np.array([[c['r']['mean_' + m] for c in G] for m in MK7])
pmat = np.array([[c['r']['pct_' + m] for c in G] for m in MK7])
z = (mat - mat.mean(1, keepdims=True)) / (mat.std(1, keepdims=True) + 1e-9)
pd.DataFrame({'column': [g['grp'] + ' ' + g['kind'] for g in G], 'cluster': [g['r'].cl for g in G],
              'n_cells': [int(g['r'].n) for g in G], 'cop_rule_passed': [g['ok'] for g in G],
              **{'mean_' + m: mat[i] for i, m in enumerate(MK7)},
              **{'pct_' + m: pmat[i] for i, m in enumerate(MK7)}}).to_csv(
    os.path.join(OUT, 'figS3g_columns.csv'), index=False)

fig, ax = plt.subplots(figsize=(9.2, 3.0))
norm = matplotlib.colors.TwoSlopeNorm(vcenter=0.0, vmin=-2.2, vmax=2.2)
for j, c in enumerate(G):
    for i, m in enumerate(MK7):
        ax.scatter([j], [len(MK7) - 1 - i], s=2.5 + 40 * (pmat[i, j] / 100.0) ** 1.15,
                   c=[plt.cm.RdBu_r(norm(z[i, j]))], edgecolors='#404040', linewidths=0.25, zorder=3)
    if c['kind'] == 'COP/NFOL':
        ax.add_patch(plt.Rectangle((j - 0.46, -0.42), 0.92, len(MK7) - 0.18, fill=False,
                                   ec='#d51f28' if c['ok'] else '#999999', lw=1.0, zorder=4))
        if not c['ok']:
            ax.text(j, len(MK7) + 0.02, 'fails C1ql1', ha='center', va='bottom',
                    fontsize=5, color='#999999')
    ax.text(j, len(MK7) - 0.55, c['kind'], ha='center', va='bottom', fontsize=5)
    ax.text(j, -0.62, 'n=%d' % int(c['r'].n), ha='center', va='top', fontsize=5,
            color=('#d51f28' if c['ok'] else '#999999') if c['kind'] == 'COP/NFOL' else 'black')
for gi, (nm, _) in enumerate(GROUPS):
    ax.text(gi * 3 + 1, -1.35, nm, ha='center', va='top', fontsize=5.5)
ax.set_xlim(-0.5, len(G) - 0.5); ax.set_ylim(-0.9, len(MK7) + 0.55)
for j in range(len(G) + 1):
    ax.axvline(j - 0.5, lw=0.9 if j % 3 == 0 else 0.35,
               color=(0.72, 0.72, 0.72) if j % 3 == 0 else (0.88, 0.88, 0.88))
ax.set_yticks(np.arange(len(MK7))[::-1]); ax.set_yticklabels(MK7, fontsize=5.5, style='italic')
ax.set_xticks([]); ax.tick_params(length=0)
for sp in ax.spines.values(): sp.set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'figS3g_marker_columns.pdf')); fig.savefig(
    os.path.join(OUT, 'figS3g_marker_columns.png'), dpi=300); plt.close(fig)
print('wrote figS3g_marker_columns.pdf | columns:', [g['grp'] + ':' + g['kind'] for g in G])
show(os.path.join(OUT, 'figS3g_marker_columns.png'))""")

md("""**Expected** — twelve columns (four libraries × COP/NFOL, OPC, MOL). The
dsciCAT COP columns and the pooled-HT COP column share the programme *Gpr17*+ /
*Bmp4*+ / *C1ql1*− / *Ptprz1* mid; the single HT library's best COP candidate
**fails the rule on C1ql1** (its grey box), which is the panel's point: pooling
the HT libraries recovers the canonical marker programme, the single library
does not separate it.""")

# ===========================================================================
md("""---

## 6. Outputs

| File | Panel | Content |
|---|---|---|
| `S3_light.RDS`, `S3_avgexpr_*.csv`, `S3_counts_*.csv`, `S3_summary.txt` | §1 | labels, RNA/ATAC embeddings, type × HVG average expression |
| `figS3ab_umaps.pdf/.png` | A, B | RNA and ATAC UMAPs of the three libraries |
| `figS3cde_xcorr.pdf/.png`, `figS3_xcorr_*.csv`, `figS3_xcorr_stats.csv` | C–E | cluster–cluster correlation + the bootstrap table |
| `S3_lineage_*.rds` | §4.1 | the four lineage subsets (RNA counts) |
| `S3_u13_clusters.csv`, `S3_marker_table.csv`, `S3_strict4029.csv` | §4.2 | every cluster of every condition and seed |
| `figS3f_resolvability.pdf/.png`, `figS3_f_conditions.csv` | F | resolvability vs depth |
| `figS3g_marker_columns.pdf/.png`, `figS3g_columns.csv` | G | marker identity of the twelve clusters |

Run end to end:

```bash
cd figS3
jupyter nbconvert --to notebook --execute --inplace \\
    --ExecutePreprocessor.timeout=-1 Figure_S3.ipynb
```

§1 and §4.1/§4.2 are guarded by their outputs; delete the named files (or set
`FIGS3_FORCE=1`) to redo them. With all caches warm the notebook takes about a
minute; cold, ≈ 12 min.""")

nb.cells = C
if __name__ == '__main__':
    path = os.path.join(OUT_DIR, 'Figure_S3.ipynb')
    nbf.write(nb, path)
    print('wrote', path, len(C), 'cells')
