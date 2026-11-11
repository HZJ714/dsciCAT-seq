#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build Figure_3.ipynb — reproduction of Figure 3 (three-library WNN atlases,
Gpr17+ COP/NFOL, oligodendrocyte-lineage states).

Design (same rules the author set for Figure 1 / Figure 2):
  * reproduction first: every panel states the input files it consumes and the
    shortest complete code that turns them into the panel;
  * styling is deliberately minimal — no journal geometry, no panel letters, no
    compositing (those live in figures/scripts/compose_fig3_v10.py and are not
    part of the science);
  * the one step every panel depends on (three per-library Seurat objects + the
    joint multiomic object -> one small derived cache) is written out in full in
    §1, because panels A-C, the insets and Figure S3 A/B all read its output;
  * provenance is stated honestly: the delivered F panel is a crop of an ArchR
    plotBrowserTrack output, so §5 re-runs that call on the same ArchR project
    instead of imitating the picture, and §6 measures the result against the
    vector files the figure was actually made from;
  * comments in English, code as short as possible while complete.

Run:  python build_fig3_notebook.py
Out:  <repo>/fig3/Figure_3.ipynb
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
md("""# Figure 3 — One identical parameter set applied to three libraries resolves the Gpr17⁺ COP/NFOL state

| Panel | What it shows | Input | Code |
|---|---|---|---|
| A–C | WNN UMAPs of HT-scCAT-seq (4,029), dsciCAT-seq-50k (17,093), dsciCAT-seq-100k (30,350), coloured by cell type; COP/NFOL outlined in B and C | 3 per-library RDS + joint multiomic RDS (§1) | R/Seurat (§2) |
| B/C insets | *Gpr17* feature plots of the re-clustered COP/NFOL subsets (60 / 133 cells) | same cache + RNA counts of those cells | R/Seurat (§2.1) |
| D | UMAP of the oligodendrocyte-lineage subset (5,977 cells, six states) | `3k_dsci_oligo_subset_umap.csv` + `_celltype.csv` | Python (§3) |
| E | 17 lineage markers × six states: dot size = % expressing, colour = average expression | `3k_dsci_oligo_subset_rna.h5ad` + `_celltype.csv` | Python (§4) |
| F | coverage of the five loci of E over the same six states, one row per state, aligned with E | ArchR project of the lineage subset (`3kdsci/archr-subset`) | R/ArchR (§5) |

```
              3k_dsci_ht_10x_gs_subset.RDS   the joint multiomic object Figure 2 ends in (70,643 nuclei)
                        │
                        │  §1.0  one subset() per library:
                        │        EL14  Cells[which(batch      == 'Cortex_EL14')    ] ->  4,029
                        │        FY13  Cells[which(orig.ident == 'dsciCAT-seq-50k') ] -> 17,093
                        │        FY14  Cells[which(orig.ident == 'dsciCAT-seq-100k')] -> 30,350
                        │        = EL14_subset.RDS / 3k_FY13_subset.RDS / 3k_FY14_subset.RDS,
                        │          which carry the same cells, counts and QC metadata plus their own
                        │          library-level PCA/UMAP/graphs/RNA_snn_res.1 and cell_type table
                        ▼
              §1  fig3_light.RDS  (labels + RNA/ATAC/WNN embeddings + counts of the inset cells)
              │        └── §1.0 check: fig3_subset_derivation.csv (cells, order, counts, metadata)
                                │
              ┌─────────────────┴──────────────────┐
              §2  panels A-C (WNN UMAP)            §2.1  Gpr17 insets (subset re-cluster)
              └──── A-C COP/NFOL outline coordinates -> fig3_copnfol_hull.csv

3k_dsci_oligo_subset_umap.csv ┐                 3k_dsci_oligo_subset_rna.h5ad ┐
3k_dsci_oligo_subset_celltype.csv ├─ §3 panel D   3k_dsci_oligo_subset_celltype.csv ├─ §4 panel E
                                  │                                               │
              3kdsci/archr-subset (arrow files + cell_type) ── §5 panel F  (same 6 rows, same order as E)

                                        §6 validation: both panels measured against their own source PDFs
```

Panels D–F are the **same 5,977 lineage cells** and the **same six states in the
same order**, so E and F line up row by row; that alignment is the point of the
panel.""")

# ===========================================================================
md("""---

## 0. Environment, inputs, outputs

| Variable | Default | Contents |
|---|---|---|
| `FIG3_DATA` | `<work>/FY` | the joint RDS object of §0.1 plus the three per-library objects it is subset into (§1.0) |
| `FIG3_OLIGO` | `<work>/FY/oligo` | the oligodendrocyte-lineage subset (RNA h5ad, ATAC h5ad, two CSVs) |
| `FIG3_OUT` | `fig3_out` | created here; every panel PDF/PNG/CSV lands here |
| `FIG3_FORCE` | unset | `1` re-runs §1 even if the cache exists |
| `FIG3_SAVE_DERIVED` | unset | `1` makes §1 also write the plain `subset()`s as `fig3_out/{EL14,FY13,FY14}_from_joint.RDS` (≈ 9 GB) |

Python ≥3.8 with `numpy pandas scipy anndata matplotlib`. §1, §2 and §2.1 need
`Rscript` with `Seurat`, `ggplot2`, `dplyr`, `tidyr`, `ggrepel`, `pals`,
`scales`, `patchwork`, `grid`; §5 needs `ArchR` (mm10). §6 additionally reads
vector PDFs and is skipped when PyMuPDF is not importable.""")

code("""
import matplotlib
matplotlib.use('Agg')
import os, sys, subprocess, time
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager as _fm

def show(*paths):
    # report the saved file only -- nothing is displayed inline (section 0.2)
    for q in paths:
        print('%s  (%.0f KB)' % (q, os.path.getsize(q) / 1024.0))

DATA  = os.environ.get('FIG3_DATA', os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')))
OLIGO = os.environ.get('FIG3_OLIGO', os.path.join(os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')), 'oligo'))
OUT   = os.environ.get('FIG3_OUT', 'fig3_out')
RSCRIPT = os.environ.get('RSCRIPT', 'Rscript')
os.makedirs(OUT, exist_ok=True)

# house style = the only styling in this notebook: Arial-metric font, vector text
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
                     'savefig.dpi': 600})

def write_r(name, src):
    p = os.path.join(OUT, name)
    with open(p, 'w') as fh: fh.write(src)
    return p

def run_r(path, tag=''):
    t0 = time.time()
    rc = subprocess.run([RSCRIPT, path], env=dict(os.environ, FIG3_DATA=DATA,
                                                  FIG3_OLIGO=OLIGO, FIG3_OUT=OUT))
    dt = time.time() - t0
    assert rc.returncode == 0, '%s failed (exit %d)' % (path, rc.returncode)
    print('%s done in %.1f min' % (tag or os.path.basename(path), dt / 60))

# the six lineage states, in delivered order (bottom -> top in E/F, legend in D)
ST  = ['COP/NFOL', 'MFOL', 'MOL2', 'MOL3', 'MOL3/6', 'MOL6']
# one palette for the six states, used by D, E rows and F rows
COL = dict(zip(ST, ['#D51F28', '#272E6A', '#208A42', '#89288F', '#F47D2B', '#EEE500']))
# the 17 markers of panel E, in delivered column order, and the 5 loci of panel F
MK  = ['Gpr17', 'Bmp4', 'Cemip2', '9630013A20Rik', 'Synpr', 'Mbp', 'Mag', 'Apod', 'Trf',
       'Anln', 'Ryr2', 'Nrxn3', 'Mef2c', 'Ptgds', 'Grm3', 'Car2', 'Il33']
LOCI = ['Bmp4', 'Synpr', 'Anln', 'Nrxn3', 'Grm3']
print('DATA =', DATA, '| OLIGO =', OLIGO, '| OUT =', OUT, '| font =', plt.rcParams['font.family'])""")

md("""### 0.1 Input files

| Used by | Path | What it is |
|---|---|---|
| §1 | `FIG3_DATA/3k_dsci_ht_10x_gs_subset.RDS` (13 GB) | **the** input: the joint multiomic object of Figure 2, 70,643 nuclei, four libraries; carries `umap.atac`, `wnn.umap`, and every cell of the three libraries below |
| §1 | `FIG3_DATA/FY13/3k_FY13_subset.RDS` (1.7 GB) | dsciCAT-seq-50k, 17,093 nuclei = `subset(joint, orig.ident == 'dsciCAT-seq-50k')` (§1.0), plus its own library-level analysis: RNA/ATAC assays, `umap`, `RNA_snn_res.1`, `cell_type` |
| §1 | `FIG3_DATA/FY14/3k_FY14_subset.RDS` (2.6 GB) | dsciCAT-seq-100k, 30,350 nuclei = the same with `orig.ident == 'dsciCAT-seq-100k'`, same layout |
| §1 | `FIG3_DATA/HT/EL14_subset.RDS` (0.6 GB) | HT-scCAT-seq, single library, 4,029 nuclei = the same with `batch == 'Cortex_EL14'`, same layout |
| §3 | `FIG3_OLIGO/3k_dsci_oligo_subset_umap.csv` | UMAP_1/UMAP_2 of the 5,977 lineage cells (row names = cell barcodes) |
| §3–§5 | `FIG3_OLIGO/3k_dsci_oligo_subset_celltype.csv` | `cell_id, cluster` — the six-state annotation used by D, E and F |
| §4 | `FIG3_OLIGO/3k_dsci_oligo_subset_rna.h5ad` | RNA of the same 5,977 cells; `.raw` = 31,263 genes × log-normalised data, `.X` = 3,000 HVGs (scaled) |
| §5 | `FIG3_OLIGO/3kdsci/archr-subset/` | ArchR project of the same 5,977 cells: `ArrowFiles/FY13.arrow` + `FY14.arrow`, `Save-ArchR-Project.rds` whose metadata carries `cell_type` (the six states) |
| §6 | `FIG3_OLIGO/3k_dsci_oligo_subset_celltype_dotplot.pdf`, `FIG3_OLIGO/3kdsci/archr-subset/Plots/oligo_subset_celltype_track.pdf` | the vector files the delivered E and F panels were cropped from (validation references) |

The lineage subset of D–F is **not** re-derived here: it is the published
6-state annotation of the dsciCAT oligodendrocyte lineage (Figure 2's `Oligo`,
`OPC` and `COP/NFOL` cells re-clustered into COP/NFOL, MFOL, MOL2, MOL3, MOL3/6,
MOL6). Its three files (coordinates, labels, matrices) share the same 5,977
barcodes, which §0.3 checks.""")

md("""### 0.2 Provenance — read this before comparing with the printed figure

* **The three per-library objects are subsets of the joint object** (§1.0): the
  same cells in the same order, with the same RNA count matrix and the same RNA
  QC metadata, cell for cell; on the ATAC side they share every peak checked but
  the two dsciCAT-seq objects carry 33 peaks more than the joint object, which is
  what moves their `nCount_ATAC` by ≤ 2.1% per cell. What they add on top of the
  subset is the library-level analysis — their
  own PCA/UMAP, neighbour graphs, `RNA_snn_res.1` and the hand-written
  `cell_type` table — which the joint object does not contain (its `cell_type`
  agrees with the library-level call for 75–80% of the cells) and which this
  notebook therefore does not re-derive. `FIG3_SAVE_DERIVED=1` writes the plain
  subsets to `fig3_out/*_from_joint.RDS`.
* **A–C** are the WNN embeddings of the joint multiomic object, restricted to the
  cells of each library, drawn with **one identical set of clustering/annotation
  parameters** for all three libraries. The cell-type labels come from each
  per-library object (`cell_type`, with COP/NFOL taken from `RNA_snn_res.1`
  cluster 29 in FY13 and cluster 31 in FY14) — that re-annotation is the whole
  point of the figure and is written out in §1.
* **The COP/NFOL marking** in B and C rings each cell of that cluster in red and
  writes the cluster's convex hull to `fig3_copnfol_hull.csv`; the delivered
  figure draws a hand-placed ellipse around the same cells.
* **E**: dot size is the percentage of expressing cells and colour is the average
  expression, both recomputed from the h5ad. §6 measures the delivered dots: the
  size channel is this `pct` column, cell for cell (Pearson r = 1.000); the
  colour channel is a per-gene rescaling of an average-expression table that the
  source files do not pin down, so §4 plots the unscaled mean and writes `avg`,
  `pct` **and** the per-gene `z` to `fig3e_marker_stats.csv`.
* **F**: the delivered rows are a crop of
  `3kdsci/archr-subset/Plots/oligo_subset_celltype_track.pdf`, the output of
  `ArchR::plotBrowserTrack` on the project of §5; §5 re-runs that call and §6
  compares the region strings page by page.""")

md("""### 0.3 Input check""")

code("""
need = [('§1 FY14',      os.path.join(DATA,  'FY14/3k_FY14_subset.RDS')),
        ('§1 FY13',      os.path.join(DATA,  'FY13/3k_FY13_subset.RDS')),
        ('§1 EL14',      os.path.join(DATA,  'HT/EL14_subset.RDS')),
        ('§1 joint',     os.path.join(DATA,  '3k_dsci_ht_10x_gs_subset.RDS')),
        ('§3 umap',      os.path.join(OLIGO, '3k_dsci_oligo_subset_umap.csv')),
        ('§3/4/5 labels',os.path.join(OLIGO, '3k_dsci_oligo_subset_celltype.csv')),
        ('§4 rna h5ad',  os.path.join(OLIGO, '3k_dsci_oligo_subset_rna.h5ad')),
        ('§5 archr proj',os.path.join(OLIGO, '3kdsci/archr-subset/Save-ArchR-Project.rds'))]
for k, p in need:
    print('%-16s %-9s %s' % (k, 'ok' if os.path.exists(p) else 'MISSING', p))
assert all(os.path.exists(p) for _, p in need), 'input missing - set FIG3_DATA / FIG3_OLIGO'

_um = pd.read_csv(need[4][1], index_col=0)
_ct = pd.read_csv(need[5][1])
lab = _ct.set_index('cell_id')['cluster']
print('lineage cells: umap %d | labels %d | shared barcodes %d'
      % (len(_um), len(_ct), len(set(_um.index) & set(lab.index))))
print('states:', lab.reindex(_um.index).value_counts().to_dict())""")

md("""**Expected** — all eight inputs present; 5,977 lineage cells in both CSVs with a
100% barcode overlap; MOL6 2,574 · MOL2 1,016 · MOL3/6 1,006 · MOL3 688 ·
MFOL 490 · COP/NFOL 203.""")

md("""### 0.4 Software""")

code("""
import anndata, scipy, matplotlib
print('python   ', sys.version.split()[0])
for m in (numpy := np, pandas := pd, scipy, anndata, matplotlib):
    print('%-9s' % m.__name__.split('.')[0], m.__version__)
subprocess.run([RSCRIPT, '--version'])""")

# ===========================================================================
md("""---

## 1. The step every panel depends on — three libraries, one annotation, one cache

**Input**: the four RDS objects of §0.1. **Output**: `fig3_out/fig3_light.RDS`
(the only thing §2 and §2.1 read) + `fig3_out/fig3_s1_summary.txt`.

Loading the 13 GB joint object is the single expensive read of the whole
notebook, so it is done once and reduced to embeddings + labels — and it is also
where the three per-library objects come from, which is §1.0 below. The rest of
the derivation is written out in full: the COP/NFOL re-annotation from
`RNA_snn_res.1`, the shared palette, the barcode intersection with the joint
object, and the RNA counts of the cells that the B/C insets re-cluster.

### 1.0 Where `EL14_subset.RDS`, `3k_FY13_subset.RDS` and `3k_FY14_subset.RDS` come from

They are not separate data. Figure 2 ends in the joint multiomic object
`3k_dsci_ht_10x_gs_subset.RDS` (70,643 nuclei, four libraries), and one library
out of it is one line of `subset()`:

```r
el14 <- subset(ht, cells = Cells(ht)[which(ht$batch       == 'Cortex_EL14')])    #  4,029
fy13 <- subset(ht, cells = Cells(ht)[which(ht$orig.ident  == 'dsciCAT-seq-50k')]) # 17,093
fy14 <- subset(ht, cells = Cells(ht)[which(ht$orig.ident  == 'dsciCAT-seq-100k')])# 30,350
```

`orig.ident` carries the dsciCAT-seq libraries but lumps all four HT-scCAT-seq
libraries into `HT-scCAT-seq-HT` (14,762 nuclei), so the single-library object
has to be picked by `batch`, which still holds `Cortex_EL13 … EL16`. §1.0 runs
those three lines and checks them against the objects actually used
(`fig3_out/fig3_subset_derivation.csv`, per-column detail in
`fig3_subset_meta_equal.csv`):

* **Inherited, cell for cell**: the same 4,029 / 17,093 / 30,350 barcodes in the
  same order, the RNA count matrix (400 genes × all cells bit-identical), the
  31,263 genes, and every RNA QC column — `nCount_RNA`, `nFeature_RNA`,
  `percent.mt`, `percent.ribo`, `outlier`, `scDblFinder.class`, `orig.ident`,
  `batch`, `rna_celltype` — at 1.0000 agreement.
* **Not exactly inherited on the ATAC side**: the 400 peaks spot-checked are
  identical, but the two dsciCAT-seq objects carry 186,354 peaks where the joint
  object carries 186,321 (EL14 has the joint's set exactly). Those 33 peaks are
  why `nCount_ATAC` / `nFeature_ATAC` differ for 58–64% of the dsciCAT-seq cells —
  always upward, by ≤ 2.1% per cell and 0.09% in total. The joint object is
  quantified on the shared peak set; the panels here never use those two columns.
* **Not inherited at all**: the library-level analysis — each object's own
  `pca`/`umap`/`umap.atac`/`harmony`/`wnn.umap`, its neighbour graphs,
  `RNA_snn_res.1` (agreement 0.16 / 0.36 / 0.52) and the hand-written
  `cell_type` table (0.76 / 0.75 / 0.80). The same pipeline was run on one
  library instead of four and the clusters were then annotated by hand, so none
  of it can be read out of the joint object. That is the one reason §1 still
  reads the three files; `FIG3_SAVE_DERIVED=1` writes the pure subsets to
  `fig3_out/*_from_joint.RDS` if you want to look at them.""")

code("""
R1 = r'''
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(pals)})
DATA <- Sys.getenv('FIG3_DATA', os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')))
OUT  <- Sys.getenv('FIG3_OUT', 'fig3_out')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SAVE_DERIVED <- Sys.getenv('FIG3_SAVE_DERIVED', '') == '1'   # write *_from_joint.RDS

sub <- c(EL14 = 'HT/EL14_subset.RDS', FY13 = 'FY13/3k_FY13_subset.RDS', FY14 = 'FY14/3k_FY14_subset.RDS')
copclust <- c(EL14 = NA_integer_, FY13 = 29L, FY14 = 31L)   # COP/NFOL cluster of each dsciCAT library
# the rule that picks a library out of the joint object (see the text above):
rule <- list(EL14 = quote(batch == 'Cortex_EL14'),
             FY13 = quote(orig.ident == 'dsciCAT-seq-50k'),
             FY14 = quote(orig.ident == 'dsciCAT-seq-100k'))
libs <- list(); log <- character(); chk <- list()
t0 <- proc.time()[3]
ht <- readRDS(file.path(DATA, '3k_dsci_ht_10x_gs_subset.RDS')); DefaultAssay(ht) <- 'RNA'
log <- c(log, sprintf('ht_10x cells=%d loaded %.1f min', ncol(ht), (proc.time()[3] - t0) / 60))
cat(tail(log, 1), '\\n')
pick <- lapply(rule, function(r) Cells(ht)[which(eval(r, ht@meta.data))])
selector <- function(nm) { r <- as.character(rule[[nm]])
                          paste0('ht$', r[2], " == '", r[3], "'") }
for (nm in names(pick)) cat(sprintf('   %s <- subset(ht, cells = Cells(ht)[which(%s)]): %d of %d cells\\n',
                                    nm, selector(nm), length(pick[[nm]]), ncol(ht)))

## one pass over the three libraries: read the delivered object, check that it is
## exactly this subset of the joint object, and keep what the panels need
RNACOLS <- c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo', 'outlier',
             'scDblFinder.class', 'orig.ident', 'batch', 'rna_celltype')
ATACOLS <- c('nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'FRiP', 'nucleosome_signal')
meta_long <- data.frame()
for (nm in names(sub)) {
  t0 <- proc.time()[3]
  o <- readRDS(file.path(DATA, sub[[nm]])); DefaultAssay(o) <- 'RNA'
  # the annotation of the figure: cell_type, but the COP/NFOL cluster by name
  ct <- as.character(o$cell_type)
  if (!is.na(copclust[[nm]])) ct[o$RNA_snn_res.1 == copclust[[nm]]] <- 'COP/NFOL'
  sel <- if (!is.na(copclust[[nm]])) which(o$RNA_snn_res.1 == copclust[[nm]]) else which(ct == 'OPC')
  libs[[nm]] <- list(n = ncol(o), ct = ct, cells = colnames(o),
                     rna = as.data.frame(Embeddings(o, 'umap')),          # panel S3A
                     inset_cells = colnames(o)[sel],
                     inset_counts = as.matrix(GetAssayData(o, assay = 'RNA', slot = 'counts')[, sel]),
                     gpr17 = as.numeric(GetAssayData(o, assay = 'RNA', slot = 'data')['Gpr17', ]))
  log <- c(log, sprintf('%s cells=%d types=%d inset=%d loaded %.1f min', nm, ncol(o),
                        length(unique(ct)), length(sel), (proc.time()[3] - t0) / 60))
  cat(tail(log, 1), '\\n')
  jj <- pick[[nm]][match(colnames(o), pick[[nm]])]        # joint rows, in this object's order
  feq <- function(x, y) identical(as.numeric(as.matrix(x)), as.numeric(as.matrix(y)))
  g  <- intersect(rownames(o[['RNA']]),  rownames(ht[['RNA']]))[1:400]
  gp <- intersect(rownames(o[['ATAC']]), rownames(ht[['ATAC']]))[1:400]
  naeq <- function(a, b) { a <- as.character(a); b <- as.character(b)
                           mean((a == b) | (is.na(a) & is.na(b))) }   # NA == NA counts as equal
  me <- sapply(intersect(colnames(o@meta.data), colnames(ht@meta.data)), function(k)
        naeq(o@meta.data[[k]], ht@meta.data[jj, k]))
  meta_long <- rbind(meta_long, data.frame(library = nm, column = names(me), equal = as.numeric(me)))
  chk[[nm]] <- data.frame(
    library = nm, selector = selector(nm),
    cells = ncol(o), joint_cells = length(pick[[nm]]),
    same_cell_set = setequal(colnames(o), pick[[nm]]), same_cell_order = identical(colnames(o), jj),
    rna_counts_equal = feq(GetAssayData(o, assay = 'RNA',  slot = 'counts')[g, ],
                           GetAssayData(ht, assay = 'RNA',  slot = 'counts')[g, jj]),
    atac_counts_equal = feq(GetAssayData(o, assay = 'ATAC', slot = 'counts')[gp, ],
                            GetAssayData(ht, assay = 'ATAC', slot = 'counts')[gp, jj]),
    rna_genes = nrow(o[['RNA']]), rna_genes_joint = nrow(ht[['RNA']]),
    atac_peaks = nrow(o[['ATAC']]), atac_peaks_joint = nrow(ht[['ATAC']]),
    rna_qc_equal = round(mean(me[intersect(RNACOLS, names(me))]), 4),
    atac_qc_equal = round(mean(me[intersect(ATACOLS, names(me))]), 4),
    cell_type_equal = round(unname(me['cell_type']), 4),
    cluster_equal = round(unname(me['RNA_snn_res.1']), 4), stringsAsFactors = FALSE)
  cat(sprintf('   %s vs subset(ht, ...): cells %s order %s RNA matrix %s ATAC matrix %s | RNA QC meta %.4f | ATAC QC meta %.4f (%d peaks here vs %d in ht) | cell_type %.3f | RNA_snn_res.1 %.3f\\n',
              nm, chk[[nm]]$same_cell_set, chk[[nm]]$same_cell_order, chk[[nm]]$rna_counts_equal,
              chk[[nm]]$atac_counts_equal, chk[[nm]]$rna_qc_equal, chk[[nm]]$atac_qc_equal,
              chk[[nm]]$atac_peaks, chk[[nm]]$atac_peaks_joint, chk[[nm]]$cell_type_equal,
              chk[[nm]]$cluster_equal))
  if (SAVE_DERIVED) {                       # the standalone subset, if you want the file
    d <- subset(ht, cells = pick[[nm]])
    saveRDS(d, file.path(OUT, paste0(nm, '_from_joint.RDS')))
    cat(sprintf('   wrote %s (%d cells, %d genes) - library-level analysis NOT included\\n',
                paste0(nm, '_from_joint.RDS'), ncol(d), nrow(d)))
    rm(d); gc(verbose = FALSE)
  }
  rm(o); gc(verbose = FALSE)
}
chk <- do.call(rbind, chk)
write.csv(chk, file.path(OUT, 'fig3_subset_derivation.csv'), row.names = FALSE)
write.csv(meta_long[order(meta_long$column, meta_long$library), ],
          file.path(OUT, 'fig3_subset_meta_equal.csv'), row.names = FALSE)
log <- c(log, sprintf('derivation check written to %s', file.path(OUT, 'fig3_subset_derivation.csv')))
E <- list(atac = Embeddings(ht, 'umap.atac'), wnn = Embeddings(ht, 'wnn.umap'))
for (nm in names(libs)) {                      # same cells as the RNA panels: barcode match
  idx <- match(libs[[nm]]$cells, colnames(ht))
  libs[[nm]]$atac <- as.data.frame(E$atac[idx, , drop = FALSE])   # panel S3B
  libs[[nm]]$wnn  <- as.data.frame(E$wnn[idx,  , drop = FALSE])   # panels A-C
  log <- c(log, sprintf('%s ht_10x barcode match: %d / %d (%.1f%%)', nm, sum(!is.na(idx)),
                        length(idx), 100 * mean(!is.na(idx))))
  cat(tail(log, 1), '\\n')
}
rm(ht, E); gc(verbose = FALSE)
# one colour per cell type across the three libraries
types <- sort(unique(unlist(lapply(libs, function(s) s$ct))))
cols <- unique(c(as.vector(pals::glasbey(32)), as.vector(pals::alphabet(26))))[1:length(types)]
names(cols) <- types
log <- c(log, sprintf('cell types (union of 3 libraries): %d', length(types)))
saveRDS(list(libs = libs, cols = cols, types = types), file.path(OUT, 'fig3_light.RDS'))
writeLines(log, file.path(OUT, 'fig3_s1_summary.txt'))
cat('wrote', file.path(OUT, 'fig3_light.RDS'), '\\n')
'''
r1 = write_r('fig3_s1_object.R', R1)
cache = os.path.join(OUT, 'fig3_light.RDS')
if os.path.exists(cache) and os.environ.get('FIG3_FORCE') != '1':
    print('cache exists, skipped §1:', cache, '(%.1f MB) - delete it or set FIG3_FORCE=1'
          % (os.path.getsize(cache) / 1e6))
    print(open(os.path.join(OUT, 'fig3_s1_summary.txt')).read())
else:
    run_r(r1, '§1')""")

md("""**Expected** — the §1.0 check to come out `TRUE` for `same_cell_set`,
`same_cell_order`, `rna_counts_equal`, `atac_counts_equal` and `rna_qc_equal =
1.0` on all three libraries, with `atac_qc_equal` 1.0 for EL14 and 0.66 / 0.68 for
FY13 / FY14 (the 33 extra peaks), `cell_type` 0.76 / 0.75 / 0.80 and
`RNA_snn_res.1` 0.16 / 0.36 / 0.52 because the library-level annotation is the
one thing a subset cannot carry. Then EL14 4,029 cells / 22 types / 76 OPC
cells for the inset control; FY13 17,093 / 32 / 60; FY14 30,350 / 33 / 133; joint
object 70,643 cells; **100% barcode match** for all three libraries (so A–C, S3A
and S3B are literally the same cells); 36 cell types in the union. Measured on
this machine, 2026-10-22 (with the derivation checks on): the 13 GB joint object loads in
**7.2 min**, the three library objects in 2.6 min total, whole step **10.3 min**.""")

# ===========================================================================
md("""---

## 2. Panels A–C — WNN UMAPs of the three libraries

**Input**: `fig3_light.RDS` (§1). **Output**:
`fig3_out/fig3abc_wnn_umaps.pdf` (+ `.png`) and `fig3_out/fig3_copnfol_hull.csv`.

Each library's cells are placed by the **joint WNN embedding**, coloured by that
library's own cell-type call, with the shared palette. Point size scales with
cell number (0.07 / 0.055 / 0.05 for 4k / 17k / 30k cells) so that the three
panels have the same ink density; labels are repelled to the median of each type.
The COP/NFOL cluster of B and C is ringed in red, cell by cell, and its outline
(convex hull of the cluster's WNN coordinates) plus the barcodes that generated
it are written to `fig3_copnfol_hull.csv` (EL14 has no COP/NFOL cells, so A has
no ring).""")

code("""
R2 = r'''
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)
  library(tidyr); library(ggrepel); library(grid)})
OUT <- Sys.getenv('FIG3_OUT', 'fig3_out')
S   <- readRDS(file.path(OUT, 'fig3_light.RDS')); libs <- S$libs; cols <- S$cols
titles <- c(EL14 = 'HT-scCAT-seq (n = 4,029) WNN', FY13 = 'dsciCAT-seq-50k (n = 17,093) WNN',
            FY14 = 'dsciCAT-seq-100k (n = 30,350) WNN')
pt    <- c(EL14 = 0.07, FY13 = 0.055, FY14 = 0.05)

# coordinates are drawn on a square canvas: normalise each library to the same span
norm2 <- function(coords, span = 2.4) {
  x <- coords[, 1]; y <- coords[, 2]; s <- max(diff(range(x)), diff(range(y)))
  data.frame(x = (x - mean(range(x))) / s * span, y = (y - mean(range(y))) / s * span)
}
make_panel <- function(df, letter, ttl, mark = NULL) {
  lab <- df %>% group_by(cell_type) %>% summarise(x = median(x), y = median(y), .groups = 'drop')
  ggplot(df, aes(x, y)) +
    geom_point(aes(color = cell_type), size = pt[[df$lib[1]]], stroke = 0, alpha = 0.95) +
    scale_color_manual(values = cols, guide = 'none') +
    geom_text_repel(data = lab, aes(x, y, label = cell_type, color = cell_type),
                    size = 5 / 2.845276, fontface = 'bold', box.padding = unit(0.25, 'lines'),
                    point.padding = unit(0.2, 'lines'), min.segment.length = 0.1,
                    segment.size = 0.3, segment.alpha = 0.7, max.overlaps = Inf,
                    force = 4, force_pull = 1.2, bg.color = 'white', bg.r = 0.15, seed = 42) +
    coord_fixed() + xlim(-1.7, 1.7) + ylim(-1.7, 1.7) + theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0)) +
    annotate('text', x = -1.6, y = 1.62, label = letter, size = 6 / 2.845276, hjust = 0) +
    annotate('text', x = -1.42, y = 1.62, label = ttl, size = 5 / 2.845276, hjust = 0) +
    if (!is.null(mark)) geom_point(data = mark, aes(x, y), inherit.aes = FALSE, shape = 21,
                                   colour = 'red', fill = NA, size = 1.6, linewidth = 0.25)
}
ps <- list(); hull <- NULL
for (i in seq_along(libs)) {
  nm <- names(libs)[i]
  df <- norm2(libs[[nm]]$wnn); df$cell_type <- libs[[nm]]$ct; df$lib <- nm
  g <- df[df$cell_type == 'COP/NFOL', ]            # outline of the COP/NFOL cluster
  k <- if (nrow(g) > 3) grDevices::chull(g$x, g$y) else NULL
  if (!is.null(k)) hull <- rbind(hull, data.frame(lib = nm, x = g$x[k], y = g$y[k],
                                                  cell = libs[[nm]]$cells[df$cell_type == 'COP/NFOL'][k],
                                                  n = nrow(g)))
  ps[[i]] <- make_panel(df, LETTERS[i], titles[[nm]], if (is.null(k)) NULL else g)
}
write.csv(hull, file.path(OUT, 'fig3_copnfol_hull.csv'), row.names = FALSE)
cat('COP/NFOL cells: FY13', sum(libs$FY13$ct == 'COP/NFOL'), '| FY14',
    sum(libs$FY14$ct == 'COP/NFOL'), '| EL14', sum(libs$EL14$ct == 'COP/NFOL'), '\\n')
draw <- function(g, x, y, w, h) { pushViewport(viewport(x = x, y = y, width = w, height = h,
  just = c('left', 'bottom'), default.units = 'cm')); grid.draw(g); popViewport() }
pdf(file.path(OUT, 'fig3abc_wnn_umaps.pdf'), width = 15 / 2.54, height = 5 / 2.54)
grid.newpage()
for (i in 1:3) draw(ggplotGrob(ps[[i]]), (i - 1) * 5, 0, 5, 5)
dev.off()
png(file.path(OUT, 'fig3abc_wnn_umaps.png'), width = 15, height = 5, units = 'cm', res = 300)
grid.newpage()
for (i in 1:3) draw(ggplotGrob(ps[[i]]), (i - 1) * 5, 0, 5, 5)
dev.off()
cat('wrote fig3abc_wnn_umaps.pdf\\n')
'''
run_r(write_r('fig3_s2_umaps.R', R2), '§2')""")

code("""png = os.path.join(OUT, 'fig3abc_wnn_umaps.png')
if os.path.exists(png): show(png)
h = pd.read_csv(os.path.join(OUT, 'fig3_copnfol_hull.csv'))
print('COP/NFOL cells ringed on the panels, hull written here:')
print(h.groupby('lib')['n'].first().to_string())""")

md("""**Expected** — three 5 × 5 cm WNN panels titled by library; COP/NFOL circled in
FY13 (60 cells) and FY14 (133 cells) and **absent from EL14** (0 cells) — that
absence at identical parameters is the claim of the figure. Measured here:
`FY13 60 | FY14 133 | EL14 0`.""")

# ===========================================================================
md("""---

## 2.1 Insets of B and C — *Gpr17* in the re-clustered COP/NFOL subsets

**Input**: `inset_counts` of §1 (60 / 133 cells, plus the 76 EL14 OPC cells as a
control). **Output**: `fig3_out/fig3bc_gpr17_insets.pdf` (+ `.png`).

The COP/NFOL cells of each library are taken out and clustered again with the
same default parameters used everywhere in this project (3,000 HVG → PCA 30 →
10 neighbours → SLM `algorithm = 3`, resolutions 0.1/0.3/0.5 → UMAP on PCA 1:10,
`seed.use = 42`), then coloured by *Gpr17* log-normalised expression on a
**common** grey→red scale shared by all three panels, so the three are directly
comparable. The delivered figure shows the two dsciCAT insets (n = 60, n = 133);
the EL14 OPC cells (76) are put through the identical call here as a control.""")

code("""
R3 = r'''
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)})
OUT <- Sys.getenv('FIG3_OUT', 'fig3_out')
S   <- readRDS(file.path(OUT, 'fig3_light.RDS')); libs <- S$libs
# common colour limit: max Gpr17 over the inset cells of all three libraries, rounded up
max_gpr <- ceiling(max(sapply(names(libs), function(k) max(libs[[k]]$gpr17[match(libs[[k]]$inset_cells,
                                                                                 libs[[k]]$cells)]))) * 10) / 10
cat('common Gpr17 scale: 0 -', max_gpr, '\\n')
ps <- list(); info <- character()
for (nm in names(libs)) {
  s <- libs[[nm]]
  sub <- CreateSeuratObject(s$inset_counts, project = nm); DefaultAssay(sub) <- 'RNA'
  sub <- NormalizeData(sub, verbose = FALSE)
  sub <- suppressWarnings(FindVariableFeatures(sub, nfeatures = min(3000, nrow(sub)), verbose = FALSE))
  sub <- ScaleData(sub, verbose = FALSE)
  sub <- RunPCA(sub, npcs = 30, verbose = FALSE)
  sub <- FindNeighbors(sub, reduction = 'pca', dims = 1:10, verbose = FALSE)
  sub <- FindClusters(sub, algorithm = 3, resolution = c(0.1, 0.3, 0.5), verbose = FALSE)
  sub <- RunUMAP(sub, reduction = 'pca', dims = 1:10, seed.use = 42, verbose = FALSE)
  tb <- sort(table(sub$RNA_snn_res.0.5), decreasing = TRUE)
  info <- c(info, sprintf('%s n=%d -> %d cluster(s) at res 0.5 (%s), Gpr17 mean %.2f', nm,
           ncol(sub), length(tb), paste(tb, collapse = ','),
           mean(GetAssayData(sub, assay = 'RNA', slot = 'data')['Gpr17', ])))
  df <- data.frame(x = Embeddings(sub, 'umap')[, 1], y = Embeddings(sub, 'umap')[, 2],
                   gpr = as.numeric(GetAssayData(sub, assay = 'RNA', slot = 'data')['Gpr17', ]))
  ps[[nm]] <- ggplot(df, aes(x, y)) +
    geom_point(aes(color = gpr), size = 1.6, alpha = 0.95) +
    scale_color_gradientn(colours = c('grey80', '#ffcccc', '#ff0000', '#b30000'),
                          limits = c(0, max_gpr), oob = scales::squish, name = 'Gpr17 (log-norm)') +
    labs(title = sprintf('%s (n = %d, re-clustered)', nm, nrow(df))) +
    coord_fixed() + theme_void() +
    theme(plot.title = element_text(size = 6, hjust = 0.5),
          legend.key.size = unit(2.2, 'mm'), legend.text = element_text(size = 5),
          legend.title = element_text(size = 5.5))
}
writeLines(info, file.path(OUT, 'fig3_insets_summary.txt')); cat(paste(info, collapse = '\\n'), '\\n')
fig <- patchwork::wrap_plots(ps, ncol = 3)
ggsave(file.path(OUT, 'fig3bc_gpr17_insets.pdf'), fig, width = 15, height = 4.2, units = 'cm')
ggsave(file.path(OUT, 'fig3bc_gpr17_insets.png'), fig, width = 15, height = 4.2, units = 'cm', dpi = 300)
cat('wrote fig3bc_gpr17_insets.pdf\\n')
'''
run_r(write_r('fig3_s21_insets.R', R3), '§2.1')""")

code("""png = os.path.join(OUT, 'fig3bc_gpr17_insets.png')
if os.path.exists(png): show(png)
print(open(os.path.join(OUT, 'fig3_insets_summary.txt')).read())""")

md("""**Expected** — FY13's 60 COP/NFOL cells stay one cluster at res 0.5, FY14's 133
split into two (81 / 52), EL14's 76 OPC cells stay one cluster. Measured *Gpr17*
means on the common 0-3.4 scale: EL14 0.67, FY13 0.71, FY14 0.93 — i.e. the inset
shows how *Gpr17* is distributed inside each re-clustered subset, and the EL14
OPC control is **not** negative for it; what separates the libraries is that
EL14's identical parameters never resolve a COP/NFOL cluster in the first place
(panel A). The delivered figure shows the two dsciCAT insets only.""")

# ===========================================================================
md("""---

## 3. Panel D — UMAP of the oligodendrocyte-lineage subset

**Input**: `3k_dsci_oligo_subset_umap.csv` + `3k_dsci_oligo_subset_celltype.csv`.
**Output**: `fig3_out/fig3d_lineage_umap.pdf` (+ `.png`).""")

code("""
_um = pd.read_csv(os.path.join(OLIGO, '3k_dsci_oligo_subset_umap.csv'), index_col=0)
_ct = pd.read_csv(os.path.join(OLIGO, '3k_dsci_oligo_subset_celltype.csv')).set_index('cell_id')['cluster']
df = _um.loc[_um.index.intersection(_ct.index)].join(_ct.rename('state'))
df['state'] = pd.Categorical(df['state'], ST, ordered=True)
fig, ax = plt.subplots(figsize=(3.4, 3.4))
for s in ST:                                            # one colour per state, no re-clustering here
    g = df[df.state == s]
    ax.scatter(g.UMAP_1, g.UMAP_2, s=0.45, c=COL[s], linewidths=0, alpha=0.9,
               label='%s (%d)' % (s, len(g)))
ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values(): sp.set_visible(False)
ax.legend(fontsize=5, markerscale=8, frameon=False, handletextpad=0.2, columnspacing=0.4, loc='center left',
          bbox_to_anchor=(1.0, 0.5))
ax.set_title('dsciCAT-seq oligodendrocyte lineage (n = %d)' % len(df), fontsize=6.5)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fig3d_lineage_umap.pdf')); fig.savefig(
    os.path.join(OUT, 'fig3d_lineage_umap.png'), dpi=300); plt.close(fig)
df.reset_index().rename(columns={'index': 'cell_id'}).to_csv(os.path.join(OUT, 'fig3d_umap_states.csv'),
                                                             index=False)
print('wrote fig3d_lineage_umap.pdf |', len(df), 'cells |', df.state.value_counts().to_dict())""")

code("""show(os.path.join(OUT, 'fig3d_lineage_umap.png'))""")

md("""**Expected** — 5,977 cells in six states (COP/NFOL 203 · MFOL 490 · MOL2 1,016 ·
MOL3 688 · MOL3/6 1,006 · MOL6 2,574); a single continuous sheet with the
COP/NFOL end narrow and the mature MOL6 end broad.""")

# ===========================================================================
md("""---

## 4. Panel E — 17 lineage markers × six states

**Input**: `3k_dsci_oligo_subset_rna.h5ad` (`.raw` = log-normalised data for
31,263 genes) + the same label file. **Output**:
`fig3_out/fig3e_marker_dotplot.pdf` (+ `.png`) and `fig3e_marker_stats.csv`.

Two numbers per cell of the grid, both recomputed here:
`pct` = % of the state's cells with the gene detected, `avg` = mean expression of
the state. Dot size uses `pct` (the channel that can be checked against the
delivered figure), colour uses `avg`.""")

code("""
import anndata as ad, scipy.sparse as sprs
a = ad.read_h5ad(os.path.join(OLIGO, '3k_dsci_oligo_subset_rna.h5ad'))
raw = a.raw.to_adata()
lab = pd.read_csv(os.path.join(OLIGO, '3k_dsci_oligo_subset_celltype.csv')).set_index('cell_id')['cluster']
lab = lab.reindex(raw.obs_names)
missing = [g for g in MK if g not in set(raw.var_names)]
assert not missing, 'markers absent from the RNA matrix: %s' % missing
R = raw[:, MK].X
R = np.asarray(R.todense()) if sprs.issparse(R) else np.asarray(R)
rows = []
for s in ST:
    X = R[(lab == s).values]
    rows.append(pd.DataFrame({'state': s, 'gene': MK, 'n': len(X),
                              'avg': X.mean(0), 'pct': 100 * (X > 0).mean(0)}))
E = pd.concat(rows, ignore_index=True)
E['z'] = E.groupby('gene')['avg'].transform(lambda v: (v - v.mean()) / (v.std(ddof=0) or 1))
E.to_csv(os.path.join(OUT, 'fig3e_marker_stats.csv'), index=False)

fig, ax = plt.subplots(figsize=(4.6, 2.6))
xs, ys = np.arange(len(MK)), np.arange(len(ST))
sz = E.pct.values.reshape(len(ST), len(MK)) ** 0.6 * 22      # size <- percent expressing
cl = E.avg.values.reshape(len(ST), len(MK))
sc = ax.scatter(np.tile(xs, len(ST)), np.repeat(ys, len(MK)), s=sz.ravel(),
                c=cl.ravel(), cmap='Reds', vmin=0, vmax=np.nanmax(cl), edgecolor='none')
ax.set_xticks(xs); ax.set_xticklabels(MK, rotation=90, fontsize=5.5)
ax.set_yticks(ys); ax.set_yticklabels(ST, fontsize=5.5)
ax.set_xlim(-0.7, len(MK) - 0.3); ax.set_ylim(-0.7, len(ST) - 0.3); ax.invert_yaxis()
for sp in ax.spines.values(): sp.set_visible(False)
ax.tick_params(length=0)
cb = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.02); cb.ax.tick_params(labelsize=5)
cb.set_label('average expression', fontsize=5.5)
ax.set_title('17 lineage markers, mean log-normalised RNA and % positive cells', fontsize=6.5)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fig3e_marker_dotplot.pdf')); fig.savefig(
    os.path.join(OUT, 'fig3e_marker_dotplot.png'), dpi=300); plt.close(fig)
print('wrote fig3e_marker_dotplot.pdf')
piv = E.pivot(index='state', columns='gene', values='pct').loc[ST, MK]
print('percent expressing:'); print(piv.round(0).to_string())""")

code("""show(os.path.join(OUT, 'fig3e_marker_dotplot.png'))""")

md("""**Expected** — the marker block of each state comes out where it should:
COP/NFOL *Gpr17* 50% / *Bmp4* 28% / *Cemip2* 35% / *9630013A20Rik* 87% with the
mature genes low there (*Grm3* 38%, *Il33* 9%); MOL2 *Anln* 72% / *Apod* 81% /
*Trf* 88%; MOL3 *Nrxn3* 92% / *Ryr2* 61%; MOL6 *Ptgds* 92% / *Grm3* 85% /
*Il33* 43%; *Mbp* ≥ 79% in every state (myelin genes are already on across this
subset).

**Verification.** §6 measures the 102 dots of the vector file the delivered panel
was cropped from: their radii are `r(cm) = 0.00170 × pct + 0.0060` with **Pearson
r = 1.000** against the `pct` column above, i.e. the printed size channel is
exactly this table. The printed colour channel is a per-gene rescaling of an
average-expression table (it is not a monotone function of the `avg` column —
*Mbp* is at the colour minimum in MOL2 although its mean here is the second
highest of the six states), so this notebook plots the unscaled mean and leaves
the scaling to whoever composes the figure; `fig3e_marker_stats.csv` carries
`avg`, `pct` and the per-gene `z` side by side.""")

# ===========================================================================
md("""---

## 5. Panel F — ATAC coverage of the five loci of E, one row per state

**Input**: the ArchR project of the same 5,977 lineage cells —
`FIG3_OLIGO/3kdsci/archr-subset` (`ArrowFiles/FY13.arrow`, `FY14.arrow`,
`Save-ArchR-Project.rds`, whose `cell_type` column carries the six states).
**Output**: `fig3_out/fig3f_archr_tracks.pdf` (one page per locus) plus one PNG
per page.

The delivered F panel is a crop of
`FIG3_OLIGO/3kdsci/archr-subset/Plots/oligo_subset_celltype_track.pdf`, which is
itself the direct output of the call below: pseudo-bulk coverage of the cells of
each state, normalised by reads-in-TSS, drawn over the gene model at TSS ± 50 kb.
So this section re-runs the original computation rather than imitating its
picture: the same project, the same grouping, the same window, the same
normalisation. §6 checks the regenerated pages against the delivered source page
by page.""")

code("""
R5 = r'''
suppressPackageStartupMessages({library(ArchR); library(grid)})
addArchRGenome('mm10'); addArchRThreads(threads = 8)
OLIGO <- Sys.getenv('FIG3_OLIGO', os.path.join(os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')), 'oligo'))
OUT   <- Sys.getenv('FIG3_OUT', 'fig3_out')
proj  <- loadArchRProject(file.path(OLIGO, '3kdsci/archr-subset'), showLogo = FALSE)
ST6 <- c('COP/NFOL', 'MFOL', 'MOL2', 'MOL3', 'MOL3/6', 'MOL6')
stopifnot(all(ST6 %in% unique(proj$cell_type)))
cat('project cells:', nCells(proj), '| groups:',
    paste(names(table(proj$cell_type)), table(proj$cell_type), sep = '=', collapse = ' '), '\\n')
p <- plotBrowserTrack(ArchRProj = proj, groupBy = 'cell_type',
                      geneSymbol = c('Bmp4', 'Synpr', 'Anln', 'Nrxn3', 'Grm3'),
                      upstream = 50000, downstream = 50000,
                      plotSummary = c('bulkTrack', 'geneTrack'),
                      features = NULL, loops = NULL, verbose = FALSE)
pdf(file.path(OUT, 'fig3f_archr_tracks.pdf'), width = 10, height = 14)
for (i in seq_along(p)) { grid.newpage(); grid.draw(p[[i]]) }
dev.off()
for (i in seq_along(p)) {
  png(file.path(OUT, sprintf('fig3f_archr_track_%d.png', i)),
      width = 10, height = 14, units = 'in', res = 150)
  grid.newpage(); grid.draw(p[[i]]); dev.off()
}
cat('wrote fig3f_archr_tracks.pdf with', length(p), 'pages\\n')
'''
run_r(write_r('fig3_s5_tracks.R', R5), '§5')""")

code("""
import re
def track_info(pdf_path):
    \"\"\"-> [(locus, region, coverage range, row labels)] read off an ArchR track PDF.\"\"\"
    try:
        import fitz
    except ImportError:
        print('skipped: PyMuPDF not importable (run with PYTHONPATH=<work>/pylibs)')
        return None
    out = []
    doc = fitz.open(pdf_path)
    for pg in doc:
        t = pg.get_text().split('\\n')
        reg = next((x.strip() for x in t if re.match(r'^chr[0-9XY]+:', x)), None)
        rng = next((x for x in t if 'Signal Range' in x), '')
        rows = [x for x in t if x in ST]
        out.append((reg, rng, rows))
    return out

info = track_info(os.path.join(OUT, 'fig3f_archr_tracks.pdf'))
if info:
    for i, (reg, rng, rows) in enumerate(info, 1):
        print('page %d  %-28s rows=%d  %s' % (i, reg, len(rows), rng.strip()))
    show(os.path.join(OUT, 'fig3f_archr_track_1.png'))""")

md("""**Expected** — 5 pages, one per locus, each with the six state rows in `ST` order
and a gene-model row; regions `chr14:46340668−46440669` (Bmp4),
`chr14:13234779−13334780` (Synpr), `chr9:22339205−22439206` (Anln),
`chr12:88672875−88772876` (Nrxn3), `chr5:9675351−9775352` (Grm3) — i.e. TSS ± 50 kb
of each gene, exactly the windows of the delivered crop. Measured here: coverage
ranges 0–2.09, 0–1.34, 0–1.65, 0–0.28, 0–1.58 against 0–2.07, 0–1.34, 0–1.57,
0–0.28, 0–1.53 in the delivered source (ArchR's per-group normalisation moves by
a few percent between runs; the tracks themselves are identical).""")

# ===========================================================================
md("""---

## 6. Validation against the files the figure was made from

Two checks, both against vector files, both optional (they print a skip message
when the reference file or PyMuPDF is absent):

* **E** — the figure's own source dot plot,
  `FIG3_OLIGO/3k_dsci_oligo_subset_celltype_dotplot.pdf`, holds 102 dots. Their
  radii are measured and compared with the `pct` column of
  `fig3e_marker_stats.csv`; their fills are decoded through the panel's colour bar
  and compared with the `avg` column.
* **F** — the region strings and row labels of §5 are compared with the delivered
  source track PDF, page by page.""")

code("""
SRC_E = os.path.join(OLIGO, '3k_dsci_oligo_subset_celltype_dotplot.pdf')
SRC_F = os.path.join(OLIGO, '3kdsci/archr-subset/Plots/oligo_subset_celltype_track.pdf')
try:
    import fitz
except ImportError:
    print('validation skipped: PyMuPDF not importable (PYTHONPATH=<work>/pylibs)')
else:
    from scipy.stats import pearsonr

    # ---- E: measure the 102 dots of the source dot plot -------------------
    if os.path.exists(SRC_E):
        pg = fitz.open(SRC_E)[0]
        CM = 72 / 2.54
        dots = []
        for g in pg.get_drawings():
            for it in g['items']:
                if it[0] == 'c' and g.get('fill'):
                    p = np.array([[q.x, q.y] for q in it[1:5]]) / CM
                    dots.append((p[:, 0].mean(), p[:, 1].mean(),
                                 (p[:, 0].max() - p[:, 0].min()) / 2,
                                 np.array(g['fill'][:3], float) * 255))
        dots = [d for d in dots if d[0] < 33]                 # drop the legend dots
        def cluster(vals, tol=1.0):
            out = []
            for v in sorted(vals):
                if out and v - np.mean(out[-1]) <= tol: out[-1].append(v)
                else: out.append([v])
            return [float(np.mean(o)) for o in out]
        CX, CY = cluster([d[0] for d in dots]), cluster([d[1] for d in dots])
        assert (len(CX), len(CY)) == (len(MK), len(ST)), (len(CX), len(CY))
        RAD = np.full((len(ST), len(MK)), np.nan)
        for x, y, r, _ in dots:
            i, j = int(np.argmin([abs(y - v) for v in CY])), int(np.argmin([abs(x - v) for v in CX]))
            RAD[i, j] = r
        st = pd.read_csv(os.path.join(OUT, 'fig3e_marker_stats.csv'))
        pct = st.pivot(index='state', columns='gene', values='pct').loc[ST, MK].values
        rr, kk = np.polyfit(pct.ravel(), RAD.ravel(), 1)
        print('E  size channel: source dot radius (cm) = %.5f x pct + %.4f | Pearson r = %.4f'
              % (rr, kk, pearsonr(pct.ravel(), RAD.ravel()).statistic))
        print('   -> the printed size channel is this table pct, cell for cell')
        # colour channel: decode the fills through the embedded colour bar
        img = None
        try:
            from PIL import Image as PILImage
            import io
            raw = fitz.open(SRC_E).extract_image(pg.get_images(full=True)[0][0])
            img = np.array(PILImage.open(io.BytesIO(raw['image'])).convert('RGB')).astype(float)
        except Exception as e:
            print('   colour bar not decodable (%s)' % e)
        if img is not None:
            bb = pg.get_image_info()[0]['bbox']
            ticks = sorted({(round((s['bbox'][1] + s['bbox'][3]) / 2 / CM, 2), s['text'].strip())
                            for bl in pg.get_text('dict')['blocks'] for ln in bl.get('lines', [])
                            for s in ln['spans']
                            if s['text'].strip() in ('0.0', '0.5', '1.0', '1.5', '2.0')
                            and s['bbox'][0] / CM > 33})
            a, b = np.polyfit([t[0] for t in ticks], [float(t[1]) for t in ticks], 1)
            ytop, ybot = bb[1] / CM, bb[3] / CM
            lutv = a * (ytop + (ybot - ytop) * np.arange(img.shape[0]) / (img.shape[0] - 1)) + b
            COLR = np.zeros((len(ST), len(MK), 3))
            for x, y, r, f in dots:
                i, j = int(np.argmin([abs(y - v) for v in CY])), int(np.argmin([abs(x - v) for v in CX]))
                COLR[i, j] = f
            V = np.array([[lutv[int(np.argmin(((img[:, 0, :] - c) ** 2).sum(1)))]
                           for c in row] for row in COLR])
            avg = st.pivot(index='state', columns='gene', values='avg').loc[ST, MK].values
            z = st.pivot(index='state', columns='gene', values='z').loc[ST, MK].values
            print('E  colour channel: decoded colour-bar value vs recomputed statistics')
            print('   vs avg      Pearson %.3f | vs per-gene z  Pearson %.3f'
                  % (pearsonr(V.ravel(), avg.ravel()).statistic,
                     pearsonr(V.ravel(), z.ravel()).statistic))
            print('   -> the printed colour channel is NOT a monotone function of this')
            print('      avg table (Mbp sits at the colour minimum in MOL2 although its')
            print('      mean log-normalised RNA is the 2nd highest of the six states);')
            print('      the delivered figure applies a per-gene rescaling of a related')
            print('      table that the source files do not pin down. Size is exact,')
            print('      colour is the recomputed mean - see fig3e_marker_stats.csv.')
    else:
        print('E source dot plot not found, skipped:', SRC_E)

    # ---- F: region strings, page by page ---------------------------------
    if info and os.path.exists(SRC_F):
        ref = {r: (rng, rows) for r, rng, rows in track_info(SRC_F) if r}
        hit = sum(1 for r, _, _ in info if r in ref)
        print('F  regions: %d of %d regenerated pages found in the delivered source (%d pages)'
              % (hit, len(info), len(ref)))
        for r, rng, rows in info:
            print('   %-28s %s | rows %d' % (r, 'same region' if r in ref else 'NOT IN SOURCE', len(rows)))
    else:
        print('F source track PDF not found, skipped:', SRC_F)""")

# ===========================================================================
md("""---

## 7. Outputs

| File | Panel | Content |
|---|---|---|
| `fig3_light.RDS`, `fig3_s1_summary.txt` | §1 | labels, three embeddings per library, counts of the inset cells |
| `fig3_subset_derivation.csv`, `fig3_subset_meta_equal.csv` | §1.0 | per library: the `subset()` selector, cell counts, and the equality checks against the joint object (cell set, order, RNA/ATAC matrices, per-column metadata agreement) |
| `{EL14,FY13,FY14}_from_joint.RDS` | §1.0 | only with `FIG3_SAVE_DERIVED=1`: the plain subsets, without the library-level analysis |
| `fig3abc_wnn_umaps.pdf/.png`, `fig3_copnfol_hull.csv` | A–C | three WNN UMAPs + the COP/NFOL outline coordinates |
| `fig3bc_gpr17_insets.pdf/.png`, `fig3_insets_summary.txt` | B/C insets | *Gpr17* on the re-clustered subsets |
| `fig3d_lineage_umap.pdf/.png`, `fig3d_umap_states.csv` | D | lineage UMAP + the coordinates/labels used |
| `fig3e_marker_dotplot.pdf/.png`, `fig3e_marker_stats.csv` | E | dot plot + `avg`, `pct`, per-gene `z` |
| `fig3f_archr_tracks.pdf`, `fig3f_archr_track_1..5.png` | F | coverage tracks, one page per locus |
| `fig3_s*.R` | all | the R sources exactly as run |

Run end to end:

```bash
cd fig3
PYTHONPATH=<work>/pylibs jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=-1 Figure_3.ipynb
```

(`PYTHONPATH` only feeds §6 the PyMuPDF build in `<work>/pylibs`; without it
§1–§5 still run and §6 prints a skip message.)

§1 is cached: delete `fig3_out/fig3_light.RDS` (or set `FIG3_FORCE=1`) to redo
the 13 GB read. §1 is also the only step that reads the per-library objects: it
first checks them against `subset(joint, ...)` (§1.0), so the three RDS files are
a convenience copy of the library-level analysis, not extra data. Everything after §1 reads only the 9 MB cache, the three lineage
files and the ArchR project, and takes a few minutes.""")

nb.cells = C
path = os.path.join(OUT_DIR, 'Figure_3.ipynb')
nbf.write(nb, path)
print('wrote', path, len(C), 'cells')
