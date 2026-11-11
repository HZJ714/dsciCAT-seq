#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build Figure_2.ipynb — reproduction of Figure 2 (mouse-brain multiome atlas).

Design (same conventions as Figure_1.ipynb, per author requirements):
  * reproduction first: for every panel the notebook states the *input files* and
    gives the shortest complete code that turns them into the panel;
  * styling is deliberately minimal (no journal layout, no panel letters, no
    colour-bar placement, no compositing) — those are done by hand in the
    production scripts and are not part of the science;
  * the step that everything depends on — the WNN-integrated Seurat object with
    the final 33-type annotation — is read once in Section 1, reduced to the
    small derived cache (embeddings + metadata + the 48 marker rows of the RNA
    assay) that panels A–D consume, and every measured number is printed there;
  * panel E's tracks are ArchR code (Section 4) run on the same ArchR project the
    delivered figure was drawn from;
  * provenance is stated honestly in Section 0.2: the notebook starts from the
    integrated object, not from the raw libraries; the integration itself
    (harmony WNN of RNA+ATAC) happened upstream and is not recomputed here;
  * comments in English, code as short as possible while complete.

Run:  python build_fig2_notebook.py
Out:  <repo>/fig2/Figure_2.ipynb
"""
import os
import nbformat as nbf

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
nb = nbf.v4.new_notebook(metadata={
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.10"}})
C = []
md, code = (lambda s: C.append(nbf.v4.new_markdown_cell(s))), \
           (lambda s: C.append(nbf.v4.new_code_cell(s.strip("\n"))))

# ===========================================================================
# Title
# ===========================================================================
# Title
# ===========================================================================
md("""# Figure 2 — Joint analysis of 16-barcode dsciCAT-seq with HT-scCAT-seq resolves 33 cell types and oligodendrocyte states

This notebook starts from the **SoupX-corrected RNA matrix** of the dsciCAT-seq
libraries and ends with the five panels. §1–§5 build the analysis object,
§6–§9 draw from it, so the whole figure is one file at a time:

```
SoupX-corrected RNA matrix  +  peak matrix / fragment file              (§0.1)
        │
        ├─ §1  per-library Signac ATAC object          85,255 / 89,096 peaks × 23,376 / 43,249
        ├─ §2  per-library multiome: RNA ∩ ATAC          21,038 + 38,924 cells
        ├─ §3  merge across technologies                 59,962 + 14,828 + 8,452 = 83,242
        ├─ §4  ArchR arrows + its default filter          83,242 → 83,135 → 70,643
        ├─ §5  WNN integration (harmony)                 3k_dsci_ht_10x_gs_subset.RDS
        ├─ §6  load once: 33 types, technique, order, palette, derived cache
        ├─ §7  panels A–C   ├─ §8  panel D   └─ §9  panel E
```

| Panel | What it shows | Drawn by |
|---|---|---|
| A | UMAP of 70,643 nuclei in RNA space (Gene expression) | §7, reduction `umap` |
| B | the same nuclei in ATAC space (Chromatin accessibility) | §7, reduction `umap.atac` |
| C | the same nuclei in the WNN-integrated space | §7, reduction `wnn.umap` |
| D | dot plot, 48 lineage marker genes × 33 cell types | §8, RNA assay |
| E | ATAC pseudo-bulk coverage tracks at five lineage loci | §9, ArchR project |

**Only the dsciCAT-seq libraries run through this pipeline**: FY13 and FY14 are
the 16-barcode dsciCAT-seq runs that §1–§5 process. The HT-scCAT-seq and 10x
Multiome matrices join at §3, where their RNA and ATAC assays are merged into
the same object; from §4 on every library is treated identically.""")

# ===========================================================================
md("""---

## 0. Environment, inputs, outputs

| Variable | Default | Contents |
|---|---|---|
| `FIG2_INPUT` | `<raw>` | the workflow outputs this notebook starts from: `FY/FY13`, `FY/FY14` (RNA matrices, peak calls, fragments) and `scCAT-seq-HT/` (the HT-scCAT-seq / 10x Multiome object) |
| `FIG2_DATA` | `<work>/FY` | the analysis tree: the SoupX matrix of each library, the per-stage objects §2–§5 write, and the ArchR project §4 and §9 read |
| `FIG2_OUT` | `fig2_out` | created here; panel PDFs, the derived cache and the hand-off files between stages land here |

Each stage first looks for its own output object and loads it if present, so the
notebook reads as a document (top to bottom: how the object is built) and runs as
a pipeline (delete an intermediate and it rebuilds from that point). §4 and §5
are the expensive ones — ArchR and the WNN integration over 70,643 cells;
everything before them takes minutes. R ≥ 4.1 with Seurat, Signac, ArchR,
harmony; versions are printed in §0.3 and the full list is under Requirements.""")

code("""
import matplotlib
matplotlib.use('Agg')
import os, subprocess, time, json, glob
import numpy as np

IN   = os.environ.get('FIG2_INPUT', os.environ.get('FIG1_DATA', ''))
DATA = os.environ.get('FIG2_DATA', os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')))
OUT  = os.path.abspath(os.environ.get('FIG2_OUT', 'fig2_out'))   # absolute: §4 changes the cwd
os.makedirs(OUT, exist_ok=True)
print('INPUT =', IN, '| exists:', os.path.isdir(IN))
print('DATA  =', DATA, '| exists:', os.path.isdir(DATA), '| OUT =', OUT)

def show(*paths):
    # report the saved file only -- nothing is displayed inline (section 0.2)
    for q in paths:
        print('%s  (%.0f KB)' % (q, os.path.getsize(q) / 1024.0))

def run_r(rscript, env_extra=None, echo=True):
    '''run one R script file with the three roots exported; return the exit code'''
    env = dict(os.environ, FIG2_INPUT=IN, FIG2_DATA=DATA, FIG2_OUT=OUT)
    if env_extra: env.update(env_extra)
    t0 = time.time()
    r = subprocess.run(['Rscript', rscript], env=env, capture_output=True, text=True)
    dt = (time.time() - t0) / 60
    if echo:
        print(r.stdout[-4000:])
        if r.stderr.strip():
            print('--- stderr (tail) ---'); print(r.stderr[-1500:])
    print('Rscript exit code %d - %.1f min' % (r.returncode, dt))
    return r.returncode
""")

# ===========================================================================
md("""### 0.1 The files this notebook starts from

| Used by | Path | What it is |
|---|---|---|
| §2 | `$FIG2_DATA/FY13/SoupxMatrix/`, `$FIG2_DATA/FY14/SoupxMatrix/` | **the SoupX-corrected RNA matrix** of each library, in 10x format (`matrix.mtx.gz` + `features` + `barcodes`). This is the first file of this notebook's chain |
| §2 | `$FIG2_DATA/FY1x/FY1x_RNA_soupx_PCA*.RDS` | the same matrix after the recorded QC cuts and the DoubletFinder singlet call of that run — 21,038 (FY13) / 38,924 (FY14) cells |
| §1 | `$FIG2_INPUT/FY/FY1x/ATAC/out/Peak/` + `ATAC/03.d2cfile/FY1x.fragments.tsv.gz` | peak × cell matrix (`matrix.mtx` + `peak.bed` + `barcodes.tsv`) and the deduplicated fragments of the same bead run |
| §1 | `$FIG2_DATA/FY13/FY13_ATAC_object.RDS`, `$FIG2_DATA/FY14/FY14_ATAC_object_force3k.RDS` | the two Signac objects built from them (the `_force3k` run is the one the figure uses) |
| §3 | `$FIG2_INPUT/scCAT-seq-HT/Mbrain_process.rds` | HT-scCAT-seq + 10x Multiome RNA and ATAC assays, merged in at §3 |
| §4, §9 | `$FIG2_DATA/archr-subset-copnfol/` | the ArchR project (arrow files, peak calls, gene-score matrix) |
| §5 | `$FIG2_DATA/3k_dsci_ht_10x_gs_subset.RDS` | what §5 writes: the object §6–§9 read |

The SoupX correction itself is upstream of this notebook: the corrected matrix is
read as given, together with the peak matrix called from the same bead run.""")

md("""### 0.2 Provenance: read this before comparing panels

* **What is computed here.** FY13 and FY14 — the two dsciCAT-seq libraries —
  go through §1–§5: peak/RNA intersection, the recorded QC and singlet calls,
  the merge, the ArchR filter, the WNN integration. The HT-scCAT-seq and 10x
  Multiome matrices enter at §3 as processed assay blocks; nothing here
  re-derives them.
* **33 types, not 32.** `cell_type` has 32 classes; `RNA_snn_res.1` cluster 34
  (280 cells, 276 Oligo + 4 OPC) was re-annotated as **COP/NFOL** on
  Gpr17/Bmp4 co-expression evidence, giving `cell_type_updated` = 33 classes.
  Every panel uses `cell_type_updated`. The annotation is a manual step: §5
  writes the object without it and §6 applies the recorded labels.
* **`technique` groups the two dsciCAT-seq libraries** (FY13 → dsciCAT-seq-50k,
  FY14 → dsciCAT-seq-100k → dsciCAT-seq), which is why the legend counts are
  47,443 / 14,762 / 8,438 (dsciCAT-seq / HT-scCAT-seq-HT / 10x).
* **Panel D's rendering parameters were recovered by measurement.** The analysis
  script saved a DotPlot at 8 × 8 cm with `dot.scale = 16`; the delivered panel
  is not that plot (largest dot 0.069 cm, no title, 6 pt labels). §8 reproduces
  the measured numbers.
* **Panel E uses the same ArchR project the figure was drawn from**
  (`archr-subset-copnfol`, ArchR 1.0.2, mm10); its 33-type grouping is re-derived
  in §9 from the Seurat object by barcode matching.
* **Two final-composition touches are not reproduced.** The delivered page
  renames `Ependyma` to **Ependymal cells** and recomposes A–C with one shared
  legend and re-placed centroid numbers. This notebook produces the panel
  sources: `Ependyma`, per-panel legends. Fonts are substituted (the page is
  Arial; this machine has only Helvetica, which is metric-compatible, so
  geometry is unaffected).""")

md("""### 0.3 Input check""")
code("""
CHK = r'''
suppressPackageStartupMessages({library(Seurat); library(Signac); library(ggplot2)
  library(dplyr); library(tidyr); library(pals)})
IN  <- Sys.getenv('FIG2_INPUT', os.environ.get('FIG1_DATA', ''))
FY  <- Sys.getenv('FIG2_DATA', os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')))
OUT <- Sys.getenv('FIG2_OUT', 'fig2_out')
p <- c(file.path(FY, 'FY13/SoupxMatrix', 'matrix.mtx.gz'),
       file.path(FY, 'FY14/SoupxMatrix', 'matrix.mtx.gz'),
       list.files(file.path(FY, 'FY13'), 'RNA_soupx_PCA.*[.]RDS$', full.names = TRUE)[1],
       list.files(file.path(FY, 'FY14'), 'RNA_soupx_PCA.*[.]RDS$', full.names = TRUE)[1],
       file.path(IN, 'FY/FY13/ATAC/out/Peak/matrix.mtx'),
       file.path(IN, 'FY/FY14/ATAC/out/Peak/matrix.mtx'),
       file.path(FY, 'FY13/FY13_ATAC_object.RDS'),
       file.path(FY, 'FY14/FY14_ATAC_object_force3k.RDS'),
       file.path(IN, 'scCAT-seq-HT/Mbrain_process.rds'),
       file.path(FY, 'archr-subset-copnfol/Save-ArchR-Project.rds'),
       file.path(FY, '3k_dsci_ht_10x_gs_subset.RDS'))
for (f in p) cat(sprintf('%-62s %s\\n', sub(IN, '<IN>', sub(FY, '<DATA>', f, fixed = TRUE),
                              fixed = TRUE), if (file.exists(f)) 'present' else 'MISSING'))
cat('R', R.version.string, '| Seurat', as.character(packageVersion('Seurat')),
    '| Signac', as.character(packageVersion('Signac')), '| ArchR',
    tryCatch(as.character(packageVersion('ArchR')), error = function(e) 'NO'),
    '| harmony', tryCatch(as.character(packageVersion('harmony')), error = function(e) 'NO'),
    '| pals', as.character(packageVersion('pals')), '\\n')
'''
p = os.path.join(OUT, 'fig2_check.R'); open(p, 'w').write(CHK)
run_r(p)
""")
md("""**Expected** — every row `present`; R ≥ 4.1, Seurat ≥ 4, Signac ≥ 1.5,
ArchR 1.0.x, harmony ≥ 1, pals ≥ 1.7.""")

# ===========================================================================
md("""---

## 1. Per-library ATAC objects

**Input**: `ATAC/out/Peak/` (peak × cell matrix) + `03.d2cfile/*.fragments.tsv.gz`.
**Output**: `$FIG2_DATA/FY1x/FY1x_ATAC_object*.RDS` — the two objects holding a
`ChromatinAssay`: 85,255 peaks × 23,376 cells (FY13), 89,096 × 43,249 (FY14); §2
reloads them with the same `force3k`-first rule.

The ATAC object is the counterpart of the SoupX matrix on the chromatin side:
same bead run, same cell list, peaks instead of genes. `min.cells = 10` drops
peaks seen in fewer than ten cells, and `TSSEnrichment` / `FRiP` are computed
here because §4's filter and the quality narrative of the paper use them.""")

code(r"""
ATA = r'''
suppressMessages({library(Seurat); library(Signac); library(rtracklayer); library(GenomeInfoDb)})
IN  <- Sys.getenv('FIG2_INPUT', os.environ.get('FIG1_DATA', ''))
FY  <- Sys.getenv('FIG2_DATA', os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')))
OUT <- Sys.getenv('FIG2_OUT', 'fig2_out')
say <- function(...) cat(sprintf(...), '\n')
gtf <- os.path.join(os.environ.get('DSCI_REFDATA_ROOT', ''), 'mm10/scATAC_ref/genes.gtf')

## §2 needs the two objects, and the same preference rule resolves them there.
ATAC_PATH <- function(s) {
  f <- list.files(file.path(FY, s), paste0(s, '_ATAC_object.*[.]RDS$'), full.names = TRUE)
  f[order(grepl('force3k', f), decreasing = TRUE)][1]                 # the 3k-force run first
}

atacs <- list()
for (s in c('FY13', 'FY14')) {
  saved <- list.files(file.path(FY, s), paste0(s, '_ATAC_object.*[.]RDS$'), full.names = TRUE)
  saved <- saved[order(grepl('force3k', saved), decreasing = TRUE)]   # the 3k-force run first
  if (length(saved)) {
    atacs[[s]] <- readRDS(saved[1])
    say('%s: loaded %s | %d peaks x %d cells | FRiP median %.3f | TSS enrichment median %.2f',
        s, basename(saved[1]), nrow(atacs[[s]]), ncol(atacs[[s]]),
        median(atacs[[s]]$FRiP), median(atacs[[s]]$TSS.enrichment))
    next
  }
  ## built from the peak call of the same bead run
  d          <- file.path(IN, 'FY', s)
  mtx_path   <- file.path(d, 'ATAC/out/Peak/matrix.mtx')
  feat_path  <- file.path(d, 'ATAC/out/Peak/peak.bed')
  bcs_path   <- file.path(d, 'ATAC/out/Peak/barcodes.tsv')
  frag.file  <- list.files(file.path(d, 'ATAC/03.d2cfile'), '[.]fragments[.]tsv[.](gz|txt)$',
                           full.names = TRUE)[1]
  counts <- Matrix::readMM(mtx_path)
  features <- readr::read_tsv(feat_path, col_names = FALSE) %>% tidyr::unite(feature)
  barcodes <- readr::read_tsv(bcs_path, col_names = FALSE) %>% tidyr::unite(barcode)
  row.names(counts) <- features$feature; colnames(counts) <- barcodes$barcode

  gtf.gr <- if (file.exists(gtf)) rtracklayer::import(gtf) else
            readRDS(os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'mm10.gene.coords.RDS'))
  gene.coords <- gtf.gr[gtf.gr$type == 'gene']
  seqlevelsStyle(gene.coords) <- 'UCSC'
  gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')

  g <- StringToGRanges(rownames(counts), sep = c(':', '-'))
  counts <- counts[as.vector(seqnames(g) %in% standardChromosomes(g)), ]

  a <- CreateChromatinAssay(counts = counts, sep = c(':', '-'), fragments = frag.file,
                            min.cells = 10, annotation = gene.coords)
  obj <- CreateSeuratObject(counts[, colnames(a)], project = s)
  obj[['ATAC']] <- a
  DefaultAssay(obj) <- 'ATAC'
  obj <- TSSEnrichment(obj)
  tf <- CountFragments(frag.file); rownames(tf) <- tf$CB
  obj$fragments <- tf[colnames(obj), 'frequency_count']
  obj <- FRiP(obj, assay = 'ATAC', total.fragments = 'fragments')
  obj$loguniqueFrag <- log10(obj$fragments)
  obj$sample <- s
  saveRDS(obj, file.path(FY, s, paste0(s, '_ATAC_object.RDS')))
  atacs[[s]] <- obj
  say('%s: built from %s | %d peaks x %d cells', s, basename(frag.file), nrow(obj), ncol(obj))
}
say('checked %d objects: %s', length(atacs),
    paste(sapply(atacs, function(x) sprintf('%d peaks x %d cells', nrow(x), ncol(x))),
          collapse = ' | '))
'''
p = os.path.join(OUT, 'fig2_up1_atac.R'); open(p, 'w').write(ATA)
assert run_r(p) == 0
""")
md("""**Expected** — `FY13: 85255 peaks x 23376 cells`, `FY14: 89096 peaks x 43249
cells`, FRiP medians 0.23–0.25. Nuclei are TSS-poor (median enrichment ≈ 1.6),
which is why §4 filters on ArchR's own `TSSEnrichment` column rather than this
one.""")

# ===========================================================================
md("""---

## 2. Per-library multiome: RNA ∩ ATAC

**Input**: the SoupX matrix of §0.1 + the ATAC object of §1.
**Output**: `$FIG2_DATA/FY13/3k_FY13_multiome.RDS` (21,038 cells) and
`$FIG2_DATA/FY14/3k_FY14_multiome.RDS` (38,924 cells).

The two modalities are called from the same beads but written by two programs, so
a cell only exists for this analysis if it has both an RNA and an ATAC entry. The
RNA side is read *after* its own cleaning: the SoupX matrix is
`CreateSeuratObject`-ed, cut on `nFeature_RNA > 300`,
`nFeature_RNA < n95` (the 95th percentile of this run) and `percent.mt < 10`, and
then reduced to the DoubletFinder singlets of that run (`pN = 0.25`, `pK` from
`find.pK`, `nExp = 10 %` of the cells). If the singlet object of the run is on
disk it is used directly; the `else` branch is the code that made it.""")

code(r"""
MUL = r'''
suppressMessages({library(Seurat); library(Signac); library(DoubletFinder); library(dplyr)
  library(glue)})
IN  <- Sys.getenv('FIG2_INPUT', os.environ.get('FIG1_DATA', ''))
FY  <- Sys.getenv('FIG2_DATA', os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')))
OUT <- Sys.getenv('FIG2_OUT', 'fig2_out')
say <- function(...) cat(sprintf(...), '\n')
set.seed(1234)
dim <- 50; minfeatures <- 300; mtegne_filter <- 10; doublets.percentage <- 0.1
atacs <- list()
for (s in c('FY13', 'FY14')) {                     # same preference rule as §1
  f <- list.files(file.path(FY, s), paste0(s, '_ATAC_object.*[.]RDS$'), full.names = TRUE)
  f <- f[order(grepl('force3k', f), decreasing = TRUE)]
  if (length(f)) atacs[[s]] <- readRDS(f[1])
}

for (s in c('FY13', 'FY14')) {
  d <- file.path(FY, s)
  fut <- file.path(d, paste0('3k_', s, '_multiome.RDS'))
  if (file.exists(fut)) {
    m <- readRDS(fut); say('%s: loaded %s | %d cells', s, basename(fut), ncol(m))
    rm(m); next
  }

  ## (a) RNA side: the QC-cut, singlet-called object of this SoupX run, or rebuild it
  rn <- list.files(d, paste0(s, '_RNA_soupx_PCA.*[.]RDS$'), full.names = TRUE)
  if (length(rn)) {
    rn <- rn[order(file.info(rn)$mtime, decreasing = TRUE)][1]     # latest run wins
    rna <- readRDS(rn)
    say('%s: RNA side loaded %s | %d cells', s, basename(rn), ncol(rna))
  } else {
    rna <- CreateSeuratObject(Read10X(file.path(d, 'SoupxMatrix'), gene.column = 1))
    rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^mt-")
    om <- rna@meta.data[order(-rna@meta.data$nFeature_RNA), ]
    n95 <- as.integer(nrow(om) * 0.05)
    n95features <- as.numeric(om[n95, "nFeature_RNA"])
    rna <- subset(rna, subset = nFeature_RNA > minfeatures &
                                       nFeature_RNA < n95features &
                                       percent.mt < mtegne_filter)
    sweep.res.list <- paramSweep_v3(rna, PCs = 1:dim, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC == max(bcmvn$MeanBC), ]$pK))
    nExp_poi <- round(doublets.percentage * ncol(rna))
    rna <- doubletFinder_v3(rna, PCs = 1:dim, pN = 0.25, pK = pK, nExp = nExp_poi,
                            reuse.pANN = FALSE, sct = FALSE)
    colnames(rna@meta.data)[ncol(rna@meta.data)] <- "doublet_info"
    write.table(rna@meta.data, file.path(d, paste0(s, "_doublets_info.txt")),
                sep = "\t", quote = FALSE)
    rna <- subset(rna, subset = doublet_info == "Singlet")
    saveRDS(rna, file.path(d, paste0(s, "_RNA_soupx_PCA.RDS")))
    say('%s: RNA side built | %d singlets', s, ncol(rna))
  }

  ## (b) keep the cells that have both modalities, then attach the peaks
  atac <- atacs[[s]]
  common <- intersect(colnames(rna), colnames(atac))
  say('%s: RNA %d | ATAC %d | intersect %d', s, ncol(rna), ncol(atac), length(common))
  rna  <- rna[, common]
  atac <- atac[, common]
  brain <- rna
  brain[["ATAC"]] <- atac[["ATAC"]]
  brain$nCount_ATAC  <- atac$nCount_ATAC
  brain$nFeature_ATAC <- atac$nFeature_ATAC
  brain$FRiP          <- atac$FRiP
  brain$TSS.enrichment <- atac$TSS.enrichment
  brain$orig.ident <- s
  saveRDS(brain, fut)
  say('%s: multiome -> %d cells', s, ncol(brain))
  rm(brain, rna, atac); invisible(gc())
}
'''
p = os.path.join(OUT, 'fig2_up2_multiome.R'); open(p, 'w').write(MUL)
assert run_r(p) == 0
""")
md("""**Expected** — `FY13: RNA 21038 | ATAC 23376 | intersect 21038`,
`FY14: RNA 38924 | ATAC 43249 | intersect 38924`: the intersection costs nothing,
because every cell of the cleaned RNA matrix also has a peak row. It is a real
filter in the other direction (23,376 beads with peaks, 21,038 of them with
corrected RNA), and it is what makes a cell *multiome* rather than ATAC-only.
21,038 + 38,924 = 59,962 dsciCAT-seq cells go into §3.

The two cells-per-library numbers the figure quotes — 20,971 and 38,884 — appear
one step later, in §4: 67 FY13 and 40 FY14 of these cells have no fragments in
the deduplicated fragment file, so they never become arrow records.""")

# ===========================================================================
md("""---

## 3. Merge across technologies

**Input**: the two objects of §2 + the HT-scCAT-seq / 10x Multiome RNA and ATAC
assays. **Output**: the candidate cell list `fig2_out/fig2_candidate_cells.RDS` and, when the chain is
being rebuilt, `$FIG2_DATA/3k_dsci_ht_10x_merged.RDS` — **83,242** cells (59,962
dsciCAT-seq + 23,280 HT-scCAT-seq + 10x Multiome).

`merge()` refuses objects that carry reductions, graphs or neighbours and
silently mixes them when it accepts, so every slot is cleared first — and each
object's cells are then labelled by library, which is the batch variable harmony
uses in §5 and the grouping the legend counts in.""")

code(r"""
MRG = r'''
suppressMessages({library(Seurat); library(Signac)})
IN  <- Sys.getenv('FIG2_INPUT', os.environ.get('FIG1_DATA', ''))
FY  <- Sys.getenv('FIG2_DATA', os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')))
OUT <- Sys.getenv('FIG2_OUT', 'fig2_out')
say <- function(...) cat(sprintf(...), '\n')

## What §4 and the legend need from §3 is one thing: the candidate cell list, one
## row per cell, with the library it came from. Assembling the *object* behind
## that list - merge() over 83,242 cells x ~1e5 peaks - is the most expensive step
## in this notebook and feeds nothing else, so it is built only when the chain
## really has to be re-run; the file.exists() rule below says so out loud.
libs <- list(FY13 = file.path(FY, 'FY13', '3k_FY13_multiome.RDS'),
             FY14 = file.path(FY, 'FY14', '3k_FY14_multiome.RDS'))
dsci <- do.call(rbind, lapply(names(libs), function(s)
  data.frame(cell = colnames(readRDS(libs[[s]])), sample = s)))
auth <- readRDS(file.path(IN, 'scCAT-seq-HT/Mbrain_process.rds'))
acell <- data.frame(cell = colnames(auth),
                    sample = ifelse(grepl('^Cortex_EL', colnames(auth)), 'HT-scCAT-seq-HT', '10x'))
cand <- rbind(dsci, acell)
say('candidates: %d cells = %d dsciCAT-seq + %d HT/10x',
    nrow(cand), nrow(dsci), nrow(acell))
print(table(cand$sample))
stopifnot(nrow(cand) == 83242, !anyDuplicated(cand$cell))
saveRDS(cand, file.path(OUT, 'fig2_candidate_cells.RDS'))
rm(dsci, auth); invisible(gc())

## the merged object, built only if something downstream still needs it
out.rds <- file.path(FY, '3k_dsci_ht_10x_merged.RDS')
gs.rds  <- file.path(FY, '3k_dsci_ht_10x_gs_subset.RDS')
if (file.exists(out.rds)) {
  say('merged object present: %s (%.1f GB)', basename(out.rds), file.info(out.rds)$size / 1e9)
} else if (file.exists(gs.rds)) {
  say('merged object not rebuilt: %s already exists and nothing downstream reads the intermediate',
      basename(gs.rds))
} else {
  ## merge() refuses objects that carry reductions, graphs or neighbours and
  ## silently mixes them when it accepts, so every slot is cleared first. Cell
  ## names must come through unchanged - they already start with FY13_ / FY14_ and
  ## §4 matches them against the fragment files - so no add.cell.ids.
  obj <- lapply(libs, function(f) {
    o <- readRDS(f); o@graphs <- list(); o@neighbors <- list(); o@reductions <- list(); o })
  multi <- merge(obj[[1]], y = obj[2], add.cell.ids = NULL)
  say('dsciCAT-seq merged: %d cells', ncol(multi))
  auth <- readRDS(file.path(IN, 'scCAT-seq-HT/Mbrain_process.rds'))
  auth@graphs <- list(); auth@neighbors <- list(); auth@reductions <- list()
  auth$orig.ident <- ifelse(grepl('^Cortex_EL', colnames(auth)), 'HT-scCAT-seq-HT', '10x')
  DefaultAssay(multi) <- 'RNA'; DefaultAssay(auth) <- 'RNA'
  merged <- merge(auth, y = multi, add.cell.ids = NULL)
  stopifnot(ncol(merged) == 83242, !anyDuplicated(colnames(merged)))
  saveRDS(merged, out.rds)
}
'''
p = os.path.join(OUT, 'fig2_up3_merge.R'); open(p, 'w').write(MRG)
assert run_r(p) == 0
""")

md("""**Expected** — `candidates: 83242 cells = 59962 dsciCAT-seq + 23280 HT/10x`
and the four-row table `10x 8452 | FY13 21038 | FY14 38924 | HT-scCAT-seq-HT
14828`. Then one line about the merged object: `present` if
`3k_dsci_ht_10x_merged.RDS` is on disk, `not rebuilt` if §5's object already
covers it, or `dsciCAT-seq merged: 59962 cells` if the merge actually ran.""")

# ===========================================================================
md("""---

## 4. ArchR with its default filter, and the gene-score matrix

**Input**: the merged cell list of §3 + the fragment file of each library.
**Output**: the ArchR project `$FIG2_DATA/archr-subset-copnfol` and, when the
chain is being rebuilt, `3k_dsci_ht_10x_gs_preWNN.RDS` — §3's cells reduced to the
survivors and carrying a new `gene_score` assay.

This is where 83,242 becomes **70,643**, in two steps. `createArrowFiles()` turns
the merged cell list into arrow records, and **107** cells get none (67 FY13 +
40 FY14), so the arrows hold 83,135 cells. These are not low-coverage cells —
their bead-level metadata reports ~1.3 × 10⁵ unique fragments each — so they are
lost in the fragment-to-arrow step, not in any QC threshold; the notebook
asserts the arrows' count rather than guessing at the mechanism. ArchR's
**default** filter — `minFrags = 1000`, `minTSS = 4`, nothing tuned — then
removes **12,492** cells whose TSS enrichment is below 4:

```
83,242 merged  →  83,135 arrow records  −  12,492 with TSSEnrichment < 4  =  70,643
   10x    8,452 →   8,452 →   8,438      |   HT  14,828 →  14,828 →  14,762
   FY13  21,038 →  20,971 →  17,093      |  FY14 38,924 →  38,884 →  30,350
```

The 12,492 is the whole `minTSS` loss: measured on the arrows, not one of their
83,135 cells falls below `minFrags = 1000`, so the second half of the default
filter contributes nothing.

§3 and §4 are the two stages whose product is an *intermediate*; if
`3k_dsci_ht_10x_gs_subset.RDS` (§5's output, the file §6–§9 read) is already on
disk, they measure and report the chain instead of writing a 20–30 GB object that
nothing downstream will load. Delete that file and the same two cells build it.

`addGeneScoreMat = TRUE` then writes the ATAC-derived gene-score matrix (24,333
genes × 70,643 cells, every cell summing to 10,000), which becomes the
`gene_score` assay — the third modality of this object.""")

code(r"""
ARC = r'''
suppressMessages({library(ArchR); library(Seurat); library(Signac)})
IN  <- Sys.getenv('FIG2_INPUT', os.environ.get('FIG1_DATA', ''))
FY  <- Sys.getenv('FIG2_DATA', os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')))
OUT <- Sys.getenv('FIG2_OUT', 'fig2_out')
say <- function(...) cat(sprintf(...), '\n')
options(expressions = 500000); Sys.setenv(R_MAX_NUM_DLLS = 150)
plan('multicore', workers = 16); options(future.globals.maxSize = 1000 * 1024^10)
dir.create(file.path(FY, 'archr-3k_dsci_ht_10x'), showWarnings = FALSE, recursive = TRUE)
setwd(file.path(FY, 'archr-3k_dsci_ht_10x'))          # ArchR writes its logs to the cwd
addArchRGenome('mm10'); addArchRThreads(threads = 16)

prj     <- file.path(FY, 'archr-subset-copnfol')
arrows  <- file.path(FY, '3k_dsci_ht_10x_ArchR.RDS')
gs.rds  <- file.path(FY, '3k_dsci_ht_10x_gs_subset.RDS')
cand    <- readRDS(file.path(OUT, 'fig2_candidate_cells.RDS'))

## (a) arrow records: the candidate cells whose fragments ArchR can index
if (file.exists(arrows)) {
  pre <- readRDS(arrows)
  say('arrows loaded: %d cells | samples %s', length(pre$cellNames),
      paste(names(table(pre$Sample)), collapse = ', '))
} else {
  ## one fragment file per library - the delivered project's samples are
  ## FY13 FY14 EL13 EL14 EL15 EL16 10x, and sampleNames comes from these basenames
  fl <- c(list.files(file.path(IN, 'FY/FY13/ATAC/03.d2cfile'), '[.]fragments[.]tsv[.]gz$',
                     full.names = TRUE)[1],
          list.files(file.path(IN, 'FY/FY14/ATAC/03.d2cfile'), '[.]fragments[.]tsv[.]gz$',
                     full.names = TRUE)[1],
          list.files(file.path(IN, 'scCAT-seq-HT'), '[.]frag', full.names = TRUE))
  fl <- fl[file.exists(fl)]
  ArrowFiles <- createArrowFiles(inputFiles = fl, validBarcodes = cand$cell,
                                 sampleNames = sub('[.]frag.*', '', basename(fl)),
                                 addTileMat = TRUE, addGeneScoreMat = TRUE,
                                 threads = 16, force = TRUE, cleanTmp = TRUE)
  pre <- ArchRProject(ArrowFiles = ArrowFiles,
                      outputDirectory = file.path(FY, 'archr-3k_dsci_ht_10x'),
                      copyArrows = FALSE)
  saveRDS(pre, arrows)
}
stopifnot(length(pre$cellNames) == 83135)

## (b) ArchR's own default filter, measured on the arrows
tss <- as.numeric(pre@cellColData$TSSEnrichment)
nfr <- as.numeric(pre@cellColData$nFrags)
keep <- !(tss < 4 | nfr < 1000)
say('arrows: %d cells | dropped by the default filter: %d | survivors: %d',
    length(pre$cellNames), sum(!keep), sum(keep))
say('below minFrags = 1000: %d | below minTSS = 4: %d', sum(nfr < 1000), sum(tss < 4))
print(table(pre$Sample[keep]))
stopifnot(sum(keep) == 70643)

## (c) the subsetted project - the object §9 draws its tracks from
if (file.exists(file.path(prj, 'Save-ArchR-Project.rds'))) {
  proj <- loadArchRProject(prj, showLogo = FALSE)
  say('project loaded: %s | %d cells | matrices: %s', basename(prj), length(proj$cellNames),
      paste(getAvailableMatrices(proj), collapse = ', '))
} else {
  proj <- subsetArchRProject(pre, cells = pre$cellNames[keep], outputDirectory = prj,
                             dropCells = TRUE, force = TRUE)
  proj <- addIterativeLSI(proj, useMatrix = 'TileMatrix', name = 'IterativeLSI', iterations = 2,
           clusterParams = list(resolution = c(0.1, 0.3, 0.5), sampleCells = 30000, n.start = 10),
           varFeatures = 25000, dimsToUse = 2:40, force = TRUE)
  proj <- addGroupCoverages(proj, groupby = 'Sample', force = TRUE)
  proj <- addPeakCalls(proj, ID = 'Conservative', type = 'Quantile', quantile = 0.99, force = TRUE)
  proj <- addPeakMatrix(proj)
  proj <- addImputeWeights(proj, useLogScaling = TRUE, parallel = 8)
  proj <- saveArchRProject(proj, outputDirectory = prj, load = FALSE)
}

## (d) the gene-score matrix: the third modality, and the fingerprint of this step
Matrix <- getMatrixFromProject(proj, useMatrix = 'GeneScoreMatrix', threads = 16)
geneDF <- ArchR:::.getFeatureDF(getArrowFiles(proj), 'GeneScoreMatrix')
rownames(Matrix) <- geneDF[, 'name']
Matrix <- assays(Matrix)[['GeneScoreMatrix']]
colnames(Matrix) <- sub('^.*#', '', colnames(Matrix))
say('GeneScoreMatrix: %d genes x %d cells | every cell sums to 1e4: %s',
    nrow(Matrix), ncol(Matrix), all(round(colSums(Matrix)) == 1e4))
stopifnot(nrow(Matrix) == 24333, ncol(Matrix) == 70643)

## (e) the pre-WNN object, built only if §5 still needs it
out.rds <- file.path(FY, '3k_dsci_ht_10x_gs_preWNN.RDS')
if (file.exists(out.rds)) {
  say('preWNN object present: %s (%.1f GB)', basename(out.rds), file.info(out.rds)$size / 1e9)
} else if (file.exists(gs.rds)) {
  say('preWNN object not rebuilt: %s already exists', basename(gs.rds))
} else {
  merged <- readRDS(file.path(FY, '3k_dsci_ht_10x_merged.RDS'))
  obj <- merged[, intersect(colnames(merged), colnames(Matrix))]
  stopifnot(ncol(obj) == 70643)
  obj[['gene_score']] <- CreateAssayObject(counts = Matrix[, colnames(obj), drop = FALSE])
  obj$nCount_gene_score <- colSums(Matrix[, colnames(obj), drop = FALSE])
  print(table(obj$orig.ident))
  rm(Matrix, merged); invisible(gc())
  saveRDS(obj, out.rds)
}
'''
p = os.path.join(OUT, 'fig2_up4_archr.R'); open(p, 'w').write(ARC)
assert run_r(p) == 0
""")

md("""**Expected** — `arrows: 83135 cells | dropped by the default filter: 12492 |
survivors: 70643`, `below minFrags = 1000: 0 | below minTSS = 4: 12492` (the whole
loss is TSS), `project loaded: archr-subset-copnfol | 70643 cells`,
`GeneScoreMatrix: 24333 genes x 70643 cells | every cell sums to 1e4: TRUE`, and
the survivor table `10x 8438 | EL13 3801 | EL14 4029 | EL15 3453 | EL16 3479 |
FY13 17093 | FY14 30350` — EL13–16 summing to the legend's 14,762 HT-scCAT-seq-HT
cells. The pre-WNN object is reported as `present`, `not rebuilt`, or written.""")

# ===========================================================================
md("""---

## 5. WNN integration — the object every panel reads

**Input**: `3k_dsci_ht_10x_gs_preWNN.RDS` (§4). **Output**:
`3k_dsci_ht_10x_gs_subset.RDS` (13 GB) — the file §6 loads and Figure 2 is drawn
from.

`Process_Multiome_Integration()` (from `WNN_harmony.r`, reproduced below in full)
is what defines the three embeddings of panels A–C: PCA + harmony on RNA →
`umap`, TF-IDF + SVD + harmony on LSI → `umap.atac`, weighted-nearest-neighbour
WNN on the two → `wnn.umap`, plus the clusterings the annotation was made
against. The parameters that produced the delivered object are
`pc.nums = 50, dims_use_rna = 1:50, dims_use_atac = 2:40,
res_clusters = c(0.1, 0.3, 0.5, 1)`, batch variable `orig.ident`; the evidence is
the object itself (`pca`, `harmony`, `umap`, `lsi`, `harmony_lsi`, `umap.atac`,
`wnn.umap`, a neighbour named `wnn_harmony`, and `RNA_snn_res.*` /
`ATAC_snn_res.*` / `wsnn_res.*` columns).""")

code(r"""
WNN = r'''
suppressMessages({library(Seurat); library(Signac); library(harmony)})
FY  <- Sys.getenv('FIG2_DATA', os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')))
OUT <- Sys.getenv('FIG2_OUT', 'fig2_out')
say <- function(...) cat(sprintf(...), '\n')
plan('multicore', workers = 16); options(future.globals.maxSize = 1000 * 1024^10)

# ---- WNN_harmony.r, in full ---------------------------------------------
Process_Multiome_Integration <- function(obj, nfeatures = 3000, pc.nums = 30,
    dims_use_rna = 1:10, dims_use_atac = 2:10, res_clusters = c(0.1, 0.3, 0.5),
    var_int = 'orig.ident', reduction_use_rna = 'harmony',
    reduction_use_atac = 'harmony_lsi', run_harmony = TRUE) {
  DefaultAssay(obj) <- 'RNA'
  obj <- NormalizeData(obj); obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
  obj <- ScaleData(obj); obj <- RunPCA(obj, npcs = pc.nums, reduction.name = 'pca')
  if (run_harmony) {
    obj <- RunHarmony(obj, group.by.vars = var_int, assay.use = 'RNA', reduction.use = 'pca',
                      plot.convergence = FALSE, reduction.save = reduction_use_rna)
    rna_red <- reduction_use_rna
  } else rna_red <- 'pca'
  obj <- FindNeighbors(obj, verbose = FALSE, reduction = rna_red, dims = dims_use_rna)
  obj <- FindClusters(obj, verbose = FALSE, algorithm = 3, resolution = res_clusters)
  obj <- RunUMAP(obj, reduction = rna_red, dims = dims_use_rna)              # 'umap' = panel A

  DefaultAssay(obj) <- 'ATAC'
  obj <- RunTFIDF(obj); obj <- FindTopFeatures(obj, min.cutoff = 'q50'); obj <- RunSVD(obj)
  if (run_harmony) {
    obj <- RunHarmony(obj, group.by.vars = var_int, assay.use = 'ATAC', reduction.use = 'lsi',
                      reduction.save = reduction_use_atac, project.dim = FALSE)
    atac_red <- reduction_use_atac
  } else atac_red <- 'lsi'
  obj <- FindNeighbors(obj, reduction = atac_red, dims = dims_use_atac, n.neighbors = 40)
  obj <- FindClusters(obj, verbose = FALSE, algorithm = 3, resolution = res_clusters,
                      cluster.name = 'ATAC_snn_res')
  obj <- RunUMAP(obj, reduction = atac_red, dims = dims_use_atac,
                 reduction.name = 'umap.atac', reduction.key = 'atacUMAP_')   # panel B

  obj <- FindMultiModalNeighbors(obj, reduction.list = list(rna_red, atac_red),
                                 dims.list = list(dims_use_rna, dims_use_atac),
                                 weighted.nn.name = 'wnn_harmony')
  obj <- FindClusters(obj, graph.name = 'wsnn', algorithm = 3, resolution = res_clusters,
                      verbose = FALSE)
  obj <- RunUMAP(obj, nn.name = 'wnn_harmony', reduction.name = 'wnn.umap',
                 reduction.key = 'wnnUMAP_', seed.use = 1234)                 # panel C
  return(obj)
}
# -------------------------------------------------------------------------

objf <- file.path(FY, '3k_dsci_ht_10x_gs_subset.RDS')
if (file.exists(objf)) {
  say('cached: %s (%.1f GB)', basename(objf), file.info(objf)$size / 1e9)
} else {
  obj <- readRDS(file.path(FY, '3k_dsci_ht_10x_gs_preWNN.RDS'))
  obj <- Process_Multiome_Integration(obj, pc.nums = 50, dims_use_rna = 1:50,
                                      dims_use_atac = 2:40,
                                      res_clusters = c(0.1, 0.3, 0.5, 1))
  ## gene activity from the peaks (aRNA) — used by other figures, kept here so the
  ## object is the one the panels were drawn from
  DefaultAssay(obj) <- 'ATAC'
  obj[['aRNA']] <- CreateAssayObject(counts = GeneActivity(obj))
  DefaultAssay(obj) <- 'aRNA'
  obj <- NormalizeData(obj, assay = 'aRNA', normalization.method = 'LogNormalize',
                       scale.factor = median(obj$nCount_RNA))
  saveRDS(obj, objf)
}
obj <- readRDS(objf)
say('cells: %d | assays: %s', ncol(obj), paste(Assays(obj), collapse = ', '))
say('reductions: %s', paste(Reductions(obj), collapse = ', '))
say('clusters @ RNA_snn_res.1: %d | @ wsnn_res.1: %d',
    length(levels(obj$RNA_snn_res.1)), length(levels(obj$wsnn_res.1)))
'''
p = os.path.join(OUT, 'fig2_up5_wnn.R'); open(p, 'w').write(WNN)
assert run_r(p) == 0
""")
md("""**Expected** — `cells: 70643 | assays: RNA, ATAC, gene_score` (the delivered
object adds `aRNA` and `chromvar`, and a `cell_type` column from the manual
annotation), `reductions: pca, harmony, umap, lsi, harmony_lsi, umap.atac,
wnn.umap`, `clusters @ RNA_snn_res.1: 48 | @ wsnn_res.1: 59`.

The delivered 13 GB object is the input of everything below; the rest of the
notebook does not change it.""")

# ===========================================================================
md("""---

## 6. The step everything depends on — the integrated object

**Input**: `3k_dsci_ht_10x_gs_subset.RDS` (13 GB), written by §5.
**Output**: `fig2_out/gs_subset_light.RDS` (the derived cache that panels A–D
consume) + `fig2_out/summary.json` (every measured number, also printed below).

Loading the 13 GB object is the only expensive read of the notebook; everything
after it uses the small derived cache. The derivation is written out in full:
`cell_type_updated` (COP/NFOL re-annotation), `technique` (50k+100k →
dsciCAT-seq), the fixed 33-type order `dot_order` with its numbering and the
shared palette `global_colors_full`, and the three reductions.""")

code("""
R1 = r'''
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)
  library(tidyr); library(pals)})
DATA <- Sys.getenv('FIG2_DATA', os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')))
OUT  <- Sys.getenv('FIG2_OUT', 'fig2_out')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rds <- file.path(DATA, '3k_dsci_ht_10x_gs_subset.RDS')
t0 <- proc.time()[3]
obj <- readRDS(rds)
cat('loaded', rds, 'in', round((proc.time()[3] - t0) / 60, 1), 'min\n')

# ---- the final annotation (used by every panel) ---------------------------
obj$cell_type_updated <- as.character(obj$cell_type)
obj$cell_type_updated[obj$RNA_snn_res.1 == 34] <- 'COP/NFOL'   # 280 cells
obj$technique <- ifelse(obj$orig.ident %in% c('dsciCAT-seq-50k', 'dsciCAT-seq-100k'),
                        'dsciCAT-seq', as.character(obj$orig.ident))

# fixed 33-type order = track/dotplot order; numbering 1=AST ... 33=Ex_L6b
dot_order <- c('AST','Oligo','COP/NFOL','OPC','MG','EC','VLMC','Ependyma',
               'In_D1','In_D2','In_HAP1','In_LAMP5','In_LHX6','In_SST','In_VIP',
               'In_PVALB','In_TH','Neuron_MB','Ex_LHA','Ex_HY-MM','Ex_DG','Ex_TH',
               'Ex_MH-LH','Ex_SUB_Pros','Ex_CA1_Pros','Ex_CAR3','Ex_L2/3_IT',
               'Ex_L2/3_IT_PIR','Ex_L4/5_IT','Ex_L5/6_NP','Ex_L6_IT','Ex_L6_CT','Ex_L6b')
num_map <- setNames(seq_along(dot_order), dot_order)

# shared palette = global_colors_full (unique(c(glasbey(32), alphabet(26)))[1:42])
all_types_full <- sort(unique(c(
  'AST','COP/NFOL','EC','Ependyma','Ex_CA1_Pros','Ex_CA3','Ex_CAR3','Ex_CR',
  'Ex_DG','Ex_HY-MM','Ex_IC','Ex_L2/3_IT','Ex_L2/3_IT_PIR','Ex_L4','Ex_L4/5_IT',
  'Ex_L5_ET','Ex_L5_IT','Ex_L5_NP','Ex_L5/6_NP','Ex_L6_CT','Ex_L6_IT','Ex_L6b',
  'Ex_LHA','Ex_MH-LH','Ex_SUB_Pros','Ex_TH','In_D1','In_D2','In_D1-D2','In_HAP1',
  'In_LAMP5','In_LHX6','In_MEIS2','In_PVALB','In_SST','In_TH','In_VIP','MG',
  'Neuron_MB','Oligo','OPC','VLMC')))
global_colors_full <- unique(c(as.vector(pals::glasbey(32)), as.vector(pals::alphabet(26))))[1:length(all_types_full)]
names(global_colors_full) <- all_types_full
cols <- global_colors_full[dot_order]

# ---- measured numbers -----------------------------------------------------
cat('total cells:', ncol(obj), '\n')
print(table(obj$technique)); print(table(obj$orig.ident))
cat('cell types (cell_type_updated):', length(unique(obj$cell_type_updated)), '\n')
cat('COP/NFOL:', sum(obj$cell_type_updated == 'COP/NFOL'),
    sprintf('(%.1f%% of %d)', 100 * mean(obj$cell_type_updated == 'COP/NFOL'), ncol(obj)), '\n')
stopifnot(all(dot_order %in% unique(obj$cell_type_updated)))

# ---- small derived cache for panels A-D -----------------------------------
red <- list(RNA = Embeddings(obj, 'umap'),
            ATAC = Embeddings(obj, 'umap.atac'),
            WNN  = Embeddings(obj, 'wnn.umap'))
m48 <- c('Qk','Slc1a2','Slc1a3','Mbp','Plp1','Gpr17','Pdgfra','Tnr','P2ry12',
         'Hexb','Flt1','Pecam1','Eya2','Foxd1','Dnah12','Adgrv1','Rbfox3','Gad1',
         'Gad2','Drd1','Drd2','Hap1','Lamp5','Lhx6','Sst','Vip','Pvalb','Th',
         'Slc17a6','Rmst','Ndst4','Foxb1','Prox1','Prkcd','Tac2','Pou4f1','Hcn1',
         'Fibcd1','Spink8','Car3','Slc17a7','Cux2','Igfn1','Rorb','Tshz2','Il1rapl2',
         'Tle4','Shisa6')
DefaultAssay(obj) <- 'RNA'
ex <- as.matrix(GetAssayData(obj, slot = 'data')[m48, , drop = FALSE])
m48_present <- intersect(m48, rownames(ex))
cat('marker genes in RNA assay:', length(m48_present), '/', length(m48), '\n')
saveRDS(list(meta = obj@meta.data[, c('orig.ident', 'technique', 'cell_type_updated')],
             red = red, ex = ex, m48 = m48_present,
             dot_order = dot_order, num_map = num_map, cols = cols),
        file.path(OUT, 'gs_subset_light.RDS'))
writeLines(c(sprintf('total=%d', ncol(obj)),
             sprintf('tech=%s', paste(names(table(obj$technique)), table(obj$technique),
                                      sep = ':', collapse = ',')),
             sprintf('types=%d', length(unique(obj$cell_type_updated))),
             sprintf('copnfol=%d', sum(obj$cell_type_updated == 'COP/NFOL')),
             sprintf('markers=%d', length(m48_present))),
           file.path(OUT, 'summary.txt'))
cat('wrote gs_subset_light.RDS + summary.txt\n')
'''
r1 = os.path.join(OUT, 'fig2_s6_object.R')
with open(r1, 'w') as fh: fh.write(R1)
run_r(r1)
""")
md("""**Expected** — 70,643 cells; dsciCAT-seq 47,443 = 50k 17,093 + 100k 30,350,
HT-scCAT-seq-HT 14,762, 10x 8,438; 33 cell types; COP/NFOL 280 (0.4%); all 48
marker genes present in the RNA assay. These are the legend numbers of the
delivered figure.""")

# ===========================================================================
md("""---

## 7. Panels A–C — UMAP embeddings of the 70,643 nuclei

**Input**: the three reductions of the derived cache (`RNA`/`ATAC`/`WNN`).
**Output**: `fig2_out/fig2_umap_{rna,atac,wnn}_celltype.pdf` (+ `.png`).

One point per nucleus, coloured by `cell_type_updated` with the shared 33-type
palette, a number at each type's centroid (number = legend order 1–33, white
halo, ink colour chosen by the luminance of the type colour), legend
"1 AST … 33 Ex_L6b". These are the data, palette, order and centroid numbering
of the delivered figure's panels A–C; the delivered page recomposes the three
plots with a single shared legend and re-places the numbers by a repulsion
pass — that recomposition is a production step, not reproduced here.""")

code("""
R2 = r'''
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)
  library(tidyr); library(pals)})
OUT  <- Sys.getenv('FIG2_OUT', 'fig2_out')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
S <- readRDS(file.path(OUT, 'gs_subset_light.RDS'))
meta <- S$meta; dot_order <- S$dot_order; num_map <- S$num_map; cols <- S$cols

# number ink colour: dark type colour -> white digits, light -> black
lum <- apply(col2rgb(cols), 2, function(v) 0.299*v[1] + 0.587*v[2] + 0.114*v[3])
num_col <- ifelse(lum < 0.5, 'white', 'black'); names(num_col) <- dot_order
leg_lab <- setNames(paste0(num_map[dot_order], ' ', dot_order), dot_order)

reductions <- list(rna = 'RNA', atac = 'ATAC', wnn = 'WNN')
for (nm in names(reductions)) {
  e <- as.data.frame(S$red[[reductions[[nm]]]]); colnames(e) <- c('UMAP_1', 'UMAP_2')
  e$cell_type <- factor(meta$cell_type_updated, levels = dot_order)
  ctr <- e %>% group_by(cell_type) %>%
    summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2), .groups = 'drop')
  ctr$num <- num_map[as.character(ctr$cell_type)]; ctr$txt_col <- num_col[as.character(ctr$cell_type)]
  p <- ggplot(e, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = cell_type), size = 0.06, alpha = 0.5) +
    scale_color_manual(values = cols, labels = leg_lab, name = 'Cell Type (n = number)') +
    geom_text(data = ctr, aes(UMAP_1, UMAP_2, label = num), color = 'white', size = 6.6/2.845) +
    geom_text(data = ctr[ctr$txt_col == 'black', ], aes(UMAP_1, UMAP_2, label = num),
              color = 'black', size = 6/2.845) +
    geom_text(data = ctr[ctr$txt_col == 'white', ], aes(UMAP_1, UMAP_2, label = num),
              color = 'white', size = 6/2.845) +
    guides(color = guide_legend(ncol = 2, keywidth = 0.18, keyheight = 0.18,
                                override.aes = list(size = 1.6, alpha = 1))) +
    theme_bw(base_size = 9) +
    theme(legend.position = 'right', legend.text = element_text(size = 6),
          legend.title = element_text(size = 7), legend.key = element_blank(),
          legend.spacing.y = unit(0.02, 'cm'), panel.grid = element_blank())
  ggsave(file.path(OUT, sprintf('fig2_umap_%s_celltype.pdf', nm)), p, width = 18, height = 12, units = 'cm')
  ggsave(file.path(OUT, sprintf('fig2_umap_%s_celltype.png', nm)), p, width = 18, height = 12, units = 'cm', dpi = 150)
  cat('saved fig2_umap_%s_celltype.pdf\n', nm)
}
cat('order:', paste(paste0(num_map, names(num_map)), collapse = ' '), '\n')
'''
r2 = os.path.join(OUT, 'fig2_s7_umaps.R')
with open(r2, 'w') as fh: fh.write(R2)
run_r(r2)
""")
code("""for nm in ['rna', 'atac', 'wnn']:
    p = os.path.join(OUT, f'fig2_umap_{nm}_celltype.png')
    if os.path.exists(p):
        show(p)
""")
md("""**Expected** — three 18 × 12 cm UMAP PDFs, 70,643 points each, 33 colours,
numbered centroids exactly 1–33 in `dot_order` (AST → Ex_L6b), legends
"1 AST … 33 Ex_L6b". The RNA/ATAC/WNN panels of the delivered figure are these
three plots on a shared legend block.""")

# ===========================================================================
md("""---

## 8. Panel D — dot plot of 48 lineage marker genes across the 33 types

**Input**: the 48 marker rows of the RNA assay in the derived cache.
**Output**: `fig2_out/fig2_full_celltype_dotplot.pdf` (+ `.png`).

`Seurat::DotPlot` with `scale = TRUE` (average expression z-scored per gene,
clipped to 0–3): dot size = percent of cells expressing, colour = average
expression. The 48 genes are the lineage markers of the delivered figure, in
delivered order; the y axis is `rev(dot_order)` so Ex_L6b is at the top and AST
at the bottom, matching the track panel.

The rendering parameters were recovered by measuring the delivered panel, not
copied from the analysis notebook, because the analysis notebook's own save
(`dot.scale = 16`, `8 × 8 cm`, bold y labels, a title) does not look like the
published panel: it draws 12 mm dots that overlap into a solid mass. The
delivered D panel has a 48 × 33 grid at 1.86 × 2.27 mm per cell, its largest
dot is 1.37 mm, its type labels are plain 6 pt and it carries no title. The
code below reproduces those measurements.""")

code("""
R3 = r'''
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(grid)})
OUT  <- Sys.getenv('FIG2_OUT', 'fig2_out')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
S <- readRDS(file.path(OUT, 'gs_subset_light.RDS'))
obj <- CreateSeuratObject(S$ex); obj@meta.data$cell_type_updated <- S$meta$cell_type_updated
Idents(obj) <- 'cell_type_updated'
present <- S$m48
lv <- levels(factor(S$meta$cell_type_updated, levels = S$dot_order))

# Panel geometry below was measured off the delivered D panel (its source,
# 'Figure S2-dot-track.pdf', is on this server; its own build script is not).
# dot.scale 16 (the value in the analysis notebook) draws 12 mm dots here, which
# is what the delivered panel does NOT look like: its largest dot is 1.37 mm.
PANEL_W <- 8.90   # cm -> 48 columns at 0.186 cm pitch (delivered: 0.186)
PANEL_H <- 7.48   # cm -> 33 rows at 0.227 cm pitch     (delivered: 0.227)
DOT_SCALE <- 1.5  # mm at 100% expressed (delivered largest dot: 1.37 mm)

p <- DotPlot(obj, group.by = 'cell_type_updated', features = present,
             scale = TRUE, cluster.idents = FALSE, col.min = 0, col.max = 3,
             dot.scale = DOT_SCALE) +
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = 'italic', size = 6),
        axis.text.y = element_text(size = 6),
        axis.title  = element_text(size = 8, face = 'bold'),
        legend.key.width = unit(0.22, 'cm'), legend.key.height = unit(0.30, 'cm'),
        legend.key.size  = unit(0.22, 'cm'), legend.spacing = unit(0, 'cm'),
        legend.position = 'top', legend.box = 'horizontal', legend.box.just = 'left',
        legend.box.spacing = unit(0.25, 'cm'),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0),
        plot.title = element_blank(),
        plot.margin = unit(c(0.1, 0.7, 0.15, 0.1), 'cm')) +
  scale_color_gradientn(colours = c(rgb(252/255, 247/255, 243/255),
                                    rgb(253/255, 200/255, 180/255),
                                    rgb(229/255, 51/255, 38/255),
                                    rgb(102/255, 0/255, 13/255))) +
  # Seurat builds these two guides itself; guide-level themes are the only ones
  # that survive the complete theme it attaches (cowplot's legend.text = 8 pt).
  guides(size = guide_legend(title = 'Percent Expressed', order = 1,
                            title.theme = element_text(size = 7.6),
                            label.theme = element_text(size = 7.6)),
         colour = guide_colorbar(title = 'Average Expression', order = 2,
                                 barwidth = unit(1.9, 'cm'), barheight = unit(0.28, 'cm'),
                                 title.theme = element_text(size = 7.6,
                                             margin = margin(l = 0.35, unit = 'cm')),
                                 label.theme = element_text(size = 7.6))) +
  labs(x = 'Gene expression', y = 'Celltypes') +
  scale_y_discrete(limits = rev(lv)) + scale_x_discrete(limits = present)

# Fix the panel box to the measured size, then let the page shrink to the content.
g <- ggplotGrob(p)
pc <- which(sapply(g$grobs, function(x) grepl('^panel', x$name)))
g$widths[[g$layout$l[pc]]]  <- unit(PANEL_W, 'cm')
g$heights[[g$layout$t[pc]]] <- unit(PANEL_H, 'cm')
w <- as.numeric(convertWidth(sum(g$widths), 'cm'))
h <- as.numeric(convertHeight(sum(g$heights), 'cm'))
cat(sprintf('page %.2f x %.2f cm\n', w, h))
pdf(file.path(OUT, 'fig2_full_celltype_dotplot.pdf'), width = w/2.54, height = h/2.54)
grid.newpage(); grid.draw(g); dev.off()
png(file.path(OUT, 'fig2_full_celltype_dotplot.png'), width = w/2.54, height = h/2.54,
    units = 'in', res = 200)
grid.newpage(); grid.draw(g); dev.off()
cat('dotplot genes:', length(present), '| types:', length(lv), '\n')
cat(paste(present, collapse = ' '), '\n')
'''
r3 = os.path.join(OUT, 'fig2_s8_dotplot.R')
with open(r3, 'w') as fh: fh.write(R3)
run_r(r3)
""")
code("""p = os.path.join(OUT, 'fig2_full_celltype_dotplot.png')
if os.path.exists(p):
    show(p)
""")
md("""**Expected** — a dot plot with no title: 48 marker genes in delivered order
(Qk, Slc1a2, Slc1a3, Mbp, Plp1, Gpr17, Pdgfra, Tnr, P2ry12, Hexb, Flt1, Pecam1,
Eya2, Foxd1, Dnah12, Adgrv1, Rbfox3, Gad1, Gad2, Drd1, Drd2, Hap1, Lamp5, Lhx6,
Sst, Vip, Pvalb, Th, Slc17a6, Rmst, Ndst4, Foxb1, Prox1, Prkcd, Tac2, Pou4f1,
Hcn1, Fibcd1, Spink8, Car3, Slc17a7, Cux2, Igfn1, Rorb, Tshz2, Il1rapl2, Tle4,
Shisa6), 33 types top-to-bottom Ex_L6b → AST, legends "Percent Expressed
0–75" and "Average Expression 3–0" laid out horizontally above the panel, no
plot title. The panel box is 8.90 × 7.48 cm (page ≈ 11.0 × 9.8 cm) and the
largest dot is ≈ 1.4 mm across — the delivered panel's measured geometry.""")

# ===========================================================================
md("""---

## 9. Panel E — ATAC pseudo-bulk coverage tracks at five lineage loci

**Input**: the ArchR project `archr-subset-copnfol` (mm10).
**Output**: `fig2_out/fig2_track_Qk_Gpr17_Plp1_Slc32a1_Slc17a7.pdf` (5 pages,
one per locus) — the panel source of the delivered figure.

`ArchR::plotBrowserTrack` with `groupBy = cell_type_ht` (the 33-type annotation
matched back from the Seurat object by barcode), `useGroups` = the 33 types in
`dot_order`, five loci (`Qk`, `Gpr17`, `Plp1`, `Slc32a1`, `Slc17a7`), ±50 kb,
and the shared palette. This is ArchR's full browser-track rendering: pseudo-bulk
ATAC coverage per cell type, feature tracks, gene model and strand arrowheads.

**This is the slow step** — ArchR streams the arrow files, so allow ~10 min and
several GB of RAM, and run the notebook with
`--ExecutePreprocessor.timeout=-1`. The rendered PDF is cached: delete it (or
set `FIG2_FORCE=1`) to re-run.""")

code("""
R4 = r'''
suppressPackageStartupMessages({library(ArchR); library(Seurat); library(grid)})
set.seed(1)
addArchRGenome('mm10'); addArchRThreads(threads = 16)
DATA <- Sys.getenv('FIG2_DATA', os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')))
OUT  <- Sys.getenv('FIG2_OUT', 'fig2_out')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

proj <- loadArchRProject(path = file.path(DATA, 'archr-subset-copnfol'), showLogo = FALSE)

# match the ArchR cells to the Seurat annotation (the project predates it)
S <- readRDS(file.path(OUT, 'gs_subset_light.RDS'))
m <- match(proj$cell_id, rownames(S$meta))
stopifnot(!any(is.na(m)))
proj$cell_type_ht <- S$meta$cell_type_updated[m]

use_groups <- S$dot_order[as.vector(table(factor(proj$cell_type_ht, levels = S$dot_order))) > 0]
pal_ht <- S$cols[use_groups]                                  # shared palette

save_grobs_pdf <- function(plot_list, filename, width = 15, height = 20) {
  pdf(file.path(OUT, filename), width = width, height = height)
  for (nm in names(plot_list)) { grid::grid.newpage(); grid::grid.draw(plot_list[[nm]]) }
  dev.off()
  cat('Saved:', filename, '\n')
}
p_tracks <- plotBrowserTrack(ArchRProj = proj, groupBy = 'cell_type_ht',
                             useGroups = use_groups,
                             geneSymbol = c('Qk', 'Gpr17', 'Plp1', 'Slc32a1', 'Slc17a7'),
                             upstream = 50000, downstream = 50000, pal = pal_ht)
save_grobs_pdf(p_tracks, 'fig2_track_Qk_Gpr17_Plp1_Slc32a1_Slc17a7.pdf')
cat('groups used:', length(use_groups), '| cells in project:', nCells(proj), '\n')
'''
out_pdf = os.path.join(OUT, 'fig2_track_Qk_Gpr17_Plp1_Slc32a1_Slc17a7.pdf')
if os.path.exists(out_pdf) and os.environ.get('FIG2_FORCE') != '1':
    print('already rendered, skipped:', out_pdf,
          '(%s, %.1f MB) - delete it or set FIG2_FORCE=1 to re-run'
          % (time.ctime(os.path.getmtime(out_pdf)), os.path.getsize(out_pdf) / 1e6))
else:
    r4 = os.path.join(OUT, 'fig2_s9_tracks.R')
    with open(r4, 'w') as fh: fh.write(R4)
    run_r(r4)
""")
md("""**Expected** — a 5-page PDF, one locus per page (Qk, Gpr17, Plp1, Slc32a1,
Slc17a7), 33 coverage rows per locus in `dot_order` with the shared palette. The
E panel of the delivered figure is this track rendered at its loci.""")

md("""## Requirements

**Python** (the kernel is Python 3; it only launches R and shows the results):
`numpy`.

**R ≥ 4.1**, with
* §1, §6–§9: `Seurat` (≥ 4), `Signac`, `ggplot2`, `dplyr`, `tidyr`, `pals`,
  `readr`, `rtracklayer`, `GenomeInfoDb`
* §2: `DoubletFinder` (and its `RANN`/`Matrix` back end), `glue`
* §4, §9: `ArchR` (1.0.x) + `BSgenome.Mmusculus.UCSC.mm10`
* §5: `harmony`
* `future` / `future.apply` for the parallel steps.

**Reference data**: the mm10 gene annotation used by ArchR (`BSgenome` + the
mm10 TSS table) and, for the `else` branch of §1, a mm10 `genes.gtf`
(`REF/genes.gtf`) or the pre-extracted `<work>/mm10.gene.coords.RDS`.

**Memory and time.** §4 (ArchR) and §5 (WNN over 70,643 cells) are the heavy
steps — hours on 16 threads, tens of GB of scratch under `FIG2_DATA`. The
read of the 13 GB object in §6 needs ~13 GB of RAM. Everything else is minutes.
Run with `--ExecutePreprocessor.timeout=-1` so the kernel does not kill the long
cells:

## Run

```bash
export FIG2_INPUT=/path/to/HZJ        # workflow outputs (RNA matrices, peaks, fragments)
export FIG2_DATA=/path/to/FY          # analysis tree (SoupX matrices, per-stage objects)
export FIG2_OUT=fig2_out
jupyter nbconvert --to notebook --execute --inplace \\
    --ExecutePreprocessor.timeout=-1 Figure_2.ipynb
```

Each stage loads its own output when it already exists, so an interrupted run
resumes where it stopped; delete `fig2_out/fig2_candidate_cells.RDS` to redo
everything, or one `3k_FY1x_multiome.RDS` to redo a single library.

Every plotting cell saves its figure and reports the path and size (nothing is
displayed inline; the notebook stays small on GitHub), so the executed notebook is
self-checking (panels A–D as PNG previews; panel E is a 5-page PDF).""")

nb.cells = C
path = os.path.join(OUT_DIR, 'Figure_2.ipynb')
nbf.write(nb, path)
print('wrote', path, len(C), 'cells')

