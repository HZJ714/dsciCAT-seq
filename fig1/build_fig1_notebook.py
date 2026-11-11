#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build Figure_1.ipynb — reproduction of Figure 1 (cell-line benchmarking).

Design (per author requirements):
  * reproduction first: for every panel the notebook states the *input files* and
    gives the shortest complete code that turns them into the panel;
  * styling is deliberately minimal (no journal layout, no panel letters, no
    colour-bar placement, no compositing) — those are done by hand in the
    production scripts and are not part of the science;
  * the upstream step that the whole figure depends on — turning a mixed-species
    (GRCh38 + mm10) count matrix into per-species matrices and per-cell metrics —
    is written out in full (Section 1), because panels B/D/E/F all consume its
    output;
  * provenance is stated honestly: three technologies were re-processed through
    our own pipeline with the reads down-sampled (uniform upstream), the other
    four are published matrices used as-is (EmptyDrops_CR-filtered by their
    authors, no depth matching) — see Section 0.2;
  * panel C's ATAC rows are Signac code (Section 5);
  * panel A is a vector illustration, not a plot;
  * comments in English, code as short as possible while complete.

Run:  python build_fig1_notebook.py
Out:  <repo>/fig1/Figure_1.ipynb
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
md("""# Figure 1 — dsciCAT-seq: library design and benchmarking against mainstream single-cell multiomic methods in two cell lines

Every panel here is a *measurement on files that exist before the figure is drawn*.
This notebook shows, for each panel, which files go in and the shortest complete
code that produces the numbers and the plot. Journal layout (panel letters,
colour keys, exact cm geometry, compositing) is done separately and is not
reproduced here.

| Panel | What it shows | Input | Code |
|---|---|---|---|
| A | combined experimental + computational workflow | — | none: vector illustration |
| B | species collision in both channels of library HG14 | `HG14/sub/ISSAAC_joint_metrics_sub.csv` | Python (§2) |
| C | genome tracks at one locus per cell line | `*_fragments.tsv.gz` + per-species peak matrices | R/Signac (§5) |
| D | RNA capture: lg(UMIs), lg(detected genes) | 26 × `metrics.csv` | Python (§3) |
| E | ATAC capture: lg(unique fragments) | the same 26 × `metrics.csv` | Python (§3) |
| F | RNA similarity between platforms | per-species RNA `matrix.mtx` + `features.tsv` | Python (§4) |

```
WDL pipeline (UPSTREAM.md): FASTQ -> mixed-species (GRCh38 + mm10) matrices
        |                              + per-bead/per-cell species counts
        v
§1  bead->cell merge -> species call (>0.8) -> split matrix -> metrics.csv
        |                                          |              |
        |                                          v              v
        |                              §4 panel F (pseudobulk   §3 panels D,E
        |                                  -> CPM -> Pearson)   (26 metrics.csv)
        v
   joint per-cell hs/mm table -> §2 panel B (collision)
                              -> §5 panel C (ATAC tracks, R/Signac)
```

§1 is the step that has to be understood to read the rest: **all six
benchmarking platforms are reduced to the same three-column cell table**
(`cell, nCounts, nFeatures`) and the same 10x-style matrix directory
(`matrix.mtx`, `features.tsv`, `barcodes.tsv`).""")

# ===========================================================================
md("""---

## 0. Environment, inputs, outputs

| Variable | Default | Contents |
|---|---|---|
| `FIG1_DATA` | `<raw>` | per-platform results tree, layout in §0.1 |
| `FIG1_OUT` | `fig1_out` | created here; every panel PDF / CSV lands here |

Python ≥3.8 with `numpy pandas scipy matplotlib seaborn`. Section 5 additionally
needs `Rscript` with `Seurat`, `Signac`, `Matrix`, `readr`, `tidyr`, `magrittr`,
`ggplot2`, `patchwork` and `EnsDb.Hsapiens.v86` / `EnsDb.Mmusculus.v79`.""")

code("""
import matplotlib
matplotlib.use('Agg')
import os, sys, subprocess, time
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def show(*paths):
    # report the saved file only -- nothing is displayed inline (section 0.2)
    for q in paths:
        print('%s  (%.0f KB)' % (q, os.path.getsize(q) / 1024.0))
from scipy.io import mmread
from scipy.cluster.hierarchy import linkage, leaves_list

HZJ = os.environ.get('FIG1_DATA', '')
OUT = os.environ.get('FIG1_OUT', 'fig1_out')
os.makedirs(OUT, exist_ok=True)

# house style = the only styling in this notebook: Arial-metric font, vector text
from matplotlib import font_manager as _fm
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
sns.set_style('ticks')
print('HZJ =', HZJ, '| exists:', os.path.isdir(HZJ), '| OUT =', OUT)""")

md("""### 0.1 Input files

Paths are relative to `FIG1_DATA`. `Human` / `Mouse` are the per-species
directories produced by §1 (HEK293T = human reads, NIH3T3 = mouse reads).

| Used by | File (relative to `FIG1_DATA`) | What it is |
|---|---|---|
| §1, §2 | `HG/HG14/sub/outs/per_cell_barcode_total_read_count_GRCh38.tsv`, `..._mm10.tsv`, `.../per_cell_barcode_total_fragment_count_GRCh38.tsv`, `..._mm10.tsv` | STAR-Solo / chromap per-bead counts against each genome separately (mixed-species reference) |
| §1 | `HG/HG14/HG14.barcodeMerge.tsv` (same layout in every `HG*/sub`) | bead barcode → cell barcode, from the ATAC bead deconvolution (`d2c`, see §0.2) |
| §1 | `HG/HG14/sub/filtered_mtx/`, `HG/HG14/FilterMatrix/`, `HG/HG14/sub/ATAC/{Human,Mouse}/` | raw RNA gene × cell and ATAC peak × cell matrices |
| §2 | `HG/HG14/sub/ISSAAC_joint_metrics_sub.csv` | §1 output: one row per *paired* cell, `count_hs/mm` for both channels + `species_rna/atac/joint` |
| §3 | 26 × `metrics.csv` (one per platform × cell line × modality) | `cell, nCounts, nFeatures` — RNA: UMIs / detected genes; ATAC: unique fragments / detected peaks |
| §4 | `HG/HG12/sub/RNA/{Human,Mouse}/`, `HG/HG14/sub/RNA/{Human,Mouse}/`, `scCAT-seq-HT/IL9/RNA-sub/{Human,Mouse}/`, `ISAAC/sub/RNA-sub/{Human,Mouse}/`, `10x/RNA-sub0.74/outs/Solo.out/Gene/filtered/`, `sciCAR/{HEK293T,NIH3T3}/RNA/`, `SUM-seq/RNA/3t3filtered/corrected_matrix/`, `paired-seq/Paired_seq_Cell_Mix_matrix/HEK293T/RNA/` | §1 output (first four) or published matrix (last four), 10x-style triplets |
| §5 | `*/sub/*_fragments.tsv.gz` (+ `.tbi`) and the matching `ATAC/{Human,Mouse}/` matrices | fragment records for Signac `CoveragePlot` / `TilePlot` |
| reference | `<work>/hg38.gene.coords.RDS`, `mm10.gene.coords.RDS` | gene models for the track annotation |

The 26 `metrics.csv` paths are listed verbatim in §3; §0.3 checks that all of them exist.""")

md("""### 0.2 Provenance: what was re-processed and what was taken as a matrix

This is the one thing a reader has to know to interpret panels D–F. The
benchmarking set is **not** homogeneous:

| Platform | Cell line | Upstream | Sequencing depth |
|---|---|---|---|
| dsciCAT-seq-50k (HG12) | both | **our WDL, re-run on down-sampled FASTQ** (`seqtk sample`, then re-aligned) | matched by re-processing (`sub/`) |
| dsciCAT-seq-100k (HG14) | both | our WDL, re-run on down-sampled FASTQ | matched (`sub/`) |
| ISSAAC-seq (Droplet) | both | our WDL (`ISSAAC-full-10xformat`), re-run on down-sampled FASTQ | matched (`ISAAC/sub/`, cf. `ISAAC/RNA-full/`) |
| HT-scCAT-seq (IL9) | both | **our WDL re-run on the vendor's FASTQ** (`ISSAAC-full-10xformat` / `ISSAAC-ATAC`), both a full and a down-sampled run | matched: the figure uses `IL9/RNA-sub`, `IL9/ATAC-sub` (cf. `*-full`) |
| 10x Multiome | NIH3T3 | **our WDL re-run on the vendor FASTQ** (STAR-Solo + chromap/macs2, see `10x/RNA-sub0.74/outs/Log.final.out`) | matched by down-sampling the reads (RNA 0.74, ATAC 0.62 of the library) |
| sciCAR-seq | both | published matrix (GEO GSE117089), used as-is | **not** depth-matched |
| SUM-seq | NIH3T3 | published matrix (GEO), used as-is | not depth-matched |
| Paired-seq | HEK293T | published matrix (GEO GSE130399), used as-is | not depth-matched |

**One upstream detail behind the two `dsciCAT-seq` rows.** The `HG*` runs went through
two pipelines, not one: the ISSAAC-style pipeline that produced the matrices and
`metrics.csv` above, and our own dsciCAT-seq pipeline — which is here only for the ATAC
barcode merge. The dsciCAT bead carries its barcode in a structure that neither the ISSAAC
nor the 10x chemistry has, so the ATAC FASTQ is aligned with `chromap` and the bead
barcodes are re-assembled by `d2c`, whose output is the `*.barcodeMerge.tsv` file §1 uses
to turn beads into cells. Every library was merged under the same `forcefrag2000`
fragment-length cap, so the merge is a constant of the pipeline and not a per-library
parameter.

Consequences, both of which are visible in the figure:

1. For the first five platforms the cell calling is *ours* (species split, §1),
   so their `metrics.csv` and their matrix describe the same cells.
   For sciCAR / SUM / Paired the matrix and the QC table are both the authors'
   own products; in two of them they are not even the same cell set
   (`sciCAR/HEK293T/RNA`: 465 rows in `metrics.csv` vs 812 barcodes;
   `paired-seq/.../RNA`: 1,283 vs 593). Panels D–E therefore plot *the cells each
   platform's own QC table describes*.
2. Panel F normalises by depth *inside* the comparison (CPM on pseudobulk, §4),
   which is the only place where the missing depth matching of the last three
   rows is compensated. It is not compensated in panels D and E — those are
   deliberately shown at each platform's own depth, which is part of what
   "capture efficiency at the sequenced depth" means.

The measured effect of the down-sampling (median UMIs per cell, HEK293T RNA,
full vs `sub`) — printed by the next cell:

| | cells (full) | med UMIs | cells (sub) | med UMIs |
|---|---|---|---|---|
| dsciCAT-seq-50k | 12,334 | 3,126 | 12,157 | 1,897 |
| dsciCAT-seq-100k | 20,116 | 2,168 | 19,790 | 1,587 |
| ISSAAC-seq | 2,664 | 21,468 | 2,571 | 11,912 |
| 10x Multiome | 5,648 | 42,662 | 6,603 | 11,412 |
| HT-scCAT-seq (IL9) | 1,953 | 5,105 | 1,885 | 3,056 |""")

md("""### 0.3 Input check

This cell reads only the headers of every file the notebook needs, so it fails
fast on a machine where the data live elsewhere.""")

code("""
# 0.3 existence + shape check for every input of this notebook (reads headers only)
METRICS = {  # (platform, cell line, modality) -> metrics.csv, relative to FIG1_DATA
 ('ISSAAC-seq(Droplet)','HEK293T','RNA'):  'ISAAC/sub/RNA-sub/human_metrics.csv',
 ('ISSAAC-seq(Droplet)','HEK293T','ATAC'): 'ISAAC/sub/ATAC-sub/human_metrics.csv',
 ('ISSAAC-seq(Droplet)','NIH3T3','RNA'):   'ISAAC/sub/RNA-sub/mouse_metrics.csv',
 ('ISSAAC-seq(Droplet)','NIH3T3','ATAC'):  'ISAAC/sub/ATAC-sub/mouse_metrics.csv',
 ('10xMultiome','NIH3T3','RNA'):           '10x/RNA-sub0.74/metrics.csv',
 ('10xMultiome','NIH3T3','ATAC'):          '10x/ATAC_sub0.62/ATAC-sub/metrics.csv',
 ('HT-scCAT-seq','HEK293T','RNA'):         'scCAT-seq-HT/IL9/RNA-sub/Human/metrics.csv',
 ('HT-scCAT-seq','HEK293T','ATAC'):        'scCAT-seq-HT/IL9/ATAC-sub/Human/metrics.csv',
 ('HT-scCAT-seq','NIH3T3','RNA'):          'scCAT-seq-HT/IL9/RNA-sub/Mouse/metrics.csv',
 ('HT-scCAT-seq','NIH3T3','ATAC'):         'scCAT-seq-HT/IL9/ATAC-sub/Mouse/metrics.csv',
 ('dsciCAT-seq-50k','HEK293T','ATAC'):     'HG/HG12/sub/ATAC/Human/metrics.csv',
 ('dsciCAT-seq-50k','NIH3T3','RNA'):       'HG/HG12/sub/RNA/mouse_metrics.csv',
 ('dsciCAT-seq-50k','NIH3T3','ATAC'):      'HG/HG12/sub/ATAC/Mouse/metrics.csv',
 ('dsciCAT-seq-50k','HEK293T','RNA'):      'HG/HG12/sub/RNA/human_metrics.csv',
 ('dsciCAT-seq-100k','HEK293T','ATAC'):    'HG/HG14/sub/ATAC/Human/metrics.csv',
 ('dsciCAT-seq-100k','NIH3T3','RNA'):      'HG/HG14/sub/RNA/mouse_metrics.csv',
 ('dsciCAT-seq-100k','NIH3T3','ATAC'):     'HG/HG14/sub/ATAC/Mouse/metrics.csv',
 ('dsciCAT-seq-100k','HEK293T','RNA'):     'HG/HG14/sub/RNA/human_metrics.csv',
 ('Paired-seq','HEK293T','RNA'):           'paired-seq/Paired_seq_Cell_Mix_matrix/HEK293T/RNA/metrics.csv',
 ('Paired-seq','HEK293T','ATAC'):          'paired-seq/Paired_seq_Cell_Mix_matrix/HEK293T/ATAC/metrics.csv',
 ('sciCAR-seq','HEK293T','RNA'):           'sciCAR/HEK293T/RNA/metrics.csv',
 ('sciCAR-seq','HEK293T','ATAC'):          'sciCAR/HEK293T/ATAC/metrics.csv',
 ('sciCAR-seq','NIH3T3','RNA'):            'sciCAR/NIH3T3/RNA/metrics.csv',
 ('sciCAR-seq','NIH3T3','ATAC'):           'sciCAR/NIH3T3/ATAC/metrics.csv',
 ('SUM-seq','NIH3T3','RNA'):               'SUM-seq/RNA/3t3filtered/corrected_matrix/metrics.csv',
 ('SUM-seq','NIH3T3','ATAC'):              'SUM-seq/ATAC/3t3filtered/metrics.csv',
}
JOINT = 'HG/HG14/sub/ISSAAC_joint_metrics_sub.csv'
RNA_MTX = {  # platform -> (matrix dir relative to FIG1_DATA, matrix is genes x cells?)
 'dsciCAT-seq-50k':     ('HG/HG12/sub/RNA/{sp}', False),
 'dsciCAT-seq-100k':    ('HG/HG14/sub/RNA/{sp}', False),
 'HT-scCAT-seq':        ('scCAT-seq-HT/IL9/RNA-sub/{sp}', False),
 'ISSAAC-seq(Droplet)': ('ISAAC/sub/RNA-sub/{sp}', False),
 '10xMultiome':         ('10x/RNA-sub0.74/outs/Solo.out/Gene/filtered', True),
 'sciCAR-seq':          ('sciCAR/{line}/RNA', True),
 'SUM-seq':             ('SUM-seq/RNA/3t3filtered/corrected_matrix', False),
 'Paired-seq':          ('paired-seq/Paired_seq_Cell_Mix_matrix/HEK293T/RNA', False),
}
missing = [p for p in list(METRICS.values()) + [JOINT] if not os.path.exists(os.path.join(HZJ, p))]
for _p, _ in RNA_MTX.values():
    for _sp, _line in (('Human', 'HEK293T'), ('Mouse', 'NIH3T3')):
        if not os.path.exists(os.path.join(HZJ, _p.format(sp=_sp, line=_line), 'matrix.mtx')):
            missing.append(_p.format(sp=_sp, line=_line) + '/matrix.mtx')
print('%d metrics + %d matrix dirs + 1 joint table checked; missing: %s'
      % (len(METRICS), len(RNA_MTX), missing or 'none'))

for tag, a, b in [('dsciCAT-seq-50k', 'HG/HG12/Solo/human_metrics.csv', 'HG/HG12/sub/RNA/human_metrics.csv'),
                  ('dsciCAT-seq-100k', 'HG/HG14/Solo/human_metrics.csv', 'HG/HG14/sub/RNA/human_metrics.csv'),
                  ('ISSAAC-seq', 'ISAAC/RNA-full/human_metrics.csv', 'ISAAC/sub/RNA-sub/human_metrics.csv'),
                  ('10xMultiome', '10x/outs/metrics.csv', '10x/RNA-sub0.74/metrics.csv'),
                  ('HT-scCAT-seq', 'scCAT-seq-HT/IL9/RNA-full/outs/Human/metrics.csv',
                                   'scCAT-seq-HT/IL9/RNA-sub/Human/metrics.csv')]:
    f = pd.read_csv(os.path.join(HZJ, a)); s = pd.read_csv(os.path.join(HZJ, b))
    print('  %-16s full %5d cells / %6.0f UMIs   sub %5d cells / %6.0f UMIs'
          % (tag, len(f), f.nCounts.median(), len(s), s.nCounts.median()))""")

md("""**Expected** — nothing missing; 26 metrics tables, 8 matrix directories, 1 joint
table; the five down-sampled platforms show the medians quoted in §0.2
(3,126→1,897 · 2,168→1,587 · 21,468→11,912 · 42,662→11,412 · 5,105→3,056).""")

md("""### 0.4 Software

Everything below was executed with these versions; the cell prints what is in use
on your machine.""")

code("""import importlib.metadata as im
print('python', sys.version.split()[0])
for p in ['numpy', 'pandas', 'scipy', 'matplotlib', 'seaborn']:
    try: print('  %-12s %s' % (p, im.version(p)))
    except Exception: print('  %-12s -' % p)
try:
    print(subprocess.run(['Rscript', '-e',
        'cat(R.version.string, "|", paste0("Seurat ", packageVersion("Seurat")), "|", paste0("Signac ", packageVersion("Signac")))'],
        capture_output=True, text=True, timeout=300).stdout.strip())
except Exception as e:
    print('  Rscript not available (only needed for §5):', e)""")

# ===========================================================================
md("""---

## 1. Mixed-species matrix → per-species matrices and per-cell metrics

**Input** (from the WDL runs, see `UPSTREAM.md`): the libraries are human + mouse
nuclei co-sequenced against a combined GRCh38 + mm10 reference, so the aligner
reports, for every barcode, how many reads (RNA) or fragments (ATAC) landed on
each genome, and the count matrix carries both genomes' features with a
`GRCh38_` / `mm10___` prefix. Four steps:

1. **bead → cell merge** (dsciCAT-seq only): the ATAC bead deconvolution
   (`d2c merge`) assigns several bead barcodes to one cell barcode; the
   per-bead species counts are summed per cell (`barcodeMerge.tsv`).
   For 10x/ISSAAC the barcode *is* the cell and this step is skipped.
2. **species call**, threshold 0.8 on the human (resp. mouse) fraction:
   `Human` if `count_hs/total > 0.8`, `Mouse` if `count_mm/total > 0.8`, else
   `Mix`. The threshold is applied identically to both channels; the measured
   boundary in the delivered table is 0.8001 (RNA) / 0.8004 (ATAC), i.e. >0.8.
   `Mix` cells are *not* used for any downstream matrix — they are only plotted
   in panel B.
3. **matrix split**: keep the rows (cells) of one species, drop the other
   species' features by prefix, strip the prefix, write a 10x-style directory.
4. **metrics**: `nCounts` = column sum, `nFeatures` = number of non-zeros, per cell.

Everything below is pure sparse arithmetic: a 318k × 30k peak matrix never
becomes dense.""")

code("""
# ---- 1a. per-bead species counts -> per-cell (dsciCAT-seq needs the bead merge) ----
def read_counts(path, value_col=0, key_col=1):
    \"\"\"Two-column TSV of the aligner's per-barcode counts against ONE genome.
    RNA (`..._read_count_*.tsv`, `rna_*.tsv`) and ATAC (`..._fragment_count_*.tsv`,
    `fra_*.tsv`) differ in column order, hence the two indices.\"\"\"
    df = pd.read_csv(path, sep='\\t', header=None)
    return df.set_index(df.columns[key_col])[df.columns[value_col]].rename('n')

def bead_to_cell(per_bead, merge_tsv):
    \"\"\"Sum bead-level counts into cell barcodes using the d2c bead->cell table.\"\"\"
    m = pd.read_csv(merge_tsv, sep='\\t', header=None, index_col=0).iloc[:, 0].rename('cell')
    return (pd.concat([per_bead, m], axis=1, join='inner')
              .groupby('cell')['n'].sum())

def species_call(hs, mm, thresh=0.8):
    \"\"\"hs / mm: per-cell counts on each genome (Series). Returns a DataFrame with
    count_hs, count_mm, total, species in {Human, Mouse, Mix}.\"\"\"
    df = pd.concat([hs.rename('count_hs'), mm.rename('count_mm')], axis=1, join='inner')
    df['total'] = df.count_hs + df.count_mm
    f = df.count_hs / df.total
    df['species'] = np.where(f > thresh, 'Human', np.where(1 - f > thresh, 'Mouse', 'Mix'))
    return df

# Example: library HG14 (dsciCAT-seq, 100k cells, 16 barcodes) -> panel B's table.
# Run once; ~10 min, dominated by the per-bead tables (hundreds of millions of rows).
def joint_table(hzj, lib='HG14', thresh=0.8, write=True):
    sub = os.path.join(hzj, 'HG', lib, 'sub')
    merge = os.path.join(hzj, 'HG', lib, lib + '.barcodeMerge.tsv')
    out = {}
    for chan, pat in (('rna', 'per_cell_barcode_total_read_count_%s.tsv'),
                      ('atac', 'per_cell_barcode_total_fragment_count_%s.tsv')):
        for sp, col in (('hs', 'GRCh38'), ('mm', 'mm10')):
            s = bead_to_cell(read_counts(os.path.join(sub, 'outs', pat % col), 0, 1), merge)
            out[(chan, sp)] = s
    j = pd.concat({'rna':  species_call(out[('rna', 'hs')],  out[('rna', 'mm')], thresh),
                   'atac': species_call(out[('atac', 'hs')], out[('atac', 'mm')], thresh)},
                  axis=1, keys=['rna', 'atac'])
    j.columns = ['%s_%s' % (b, a) for a, b in j.columns]      # count_hs_rna, species_atac, ...
    j['species_joint'] = np.select(
        [(j.species_rna == j.species_atac) & (j.species_rna != 'Mix'),
         (j.species_rna == 'Mix') | (j.species_atac == 'Mix')],
        [j.species_rna, 'Mix'], default='Mix')
    return j

print('1a helpers defined')""")

code("""
# ---- 1b. split one mixed-species matrix into Human/ + Mouse/ ----
import scipy.sparse
from scipy.io import mmwrite

def load_mtx(prefix):
    \"\"\"10x-style triplet -> (sparse cells x features, feature names, barcodes).\"\"\"
    X = mmread(os.path.join(prefix, 'matrix.mtx')).tocsr()
    feat = pd.read_csv(os.path.join(prefix, 'features.tsv'), sep='\\t', header=None).iloc[:, 0].values
    bc   = pd.read_csv(os.path.join(prefix, 'barcodes.tsv'), sep='\\t', header=None).iloc[:, 0].values
    if X.shape[0] == len(feat):            # genes x cells on disk -> transpose to cells x genes
        X = X.T.tocsr()
    return X, feat, bc

def split_matrix(prefix, df_species, dest_root):
    \"\"\"Write dest_root/Human and dest_root/Mouse: same cells as df_species,
    only this species' features, prefix stripped. Returns the row counts.\"\"\"
    X, feat, bc = load_mtx(prefix)
    keep_feat = {'Human': [f for f in feat if 'mm10' not in f],
                 'Mouse': [f for f in feat if 'GRCh38' not in f]}
    strip = {'Human': lambda f: f.replace('GRCh38_', ''), 'Mouse': lambda f: f.replace('mm10___', '')}
    n = {}
    for sp in ('Human', 'Mouse'):
        cells = df_species.index[df_species.species == sp]
        idx = np.isin(bc, cells)
        cols = np.where(np.isin(feat, keep_feat[sp]))[0]
        sub = X[idx][:, cols].tocoo()
        d = os.path.join(dest_root, sp); os.makedirs(d, exist_ok=True)
        pd.Series([strip[sp](feat[c]) for c in cols]).to_csv(os.path.join(d, 'features.tsv'),
                                                             sep='\\t', index=False, header=False)
        pd.Series(bc[idx]).to_csv(os.path.join(d, 'barcodes.tsv'), sep='\\t', index=False, header=False)
        mmwrite(os.path.join(d, 'matrix.mtx'),
                scipy.sparse.coo_matrix((sub.data, (sub.col, sub.row)), shape=(len(cols), idx.sum())))
        n[sp] = int(idx.sum())
    return n

print('1b helpers defined')""")

code("""
# ---- 1c. metrics.csv from a split matrix (this is what panels D/E consume) ----
def write_metrics(dest_root, sp):
    \"\"\"nCounts = UMIs/fragments per cell, nFeatures = detected genes/peaks per cell.
    Column-oriented on disk (cells x genes), so a column sum is a cell.\"\"\"
    X, _, _ = load_mtx(os.path.join(dest_root, sp))
    X = X.T.tocsc()                                   # features x cells
    out = pd.DataFrame({'nCounts': np.asarray(X.sum(axis=0)).ravel(),
                        'nFeatures': np.asarray((X != 0).sum(axis=0)).ravel()})
    out.to_csv(os.path.join(dest_root, sp, 'metrics.csv'), index_label='cell')
    return out

# 1d. the joint table of panel B, from the delivered per-species matrices:
#     (joint_table() above regenerates it from the per-bead files; the pipeline
#      output below is the file the figure was drawn from.)
df_joint = pd.read_csv(os.path.join(HZJ, JOINT))
print('joint table %d rows' % len(df_joint))
for c in ['species_rna', 'species_atac', 'species_joint']:
    v = df_joint[c].value_counts()
    print('  %-13s %s' % (c, {k: (int(n), round(n / len(df_joint) * 100, 2)) for k, n in v.items()}))
for c in ['rna', 'atac']:                              # the 0.8 rule, measured
    f = df_joint['count_hs_' + c] / df_joint['total_' + c]
    print('  %-4s Human min frac %.4f | Mix max frac %.4f  -> threshold 0.8'
          % (c, f[df_joint['species_' + c] == 'Human'].min(), f[df_joint['species_' + c] == 'Mix'].max()))
df_joint.to_csv(os.path.join(OUT, 'joint_metrics_sub.csv'), index=False)""")

md("""**Expected** — 42,816 paired cells;
RNA Human 19,766 (46.16%) / Mouse 20,097 (46.94%) / Mix 2,953 (6.90%);
ATAC Human 19,639 (45.87%) / Mouse 20,616 (48.15%) / Mix 2,561 (5.98%);
joint Human 19,397 / Mouse 19,597 / Mix 3,822.
The measured boundary (0.8001 / 0.8004 as the lowest human fraction still called
`Human`, 0.7998 / 0.7999 as the highest still called `Mix`) confirms threshold 0.8
for **both** channels.""")

# ===========================================================================
md("""---

## 2. Panel B — species collision in both channels (library HG14)

**Input**: the joint table of §1 (`count_hs/mm` for RNA and ATAC + `species_*`).
**Output**: `fig1b_contamination.pdf` — two square scatter plots, RNA on top,
ATAC below, one point per paired cell, coloured by the call of *that channel*.
The counts/percentages in the legends are the numbers quoted in the legend text.

The point of the panel is the *shape* of the `Mix` class: it is the band on the
diagonal between the two populations — cells with substantial counts on **both**
genomes, so neither fraction clears 0.8. They are not empty droplets: their median
total is 0.54× (RNA) / 0.71× (ATAC) of the `Human` population's median, and their
median human fraction is 0.38 (RNA) / 0.43 (ATAC), i.e. close to an even mixture.""")

code("""sp_order, sp_color = ['Human', 'Mouse', 'Mix'], {'Human': 'red', 'Mouse': 'blue', 'Mix': 'grey'}
spec = [('species_rna', 'count_hs_rna', 'count_mm_rna', 'GRCh38 RNA reads', 'mm10 RNA reads',
         (-5000, 100000), [0, 20000, 40000, 60000, 80000, 100000]),
        ('species_atac', 'count_hs_atac', 'count_mm_atac', 'GRCh38 ATAC fragments', 'mm10 ATAC fragments',
         (-10000, 150000), [0, 50000, 100000, 150000])]
fig, axs = plt.subplots(2, 1, figsize=(3.82 / 2.54, 6.53 / 2.54))
for ax, (spc, xc, yc, xl, yl, lim, tk) in zip(axs, spec):
    cnt = df_joint[spc].value_counts(); rat = cnt / cnt.sum() * 100
    for sp in sp_order:
        s = df_joint[df_joint[spc] == sp]
        ax.scatter(s[xc], s[yc], c=sp_color[sp], s=1.2, alpha=0.18, linewidths=0,
                   label='%s: %d (%.2f%%)' % (sp, cnt[sp], rat[sp]))
    ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_xticks(tk); ax.set_yticks(tk)
    ax.set_xticklabels([f'{v / 1e4:.0f}' for v in tk]); ax.set_yticklabels([f'{v / 1e4:.0f}' for v in tk])
    ax.set_xlabel(xl + r' ($\\times 10^4$)'); ax.set_ylabel(yl + r' ($\\times 10^4$)')
    ax.legend(loc='upper right', frameon=False, markerscale=2.6)
fig.tight_layout(); fig.savefig(f'{OUT}/fig1b_contamination.pdf'); show(f'{OUT}/fig1b_contamination.pdf')
print('fig1b_contamination.pdf')
print('  Mix vs Human, per channel:')
for c in ['rna', 'atac']:
    t = 'total_' + c
    hum, mix = df_joint.loc[df_joint['species_' + c] == 'Human', t], df_joint.loc[df_joint['species_' + c] == 'Mix', t]
    print('    %-4s median total Human %8.0f  Mix %8.0f  (ratio %.2f) | Mix median human fraction %.3f'
          % (c, hum.median(), mix.median(), mix.median() / hum.median(),
             (df_joint.loc[df_joint['species_' + c] == 'Mix', 'count_hs_' + c] /
              df_joint.loc[df_joint['species_' + c] == 'Mix', t]).median()))""")

md("""**Expected** — legend counts exactly as §1. The `Mix` class:
median total 6,637 (RNA) / 11,630 (ATAC) against 12,246 / 16,342 for `Human`
(ratios 0.54 / 0.71); median human fraction 0.377 (RNA) / 0.431 (ATAC).
34.1% (RNA) and 22.7% (ATAC) of `Mix` cells fall below the 10th percentile of the
`Human` distribution — a low-to-intermediate signal population, not a set of
essentially empty droplets.""")

# ===========================================================================
md("""---

## 3. Panels D and E — capture efficiency per platform

**Input**: 26 `metrics.csv` tables (the dict `METRICS` of §0.3). Each is
`cell, nCounts, nFeatures`; `nCounts` is UMIs for RNA and unique fragments for
ATAC, `nFeatures` is detected genes resp. peaks.
**Output**: `fig1d_rna_qc.pdf` (2 × 2 box plots: lg(UMIs) over lg(detected
genes), HEK293T over NIH3T3) and `fig1e_atac_qc.pdf` (1 × 2 violins of
lg(unique fragments)).

Two platforms are excluded from the figure although their `metrics.csv` exist in
the same tree — `scifiRNA-seq` (RNA only, no paired ATAC) and `dsciATAC-seq`
(ATAC only, single species): they were dropped from the benchmarking set, so the
26 tables above are the complete set.

Method order and colours are fixed (the same key is shared by D, E and F):""")

code("""ORDER = {'HEK293T': ['dsciCAT-seq-50k', 'dsciCAT-seq-100k', 'HT-scCAT-seq',
                     'ISSAAC-seq(Droplet)', 'Paired-seq', 'sciCAR-seq'],
         'NIH3T3':  ['dsciCAT-seq-50k', 'dsciCAT-seq-100k', 'SUM-seq', 'HT-scCAT-seq',
                     'ISSAAC-seq(Droplet)', '10xMultiome', 'sciCAR-seq']}
PALETTE = {'dsciCAT-seq-50k': '#AF2337', 'dsciCAT-seq-100k': '#AF2337', 'SUM-seq': '#FF6B35',
           'HT-scCAT-seq': '#2F3C28', 'ISSAAC-seq(Droplet)': '#3A7CA5', 'Paired-seq': '#EF476F',
           'sciCAR-seq': '#5F4BB6', '10xMultiome': '#ECC342'}

rows = []
for (m, cl, mod), rel in METRICS.items():
    d = pd.read_csv(os.path.join(HZJ, rel))
    d['Method'], d['Cell_line'], d['Modality'] = m, cl, mod
    rows.append(d)
cell_lines = pd.concat(rows, ignore_index=True)
cell_lines['logC'] = np.log10(cell_lines.nCounts + 1)     # lg(UMIs) / lg(fragments)
cell_lines['logF'] = np.log10(cell_lines.nFeatures + 1)   # lg(genes)  / lg(peaks)
print('%d cells in %d tables' % (len(cell_lines), len(METRICS)))
cell_lines.groupby(['Modality', 'Cell_line', 'Method'])[['nCounts', 'nFeatures']].median()""")

code("""# ---- Panel D: RNA box plots ----
fig, axs = plt.subplots(2, 2, figsize=(8.45 / 2.54, 10.16 / 2.54))
for j, cl in enumerate(['HEK293T', 'NIH3T3']):
    for i, (y, lab) in enumerate([('logC', 'lg(UMIs)'), ('logF', 'lg(detected genes)')]):
        ax = axs[i, j]
        sns.boxplot(data=cell_lines.query('Modality=="RNA" & Cell_line==@cl'), x='Method', y=y,
                    ax=ax, order=ORDER[cl], palette=[PALETTE[m] for m in ORDER[cl]],
                    width=0.5, saturation=1, fliersize=0.4)
        ax.set_ylabel(lab if j == 0 else ''); ax.set_xlabel('')
        ax.tick_params(axis='x', labelrotation=90, labelsize=6)
        if i == 0: ax.set_xticks([])
        if j == 1: ax.tick_params(axis='y', labelleft=False, left=False)
        if i == 1: ax.text(0.02, 0.02, cl, transform=ax.transAxes, fontweight='bold')
fig.tight_layout(); fig.savefig(f'{OUT}/fig1d_rna_qc.pdf'); show(f'{OUT}/fig1d_rna_qc.pdf')

# ---- Panel E: ATAC violins ----
fig, axs = plt.subplots(1, 2, figsize=(8.35 / 2.54, 3.66 / 2.54))
for j, cl in enumerate(['HEK293T', 'NIH3T3']):
    ax = axs[j]
    sns.violinplot(data=cell_lines.query('Modality=="ATAC" & Cell_line==@cl'), x='Method', y='logC',
                   ax=ax, order=ORDER[cl], palette=[PALETTE[m] for m in ORDER[cl]],
                   width=0.6, inner='box', linewidth=0.4)   # scale='area' is seaborn's default
    ax.set_ylabel('lg(unique fragments)' if j == 0 else ''); ax.set_xlabel(''); ax.set_xticks([])
    ax.text(0.02, 0.02, cl, transform=ax.transAxes, fontweight='bold')
    if j == 1: ax.tick_params(axis='y', labelleft=False, left=False)
fig.tight_layout(); fig.savefig(f'{OUT}/fig1e_atac_qc.pdf'); show(f'{OUT}/fig1e_atac_qc.pdf')
print('fig1d_rna_qc.pdf, fig1e_atac_qc.pdf')""")

md("""**Expected** — medians that the legend quotes (all reproduced from these tables;
`nCounts` / `nFeatures`):

| | ISSAAC | dsciCAT-50k | dsciCAT-100k | HT-scCAT | 10x | SUM | Paired | sciCAR |
|---|---|---|---|---|---|---|---|---|
| D, HEK293T UMIs / genes | 11,912 / 4,702 | 1,897 / 1,199 | 1,587 / 1,092 | 3,056 / 1,560 | — | — | 575 / 437 | 2,953 / 1,672 |
| D, NIH3T3 UMIs / genes | 12,401 / 4,100 | 2,169 / 936 | 1,808 / 857 | 4,340 / 1,672 | 11,412 / 3,329 | 1,840 / 1,343 | — | 2,709 / 1,248 |
| E, HEK293T fragments / peaks | 14,828 / 7,056 | 22,384 / 7,779 | 19,188 / 7,204 | 19,261 / 2,718 | — | — | 1,083 / 936 | 412 / 187 |
| E, NIH3T3 fragments / peaks | 15,680 / 6,663 | 20,319 / 6,668 | 17,585 / 6,182 | 16,827 / 5,481 | 30,124 / 12,070 | 6,953 / 4,312 | — | **6 / 3** |

Three readings survive from this panel pair: doubling dsciCAT-seq throughput
(50k → 100k) costs **less than 0.1 log** in every metric; the ATAC channel of
dsciCAT-seq is ~10× its own RNA channel; and sciCAR-seq's ATAC channel is not
measurable in NIH3T3 at all (median **6** fragments / 3 peaks per nucleus, against
412 / 187 in HEK293T) — that column is drawn for completeness. Note also that the
sciCAR ATAC matrices are not part of this tree: only their QC tables are
(`sciCAR/{HEK293T,NIH3T3}/ATAC/metrics.csv`), so panel E's sciCAR columns cannot be
recomputed from data, only transcribed.""")

code("""# the table above, computed
t = cell_lines.groupby(['Modality', 'Cell_line', 'Method'])[['nCounts', 'nFeatures']].median().round(0)
for mod in ['RNA', 'ATAC']:
    for cl in ['HEK293T', 'NIH3T3']:
        print(mod, cl, {m: (int(t.loc[(mod, cl, m)].nCounts), int(t.loc[(mod, cl, m)].nFeatures))
                        for m in ORDER[cl] if (mod, cl, m) in t.index})""")

# ===========================================================================
md("""---

## 4. Panel F — RNA similarity between platforms

**Input**: the per-species RNA matrices of §0.3 (`RNA_MTX`) — for the five
re-processed platforms these are the §1 output, for sciCAR / SUM / Paired the
published matrices.
**Output**: `fig1f_rna_similarity.pdf` plus `corr_HEK293T.csv` / `corr_NIH3T3.csv`
(the numeric matrices, so the panel can be re-drawn without re-reading 1 GB of
MatrixMarket).

The computation is four lines and is the whole panel:

1. **pseudobulk**: sum each matrix over its cells → one total-count vector per
   platform (all cells of that platform's species class, no filtering);
2. intersect the feature IDs across platforms (26,041 / 29,621 genes);
3. `log2(count + 1)` → **CPM** using each platform's own total → `log2(+1)` again;
4. Pearson correlation between platforms, rows/columns ordered by the leaves of
   an average-linkage clustering of that correlation matrix.

Step 3 is where the depth difference of §0.2 is removed: dividing by the
platform's own pseudobulk total makes the profile a composition, so a platform
sequenced 10× deeper is comparable. Note what that does *not* fix — a platform
with fewer captured genes has a different composition, which is the signal the
panel is meant to show.""")

code("""def gene_sums(prefix, genes_x_cells=False):
    \"\"\"Pseudobulk: total counts per gene. `genes_x_cells` marks the two published
    matrices that are stored transposed (sciCAR, 10x Solo).\"\"\"
    X = mmread(os.path.join(prefix, 'matrix.mtx')).tocsr()
    feat = pd.read_csv(os.path.join(prefix, 'features.tsv'), sep='\\t', header=None).iloc[:, 0]
    v = np.asarray(X.sum(axis=1 if genes_x_cells else 0)).ravel()
    return pd.Series(v, index=feat.values)

corr = {}
for cl in ['HEK293T', 'NIH3T3']:
    sp = 'Human' if cl == 'HEK293T' else 'Mouse'
    per = [gene_sums(os.path.join(HZJ, RNA_MTX[m][0].format(sp=sp, line=cl)), RNA_MTX[m][1])
               .rename(m) for m in ORDER[cl]]
    common = sorted(set.intersection(*[set(s.index) for s in per]))
    counts = np.log2(pd.concat([s[common] for s in per], axis=1) + 1)
    corr[cl] = np.log2(counts / counts.sum(axis=0) * 1e6 + 1).corr()
    corr[cl].to_csv(os.path.join(OUT, 'corr_%s.csv' % cl))
    print('%s: %d methods, %d shared genes' % (cl, len(ORDER[cl]), len(common)))""")

md("""**Expected** — `HEK293T: 6 methods, 26041 shared genes`,
`NIH3T3: 7 methods, 29621 shared genes`.
HEK293T: max r = **0.882** (dsciCAT-50k vs -100k), min = **0.377** (Paired-seq vs
dsciCAT-50k). NIH3T3: max = **0.883** (the same dsciCAT batch pair), min =
**0.718** (sciCAR vs dsciCAT). 10x Multiome is closest to ISSAAC-seq
(**0.8826** — within 0.0007 of the dsciCAT batch pair) and farthest from
sciCAR-seq (**0.7225**).""")

code("""for cl in ['HEK293T', 'NIH3T3']:
    c = corr[cl].copy(); np.fill_diagonal(c.values, np.nan)
    print(cl, 'max %.3f (%s)' % (c.max().max(), c.stack().idxmax()),
          '| min %.3f (%s)' % (c.min().min(), c.stack().idxmin()))
print('10x vs ISSAAC %.4f | 10x vs sciCAR %.4f'
      % (corr['NIH3T3'].loc['10xMultiome', 'ISSAAC-seq(Droplet)'],
         corr['NIH3T3'].loc['10xMultiome', 'sciCAR-seq']))""")

code("""fig, axs = plt.subplots(1, 2, figsize=(9.0 / 2.54, 6.06 / 2.54))
for ax, cl in zip(axs, ['HEK293T', 'NIH3T3']):
    order = leaves_list(linkage(corr[cl].values, method='average'))
    names = [ORDER[cl][i] for i in order]
    sns.heatmap(corr[cl].iloc[order, order], ax=ax, cmap='OrRd', vmin=0.2, vmax=1.0,
                xticklabels=names, yticklabels=False, cbar_kws={'label': 'Pearson r'})
    ax.set_title(cl, fontweight='bold'); ax.tick_params(axis='x', labelrotation=90, labelsize=6)
fig.tight_layout(); fig.savefig(f'{OUT}/fig1f_rna_similarity.pdf'); show(f'{OUT}/fig1f_rna_similarity.pdf')
print('fig1f_rna_similarity.pdf')""")

# ===========================================================================
md("""---

## 5. Panel C — genome tracks (ATAC rows, Signac/R)

Panel C stacks, for one locus per cell line, the ATAC coverage of the three
technologies (dsciCAT-seq-50k red, -100k green, HT-scCAT-seq blue), the same
coverage as single-cell tiles (300 cells), and the gene model.

**Input**: the per-species peak matrix + `fragments.tsv.gz` of each library
(§0.1), plus a gene-model RDS for the annotation row. The three libraries are
loaded as separate Seurat objects, merged on the shared peak set, and plotted
with `CoveragePlot` (grouped by `orig.ident`) over `TilePlot(tile.cells = 300)`.
Region: `chr13:91,340,000-91,370,000` (hg38, *Mir17hg*) for HEK293T and
`chr11:69,090,000-69,130,000` (mm10, *Vamp2*) for NIH3T3; the other loci of the
panel are the same call with a different region.""")

md("""> The NIH3T3 half is the same script with the mm10 paths
> (`HG/HG12/sub/mouse_fragments.tsv.gz`, `ATAC/Mouse`, `mm10.gene.coords.RDS`,
> `EnsDb.Mmusculus.v79`) and `region = "chr11-69090000-69130000"`.
>
> The three libraries are plotted **separately** and stacked with `patchwork`.
> They are deliberately *not* `merge`d into one object: `CoveragePlot` and
> `TilePlot` read the `fragments.tsv.gz` file, not the assay, so a merge would pay
> for a shared peak set (`BindMatrices` over ~780k peak names — measured: still
> running after 40 min, killed) that no row of this panel uses. Loading them
> separately also avoids intersecting the peak sets, which the delivered panel did
> not need either.
>
> Measured on this machine, 2026-10-21: §5 takes **10-15 min** (the notebook run
> was 15 min, a standalone run of the same script 9.6 min) and peaks at ~14 GB RAM;
> peak/cell counts
> dsciCAT-seq-50k 318,130 peaks × 16,165 cells · dsciCAT-seq-100k 345,484 × 27,760 ·
> HT-scCAT-seq 84,125 × 1,970 (after `min.cells = 10`), peak widths already
> 201-4614 bp so the 100-5000 bp filter is a no-op. RAM peaks at ~14 GB.
>
> Two Signac 1.9.0 quirks had to be worked around:
>
> * `CoveragePlot` returns a *nested patchwork* whose row titles ("Normalized
>   signal (range 0 - 130)", "Genes") are `annotate()` text, not axis titles - so a
>   plain `+theme(axis.title.y = element_blank())` does nothing and stacking three
>   rows makes those titles overlap. The `clean()` helper walks the tree (testing
>   `patchwork` before `ggplot`, since a patchwork *is* a ggplot) and drops both
>   the text layers and the inner axes. The locus is then stated once, as the
>   figure title.
> * `TilePlot` labels its x axis with bin indices (a 30 kb window with
>   `tile.cells = 300` reads 91,340,001-91,340,300) although the tiles sit over the
>   right window. The delivered panel has no x axis on that row either, so the
>   labels are blanked here.""")


code("""R_TRACKS = r\'\'\'
suppressPackageStartupMessages({library(Seurat); library(Signac); library(Matrix)
  library(readr); library(magrittr); library(ggplot2); library(patchwork)
  library(GenomicRanges); library(EnsDb.Hsapiens.v86)})
DATA   <- Sys.getenv("FIG1_DATA", os.environ.get('FIG1_DATA', ''))
OUT    <- Sys.getenv("FIG1_OUT", "fig1_out")
REGION <- "chr13-91340000-91370000"   # hg38, Mir17hg. NIH3T3: "chr11-69090000-69130000" (mm10)
libs <- list(                       # library -> (fragments.tsv.gz, per-species peak matrix dir)
  "dsciCAT-seq-50k"  = c("HG/HG12/sub/human_fragments.tsv.gz",  "HG/HG12/sub/ATAC/Human"),
  "dsciCAT-seq-100k" = c("HG/HG14/sub/human_fragments.tsv.gz",  "HG/HG14/sub/ATAC/Human"),
  "HT-scCAT-seq"     = c("scCAT-seq-HT/IL9/ATAC-sub/human_fragments.tsv.gz",
                         "scCAT-seq-HT/IL9/ATAC-sub/Human"))
col <- c("dsciCAT-seq-50k" = "#E15759", "dsciCAT-seq-100k" = "#59A14F", "HT-scCAT-seq" = "#4E79A7")

ann <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)   # gene models for the annotation row
seqlevels(ann) <- paste0("chr", seqlevels(ann)); genome(ann) <- "hg38"

load_one <- function(id, frags, mtx) {                  # peak matrix -> ChromatinAssay object
  f <- as.character(read_tsv(file.path(mtx, "features.tsv"), col_names = FALSE, show_col_types = FALSE)[[1]])
  b <- as.character(read_tsv(file.path(mtx, "barcodes.tsv"), col_names = FALSE, show_col_types = FALSE)[[1]])
  X <- readMM(file.path(mtx, "matrix.mtx"))
  if (nrow(X) == length(b)) X <- t(X)                   # stored cells x peaks -> peaks x cells
  pk  <- GenomicRanges::GRanges(f, sep = c(":", "-"))   # names are chr1:100000-100200
  sel <- which(width(pk) > 100 & width(pk) < 5000)      # the peak-width filter of the original
  ca <- CreateChromatinAssay(X[sel, , drop = FALSE] %>% set_rownames(f[sel]) %>% set_colnames(b),
        sep = c(":", "-"), fragments = frags, min.cells = 10, annotation = ann)
  o <- CreateSeuratObject(ca, assay = "ATAC"); o$orig.ident <- id
  DefaultAssay(o) <- "ATAC"; Idents(o) <- "orig.ident"
  cat(id, ":", ncol(o), "cells,", nrow(o), "peaks\n")
  o
}

# The three libraries are plotted separately, not merged: CoveragePlot/TilePlot read
# the fragments file, not the assay, so merging would only pay for a shared peak set
# (BindMatrices over ~780k peak names) that no row of this panel uses.
# Signac draws its row titles ("Normalized signal (range 0 - 130)", "Genes") as
# annotate() text and puts the real axes on inner plots of a patchwork tree, so the
# only way to get a stackable row is to walk the tree and blank them. The
# patchwork branch has to be tested first: a patchwork is also a ggplot.
kill <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              axis.title.x = element_blank(), axis.text.y = element_blank(),
              axis.ticks.y = element_blank(), axis.title.y = element_blank(),
              strip.text = element_blank(), plot.margin = margin(1, 1, 0, 1),
              legend.position = "none")
clean <- function(x) {
  if (inherits(x, "patchwork")) {
    for (i in seq_along(x)) x[[i]] <- tryCatch(clean(x[[i]]), error = function(e) x[[i]])
    return(x)
  }
  if (!inherits(x, "ggplot")) return(x)
  keep <- vapply(seq_along(x$layers), function(i) {
    d <- tryCatch(ggplot2::layer_data(x, i), error = function(e) NULL)
    if (is.null(d)) return(TRUE)
    bad <- FALSE
    for (cc in d) if (is.character(cc) && any(grepl("Normalized signal|^Genes$", cc))) bad <- TRUE
    !bad
  }, TRUE)
  x$layers <- x$layers[keep]
  x + kill
}
# CoveragePlot maps the coverage area to fill = group and then installs its own
# fill scale, so a scale_fill_manual() added here is silently replaced on the row
# that also carries the gene models. aes_params wins over the mapping, so the
# colour is set on the GeomArea layer directly.
setfill <- function(x, colour) {
  if (inherits(x, "patchwork")) {
    for (i in seq_along(x)) x[[i]] <- tryCatch(setfill(x[[i]], colour), error = function(e) x[[i]])
    return(x)
  }
  if (!inherits(x, "ggplot")) return(x)
  for (i in seq_along(x$layers))
    if (inherits(x$layers[[i]]$geom, "GeomArea")) x$layers[[i]]$aes_params$fill <- colour
  x
}
last <- tail(names(libs), 1)
cov <- list(); til <- list()
for (id in names(libs)) {
  o <- load_one(id, file.path(DATA, libs[[id]][1]), file.path(DATA, libs[[id]][2]))
  # the gene models ride on the last coverage row (annotation = TRUE)
  cov[[id]] <- setfill(clean(CoveragePlot(object = o, region = REGION, annotation = (id == last),
                                          peaks = FALSE, extend.upstream = 0,
                                          extend.downstream = 0)), col[id])
  til[[id]] <- clean(TilePlot(object = o, region = REGION, tile.cells = 300, order.by = "top") +
                     scale_fill_gradientn(colours = c("white", "blue", "darkblue"),
                                          limits = c(0, 2), breaks = c(0, 1, 2)))
}
til[[last]] <- til[[last]] + theme(legend.position = "right")   # one shared key

fig <- wrap_plots(c(cov, til), ncol = 1, heights = c(2, 2, 2.8, 2, 2, 2)) +
  plot_annotation(title = "chr13:91,340,000-91,370,000 (hg38, MIR17HG)",
                  theme = theme(plot.title = element_text(size = 6.5, hjust = 0,
                                                          margin = margin(0, 0, 2, 0))))
ggsave(file.path(OUT, "fig1c_atac_tracks_HEK293T.pdf"), fig, width = 9, height = 13, units = "cm")
ggsave(file.path(OUT, "fig1c_atac_tracks_HEK293T.png"), fig, width = 9, height = 13, units = "cm", dpi = 300)
cat("wrote", file.path(OUT, "fig1c_atac_tracks_HEK293T.pdf"), "\n")
\'\'\'
rfile = os.path.join(OUT, 'fig1c_atac_tracks.R')
with open(rfile, 'w') as fh:
    fh.write(R_TRACKS)
print('wrote', rfile)""")

code("""# this is the slow step: CountFragments/TilePlot stream the fragments.tsv.gz files
# (2.8 GB + 4.4 GB + 188 MB), so it costs tens of minutes and ~20 GB of RAM.
# The rendered figure is cached: delete it (or set FIG1_FORCE=1) to re-run.
import time
png = os.path.join(OUT, 'fig1c_atac_tracks_HEK293T.png')
if os.path.exists(png) and os.environ.get('FIG1_FORCE') != '1':
    print('already rendered, skipped:', png,
          '(%s, %.1f MB) - delete it or set FIG1_FORCE=1 to re-run'
          % (time.ctime(os.path.getmtime(png)), os.path.getsize(png) / 1e6))
else:
    t0 = time.time()
    rc = subprocess.run(['Rscript', rfile], env=dict(os.environ, FIG1_OUT=OUT, FIG1_DATA=HZJ))
    print('Rscript exit code %d - %.1f min' % (rc.returncode, (time.time() - t0) / 60))""")

code("""png = os.path.join(OUT, 'fig1c_atac_tracks_HEK293T.png')
if os.path.exists(png):
    show(png)
else:
    print('not rendered yet:', png)""")


nb.cells = C
path = os.path.join(OUT_DIR, 'Figure_1.ipynb')
nbf.write(nb, path)
print('wrote', path, len(C), 'cells')
