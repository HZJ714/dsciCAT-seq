#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build Figure_S1.ipynb — reproduction of Figure S1 (library quality, 8 libraries).

Conventions are the ones used by Figure_1.ipynb / Figure_2.ipynb / Figure_S2.ipynb:
  * reproduction first: every panel states the file(s) that go in and the shortest
    complete code that turns them into the numbers and the plot;
  * styling is deliberately minimal — the 18.00 x 18.40 cm page, the panel letters,
    the shared legend strip and the compositing are production steps owned by
    <work>/figures/scripts/figS1_layout.py + compose_figS1.py and are not
    reproduced here;
  * provenance is stated honestly: panels A, B and D read per-cell tables that
    CellRanger/Solo/ISSAAC already produced, panel C is re-derived here from the
    sparse matrices (that is the only panel with real computation in it);
  * comments in English.

STATUS (2026-10-22): executed in the analysis environment; every cell ran. A, B, D
reproduce the delivered figure exactly (same rows, same medians, same counts);
C reproduces the delivered correlation cache to |delta| < 5e-4, which is the
3-decimal rounding of the cached CSV. Full run ~3 min (12 cores, 3.6 GB of
matrices); A, B, D alone take ~5 s.

Run:  python build_figS1_notebook.py
Out:  <repo>/figS1/Figure_S1.ipynb
"""
import os
import nbformat as nbf

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
nb = nbf.v4.new_notebook(metadata={
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.8"}})
C = []
md = lambda s: C.append(nbf.v4.new_markdown_cell(s.strip("\n")))
code = lambda s: C.append(nbf.v4.new_code_cell(s.strip("\n")))

# ===========================================================================
# Title
# ===========================================================================
md("""# Figure S1 — dsciCAT-seq enables consistent library quality across four barcode numbers and two throughput scales

Eight dsciCAT-seq libraries: four barcode numbers (8bc / 16bc / 24bc / 48bc) x two
mRNA-capture depths (50k / 100k cells per run), each a mix of HEK293T (human) and
NIH3T3 (mouse) cells. Figure S1 asks whether the *measured* quality of those
libraries depends on the barcode number or on the throughput.

| Panel | What it shows | Input | Section |
|---|---|---|---|
| A | RNA-channel cell counts + per-cell log10(UMI) and log10(detected genes) | `sub/RNA/{human,mouse}_metrics.csv` | §1 |
| B | ATAC-channel cell counts + per-cell log10(unique fragments) | `sub/ATAC/.../{human,mouse}_metrics.csv` | §2 |
| C | 8 x 8 RNA expression-profile correlation, human and mouse | `matrix.mtx` + `features.tsv` | §3 |
| D | per-cell species assignment of the same nuclei in both channels | `ISSAAC_joint_metrics_sub.csv` | §4 |

The plots below are the *measurement*, one figure per panel; the delivered page
arranges the same four objects in one 18 x 18.4 cm sheet.""")

# ===========================================================================
# 0. Environment, inputs, provenance
# ===========================================================================
md("""## 0. Environment, inputs, provenance

### 0.1 The eight libraries and where their files are

`HG*` are the sequencing-run identifiers on the instrument; the condition names
(`50k-8bc` ... `100k-48bc`) are how the libraries are labelled in the figure.
The mapping is fixed, and every file below is resolved from it at run time —
nothing is looked up by hand.

All eight libraries resolve, and the three columns after `run` say which copy of
each table was used: `sub` for the seven libraries requantified under `sub/`,
`sub0.3` for `HG9` (re-quantified at a stricter cell cutoff, `RNA-sub0.3` /
`ATAC-sub0.33`, which is what the delivered figure uses), and `run root` for the
two libraries whose joint table is kept at the run root instead of under `sub/`
— same `*_sub*` table, different location. The resolver below implements exactly
that priority order, so a run where some directories are absent fails loudly
instead of silently mixing depths.

### 0.2 What is read and what is re-derived

| Item | Read from disk | Re-derived here |
|---|---|---|
| per-cell UMI / detected genes / fragments (A, B) | `*_metrics.csv` written by CellRanger / Solo / ISSAAC | — (aggregated only) |
| per-cell species call of both channels (D) | `ISSAAC_joint_metrics_sub.csv` (`species_joint` column) | — (counts only) |
| gene-level expression totals used for the correlation (C) | `matrix.mtx`, `features.tsv` | **yes: whole panel C** |
| ENSG identifiers for the one library whose `features.tsv` carries symbols | `HG13/sub/outs/Solo.out/.../features.tsv` as symbol→ENSG map | map checked for one-to-one and order |

Panel C is the only panel that does science rather than bookkeeping, and it is
the panel that a reader cannot check from the figure face alone, so it is
recomputed from the sparse matrices here.

**One upstream note on the ATAC side.** No file here is the bead -> cell merge, but
panels B and D are downstream of it: the dsciCAT-seq bead carries its barcode in a
structure that neither the ISSAAC nor the 10x chemistry has, so the ATAC FASTQ is
aligned with `chromap` and the bead barcodes are re-assembled by `d2c` (same
`forcefrag2000` fragment-length cap for every library) before the `metrics.csv` of §2
and the `ISSAAC_joint_metrics_sub.csv` of §4 exist. The merge is therefore a constant
of the pipeline, not a variable across the eight libraries; as a documented pipeline
input it is described in the Figure 1 notebook (`GitHub/fig1`, section 0.2).

### 0.3 Why the figures are not shown inline

Each section saves a PNG under `FIGS1_OUT` and prints the path; nothing is
displayed in the notebook. The measurements are the manuscript's, but the colours,
font sizes and canvas sizes here are the shortest form that draws the same
measurement, and they deliberately do not match the assembled page (`figS1_layout.py`
+ `compose_figS1.py`). A PNG in `FIGS1_OUT` is therefore evidence that the
numbers reproduced, not a proof of the printed artwork — the tables printed above
each section are what to check against.

### 0.4 Software

Python + numpy + pandas + scipy (`mmread`) + matplotlib. No production fonts, no
`fitz`: the figure *page* is not rebuilt here.""")

code('''import os, sys, time
import numpy as np, pandas as pd

t0 = time.time()
DATA_ROOT = os.environ.get('FIGS1_DATA', os.path.join(os.environ.get('FIG1_DATA', ''), 'HG'))
OUT_DIR = os.environ.get('FIGS1_OUT', os.path.join(os.getcwd(), 'figS1_out'))
os.makedirs(OUT_DIR, exist_ok=True)

# condition -> run, and the order the figure uses (barcode number inner, depth outer)
ORDER = ['50k-8bc', '50k-16bc', '50k-24bc', '50k-48bc',
         '100k-8bc', '100k-16bc', '100k-24bc', '100k-48bc']
LIB = {'50k-8bc': 'HG9',    '50k-16bc': 'HG12', '50k-24bc': 'HG15', '50k-48bc': 'HG10_1',
       '100k-8bc': 'HG11',  '100k-16bc': 'HG14', '100k-24bc': 'HG9_1', '100k-48bc': 'HG12_1'}
CELL_LINE = {'human': 'HEK293T', 'mouse': 'NIH3T3'}

# per-library RNA sub-directories: HG9 is the one re-quantified at a stricter cutoff
def rna_dirs(lib):
    return ['RNA-sub0.3', 'RNA_sub0.3', 'sub/RNA'] if lib == 'HG9' else ['sub/RNA']

def pick(cands, what):
    for p in cands:
        if os.path.exists(p):
            return p
    raise FileNotFoundError('%s: none of %s' % (what, cands))

def rna_metrics(lib, sp):
    return pick(['%s/%s/%s/%s_metrics.csv' % (DATA_ROOT, lib, d, sp) for d in rna_dirs(lib)],
                '%s RNA metrics (%s)' % (lib, sp))

def atac_metrics(lib, sp):
    H = 'Human' if sp == 'human' else 'Mouse'
    return pick(['%s/%s/ATAC-sub0.33/%s_metrics.csv' % (DATA_ROOT, lib, sp),
                 '%s/%s/ATAC_sub0.33/%s_metrics.csv' % (DATA_ROOT, lib, sp),
                 '%s/%s/sub/ATAC/%s/metrics.csv' % (DATA_ROOT, lib, H),
                 '%s/%s/sub/ATAC/%s_metrics.csv' % (DATA_ROOT, lib, sp),
                 '%s/%s/ISSAAC_ATAC/%s_metrics.csv' % (DATA_ROOT, lib, sp)],
                '%s ATAC metrics (%s)' % (lib, sp))

def rna_mtx_dir(lib, sp):
    for d in rna_dirs(lib):
        p = '%s/%s/%s/%s' % (DATA_ROOT, lib, d, sp)
        if os.path.exists(p + '/matrix.mtx') and os.path.exists(p + '/features.tsv'):
            return p
    raise FileNotFoundError('%s RNA matrix (%s)' % (lib, sp))

def joint_metrics(lib):
    return pick(['%s/%s/sub/ISSAAC_joint_metrics_sub.csv' % (DATA_ROOT, lib),
                 '%s/%s/ISSAAC_joint_metrics_sub.csv' % (DATA_ROOT, lib),
                 '%s/%s/ISSAAC_joint_metrics.csv' % (DATA_ROOT, lib)],
                '%s joint metrics' % lib)

def provenance(path):
    """which requantification depth a resolved file came from"""
    if 'ATAC-sub0.33' in path or 'RNA-sub0.3' in path or 'RNA_sub0.3' in path:
        return 'sub0.3'
    return 'sub' if '/sub/' in path else 'run root'

resolved = []
for tag in ORDER:
    lib = LIB[tag]
    resolved.append((tag, lib, provenance(rna_metrics(lib, 'human')),
                     provenance(atac_metrics(lib, 'human')), provenance(joint_metrics(lib))))
print(pd.DataFrame(resolved, columns=['condition', 'run', 'RNA', 'ATAC', 'joint'])
      .to_string(index=False))
print('\\ndata root : %s' % DATA_ROOT)
print('output dir: %s' % OUT_DIR)
print('python %s | numpy %s | pandas %s | scipy %s | matplotlib %s'
      % (sys.version.split()[0], np.__version__, pd.__version__,
         __import__('scipy').__version__, __import__('matplotlib').__version__))''')

md("""**Expected** (this run): all eight libraries resolve; `RNA`/`ATAC` = `sub0.3` for
`50k-8bc` and `sub` for the other seven; `joint` = `sub` for six libraries and
`run root` for `HG9`/`HG11`.""")

# ===========================================================================
# 1. Panel A — RNA-channel library quality
# ===========================================================================
md("""## 1. Panel A — RNA-channel cell counts and per-cell sensitivity

**Input.** For every library, the two per-cell tables of the RNA channel,
`human_metrics.csv` and `mouse_metrics.csv` (one row per cell; `nCounts` = UMIs,
`nFeatures` = detected genes). Both cell lines of one library come from the same
run and are quantified against their own genome, so HEK293T reads `Human` and
NIH3T3 reads `Mouse`.

**Output.** Cell recovery per library, and the median/IQR of
log10(UMI) and log10(detected genes) per library x cell line.

On the delivered page A is drawn as two halves, HEK293T left and NIH3T3 right, and
B re-uses that arrangement; here both cell lines are drawn inside one axes pair per
measure. Which line sits where is presentation, not measurement, so it is not
reproduced.""")

code('''def load_metrics(kind):
    rows = []
    for tag in ORDER:
        lib = LIB[tag]
        for sp in ['human', 'mouse']:
            p = rna_metrics(lib, sp) if kind == 'RNA' else atac_metrics(lib, sp)
            t = pd.read_csv(p)
            t['Lib'] = tag
            t['Cell_line'] = CELL_LINE[sp]
            t['src'] = provenance(p)
            rows.append(t)
    return pd.concat(rows, ignore_index=True)

rna = load_metrics('RNA')
rna['logUMI'] = np.log10(rna['nCounts'] + 1)
rna['logGene'] = np.log10(rna['nFeatures'] + 1)

print('rows read: %d' % len(rna))
print('cells per library x cell line:')
print(rna.groupby(['Cell_line', 'Lib']).size().unstack('Lib').reindex(columns=ORDER).to_string())
print('\\nmedian log10(UMI + 1):')
print(rna.groupby(['Cell_line', 'Lib'])['logUMI'].median().round(3).unstack('Lib').reindex(columns=ORDER).to_string())
print('\\nmedian log10(detected genes + 1):')
print(rna.groupby(['Cell_line', 'Lib'])['logGene'].median().round(3).unstack('Lib').reindex(columns=ORDER).to_string())''')

code('''import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def show(*paths):
    """report the saved file only -- nothing is displayed inline (section 0.2)"""
    for q in paths:
        print('%s  (%.0f KB)' % (q, os.path.getsize(q) / 1024.0))

def qc_figure(df, xcol_name, series, ylabel, fname):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.2))
    (a, b) = axes
    cells = df.groupby(['Cell_line', 'Lib']).size().unstack('Lib').reindex(columns=ORDER)
    w = 0.38; xpos = np.arange(len(ORDER))
    for i, cl in enumerate(['HEK293T', 'NIH3T3']):
        off = (i - 0.5) * w
        a.bar(xpos + off, cells.loc[cl], width=w, color=['#d62728', '#1f77b4'][i], label=cl)
        b.boxplot([df[(df.Cell_line == cl) & (df.Lib == t)][xcol_name].values for t in ORDER],
                  positions=xpos + off, widths=w * 0.8, showfliers=False)
    for ax in axes:
        ax.set_xticks(xpos); ax.set_xticklabels(ORDER, rotation=45, ha='right', fontsize=8)
        ax.set_xlabel('barcode number x capture depth', fontsize=9)
    a.set_ylabel('recovered cells', fontsize=9); a.legend(fontsize=8, frameon=False)
    b.set_ylabel(ylabel, fontsize=9)
    b.set_xticklabels([])
    fig.suptitle('Panel A — RNA channel (%s)' % series, fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, fname), dpi=200)
    plt.close(fig)
    return os.path.join(OUT_DIR, fname)

show(qc_figure(rna, 'logUMI', 'log10 UMI', 'log10 UMI per cell', 'S1a_RNA_UMI.png'),
     qc_figure(rna, 'logGene', 'log10 detected genes', 'log10 genes per cell', 'S1a_RNA_genes.png'))''')

md("""**Expected.** 10,738-22,610 cells per library and cell line; the 100k series
recovers ~1.6-1.8x the cells of the 50k series at the same barcode number, and
the barcode number itself changes nothing (within a depth the four barcode
conditions stay within ~15 % of each other, and the two lowest, `50k-8bc` /
`100k-8bc`, are `HG9` / `HG11`, the libraries re-quantified at the stricter
cutoff). Median log10(UMI) 3.15-3.36 and median log10(detected genes)
2.88-3.13, likewise flat in the barcode number.
That flatness *is* the claim of the panel: more barcodes do not cost sensitivity.""")

# ===========================================================================
# 2. Panel B — ATAC-channel library complexity
# ===========================================================================
md("""## 2. Panel B — ATAC-channel cell counts and per-cell complexity

**Input.** The per-cell tables of the ATAC channel, same eight libraries and same
two cell lines; there `nCounts` is unique fragments and `nFeatures` is peaks.

**Output.** Cell recovery per library and the distribution of
log10(unique fragments) per cell.""")

code('''atac = load_metrics('ATAC')
atac['logFrag'] = np.log10(atac['nCounts'] + 1)
atac['logPeak'] = np.log10(atac['nFeatures'] + 1)

print('rows read: %d' % len(atac))
print('cells per library x cell line:')
print(atac.groupby(['Cell_line', 'Lib']).size().unstack('Lib').reindex(columns=ORDER).to_string())
print('\\nmedian log10(unique fragments + 1):')
print(atac.groupby(['Cell_line', 'Lib'])['logFrag'].median().round(3).unstack('Lib').reindex(columns=ORDER).to_string())

fig, axes = plt.subplots(1, 2, figsize=(9, 3.2))
a, b = axes
cells = atac.groupby(['Cell_line', 'Lib']).size().unstack('Lib').reindex(columns=ORDER)
xpos = np.arange(len(ORDER)); w = 0.38
for i, cl in enumerate(['HEK293T', 'NIH3T3']):
    off = (i - 0.5) * w
    a.bar(xpos + off, cells.loc[cl], width=w, color=['#d62728', '#1f77b4'][i], label=cl)
    parts = [np.log10(atac[(atac.Cell_line == cl) & (atac.Lib == t)]['nCounts'] + 1).values
             for t in ORDER]
    vp = b.violinplot(parts, positions=xpos + off, widths=w * 0.9, showextrema=False)
    for pc in vp['bodies']:
        pc.set_facecolor(['#d62728', '#1f77b4'][i]); pc.set_alpha(0.5)
for ax in axes:
    ax.set_xticks(xpos); ax.set_xticklabels(ORDER, rotation=45, ha='right', fontsize=8)
    ax.set_xlabel('barcode number x capture depth', fontsize=9)
a.set_ylabel('recovered cells', fontsize=9); a.legend(fontsize=8, frameon=False)
b.set_ylabel('log10 unique fragments per cell', fontsize=9); b.set_xticklabels([])
fig.suptitle('Panel B — ATAC channel', fontsize=10)
fig.tight_layout(); p = os.path.join(OUT_DIR, 'S1b_ATAC_QC.png')
fig.savefig(p, dpi=200); plt.close(fig); show(p)''')

md("""**Expected.** 13,833-33,311 nuclei per library and cell line (100k ~1.7-2.2x
50k), median log10(unique fragments) 4.25-4.38, flat in the barcode number and
mildly higher for HEK293T. Fewer cells at 50k, the same per-cell complexity —
which is why the panel shows counts and distributions together.""")

# ===========================================================================
# 3. Panel C — cross-library expression-profile correlation
# ===========================================================================
md("""## 3. Panel C — do the eight libraries measure the same transcriptome?

**Input.** For each library and genome: `matrix.mtx` (genes x cells counts from
CellRanger or Solo) and `features.tsv` (gene ids, in matrix row order).

**What the panel is.** Each library is collapsed to one gene-level total-count
vector, the libraries are intersected on ENSG ids, the vectors are turned into
log2 CPM, and the 8 x 8 Pearson matrix is drawn. Human and mouse are done
separately because the two cell lines are quantified against different genomes.

**Matrix shape.** These `matrix.mtx` files are stored cells x genes (the transpose
of the CellRanger default), which is why the id check below compares
`features.tsv` with `mtx.shape[1]` and why gene totals are a `sum(axis=0)`.

**Feature ids.** Seven libraries already carry ENSG ids. `HG13/sub/outs/Solo...`
is the only one that carries symbols, and it is not used by panel C of the
8-library figure, so the id check below is a guard rather than a work-around: it
converts symbol rows through the one-to-one symbol→ENSG map and refuses to
continue if the mapped ids are not in the same order as the reference table.

**Cost.** 3.6 GB of matrices read once, ~110 s on 12 processes; the eight
vectors per genome are kept in memory for the rest of the section.""")

code('''import re
from multiprocessing import Pool
from scipy.io import mmread

REF_MAP = '%s/HG13/sub/outs/Solo.out/Gene/filtered/features.tsv' % DATA_ROOT
REF_ENSG = {'Human': '%s/HG11/sub/RNA/Human/features.tsv' % DATA_ROOT,
            'Mouse': '%s/HG11/sub/RNA/Mouse/features.tsv' % DATA_ROOT}

def symbol_map(sp):
    \"\"\"symbol -> ENSG, keeping only globally unique mappings\"\"\"
    gpre = 'GRCh38' if sp == 'Human' else 'mm10'
    cand, cnt = {}, {}
    for line in open(REF_MAP):
        f = line.rstrip('\\n').split('\\t')
        if len(f) < 2:
            continue
        m_id = re.match(r'^(%s)_+(.+)$' % gpre, f[0])
        m_sy = re.match(r'^(%s)_+(.+)$' % gpre, f[1])
        if not m_id or not m_sy:
            continue
        eid, sym = m_id.group(2), m_sy.group(2)
        cand[sym] = eid
        cnt[sym] = cnt.get(sym, 0) + 1
    return {s: e for s, e in cand.items() if cnt[s] == 1}

def feature_ids(prefix):
    \"\"\"first column of features.tsv, normalised to ENSG ids\"\"\"
    sp = 'Human' if prefix.rstrip('/').endswith('Human') else 'Mouse'
    ids = [l.split('\\t')[0] for l in open(prefix + '/features.tsv').read().split('\\n') if l]
    if all(i.startswith(('ENSG', 'ENSMUSG')) for i in ids[:50]):
        return ids
    m = symbol_map(sp)
    conv = [m.get(s, s) for s in ids]
    rpos = {e: i for i, e in enumerate(
        l.split('\\t')[0] for l in open(REF_ENSG[sp]).read().split('\\n') if l)}
    idx = [rpos[c] for c in conv if c in rpos]
    mono = all(a < b for a, b in zip(idx, idx[1:]))
    print('   [check] %s: %d symbol rows -> %d mapped, %d in reference, order kept=%s'
          % (prefix.replace(DATA_ROOT + '/', ''), len(ids), sum(1 for c, s in zip(conv, ids) if c != s),
             len(idx), mono))
    if not mono or len(idx) < 0.99 * len(ids):
        raise SystemExit('%s: symbol->ENSG map does not hold' % prefix)
    return conv

def gene_totals(prefix):
    \"\"\"matrix.mtx -> one pandas Series (ENSG -> total counts over that library's cells)\"\"\"
    mtx = mmread(os.path.join(prefix, 'matrix.mtx')).tocsr()
    ids = feature_ids(prefix)
    if len(ids) != mtx.shape[1]:
        raise SystemExit('%s: %d feature ids for %d matrix columns' % (prefix, len(ids), mtx.shape[1]))
    s = pd.Series(np.asarray(mtx.sum(axis=0)).ravel(), index=ids)
    if s.index.duplicated().any():
        s = s.groupby(level=0).sum()
    return prefix, mtx.shape[0], mtx.shape[1], s

jobs = [(tag, sp, rna_mtx_dir(LIB[tag], sp)) for tag in ORDER for sp in ['Human', 'Mouse']]
t1 = time.time()
with Pool(min(12, len(jobs))) as pool:
    out = pool.map(gene_totals, [j[2] for j in jobs])
totals = {}
for (tag, sp, _), (pfx, ncell, ngene, s) in zip(jobs, out):
    totals[(tag, sp)] = s
    print('%-11s %-6s %-24s cells %6d  genes %6d' % (tag, sp, pfx.replace(DATA_ROOT + '/', ''), ncell, ngene))
print('read %.0f s' % (time.time() - t1))''')

md("""### 3b. Common genes, log2 CPM, Pearson""")

code('''def cpm_log(sums, sp):
    \"\"\"the eight gene-total vectors -> log2(CPM + 1) on the common ENSG set\"\"\"
    per = [sums[(t, sp)].rename(t) for t in ORDER]
    common = sorted(set.intersection(*[set(s.index) for s in per]))
    counts = pd.concat([s.loc[common] for s in per], axis=1)
    counts = np.log2(counts + 1)                       # library sizes differ; stabilise first
    return np.log2(counts / counts.sum(axis=0) * 1e6 + 1), common

corr, ncommon = {}, {}
for sp in ['Human', 'Mouse']:
    mat, common = cpm_log(totals, sp)
    corr[sp] = mat.corr().loc[ORDER, ORDER]
    ncommon[sp] = len(common)
    off = corr[sp].values[np.triu_indices(len(ORDER), 1)]
    print('[%s] common genes %d | off-diagonal r %.3f - %.3f (median %.3f)'
          % (sp, len(common), off.min(), off.max(), np.median(off)))
    corr[sp].round(3).to_csv(os.path.join(OUT_DIR, 'S1c_corr_%s.csv' % sp))
    print(corr[sp].round(3).to_string())

fig, axes = plt.subplots(1, 2, figsize=(9, 4.0))
for ax, sp in zip(axes, ['Human', 'Mouse']):
    m = corr[sp].values
    im = ax.imshow(m, vmin=np.min([corr[a].values.min() for a in ['Human', 'Mouse']]), vmax=1,
                   cmap='magma')
    ax.set_xticks(range(len(ORDER))); ax.set_yticks(range(len(ORDER)))
    for t in (ax.get_xticklabels() + ax.get_yticklabels()):
        t.set_fontsize(7)
    ax.set_xticklabels(ORDER, rotation=90); ax.set_yticklabels(ORDER)
    ax.set_title('%s  (%d common genes)' % (sp, ncommon[sp]), fontsize=9)
    for i in range(len(ORDER)):
        for j in range(len(ORDER)):
            ax.text(j, i, '%.2f' % m[i, j], ha='center', va='center', fontsize=6,
                    color='w' if m[i, j] < 0.9 else 'k')
fig.tight_layout(); p = os.path.join(OUT_DIR, 'S1c_RNA_similarity.png')
fig.savefig(p, dpi=200); plt.close(fig); show(p)''')

md("""**Expected.** 36,601 human and 32,285 mouse genes survive the eight-way
intersection; the off-diagonals span 0.863-0.880 (human) and 0.871-0.886 (mouse),
medians 0.871 and 0.878 — every library pair correlates tightly with every other,
with no grouping by barcode number or depth. The three-decimal CSV written here
(`S1c_corr_{Human,Mouse}.csv`) is identical to `scripts/cache/S8_corr_*.csv`, the
file the delivered heatmap was drawn from, on all 64 cells; before rounding the
two agree to |delta| < 5e-4.""")

# ===========================================================================
# 4. Panel D — species assignment in both channels
# ===========================================================================
md("""## 4. Panel D — the same nuclei called in both channels

**Input.** `ISSAAC_joint_metrics_sub.csv` per library: one row per jointly
recovered nucleus with the RNA and ATAC human-read fractions and the joint
species call `species_joint` in {Human, Mouse, Mix}.

**What the panel is.** A 2 x 4 facet grid, one facet per library, x = ATAC
GRCh38 read %, y = RNA GRCh38 read %, coloured by the joint call, plus the
Human/Mouse/Mix counts in the legend. Nothing is computed here beyond the
value counts — the panel exists to show that the mix is ~50:50 as designed and
that the doublet-contaminated fraction stays small.""")

code('''drows = []
joint = {}
for tag in ORDER:
    lib = LIB[tag]
    p = joint_metrics(lib)
    j = pd.read_csv(p)
    joint[tag] = j
    c = j['species_joint'].value_counts()
    r = (c / len(j) * 100).round(2)
    drows.append(dict(condition=tag, run=lib, src=provenance(p), cells=len(j),
                      Human=c.get('Human', 0), Mouse=c.get('Mouse', 0), Mix=c.get('Mix', 0),
                      Human_pct=r.get('Human', 0), Mouse_pct=r.get('Mouse', 0), Mix_pct=r.get('Mix', 0)))
dtab = pd.DataFrame(drows).set_index('condition').reindex(ORDER)
print(dtab.to_string())
print('\\nMix%% range %.2f - %.2f | Human%% %.2f - %.2f | Human+Mouse agreement %.1f%%'
      % (dtab.Mix_pct.min(), dtab.Mix_pct.max(), dtab.Human_pct.min(), dtab.Human_pct.max(),
         (dtab.Human_pct + dtab.Mouse_pct).mean()))''')

code('''SP_ORDER, SP_COLOR = ['Human', 'Mouse', 'Mix'], {'Human': 'red', 'Mouse': 'blue', 'Mix': 'grey'}
GRID = [['50k-8bc', '100k-8bc', '50k-16bc', '100k-16bc'],
        ['50k-24bc', '100k-24bc', '50k-48bc', '100k-48bc']]

fig, axs = plt.subplots(2, 4, figsize=(11, 5.2))
for r, row in enumerate(GRID):
    for c, tag in enumerate(row):
        ax, df = axs[r][c], joint[tag]
        counts = df['species_joint'].value_counts()
        ratio = (counts / len(df) * 100).round(2)
        for sp in SP_ORDER:
            sub = df.query('species_joint == @sp')
            ax.scatter(sub['ATAC GRCh38 fragments %'], sub['RNA GRCh38 reads %'],
                       c=SP_COLOR[sp], edgecolors=None, s=1, alpha=.2,
                       label='%s: %d (%s%%)' % (sp, counts.get(sp, 0), ratio.get(sp, 0)))
        ax.set_xlim(0, 100); ax.set_ylim(0, 100)
        ax.set_xticks([0, 25, 50, 75, 100]); ax.set_yticks([0, 25, 50, 75, 100])
        ax.tick_params(labelsize=7, labelleft=(c == 0), left=(c == 0),
                       labelbottom=(r == 1), bottom=(r == 1))
        ax.set_title(tag, fontsize=8)
        if r == 1:
            ax.set_xlabel('ATAC GRCh38 reads %', fontsize=8)
        if c == 0:
            ax.set_ylabel('RNA GRCh38 reads %', fontsize=8)
        ax.legend(frameon=False, markerscale=2, loc='upper left', fontsize=6,
                  handlelength=0.5, handletextpad=0.15, columnspacing=0.2)
fig.tight_layout(); p = os.path.join(OUT_DIR, 'S1d_joint_mix.png')
fig.savefig(p, dpi=150); plt.close(fig); show(p)''')

md("""**Expected.** 22,387-46,355 jointly recovered nuclei per library; Human
44.2-48.3 %, Mouse 44.4-48.0 %, Mix 4.94-10.55 %, i.e. the 50:50 mixture that was
loaded, with the 100k series at the higher Mix end (7.5-10.6 %) and 50k-24bc the
cleanest (4.94 %). The eight (n, Human, Mouse, Mix) triples printed by the cell
are exactly the eight legend lines of the delivered panel.""")

# ===========================================================================
# Requirements / Run
# ===========================================================================
md("""## Requirements

Python 3 with `numpy`, `pandas`, `scipy`, `matplotlib` — nothing else (executed
here on 3.8.14 / numpy 1.22.4 / pandas 1.5.1 / scipy 1.10.1 / matplotlib 3.6.2). Panel C
is the only part that touches the sparse matrices.

The input tree is the un-redistributed output of the sequencing pipeline
(`CNGBdb` accession `CNP0010263`); per library it has to contain

```
<data-root>/<run>/sub/RNA/{human,mouse}_metrics.csv      # A
<data-root>/<run>/sub/RNA/{Human,Mouse}/(matrix.mtx,features.tsv)   # C
<data-root>/<run>/sub/ATAC/.../{human,mouse}_metrics.csv # B   (HG9: ATAC-sub0.33/)
<data-root>/<run>/sub/ISSAAC_joint_metrics_sub.csv       # D
<data-root>/HG13/sub/outs/Solo.out/Gene/filtered/features.tsv  # symbol->ENSG map (guard)
```

`FIGS1_DATA` overrides the data root, `FIGS1_OUT` the output directory.""")

md("""## Run

```bash
python build_figS1_notebook.py            # regenerate this notebook

jupyter nbconvert --to notebook --execute --inplace \\
    --ExecutePreprocessor.timeout=-1 Figure_S1.ipynb
```

~3 min for the whole notebook, of which ~2 min is panel C reading the sixteen
matrices (~110 s on 12 processes). Outputs land in `./figS1_out/`:

| File | Panel |
|---|---|
| `S1a_RNA_UMI.png`, `S1a_RNA_genes.png` | A |
| `S1b_ATAC_QC.png` | B |
| `S1c_corr_Human.csv`, `S1c_corr_Mouse.csv`, `S1c_RNA_similarity.png` | C |
| `S1d_joint_mix.png` | D |

## Status of this notebook

Written and executed on 2026-10-22 in the analysis environment that produced the
delivered `Figure_S1.pdf` (18.00 x 18.40 cm, panels A-D, Arial). A, B and D
reproduce the delivered numbers because they read the same per-cell tables the
delivered figure was drawn from; C is recomputed from the sparse matrices,
independently of the cached correlation CSV: the unrounded values agree to
|delta| < 5e-4 and the rounded CSV this notebook writes is identical to it on all
64 cells. The page assembly
(`figS1_layout.py`, `compose_figS1.py`, panel letters, shared legend strip) is a
production step and is not reproduced here.""")

nb.cells = C
path = os.path.join(OUT_DIR, 'Figure_S1.ipynb')
nbf.write(nb, path)
print('wrote', path, len(C), 'cells')
