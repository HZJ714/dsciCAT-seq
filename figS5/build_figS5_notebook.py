#!/usr/bin/env python3
"""
Build Figure_S5.ipynb (SCENIC+ support panels A-D of the delivered supplement).

Every code cell here is taken from the analysis that produced the delivered
figure (scripts listed in section 0), not rewritten. Editing this file and
running `python3 build_figS5_notebook.py` regenerates the notebook; running the
notebook executes the analysis and writes CSV/PDF into $FIGS5_OUT.

Layout/compositing (cm-exact placement, Arial, panel letters, shared legend) is
deliberately NOT here - see figS5/assembly/compose_figS5_noc3.py (now vendored into the repo).
"""
import os
import platform
import nbformat

MD, CODE = 'markdown', 'code'
CELLS = []
A = lambda kind, src: CELLS.append((kind, src.strip('\n')))

# --------------------------------------------------------------------------- 0
A(MD, r"""
# Figure S5 - SCENIC+ support: eRSS ranking, state networks, the Mef2c -> Nrxn3 locus, and downsampling stability

Companion to `Figure_5.ipynb`; same object, same selection rules, four support
panels. Delivered face: `Figure_S5.pdf` (18.00 x 46.75 cm, blocks A-D).

| Panel | What it shows | Reproduced here |
|---|---|---|
| A | eRegulon specificity score (eRSS) ranking, 6 oligodendrocyte-lineage states, top 10 labelled | yes - recomputed from the object and compared with the archived table the delivered panel was drawn from |
| B | TF-target networks of MOL2, MOL3_6 and MOL6 (Fig 5B shows MOL3) | yes - same `create_nx_tables` -> `create_nx_graph` -> `plot_networkx` path |
| C | *Nrxn3* locus: Mef2c-eRegulons at four enhancers + MEF2 motif hits + *Nrxn3* expression per state | expression and motif scan recomputed; track signal read from the decoded ArchR cache (see 0.1) |
| D | Downsampling stability of the eRegulon set | table read and re-plotted; the sweep itself is a full SCENIC+ re-run per depth (see 4.1) |

**Not reproduced here**: the cm-exact layout, Arial substitution, panel letters
and the compose step (`compose_figS5_noc3.py`), and the `retext_scenic2.py`
text-layer rebuild that the delivered PDF went through.
""")

A(MD, r"""
## 0.1 Inputs

Three environment variables drive every path below.

| Variable | Default | Must contain |
|---|---|---|
| `FIGS5_DATA` | `data` | `scplus_obj_final.pkl` (panel A, B, C-expression) - `NIHMS70326-supplement-Table_S1.xlsx` (panel B marker lists) - `nrxn3_expr.npz`, `nrxn3_tracks.npz` (panel C) - `eregulon_stability_table.csv` (panel D) - `nrxn3_mef2c_enhancers.fa` (panel C motif scan) |
| `FIGS5_OUT` | `figS5_out` | created here; every CSV/PDF product lands here |
| `FIGS5_ARCHIVE` | `<work>/pi/plots-scenic+` | the tables the delivered panels were drawn from; used only for the comparison prints |
| `FIGS5_MOTIF` | `<work>/aertslab_motif_colleciton/v10nr_clust_public` | `singletons/metacluster_179.4.cb` (MEF2 = JASPAR MA0052.4) |
| `FIGS5_GENOME` | `<work>/referrence/mm10/mm10.fa` | indexed mm10 genome - the file the delivered scan read through `pysam`; if absent, the scan falls back to the pre-cut `nrxn3_mef2c_enhancers.fa` |

Where each input comes from, as far as it is recorded on disk:

| Input | Origin |
|---|---|
| `scplus_obj_final.pkl` | the chain in 0.2 below |
| `nrxn3_expr.npz` | `Nrxn3` column of `X_EXP`, grouped by `ACC_celltype_l2` - recomputed in 3.1 and compared with the file |
| `nrxn3_tracks.npz` | ArchR **GroupCoverages** of the lineage subset (`FY/oligo/3kdsci/archr-subset/GroupCoverages/cell_type`), FY13+FY14 merged per state, RPM-normalised, 200-bp bins over chr12:88.55-89.60 Mb, decoded once by `grn_mef2c_nrxn3_track.py` (that script reads the `Coverage/chr12:88550000-89600000/Values` group and does the window arithmetic) |
| `nrxn3_mef2c_enhancers.fa` | the four Mef2c-eRegulon regions inside that window (sequence from the mm10 genome) |
| `eregulon_stability_table.csv` | `stability_downsample.py`, seed 123, one full SCENIC+ re-run per depth |
""")

A(MD, r"""
## 0.2 Where `scplus_obj_final.pkl` comes from

**Not repeated here.** The full upstream chain - Seurat/Signac RDS -> R-side export of the
binarised peak matrix -> cisTopic LDA at **60 topics** (`run_cgs_models`, `n_iter=150`,
`random_state=555`) -> three binarisations (otsu / top-3k / li) + per-state DARs -> pyCisTarget
enrichment against `v10nr_clust_public` (`run_pycistarget`, `run_without_promoters=True`) ->
`create_SCENICPLUS_object` -> `run_scenicplus` (which alone yields only 10 eRegulons) ->
`build_grn` + `make_rankings` + `score_eRegulons` + pseudobulk correlation and standard
filtering - is documented in **`GitHub/fig5/Figure_5.ipynb`, section 0**, together with the
parameters and the dated files each step left behind. Panel A and B below open the same object
that Fig 5 opens, so one provenance statement covers both.

What matters locally is only this: panels A and B read `uns['eRegulon_AUC_filtered']`
(96 filtered region-based eRegulons) and `uns['eRegulon_metadata_filtered']`, and both the
eRSS recomputation and the networks are asserted against the archived tables the delivered
panels were drawn from.

## 0.3 Why the figures are not shown inline

`matplotlib.use('Agg')` and every drawing cell only prints the path and the size
of the file it wrote. The delivered supplement went through a text-layer rebuild
and a cm-exact compose, so a screen PNG from this notebook is **not** the
manuscript figure: showing it inline would invite a comparison of fonts and sizes
that this notebook does not control. Compare numbers, not pixels.
""")

A(CODE, r"""
import os, json, subprocess, platform
import numpy as np
import pysam
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA    = os.environ.get('FIGS5_DATA', 'data')
OUT     = os.environ.get('FIGS5_OUT', 'figS5_out')
ARCHIVE = os.environ.get('FIGS5_ARCHIVE', os.environ.get('DSCI_SCENIC_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pi/plots-scenic+')))
MOTIF   = os.environ.get('FIGS5_MOTIF', os.path.join(os.environ.get('DSCI_MOTIF_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'aertslab_motif_colleciton')), 'v10nr_clust_public'))
GENOME  = os.environ.get('FIGS5_GENOME', os.path.join(os.environ.get('DSCI_REF_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'referrence')), 'mm10/mm10.fa'))
os.makedirs(OUT, exist_ok=True)

GROUP   = 'ACC_celltype_l2'
STATES  = ['COP_NFOL', 'MFOL', 'MOL2', 'MOL3', 'MOL3_6', 'MOL6']
AUC_KEY, SIG_KEY = 'eRegulon_AUC_filtered', 'Region_based'
SEED    = 555

def show(path):
    print('%s (%.0f KB)' % (path, os.path.getsize(path) / 1024.0))

import dill
obj = dill.load(open(os.path.join(DATA, 'scplus_obj_final.pkl'), 'rb'))
assert list(obj.cell_names) == list(obj.metadata_cell.index), 'cell order mismatch'
# Snapshot the state labels AND the Nrxn3 expression column now, straight after the load. The
# network step in section 2 calls the differential-gene / differential-region accessors, which
# mutate the object - metadata order and the expression view - and panel C has to be grouped
# exactly the way the delivered figure grouped it. Section 3.1 checks these copies against the
# cached arrays, which is where that claim gets settled.
OBS0 = obj.metadata_cell.loc[list(obj.cell_names), GROUP].astype(str).values
_gi = np.flatnonzero(np.asarray(obj.gene_names) == 'Nrxn3')
assert len(_gi) == 1, 'Nrxn3 maps to %d columns' % len(_gi)
_col = obj.X_EXP[:, _gi[0]]
LOGV0 = np.asarray(_col.todense()).ravel() if hasattr(_col, 'todense') else np.asarray(_col).ravel()
print('cells %d | genes %d | states %s' % (
    obj.X_EXP.shape[0], len(obj.gene_names), sorted(set(OBS0))))
""")

# --------------------------------------------------------------------------- 1
A(MD, r"""
---

## 1. Panel A - eRegulon specificity (eRSS) across the six states

**In** `scplus_obj_final.pkl` **Params** `auc_key='eRegulon_AUC_filtered'`,
`signature_keys=['Region_based']` (96 filtered eRegulons),
`variable='ACC_celltype_l2'`, `out_key_suffix='_cmp'`, top 10 labelled
**Out** `S5a_rss_scores_96.csv` (6 x 96), `S5a_rss_top10.csv`,
`figS5a_rss_grid.pdf`

`RSS` is dropped first: the object already carries an `uns['RSS']` from the
upstream run, and the delivered panel was drawn from a freshly computed one
(key `ACC_celltype_l2_cmp`). Recomputing is what makes the comparison with the
archived table meaningful.
""")

A(CODE, r"""
# ---- 1.1 eRSS, recomputed ----
from scenicplus.RSS import regulon_specificity_scores, plot_rss

assert list(obj.uns[AUC_KEY][SIG_KEY].columns).__len__() == 96, 'expected 96 filtered region-based eRegulons'
if 'RSS' in obj.uns:
    del obj.uns['RSS']
regulon_specificity_scores(obj, variable=GROUP, auc_key=AUC_KEY, signature_keys=[SIG_KEY],
                           selected_regulons=list(obj.uns[AUC_KEY][SIG_KEY].columns),
                           out_key_suffix='_cmp')
rss = obj.uns['RSS'][GROUP + '_cmp']
rss.to_csv(os.path.join(OUT, 'S5a_rss_scores_96.csv'))
top10 = pd.concat({ct: rss.loc[ct].nlargest(10) for ct in rss.index}, axis=0)
top10.index.names = ['celltype', 'eRegulon']
top10.to_csv(os.path.join(OUT, 'S5a_rss_top10.csv'))

print('RSS shape %s' % (rss.shape,))
for ct in STATES:
    n, v = rss.loc[ct].nlargest(1).index[0], rss.loc[ct].max()
    print('  %-9s top1 %-34s %.3f | RSS>0.5: %d' % (ct, n, v, int((rss.loc[ct] > 0.5).sum())))

arch = os.path.join(ARCHIVE, 'rss_scores_filtered_96.csv')
if os.path.exists(arch):
    a = pd.read_csv(arch, index_col=0).loc[rss.index, rss.columns]
    print('vs archived table the delivered panel used: max|delta| = %.2e' %
          np.abs(a.values - rss.values).max())
""")

A(MD, r"""
**Expected** - `RSS shape (6, 96)`; the six top-1 eRegulons and their scores are
the numbers the delivered panel labels:

| state | top 1 | eRSS |
|---|---|---|
| COP_NFOL | `Glis1_extended_+_(6r)` | **0.621** (by far the strongest) |
| MFOL | `Pax5_+_(9r)` | 0.343 |
| MOL2 | `Usf1_+_(11r)` | 0.378 |
| MOL3 | `Nr1h3_extended_+_(11r)` | 0.356 |
| MOL3_6 | `Creb5_extended_+_(172r)` | 0.331 |
| MOL6 | `Creb5_extended_+_(172r)` | 0.513 |

`max|delta|` against the archived table is at float-print precision (the
delivered table is the same computation, so it is expected to be ~1e-15, not a
rounded match).
""")

A(CODE, r"""
# ---- 1.2 the panel ----
f = os.path.join(OUT, 'figS5a_rss_grid.pdf')
plot_rss(obj, rss_key=GROUP + '_cmp', top_n=10, num_columns=2,
         figsize=(10, 12), fontsize=9, save=f)
show(f)
""")

# --------------------------------------------------------------------------- 2
A(MD, r"""
---

## 2. Panel B - TF-target networks of MOL2, MOL3_6 and MOL6

**In** `scplus_obj_final.pkl`, `NIHMS70326-supplement-Table_S1.xlsx` (state markers) **Params**
`eRegulon_metadata_key='eRegulon_metadata_filtered'`, `use_edge_tables=['TF2G']`, adaptive
`min_deg` 2-6 capped at 60 genes, `kamada_kawai_layout`, `scale_position_by=450`
**Out** `figS5b_grn_<tag>.pdf`, `S5b_grn_<tag>_nodes.csv`, `S5b_grn_<tag>_edges.csv`

The same generator as Fig 5B (MOL3), run on the other three states, so that the main-figure
network is not read as a singleton. Marker genes come from Table S1 of Marques et al. (2016);
`min_deg` is not fixed per state but climbs 2->6 until the gene set fits under 60, which is why
the printed `min_deg` is part of the expected output.
""")

A(CODE, r"""
# ---- 2.1 one pass per state (identical to Fig 5 2c, three states) ----
import re
from collections import Counter
from scenicplus.networks import create_nx_tables, create_nx_graph, plot_networkx

markers = {}
s1 = pd.ExcelFile(os.path.join(DATA, 'NIHMS70326-supplement-Table_S1.xlsx')).parse('Sheet1')
for c in range(1, s1.shape[1]):
    nm = s1.iloc[1, c]
    if isinstance(nm, str):
        markers[nm.strip()] = {g for g in s1.iloc[2:, c].dropna().astype(str).str.strip() if g}
print('Table S1 states:', len(markers))

FIGURES = [
    dict(tag='mol2',   var=GROUP + '_Log2FC_MOL2',
         tfs=['Nfix', 'Tcf4_extended', 'Klf7', 'Tcf12', 'Mbnl2_extended'],
         markers=markers['MOL2'] | {'Apod', 'Anln'}),
    dict(tag='mol3_6', var=GROUP + '_Log2FC_MOL3_6',
         tfs=['Sox6', 'Mllt10_extended', 'Creb5_extended'],
         markers=markers['MOL2/3'] | markers['MOL5/6']),
    dict(tag='mol6',   var=GROUP + '_Log2FC_MOL6',
         tfs=['Creb5_extended', 'Mbnl2_extended'],
         markers=markers['MOL6'] | markers['MOL5/6']),
]

meta_f = obj.uns['eRegulon_metadata_filtered']
for f in FIGURES:
    sub = meta_f[meta_f['Region_signature_name'].str.contains(
        '|'.join(re.escape(t) + '_' for t in f['tfs']))]
    linked = f['markers'] & set(sub['Gene'])
    deg = Counter(sub['Gene'])
    for min_deg in (2, 3, 4, 5, 6):
        hubs = {g for g, n in deg.items() if n >= min_deg}
        if len(linked | hubs) <= 60:
            break
    genes = sorted(linked | hubs)
    regions = sorted(set(sub[sub['Gene'].isin(genes)]['Region']))

    nt = create_nx_tables(scplus_obj=obj, eRegulon_metadata_key='eRegulon_metadata_filtered',
                          subset_eRegulons=f['tfs'], subset_regions=regions, subset_genes=genes,
                          add_differential_gene_expression=True,
                          add_differential_region_accessibility=True,
                          differential_variable=[GROUP])
    tfs = nt['Node']['TF']['TF'].dropna().unique().tolist()
    cmap = {tf: plt.get_cmap('tab10')(i % 10) for i, tf in enumerate(tfs)}
    G, pos, edge_tab, node_tab = create_nx_graph(
        nt, use_edge_tables=['TF2G'],
        color_edge_by={'TF2G': {'variable': 'TF', 'category_color': cmap}},
        transparency_edge_by={'TF2G': {'variable': 'TF2G_importance', 'min_alpha': 0.2, 'v_min': 0}},
        width_edge_by={'TF2G': {'variable': 'TF2G_importance', 'max_size': 4, 'min_size': 1.0}},
        color_node_by={'TF': {'variable': 'TF', 'category_color': cmap},
                       'Gene': {'variable': f['var'], 'continuous_color': 'bwr'}},
        transparency_node_by={'Gene': {'variable': f['var'], 'min_alpha': 0.2}},
        size_node_by={'TF': {'variable': 'fixed_size', 'fixed_size': 150},
                      'Gene': {'variable': 'fixed_size', 'fixed_size': 70}},
        shape_node_by={'TF': {'variable': 'fixed_shape', 'fixed_shape': 'circle'},
                       'Gene': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'}},
        label_size_by={'TF': {'variable': 'fixed_label_size', 'fixed_label_size': 22.0},
                       'Gene': {'variable': 'fixed_label_size', 'fixed_label_size': 11.0}},
        layout='kamada_kawai_layout', scale_position_by=450)
    plt.subplots(figsize=(12, 12)); plot_networkx(G, pos)
    out = os.path.join(OUT, 'figS5b_grn_%s.pdf' % f['tag'])
    plt.savefig(out, bbox_inches='tight'); plt.close()
    # 机器可读补表（Codex #91：“完整网络另作机器可读补表”）。
    # ⚠ 旧版这里手写两列 CSV，两个坑：(a) `G.nodes[n].get('TF')` 永远是 False —— 图上的节点属性
    #   不叫 'TF'（叫 label/color/size/font…），所以导出的 TF 列 58 行全 False；
    #   (b) `G.edges` 只有顶点对，丢了 TF2G_importance / 方向，而图里的边宽正是这个值。
    #   现在直接落 `create_nx_graph` 返回的两张原表（就是画图用的那两张，逐列自解释）。
    edge_tab.to_csv(os.path.join(OUT, 'S5b_grn_%s_edges.csv' % f['tag']), index=False)
    node_tab.to_csv(os.path.join(OUT, 'S5b_grn_%s_nodes.csv' % f['tag']), index=False)
    degmap = dict(G.degree)
    top = sorted(tfs, key=lambda t: -degmap.get(t, 0))[:3]
    print('%-7s min_deg=%d genes=%d -> %d nodes / %d edges | TF degrees: %s' % (
        f['tag'], min_deg, len(genes), G.number_of_nodes(), G.number_of_edges(),
        ', '.join('%s %d' % (t, degmap.get(t, 0)) for t in top)))
    show(out)
""")

A(MD, r"""
**Expected** - node/edge counts and hub degrees as quoted on the delivered face:

| state | nodes / edges | hubs (TF, degree) |
|---|---|---|
| MOL2 | 58 / 97 | Klf7 **34**, Mbnl2 28, Tcf4 24, then Nfix / Tcf12 |
| MOL3_6 | 49 / 65 | Mllt10 **36**, Creb5 33, Sox6 |
| MOL6 | 56 / 74 | Mbnl2 **37**, Creb5 28 |

(MOL3 itself, 44 nodes / 44 edges, is Fig 5B and is not redrawn here.) Changing `min_deg`,
the TF lists or the 60-gene cap changes these numbers; the cap in particular is what keeps the
drawn network readable.
""")

# --------------------------------------------------------------------------- 3
A(MD, r"""
---

## 3. Panel C - the Mef2c -> Nrxn3 locus

**In** `scplus_obj_final.pkl` (expression), `nrxn3_tracks.npz` (ArchR coverage),
`nrxn3_mef2c_enhancers.fa` + `metacluster_179.4.cb` (motif scan)
**Out** `S5c_nrxn3_expression_by_state.csv`, `S5c_mef2_motif_hits.csv`, `figS5c_nrxn3_locus.pdf`

Three claims are tested separately: (i) the four Mef2c-eRegulon regions in the
window, (ii) that MEF2 motifs sit in them, (iii) that *Nrxn3* is expressed higher
in MOL3 - so that enhancer opening and target expression are concordant and
state-specific.
""")

A(CODE, r"""
# ---- 3.1 Nrxn3 expression per state, recomputed from the object ----
from scipy.stats import mannwhitneyu

GENE = 'Nrxn3'
TSS = 88722433                                   # Nrxn3 5' TSS, chr12 (mm10)
ENH = [('E1', 88707581, 88708505, 3.85),
       ('E2', 88952720, 88953675, 3.96),
       ('E3', 89284141, 89284970, 4.31),
       ('E4', 89376068, 89377735, 6.86)]

logv = LOGV0                                              # snapshot taken in section 0
vals = {s: 2.0 ** logv[OBS0 == s] - 1.0 for s in STATES}   # linearised 2^x - 1, as on the face

rows = []
for s in STATES:
    med = float(np.median(vals[s])) if len(vals[s]) else float('nan')
    rows.append([s, len(vals[s]), med, float(np.mean(vals[s]))])
    print('%-9s n=%5d median %.2f' % (s, len(vals[s]), med))
pd.DataFrame(rows, columns=['state', 'n_cells', 'median', 'mean']).to_csv(
    os.path.join(OUT, 'S5c_nrxn3_expression_by_state.csv'), index=False)

others = np.concatenate([vals[s] for s in STATES if s != 'MOL3'])
u, p = mannwhitneyu(vals['MOL3'], others, alternative='greater')
print('MOL3 vs all other states: one-sided MWU p = %.2e' % p)

cachef = os.path.join(DATA, 'nrxn3_expr.npz')
if os.path.exists(cachef):
    z = np.load(cachef, allow_pickle=True)          # the cached arrays are on the log scale
    print('vs the cached array the delivered panel used (log2(x+1) scale):')
    for s in STATES:
        n_ok = len(z[s]) == int((OBS0 == s).sum())
        d = np.abs(np.sort(z[s]) - np.sort(logv[OBS0 == s])).max() if n_ok else float('nan')
        print('  %-9s n %5d/%5d %s max|delta| %.2e' % (
            s, int((OBS0 == s).sum()), len(z[s]), 'ok' if n_ok else 'MISMATCH', d))
""")

A(MD, r"""
**Expected** - cell counts 203 / 490 / 1,016 / 688 / 1,006 / 2,574 (= 5,977), *Nrxn3* median about
**7.0** in MOL3 (6.98 measured here) against **3.0-3.8** in the other five states (roughly 2x),
one-sided Mann-Whitney p < **1e-22**, and `max|delta| = 0` per state against the cached arrays.
Those arrays are on the **log2(x+1) scale the object stores** (MOL3 median 3.00), so the check runs
on that scale and only the plot is linearised - comparing a linearised vector against a log-scale
cache is what made an earlier pass look off by 8.75.
""")

A(CODE, r"""
# ---- 3.2 MEF2 motif hits in the four enhancers (PWM scan, as scan_mef2_motif.py) ----
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def parse_pwm(path, motif='jaspar__MA0052.4'):
    # L x 4 count matrix (A,C,G,T order) for the first block named 'motif' in a .cb file
    with open(path) as f:
        lines = [l.rstrip('\n') for l in f if l.strip()]
    i = 0
    while i < len(lines):
        if lines[i].startswith('>'):
            name, i, mat = lines[i][1:], i + 1, []
            while i < len(lines) and not lines[i].startswith('>'):
                mat.append([float(x) for x in lines[i].split()])
                i += 1
            if name == motif:
                return np.array(mat)
        else:
            i += 1
    raise ValueError('motif %s not found in %s' % (motif, path))

def pssm_from_counts(counts, pseudocount=1.0):
    p = (counts + pseudocount) / (counts.sum(axis=1, keepdims=True) + 4 * pseudocount)
    return np.log2(p / 0.25)                          # log2 odds, L x 4

def score_seq(seq, odds, rc=False):
    idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    s2 = ''.join(COMP[b] for b in reversed(seq)) if rc else seq
    return float(sum(odds[k, idx[s2[k]]] for k in range(len(s2))))

def read_fa(path):
    seqs, name, buf = {}, None, []
    for l in open(path):
        if l.startswith('>'):
            if name:
                seqs[name] = ''.join(buf)
            name, buf = l[1:].strip(), []
        else:
            buf.append(l.strip())
    if name:
        seqs[name] = ''.join(buf)
    return seqs

cbb  = os.path.join(MOTIF, 'singletons', 'metacluster_179.4.cb')
faf  = os.path.join(DATA, 'nrxn3_mef2c_enhancers.fa')
CUTOFF, TOP_N = 7.0, 5                                # scan_mef2_motif.py: score_cutoff / top_n
rows = []
cuts = {n.split('_')[0]: sq for n, sq in read_fa(faf).items()}   # fa names look like 'E1_chr12:88707581-88708505_...'
if os.path.exists(cbb):
    odds = pssm_from_counts(parse_pwm(cbb))
    L = odds.shape[0]
    print('PWM %d bp | max possible score %.2f' % (L, sum(odds[k].max() for k in range(L))))
    genome = pysam.FastaFile(GENOME) if os.path.exists(GENOME) else None
    for lab, st, en, nes in ENH:
        assert lab in cuts, 'no pre-cut record for %s in %s' % (lab, faf)
        seq = genome.fetch('chr12', st - 1, en).upper() if genome else cuts[lab].upper()
        if genome:
            assert seq == cuts[lab].upper(), '%s: mm10.fa window != pre-cut .fa record' % lab
        hits = []
        for k in range(0, len(seq) - L + 1):
            kmer = seq[k:k + L]
            if 'N' in kmer:
                continue
            for rc, strand in ((False, '+'), (True, '-')):
                sc = score_seq(kmer, odds, rc)
                hit = kmer if not rc else ''.join(COMP[b] for b in reversed(kmer))
                hits.append((sc, k + 1, st + k, strand, hit))
        hits.sort(reverse=True)
        top = [h for h in hits if h[0] >= CUTOFF][:TOP_N]
        for sc, rel, gpos, strand, hit in top:
            rows.append([lab, gpos, gpos + L - 1, strand, hit, round(sc, 2), rel])
        print('  %-3s chr12:%d-%d (%d bp, %s): %d hits, best %.2f' % (
            lab, st, en, len(seq), 'mm10.fa' if genome else 'pre-cut fa', len(top),
            top[0][0] if top else float('nan')))
else:
    print('motif singleton missing: %s' % cbb)

h = pd.DataFrame(rows, columns=['enhancer', 'start', 'end', 'strand', 'hit_sequence', 'score',
                               'rel_pos_in_enhancer'])
h.to_csv(os.path.join(OUT, 'S5c_mef2_motif_hits.csv'), index=False)
print('MEF2 (metacluster_179.4 = JASPAR MA0052.4), score >= %.1f, top %d per enhancer -> %d hits'
      % (CUTOFF, TOP_N, len(h)))
if len(h):
    b = h.loc[h.groupby('enhancer').score.idxmax()]
    print('best per enhancer: ' + ', '.join('%s %.2f' % (r.enhancer, r.score) for r in b.itertuples()))
    arch = os.path.join(ARCHIVE, 'nrxn3_mef2c_motif_hits.tsv')
    if os.path.exists(arch):
        a = pd.read_csv(arch, sep='\t')
        m = a.merge(h, on=['enhancer', 'start', 'strand', 'hit_sequence'], how='inner')
        print('vs the archived scan: %d of %d delivered rows reproduced, max|score delta| %.2f'
              % (len(m), len(a), (m.score_x - m.score_y).abs().max() if len(m) else float('nan')))
""")

A(MD, r"""
**Expected** - **19** hits, split 5 / 4 / 5 / 5 over E1-E4, best score per enhancer **12.94 / 15.12 /
12.88 / 16.22** (global best E4), and all 19 rows of the archived `nrxn3_mef2c_motif_hits.tsv`
reproduced on enhancer + start + strand + sequence with `max|score delta| = 0`. Those are the
parameters of `scan_mef2_motif.py`: log2-odds PSSM with pseudocount 1 from `metacluster_179.4`
(= JASPAR MA0052.4), both strands, `score >= 7.0`, top 5 per enhancer. A different count means the
cutoff or the top-N cap moved - the genome is checked, not assumed (the mm10.fa window is asserted
equal to the pre-cut record).
""")

A(CODE, r"""
# ---- 3.3 the panel: coverage tracks + expression violins ----
trackf = os.path.join(DATA, 'nrxn3_tracks.npz')
fig, axes = plt.subplots(2, 1, figsize=(12, 9), gridspec_kw={'height_ratios': [3, 1]})
ax = axes[0]
if os.path.exists(trackf):
    z = np.load(trackf, allow_pickle=True)
    win_s, win_e = int(z['WIN_S']), int(z['WIN_E'])
    for i, s in enumerate(STATES):
        y = (len(STATES) - i) * 1.0
        cov = z[s]
        x = win_s + np.arange(len(cov)) * 200.0
        ax.fill_between(x, y, y + np.asarray(cov, float) / (np.nanmax(cov) + 1e-12) * 0.85,
                        color='0.25', lw=0)
        ax.text(win_s - 4e4, y + 0.4, s, ha='right', va='center', fontsize=8)
    for lab, st, en, nes in ENH:
        ax.plot([st, en], [0.15, 0.15], color='crimson', lw=2.5)
        ax.text((st + en) / 2, 0.02, '%s NES=%.2f' % (lab, nes), color='crimson',
                fontsize=7, ha='center', va='top')
    ax.set_xlim(win_s, win_e); ax.set_ylim(-0.4, len(STATES) + 0.4)
    ticks = [TSS + k * 10 ** 5 for k in range(0, 10)]
    ax.set_xticks(ticks)
    ax.set_xticklabels(['TSS'] + ['+%d kb' % (k * 100) for k in range(1, 10)], fontsize=7)
    ax.set_yticks([]); ax.spines[['left', 'right', 'top']].set_visible(False)
    ax.set_title('%s locus: ATAC coverage per state (RPM, 200-bp bins) + Mef2c eRegulon regions' % GENE, fontsize=9)
else:
    ax.text(0.5, 0.5, 'nrxn3_tracks.npz not in FIGS5_DATA', ha='center')

ax = axes[1]
ax.violinplot([vals[s] for s in STATES], showextrema=False, widths=0.7)
ax.set_xticks(range(1, len(STATES) + 1)); ax.set_xticklabels(STATES, fontsize=8)
ax.axvline(len(STATES) - 1.5, color='none')
ax.set_ylabel('%s expression (2^x - 1)' % GENE, fontsize=8)
ax.set_title('MOL3 median %.1f vs %.1f-%.1f in the other states' % (
    np.median(vals['MOL3']), min(np.median(vals[s]) for s in STATES if s != 'MOL3'),
    max(np.median(vals[s]) for s in STATES if s != 'MOL3')), fontsize=8)
fig.tight_layout()
f = os.path.join(OUT, 'figS5c_nrxn3_locus.pdf'); fig.savefig(f, bbox_inches='tight'); plt.close(fig)
show(f)
""")

# --------------------------------------------------------------------------- 4
A(MD, r"""
---

## 4. Panel D - how many eRegulons survive fewer cells

**In** `eregulon_stability_table.csv` **Out** `figS5d_stability.pdf`

The sweep itself is *not* re-run here: each depth is a complete SCENIC+ run on a
fixed-seed subsample (`stability_downsample.py --n-cells 1000 --seed 123`, with
`build_grn` and AUCell inside), which is why the archive holds one
`rho_<n>c_seed123_filtered_sigs.pkl` per depth. The panel only reads the summary
table, so it is re-plottable on a laptop.
""")

A(CODE, r"""
# ---- 4.1 stability table ----
stabb = pd.read_csv(os.path.join(DATA, 'eregulon_stability_table.csv'))
print(stabb.to_string(index=False))
stabb.to_csv(os.path.join(OUT, 'S5d_stability_table.csv'), index=False)
""")

A(MD, r"""
**Expected** - three depths, the full set last:

| n_cells | eRegulons | mean Jaccard vs full | Top-10 TF recovery |
|---|---|---|---|
| 1,000 | **196** | 0.110 | 1.00 |
| 3,000 | **294** | 0.280 | 1.00 |
| 5,977 (full) | **292** | 1.000 | 1.00 |

The non-monotonicity is real and is what the panel is about: 3,000 cells give
*more* filtered eRegulons (294) than all 5,977 (292), because filtering is a
threshold on correlation-based statistics that move with n; the number of
raw cistromes (507 -> 746 -> 751) is monotone. The Top-10 TF recovery of 1.00 at
every depth is the sentence the figure legend rests on.

> ⚠ **Registered, not fixed**: the delivered legend (`figure_captions.md`, S5
> block A / stability text) quotes "1k: 208 -> 6k: 292" and "0.11 -> 3k: 0.26 ->
> 1.00". The figure face and this table say **196** and **0.28**. The PDF was
> drawn from this table (`figs5_stability.py`), so the legend numbers are the
> ones out of date; the manuscript text and `figures.docx` were signed off on
> 2026-10-21, so nothing here was edited silently.
""")

A(CODE, r"""
# ---- 4.2 the panel ----
fig, axes = plt.subplots(1, 3, figsize=(12, 3.2))
for ax, col, ttl, fmt in [(axes[0], 'n_eregulons', 'eRegulons kept', '%.2f'),
                          (axes[1], 'mean_jaccard_vs_full', 'Jaccard vs full', '%.2f'),
                          (axes[2], 'top10_tf_recovery', 'Top-10 TF recovery', '%.2f')]:
    ax.plot(stabb.n_cells, stabb[col], 'o-')
    for x, y in zip(stabb.n_cells, stabb[col]):
        ax.annotate(fmt % y, (x, y), textcoords='offset points', xytext=(0, 7), ha='center', fontsize=8)
    ax.set_xscale('log'); ax.set_xticks(stabb.n_cells)
    ax.set_xticklabels(['%d' % v for v in stabb.n_cells], fontsize=7)
    ax.set_xlabel('cells entering the run', fontsize=8); ax.set_title(ttl, fontsize=9)
fig.tight_layout()
f = os.path.join(OUT, 'figS5d_stability.pdf'); fig.savefig(f, bbox_inches='tight'); plt.close(fig)
show(f)
""")

# --------------------------------------------------------------------------- 5
A(MD, r"""
---

## 5. Software and the parameters that decide whether this reproduces
""")

A(CODE, r"""
import importlib.metadata as md
PKGS = ['numpy', 'pandas', 'scipy', 'matplotlib', 'scanpy', 'anndata', 'mudata',
        'ctxcore', 'pycistarget', 'pyscenic', 'scenicplus', 'dill', 'h5py', 'networkx', 'plotnine']
here = {}
for p in PKGS:
    try: here[p] = md.version(p)
    except Exception: here[p] = 'not installed'
print('python %s | %s' % (platform.python_version(), platform.platform()))
for k, v in here.items(): print('  %-14s %s' % (k, v))

with open(os.path.join(OUT, 'requirements.txt'), 'w') as f:
    f.write('\n'.join('%s==%s' % kv for kv in here.items() if kv[1] != 'not installed') + '\\n')

print()
print('What must not change to reproduce the delivered numbers:')
for line in ['eRSS: auc_key=eRegulon_AUC_filtered, signature_keys=[Region_based] (96 eRegulons), variable=ACC_celltype_l2, out_key_suffix=_cmp',
             'networks: eRegulon_metadata_filtered, use_edge_tables=[TF2G], min_deg 2/3/3 (MOL2/MOL3_6/MOL6), <=60 genes, kamada_kawai_layout, scale_position_by=450',
             'locus: TSS chr12:88,722,433; four enhancer regions as listed in 3.1; MEF2 = metacluster_179.4 (JASPAR MA0052.4), PSSM pseudocount 1, score cut 8',
             'stability: seed 123 per depth, full SCENIC+ re-run, table read as-is',
             'AUCell elsewhere: seed 555 (Fig 5C), auc_threshold 0.05']:
    print('  - ' + line)
""")

A(MD, r"""
Versions used for the delivered analysis (SCENIC+ environment, `pip list` at run
time): `numpy 1.22.4 / pandas 1.5.1 / scipy 1.10.1 / matplotlib 3.6.2 / scanpy
1.9.2 / ctxcore 0.2.0 / dill 0.3.6 / h5py 3.8.0 / plotnine 0.12.4 / scenicplus
1.0.1.dev4+ge4bdd9f / pysenic 0.12.1+6.g31d51a1`, python 3.9. The two packages
that decide the regulon inference (`scenicplus`, `pysenic`) are the same builds
as above; the rest of this environment is newer, which affects plots, not the
tables. R/Seurat versions of the export step are not recorded, so none is quoted.
""")


def build(path='Figure_S5.ipynb'):
    nb = nbformat.v4.new_notebook()
    nb.metadata = {'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
                   'language_info': {'name': 'python', 'version': platform.python_version()}}
    for kind, src in CELLS:
        nb.cells.append(nbformat.v4.new_markdown_cell(src) if kind == 'markdown'
                          else nbformat.v4.new_code_cell(src))
    nbformat.validate(nb)
    nbformat.write(nb, path)
    code = sum(1 for c in nb.cells if c.cell_type == 'code')
    for i, c in enumerate(nb.cells):
        if c.cell_type == 'code':
            compile(c.source, 'cell %d' % i, 'exec')
    print('%s | %d cells | %d code | all code compiles' % (path, len(nb.cells), code))


if __name__ == '__main__':
    build(os.environ.get('FIGS5_NB', 'Figure_S5.ipynb'))
