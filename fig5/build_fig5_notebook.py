#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build Figure_5.ipynb (SCENIC+ eRegulons) from cell definitions.

  python3 build_fig5_notebook.py        ->  Figure_5.ipynb

House rules (author's requirements, 2026-10-21):
  * one section per panel: what goes in -> shortest code -> what comes out
  * no compositing, no font/label tuning (journal layout is done by hand)
  * English comments; core parameters stated explicitly
  * every code cell prints the number the published panel shows
"""
import json
import os

MD, CODE = "markdown", "code"


_n = [0]


def cell(kind, src):
    _n[0] += 1
    lines = src.strip("\n").split("\n")
    src = [l + "\n" for l in lines[:-1]] + [lines[-1]]
    if kind == CODE:
        return {"id": "c%02d" % _n[0], "cell_type": kind, "metadata": {},
                "execution_count": None, "outputs": [], "source": src}
    return {"id": "m%02d" % _n[0], "cell_type": kind, "metadata": {}, "source": src}


C = []
A = C.append

# ------------------------------------------------------------------ front ---
A((MD, r"""
# Figure 5 — SCENIC+ eRegulons of the oligodendrocyte lineage

Three panels built from the SCENIC+ object of the dsciCAT-seq data plus two public
single-cell technologies. This notebook reproduces **the content of each panel**.
Journal assembly, Arial substitution and label-collision tuning are done by hand in
Illustrator and are deliberately not reproduced.

| Panel | Content | Input | Key parameters | Output |
|---|---|---|---|---|
| **A** | activity x expression dot plot, 20 eRegulons x 6 states | `scplus_obj_final.pkl` | gene-vs-region AUC `abs(rho) > 0.1`; `> 5` target genes; one row per TF; min-max 0-1 | `fig5a_eregulon_dotplot.pdf` |
| **B** | MOL3 eRegulon-target network | same object + Fig 3 marker list | 5 anchor TFs; markers + hubs(`deg >= 2..6`) capped at `60` genes; `TF2G` edges only; `kamada_kawai_layout` | `fig5b_grn_mol3.pdf` (mol2/mol3_6/mol6 -> Fig S5) |
| **C** | cross-technology reproducibility, two scatters | object + `10x_subset_harmony.RDS`, `HT_harmony.RDS` | gene-based eRegulons with `>= 3` targets (**n = 145**); AUCell **`seed = 555`**; cell-count-weighted lineage mean; Spearman | `fig5c_crosstech_scatter.pdf` |

```
scplus_obj_final.pkl   (5,977 cells, 181 eRegulons, 6 states)
   |
   +- std filtering (scenicplus default) .......... 96 eRegulons
   +- pseudobulk TF-cistrome correlation (seed 555) . rho / peak per eRegulon   [B]
   +- panel A selection .......................... 20 eRegulons -> dot plot     [A]
   +- eRegulon_metadata_filtered ................ 44 nodes / 44 edges -> GRN   [B]
   +- gene signatures with >= 3 targets .......... 145 signatures               [C]
         |
   10x / HT Seurat RDS --export mtx--> AUCell(seed=555) --> mean AUC per cell type
         |
   weighted lineage mean  vs  dsciCAT 6-state weighted mean  --> Spearman        [C]
```
"""))

A((MD, r"""
## 0. Environment, inputs, outputs

Two variables drive every path below.

| Variable | Default | Must contain |
|---|---|---|
| `FIG5_DATA` | `data` | `scplus_obj_final.pkl` · `10x_subset_harmony.RDS` · `HT_harmony.RDS` · `NIHMS70326-supplement-Table_S1.xlsx` (Marques 2016, 24 states x 50 up-regulated genes) |
| `FIG5_OUT` | `fig5_out` | created here; every CSV/PDF product lands here |

Panels A and B need the SCENIC+ python environment (`scenicplus`, `pyscenic`,
`ctxcore`, `plotnine`, `dill`); panel C additionally needs an **R with Seurat v4 or v5
plus Signac** (these two public RDS objects carry an ATAC/ChromatinAssay, and R will not
open them without Signac). `FIG5_R` points at the interpreter; the default is the shared
`/public/tools/r-packages/R-4.5.3/Rscript` (Seurat 5.5.1 / SeuratObject 5.4.0 / Signac
1.17.1) when it exists, else the `Rscript` on `PATH`. Everything after that export is plain
pandas/scipy/matplotlib and can be run on a laptop.

How the input object came to exist (upstream; not re-run here). **The object is not a KLP
table and no KLP was ever produced here**: this SCENIC+ build (1.0.1.dev4) has no `run_klp`
step at all - region-to-gene linking is a GBM classifier inside `run_scenicplus`, and the
TF-to-gene edges come from the gsea-based `build_grn`. The chain below is read off the
notebooks that produced it (`Pycistopic.ipynb`, `Pycistarget.ipynb`, `Scenic+.ipynb`,
`scenic_and_scarlink_prepare.ipynb`) and off the timestamps of the archived files in
`FY/oligo/3kdsci/GRN/`. It was run on the 5,977 Harmony-integrated oligodendrocyte-lineage
cells annotated `ACC_celltype_l2`. Two caveats: the reference notebooks still carry the MOL
sibling of this subset in their paths (the same scripts were run on both), and the dates below
are file mtimes, quoted as the ordering evidence they are.

| # | Step | Command / parameters | Product on disk |
|---|---|---|---|
| 1 | Origin is a **Seurat/Signac RDS** of the dsciCAT oligo subset (RNA + ATAC assays, per-cell `celltype_l2`) - not CellRanger-ARC, not ArchR | - | `FY/oligo/3k_dsci_oligo_harmony.RDS` -> lineage subset |
| 2 | R-side export for python: RNA counts + log-normalised values (all genes), **binarised peak matrix** (`atac_assay@data > 0 -> 1`), cell metadata | `scenic_and_scarlink_prepare.ipynb` | `FY/oligo/3k_dsci_oligo_subset_rna.h5ad` (09-04), `..._subset_atac.h5ad` (08-20) |
| 3 | ATAC side into cisTopic: peak names rewritten `chr-start-end` -> `chr:start-end`, integer sparse fragment matrix (dtype asserted before use) | `create_cistopic_object(fragment_matrix = adata.X.T.tocsr().astype(int), cell_names = ..., region_names = ..., tag_cells = False, split_pattern = ' ')` | `GRN/cistarget/oligo_subset_cisTopic.pkl` |
| 4 | LDA topic models | `run_cgs_models(n_topics = [5, 20, 40, 60, 80, 100, 130, 150, 180], n_cpu = 20, n_iter = 150, random_state = 555, alpha = 50, alpha_by_topic = True, eta = 0.1, eta_by_topic = False, top_topics_coh = 5, save_path = ...)` | `GRN/cistarget/oligo_subset_cisTopic_models.pkl` (08-26, 634 MB) |
| 5 | Topic number compared on per-topic coherence, 20 vs 40 vs 60, with cell-topic clustering / UMAP / tSNE and topic-vs-state annotation as QC | `find_clusters(target = 'cell', k = 10, res = [0.1, 0.5, 1], prefix = 'pycisTopic_')`, `run_umap`, `run_tsne`, `cell_topic_heatmap`, `topic_annotation(annot_var = 'celltype_l2', binarized_cell_topic = top_3k, general_topic_thr = 0.2)` -> **60 topics used** | `oligo_subset_cisTopic60.pkl` (09-02, the one that feeds step 10); `cisTopic20.pkl` / `cisTopic40.pkl` (09-09) are the comparison re-runs |
| 6 | Binarisation of the 60-topic model into region sets, three rules | `binarize_topics(method = 'otsu')`, `binarize_topics(method = 'ntop', ntop = 3_000)`, `binarize_topics(method = 'li')` | `oligo_subset_{otsu,top3k,li}.pkl` (09-03) |
| 7 | Differential regions per state = the DAR cistromes; accessibility imputed for the marker test | `find_marker_features(...)`, `impute_accessibility(scale_factor = 10**6)`, `normalize_scores(scale_factor = 10**4)` | `oligo_subset_markers_dict.pkl` (09-03) |
| 8 | Gene annotation for `species = 'custom'`: mm10 GTF restricted to standard chromosomes, joined to chromsizes | `<refdata>/mm10/ATAC/scATAC_ref/genes.gtf`, `<work>/referrence/mm10/standard.chrom.sizes` | `referrence/mm10/custom_annot_df.tsv` (09-16) |
| 9 | cisTarget database: cistrome regions -> genome FASTA (bedtools `getfasta -nameOnly`) plus a 1 kb-padded background, against the downloaded mouse motif collection | `v10nr_clust_public` (from `resources.aertslab.org/cistarget/`), annotations `motifs-v10-nr.mgi-m0.00001-o0.0.tbl` | `GRN/cistarget/1kb_bg_padding.fa` (09-01), `consensus_regions.bed`, `<prefix>.regions_vs_motifs.rankings.feather`, `<prefix>.regions_vs_motifs.scores.feather` |
| 10 | pyCisTarget enrichment of every cistrome: topics (top-3k and li) and DARs, with and without promoters | `run_pycistarget(region_sets = {'topics_top_3k', 'topics_li', 'DARs'}, species = 'custom', custom_annot = ..., ctx_db_path = '...rankings.feather', dem_db_path = '...scores.feather', path_to_motif_annotations = '.../motifs-v10-nr.mgi-m0.00001-o0.0.tbl', run_without_promoters = True, n_cpu = 40, _temp_dir = '/tmp', annotation_version = 'v10nr_clust')` | `CTX_topics_li_{All,No_promoters}`, `CTX_topics_top_3k_*`, `CTX_DARs_*`, `DEM_*` (60 `TopicN.html` per set, i.e. one per topic) and `menr.pkl` (09-03, 3.7 GB) |
| 11 | Object assembly: RNA side = full counts through `.raw`, chromatin side = the 60-topic object, cistromes = `menr` | `create_SCENICPLUS_object(GEX_anndata = adata.raw.to_adata(), cisTopic_obj = oligo_subset_cisTopic60.pkl, menr = menr.pkl)` | in-memory `scplus_obj` |
| 12 | Full SCENIC+ run | `run_scenicplus(variable = ['ACC_celltype_l2'], species = 'mmusculus', assembly = 'mm10', tf_file = '.../allTFs_mm.txt', save_path = 'GRN/scplus', upstream = [1000, 150000], downstream = [1000, 150000], calculate_TF_eGRN_correlation = True, calculate_DEGs_DARs = True, export_to_loom_file = False, export_to_UCSC_file = False, n_cpu = 70, _temp_dir = os.environ.get('TMPDIR', ''))` -> **10 eRegulons and 0 true repressors**, which is why step 13 exists | `GRN/scplus/scplus_obj.pkl` (09-08, 8.0 GB) |
| 13 | GRN rebuilt with looser thresholds, then re-scored per cell | `build_grn` -> `format_egrns` -> `get_eRegulons_as_signatures` -> `make_rankings` (region **and** gene) -> `score_eRegulons(auc_threshold = 0.05, normalize = False)` | `GRN/scplus/region_ranking.pkl`, `gene_ranking.pkl` (09-08) |
| 14 | Pseudobulk TF-cistrome correlation, standard filtering, eRSS | `generate_pseudobulks` + `TF_cistrome_correlation(use_pseudobulk = True)` -> rho / adj_p; filter `abs(rho) > 0.1` and `> 5` targets; `regulon_specificity_scores` | `GRN/scplus/scplus_obj_final.pkl` (09-09, 8.0 GB) - **the file opened below** |

Read the table this way: the chromatin side of the input is a **binarised peak matrix taken out
of an RDS, modelled by cisTopic at 60 topics, binarised three ways, and enriched against the
v10nr motif collection by pyCisTarget**; an eRegulon here = a motif-enriched topic or DAR
cistrome (its region set) plus gsea-derived TF-to-gene targets, scored per cell by AUCell. Panel
A's 20 eRegulons, panel B's MOL3 network and panel C's 145 signatures are all subsets of this
one object; panel C only re-uses the 145 filtered signatures on other technologies' matrices (§3).
"""))

A((CODE, r"""
import os, sys, subprocess
import numpy as np
import pandas as pd

DATA = os.environ.get("FIG5_DATA", "data")
OUT  = os.environ.get("FIG5_OUT", "fig5_out")
os.makedirs(OUT, exist_ok=True)

GROUP  = "ACC_celltype_l2"                                  # state annotation in the object
STATES = ["COP_NFOL", "MFOL", "MOL2", "MOL3", "MOL3_6", "MOL6"]   # developmental order
SEED   = 555                                                # AUCell + pseudobulk seed

import dill
obj = dill.load(open(os.path.join(DATA, "scplus_obj_final.pkl"), "rb"))   # ~1 min
print("cells %d | genes %d | eRegulons %d" % (
    len(obj.cell_names), len(obj.gene_names),
    obj.uns["eRegulon_metadata"]["Gene_signature_name"].nunique()))
"""))

A((MD, r"""
**Expected**

```
cells 5977 | genes 24033 | eRegulons 181
```
"""))

# ---------------------------------------------------------------- panel A ---
A((MD, r"""
## 1. Panel A — activity x expression dot plot (20 eRegulons x 6 states)

**In:** the object. **Out:** `fig5a_eregulon_dotplot.pdf`.
Dot **size** = mean AUCell activity (region-based AUC), cell **colour** = mean
target-gene expression; both min-max scaled to 0-1 over the displayed rows.

The 20 rows are *not* a top-20 ranking. They are the official `heatmap_dotplot`
subset, obtained by three rules applied to all 181 eRegulons:

| # | Rule | Why |
|---|---|---|
| 1 | gene-based and region-based AUC profiles agree: `abs(corrwith) > 0.1` | drops eRegulons whose two evidence layers disagree |
| 2 | keep `+_+` (activator) and `-+_` (true repressor); if a TF has both a direct and an `extended_` eRegulon, keep the direct one | one row per TF |
| 3 | keep eRegulons with more than 5 target genes | rows the dot plot cannot score are dropped |

Naming convention `TF[_extended]_{TF2G}_{R2G}_(N)`: position 1 is the TF->gene sign,
position 2 the region->gene sign. A surviving `*_extended_*` eRegulon is printed with
an asterisk (co-accessibility-supported links included); `(nr)` in the row label is
the cistrome size.
"""))

A((CODE, r"""
# ---- 1a. sign structure of all 181 eRegulons (the numbers quoted in the caption) ----
meta = obj.uns["eRegulon_metadata"]
sig = meta[["Gene_signature_name", "TF2G_rho", "R2G_rho"]].drop_duplicates("Gene_signature_name")
print("TF2G+ R2G+ : %d" % int(((sig.TF2G_rho > 0) & (sig.R2G_rho > 0)).sum()))
print("TF2G+ R2G- : %d" % int(((sig.TF2G_rho > 0) & (sig.R2G_rho < 0)).sum()))
print("true repressors (TF2G-): %d" % int((sig.TF2G_rho < 0).sum()))
print("std-filtered:", obj.uns["eRegulon_metadata_filtered"]["Gene_signature_name"].nunique())
"""))

A((MD, r"""
**Expected**

```
TF2G+ R2G+ : 124
TF2G+ R2G- : 57
true repressors (TF2G-): 0
std-filtered: 96
```

Nothing negative survives, so the panel has no activator/repressor split
(`split_repressor_activator=False`). This is a property of the activity-inference
step, not a claim that the lineage has no repression.
"""))

A((CODE, r"""
# ---- 1b. the three selection rules -> exactly 20 eRegulons ----
gb, rb = obj.uns["eRegulon_AUC"]["Gene_based"], obj.uns["eRegulon_AUC"]["Region_based"]
a = gb.copy(); a.columns = [x.split("_(")[0] for x in gb.columns]
b = rb.copy(); b.columns = [x.split("_(")[0] for x in rb.columns]
cor = a.corrwith(b, axis=0)                       # per eRegulon: gene- vs region-based AUC
cor = cor[abs(cor) > 0.1]                          # rule 1

keep = [x for x in cor.index if "+_+" in x] + [x for x in cor.index if "-_+" in x]
ext = [x for x in keep if "extended" in x]        # rule 2: direct form wins over extended
direct = [x for x in keep if "extended" not in x]
keep = direct + [x for x in ext if x.replace("extended_", "") not in direct]

sel_gene = [x for x in gb.columns if x.split("_(")[0] in set(keep)]
sel_gene = [x for x in sel_gene if int(x.split("_(")[1].replace("g)", "")) > 5]   # rule 3
sel_region = [x for x in rb.columns
              if x.split("_(")[0] in {y.split("_(")[0] for y in sel_gene}]
obj.uns["selected_eRegulons"] = {"Gene_based": sel_gene, "Region_based": sel_region}
assert len(sel_gene), "nothing selected: check the eRegulon naming in this object"
tfs = sorted({x.split("_")[0] for x in sel_gene})
print(len(tfs), "eRegulons:"); print(tfs)
"""))

A((MD, r"""
**Expected** — exactly 20; rows are ordered in the figure by expression phase
(early `Sox6, Tcf7l2, Lcorl` / mid `Nfix, Nfia, Mef2c` / late `Creb5, Klf7`):

```
20 eRegulons:
['Arnt2', 'Bhlhe40', 'Creb5', 'Etv6', 'Glis1', 'Klf7', 'Lcorl', 'Mbnl2', 'Mef2c',
 'Mllt10', 'Nfia', 'Nfix', 'Pparg', 'Rfx3', 'Rora', 'Sox6', 'Tcf12', 'Tcf4',
 'Tcf7l2', 'Zbtb20']
```
"""))

A((CODE, r"""
# ---- 1c. the panel ----
from scenicplus.plotting.dotplot import heatmap_dotplot

idx_order = [s for s in STATES if s in obj.metadata_cell[GROUP].unique()]
size = obj.uns.get("eRegulon_AUC_filtered", obj.uns["eRegulon_AUC"])["Region_based"]
f = os.path.join(OUT, "fig5a_eregulon_dotplot.pdf")
heatmap_dotplot(scplus_obj=obj, size_matrix=size, color_matrix=obj.to_df("EXP"),
                scale_size_matrix=True, scale_color_matrix=True,     # 0-1 min-max
                group_variable=GROUP, subset_eRegulons=obj.uns["selected_eRegulons"]["Gene_based"],
                index_order=idx_order, split_repressor_activator=False,
                figsize=(7, 10), orientation="vertical", save=f)
print("wrote", f)
"""))

A((MD, r"""
`heatmap_dotplot` does the per-state aggregation itself (mean AUC, mean expression
of the eRegulon's target genes), so the PDF *is* the panel content, row labels,
asterisks and `(nr)` suffixes included. Verified byte-identical on re-run (only the
PDF `CreationDate` differs).
"""))

# ---------------------------------------------------------------- panel B ---
A((MD, r"""
## 2. Panel B — MOL3 eRegulon-target network

**In:** `eRegulon_metadata_filtered` (96 eRegulons) + the curated marker list of
Fig 3. **Out:** `fig5b_grn_mol3.pdf` and the edge/node tables behind it.

One network per mature state, from hand-picked anchor TFs; gene nodes are fixed by
two rules so the picture stays readable:

| Step | Rule | MOL3 |
|---|---|---|
| anchor TFs | hand-picked from the eRegulons whose activity peaks in this state (table 2a) | `Mef2c` + 4 `*_extended` |
| gene nodes | curated markers that are targets of an anchor, **plus** hub genes co-targeted by `>= min_deg` anchors | 39 |
| edges | `TF2G` only | 44 |
| layout | `kamada_kawai_layout`, positions x 450 | 44 nodes |
| encodings | gene colour = state `Log2FC` (`bwr`); edge alpha/width = `TF2G_importance` | — |

Two traps worth writing down:

1. `create_nx_tables` builds the pattern `'<name>_[^a-zA-Z0-9]'` internally, so a TF
   name passed *with* its sign suffix (`Mef2c+`) becomes a regex quantifier, matches
   nothing and silently removes every TF from the graph. Pass plain names.
2. Keeping `R2G` edges adds one diamond node per regulatory region; at ~100 regions
   the graph is unreadable. The published panel uses `TF2G` only.
"""))

A((CODE, r"""
# ---- 2a. anchor evidence: pseudobulk TF-cistrome correlation + peak state ----
from scenicplus.cistromes import generate_pseudobulks, TF_cistrome_correlation

obj.uns.pop("Pseudobulk", None)                    # recompute deterministically
generate_pseudobulks(obj, variable=GROUP, auc_key="eRegulon_AUC_filtered",
                     signature_key="Region_based", nr_cells=10, nr_pseudobulks=500, seed=SEED)
TF_cistrome_correlation(obj, use_pseudobulk=True, variable=GROUP,
                        auc_key="eRegulon_AUC_filtered", signature_key="Region_based",
                        out_key="cmp_region_based")
corr = obj.uns["TF_cistrome_correlation"]["cmp_region_based"].copy()
corr["eRegulon"] = corr["Cistrome"].str.split("_(", regex=False).str[0]
rho_map = dict(zip(corr["eRegulon"], corr["Rho"]))                # names repeat: map, not reindex
p_map = dict(zip(corr["eRegulon"], corr["Adjusted_p-value"]))

auc = obj.uns["eRegulon_AUC_filtered"]["Region_based"]          # cells x 96 eRegulons
mean_auc = auc.groupby(obj.metadata_cell[GROUP].astype(str).values).mean().T   # eRegulon x state
mean_auc.index = [i.split("_(")[0] for i in mean_auc.index]      # 'Mef2c_+_(153r)' -> 'Mef2c_+'
# NB the *_filtered* matrices use the short 'TF_sign_(N)' naming, the unfiltered ones
# use 'TF_sign_sign_(N)' (panel A). Both strip consistently within themselves.
peak = mean_auc.idxmax(axis=1)                                   # state with the highest mean AUC
ngene = obj.uns["eRegulon_metadata_filtered"].groupby(
    "Region_signature_name")["Gene"].nunique()
ngene.index = [i.split("_(")[0] for i in ngene.index]

tab = pd.DataFrame({"eRegulon": mean_auc.index, "peak_celltype": peak.values,
                    "Rho": [rho_map.get(n, np.nan) for n in mean_auc.index],
                    "adj_p": [p_map.get(n, np.nan) for n in mean_auc.index],
                    "n_genes": ngene.reindex(mean_auc.index).values})
# AUC column names and Cistrome names must strip to the same key; if the object
# uses the older single-sign naming this join silently produces NaN instead.
assert tab.Rho.notna().all(), "eRegulon naming mismatch: AUC columns vs TF_cistrome_correlation"
tab["informative"] = (tab.Rho > 0.3) & (tab.n_genes > 5) & (tab.adj_p < 0.05)
tab.to_csv(os.path.join(OUT, "eregulon_peak.csv"), index=False)
print(tab[tab.eRegulon.str.startswith("Mef2c")].to_string(index=False))
print("peaks per state:", tab[tab.informative].peak_celltype.value_counts().to_dict())
"""))

A((MD, r"""
**Expected** — the number printed on the panel and in the caption:

```
 eRegulon  peak_celltype   Rho  adj_p  n_genes  informative
Mef2c_+             MOL3  0.18    0e+00       60        False
```

Read this honestly: `Mef2c` peaks in MOL3 but its pseudobulk TF-cistrome
`rho = +0.18` is *modest* — only 15 of the 96 filtered eRegulons pass the
`rho > 0.3 & > 5 genes & p < 0.05` screen, and `Mef2c` is not one of them. The MOL3
network is anchored on `Mef2c` for mechanistic reasons (its target *Nrxn3* can be
followed to four MEF2-motif enhancers at the locus level, Fig S5), and the caption
therefore reports `rho = +0.18` as an association, not as strong TF-cistrome
coupling. RSS (`regulon_specificity_scores`, same 96 filtered eRegulons) is the
specificity measure used in Fig S5A; it is not needed to draw this panel.
"""))

A((CODE, r"""
# ---- 2b. marker list (Fig 3 / Fig S1) and the per-state figure configuration ----
markers = {}
s1 = pd.ExcelFile(os.path.join(DATA, "NIHMS70326-supplement-Table_S1.xlsx")).parse("Sheet1")
for c in range(1, s1.shape[1]):
    nm = s1.iloc[1, c]
    if isinstance(nm, str):
        markers[nm.strip()] = {g for g in s1.iloc[2:, c].dropna().astype(str).str.strip() if g}
print("Table S1 states:", len(markers), "| MOL3 markers:", len(markers["MOL3"]))

FIGURES = [
    dict(tag="mol3",   mstate="MOL3",        var=GROUP + "_Log2FC_MOL3",
         tfs=["Mef2c", "Rora_extended", "Mbnl2_extended", "Pparg_extended", "Bcl11a_extended"],
         markers=markers["MOL3"] | {"Ryr2", "Nrxn3"}),
    dict(tag="mol2",   mstate="MOL2",        var=GROUP + "_Log2FC_MOL2",
         tfs=["Nfix", "Tcf4_extended", "Klf7", "Tcf12", "Mbnl2_extended"],
         markers=markers["MOL2"] | {"Apod", "Anln"}),
    dict(tag="mol6",   mstate="MOL6/MOL5-6", var=GROUP + "_Log2FC_MOL6",
         tfs=["Creb5_extended", "Mbnl2_extended"],
         markers=markers["MOL6"] | markers["MOL5/6"]),
    dict(tag="mol3_6", mstate="MOL2-3/MOL5-6", var=GROUP + "_Log2FC_MOL3_6",
         tfs=["Sox6", "Mllt10_extended", "Creb5_extended"],
         markers=markers["MOL2/3"] | markers["MOL5/6"]),
]
"""))

A((CODE, r"""
# ---- 2c. gene/region pruning + graph build + draw (one pass per state) ----
import re
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scenicplus.networks import create_nx_tables, create_nx_graph, plot_networkx

GRAPHS = {}
for f in FIGURES:
    sub = obj.uns["eRegulon_metadata_filtered"]
    sub = sub[sub["Region_signature_name"].str.contains(
        "|".join(re.escape(t) + "_" for t in f["tfs"]))]
    linked = f["markers"] & set(sub["Gene"])
    deg = Counter(sub["Gene"])
    for min_deg in (2, 3, 4, 5, 6):                       # adaptive hub threshold
        hubs = {g for g, n in deg.items() if n >= min_deg}
        if len(linked | hubs) <= 60:
            break
    genes = sorted(linked | hubs)
    regions = sorted(set(sub[sub["Gene"].isin(genes)]["Region"]))

    nt = create_nx_tables(scplus_obj=obj, eRegulon_metadata_key="eRegulon_metadata_filtered",
                          subset_eRegulons=f["tfs"], subset_regions=regions, subset_genes=genes,
                          add_differential_gene_expression=True,
                          add_differential_region_accessibility=True,
                          differential_variable=[GROUP])
    tfs = nt["Node"]["TF"]["TF"].dropna().unique().tolist()
    cmap = {tf: plt.get_cmap("tab10")(i % 10) for i, tf in enumerate(tfs)}
    G, pos, _, _ = create_nx_graph(
        nt, use_edge_tables=["TF2G"],
        color_edge_by={"TF2G": {"variable": "TF", "category_color": cmap}},
        transparency_edge_by={"TF2G": {"variable": "TF2G_importance", "min_alpha": 0.2, "v_min": 0}},
        width_edge_by={"TF2G": {"variable": "TF2G_importance", "max_size": 4, "min_size": 1.0}},
        color_node_by={"TF": {"variable": "TF", "category_color": cmap},
                       "Gene": {"variable": f["var"], "continuous_color": "bwr"}},
        transparency_node_by={"Gene": {"variable": f["var"], "min_alpha": 0.2}},
        size_node_by={"TF": {"variable": "fixed_size", "fixed_size": 150},
                      "Gene": {"variable": "fixed_size", "fixed_size": 70}},
        shape_node_by={"TF": {"variable": "fixed_shape", "fixed_shape": "circle"},
                       "Gene": {"variable": "fixed_shape", "fixed_shape": "ellipse"}},
        label_size_by={"TF": {"variable": "fixed_label_size", "fixed_label_size": 22.0},
                       "Gene": {"variable": "fixed_label_size", "fixed_label_size": 11.0}},
        layout="kamada_kawai_layout", scale_position_by=450)
    GRAPHS[f["tag"]] = (G, nt, f)
    plt.subplots(figsize=(12, 12)); plot_networkx(G, pos)
    out = os.path.join(OUT, "fig5b_grn_%s.pdf" % f["tag"])
    plt.savefig(out, bbox_inches="tight"); plt.close()
    print("%-7s min_deg=%d genes=%d -> %d nodes / %d edges"
          % (f["tag"], min_deg, len(genes), G.number_of_nodes(), G.number_of_edges()))
"""))

A((MD, r"""
**Expected**

```
mol3    genes=39 -> 44 nodes / 44 edges
mol2    ...                -> 58 nodes / 97 edges
mol6    ...                -> 49 nodes / 65 edges
mol3_6  ...                -> 56 nodes / 73 edges
```

MOL3 is the panel shown in Fig 5: it is the only state whose anchor eRegulon can be
followed to locus-level evidence (`Mef2c -> Nrxn3`, four MEF2-motif enhancers with
concordant accessibility, Fig S5). The other three networks go to Fig S5; no network
is drawn for the sparse COP/NFOL and MFOL states. `Mef2c` is the hub and five curated
MOL3 markers (`Camk2a, Gabrb2, Kalrn, Nrxn3, Rnf150`) are among its targets.
"""))

A((CODE, r"""
# ---- 2d. the same graph as tables (Cytoscape; also the source of the marker counts) ----
for tag, (G, nt, f) in GRAPHS.items():
    e = nt["Edge"]["TF2G"]
    edges = pd.DataFrame({"source": e["TF"], "target": e["Gene"], "weight": e["TF2G_rho"],
                          "importance": e["TF2G_importance"], "regulation": e["TF2G_regulation"]})
    edges = edges.sort_values("importance", ascending=False).drop_duplicates(["source", "target"])
    # which eRegulon each TF node stands for, and its rho / peak state (from 2a)
    src2reg = obj.uns["eRegulon_metadata_filtered"].groupby("TF")["Region_signature_name"].apply(
        lambda s: ";".join(sorted(set(s.str.split("_(", regex=False).str[0])))).to_dict()
    edges["eRegulon"] = edges["source"].map(src2reg)
    first = edges["eRegulon"].str.split(";").str[0]
    edges["TF_rho"] = first.map(dict(zip(tab["eRegulon"], tab["Rho"])))
    edges["TF_peak"] = first.map(dict(zip(tab["eRegulon"], tab["peak_celltype"])))
    ntf, nge = nt["Node"]["TF"], nt["Node"]["Gene"]
    nodes = pd.DataFrame({"Node_ID": list(ntf["TF"]) + list(nge["Gene"]),
                          "Type": ["TF"] * len(ntf) + ["Target"] * len(nge)})
    nodes["is_marker"] = nodes["Node_ID"].isin(f["markers"]).map({True: "yes", False: "no"})
    nodes["n_TFs_targeting"] = nodes["Node_ID"].map(Counter(edges["target"])).fillna(0).astype(int)
    nodes["degree"] = nodes["Node_ID"].map(dict(G.degree())).fillna(0).astype(int)
    edges.to_csv(os.path.join(OUT, "grn_%s_cytoscape_edges.csv" % tag), index=False)
    nodes.to_csv(os.path.join(OUT, "grn_%s_cytoscape_nodes.csv" % tag), index=False)
    print(tag, len(edges), "edges,", (nodes.Type == "TF").sum(), "TF,",
          (nodes.is_marker == "yes").sum(), "markers in graph")
"""))

# ---------------------------------------------------------------- panel C ---
A((MD, r"""
## 3. Panel C — cross-technology reproducibility of eRegulon activity

**In:** the object (source of the signatures) + one Seurat RDS per other technology.
**Out:** `fig5c_crosstech_scatter.pdf` (two panels) and `eregulon_merged_auc.csv`.

The comparison must be *symmetric*, and that is the whole design:

| | dsciCAT-seq side | other-technology side |
|---|---|---|
| cells | all 5,977 lineage cells | OPC + Oligodendrocyte only (the labels those assays provide) |
| collapse to one vector | cell-count-weighted mean over the 6 states | cell-count-weighted mean over OPC / Oligo |
| per eRegulon | one mean AUC | one mean AUC |

* Signatures come from the dsciCAT run only — **gene-based eRegulons with `>= 3`
  target genes = 145**. The other technologies are *scored*, not re-inferred: this is
  assay-level reproducibility on the same donors, not cohort-independent validation.
* **`seed = 555` is mandatory.** `pyscenic.aucell` shuffles gene columns before
  ranking; in sparse matrices many genes tie at zero, so an unfixed seed moves rho by
  0.02-0.12 between runs. With the seed fixed, two runs are bit-identical.
* Agreement is read as **rank** correlation (Spearman): absolute AUC levels differ
  slightly between the two merged cell sets. The dashed identity line is a visual
  guide only.
* An earlier build carried a third panel, HT library `EL14` scored *without* Harmony
  (`rho = 0.80`, `p = 4e-33`). Same assay and donors as the Harmony arm, so it was
  dropped from the figure and kept in the caption; the three arms agree within 0.014.
"""))

A((CODE, r"""
# ---- 3a. R: export log1p RNA matrix + labels from each Seurat object ----
R_EXPORT = r'''
suppressPackageStartupMessages({library(Seurat); library(SeuratObject); library(Matrix)})
export_tech <- function(rds, tag, dir) {        # writes genes-x-cells .mtx + gene list + meta
  obj <- readRDS(rds)                          # RNA + ATAC assays, so Signac must be loadable
  x <- tryCatch(obj[["RNA"]]$data, error = function(e) NULL)          # Seurat v4 slot API
  if (is.null(x)) x <- tryCatch(SeuratObject::LayerData(obj, assay = "RNA", layer = "data"),
                                error = function(e) NULL)             # Seurat v5 layer API
  if (is.null(x)) x <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "data"),
                                error = function(e) NULL)             # v4 accessor
  if (is.null(x)) stop("no RNA data layer found in ", rds)
  if (length(x@x) == 0) stop("RNA data layer is empty in ", rds)
  x <- as(x, "CsparseMatrix")
  mtxf <- file.path(dir, paste0(tag, "_expr.mtx.gz"))
  writeMM(x, mtxf)          # PATH, never a gzfile() connection: writeMM to a connection
  if (file.size(mtxf) < 1000) stop("writeMM produced an empty matrix file: ", mtxf)
  writeLines(rownames(x), file.path(dir, paste0(tag, "_genes.txt")))
  write.csv(obj@meta.data, file.path(dir, paste0(tag, "_meta.csv")), row.names = TRUE)
  cat(sprintf("%s: %d cells x %d genes, nnz %d\n", tag, ncol(x), nrow(x), length(x@x)))
}
args <- commandArgs(TRUE)
export_tech(args[1], "tenx",       args[3])
export_tech(args[2], "HT_harmony", args[3])
'''
open(os.path.join(OUT, "export_tech_expr.R"), "w").write(R_EXPORT)
R = os.environ.get("FIG5_R")
if R is None:                                    # this box's own R is 4.1.3/Seurat 4.2.1 without
    cand = "/public/tools/r-packages/R-4.5.3/Rscript"   # Signac, and there is no network to add it,
    R = cand if os.path.exists(cand) else "Rscript"      # so the shared R 4.5.3 tree is the default
print("Rscript used:", R)
exported = all(os.path.exists(os.path.join(OUT, f)) and os.path.getsize(os.path.join(OUT, f)) > 1000
               for f in ("tenx_meta.csv", "HT_harmony_meta.csv",
                         "tenx_expr.mtx.gz", "HT_harmony_expr.mtx.gz"))
if not exported:
    p = subprocess.run([R, os.path.join(OUT, "export_tech_expr.R"),
                        os.path.join(DATA, "10x_subset_harmony.RDS"),
                        os.path.join(DATA, "HT_harmony.RDS"),
                        os.path.join(OUT)], capture_output=True, text=True)
    print((p.stdout or '').strip() or (p.stderr or '')[-500:])
    if p.returncode != 0:
        raise RuntimeError("R export failed (rc %d): %s" % (p.returncode, (p.stderr or "")[-500:]))
for f in ("tenx_expr.mtx.gz", "tenx_genes.txt", "tenx_meta.csv",
          "HT_harmony_expr.mtx.gz", "HT_harmony_genes.txt", "HT_harmony_meta.csv"):
    print("%-24s %8.1f MB" % (f, os.path.getsize(os.path.join(OUT, f)) / 1e6))
"""))

A((MD, r"""
**Expected** (only the lineage cells are scored; the rest of each assay is other
lineages)

```
tenx:       8438 cells | {'OPC': 298, 'Oligodendrocyte': 561}    # 4 isolation batches
HT_harmony: 14762 cells | {'OPC': 230, 'Oligodendrocyte': 347}   # 4 libraries
```
"""))

A((CODE, r"""
# ---- 3b. 145 signatures, AUCell on all three datasets (seed fixed) ----
from pyscenic.aucell import aucell
from ctxcore.genesig import GeneSignature

sig = {"%s (%s)" % (e.transcription_factor, e.cistrome_name): list(e.target_genes)
       for e in obj.uns["eRegulons"]}
sig = {k: v for k, v in sig.items() if len(v) >= 3}
print("gene-based eRegulons with >= 3 targets:", len(sig))
gss = [GeneSignature(name=k, gene2weight={g: 1.0 for g in v}) for k, v in sig.items()]


def mean_auc(expr, groups):
    return aucell(expr, gss, num_workers=min(24, os.cpu_count() or 1), seed=SEED
                  ).groupby(np.asarray(groups)).mean()


M = {"dsci": mean_auc(pd.DataFrame(obj.X_EXP, index=obj.metadata_cell.index,
                                   columns=obj.gene_names),
                      obj.metadata_cell[GROUP].astype(str))}
M["dsci"].to_csv(os.path.join(OUT, "dsci_6state_auc.csv"))
for t in ("tenx", "HT_harmony"):
    import scipy.io
    import gzip
    mtxp = os.path.join(OUT, t + "_expr.mtx.gz")
    # R writeMM() writes PLAIN MatrixMarket text when handed a .gz path (a gzfile()
    # connection is worse: it silently writes nothing). Read by actual magic bytes,
    # so the .gz name still works on environments that really did compress.
    is_gz = open(mtxp, "rb").read(2) == b"\x1f\x8b"
    with (gzip.open(mtxp, "rt") if is_gz else open(mtxp, "rt")) as fh:
        x = scipy.io.mmread(fh).tocsr()                           # genes x cells
    genes = open(os.path.join(OUT, t + "_genes.txt")).read().split()
    mt = pd.read_csv(os.path.join(OUT, t + "_meta.csv"), index_col=0)
    M[t] = mean_auc(pd.DataFrame(x.T.toarray(), index=mt.index, columns=genes),
                    mt["rna_celltype"].astype(str))
    M[t].to_csv(os.path.join(OUT, t + "_celltype_auc.csv"))
print({k: v.shape for k, v in M.items()})
"""))

A((CODE, r"""
# ---- 3c. collapse each side to one mean-AUC vector, weighted by cell numbers ----
n_dsci = obj.metadata_cell[GROUP].astype(str).value_counts()          # not hardcoded
vec = {"dsci_merged": sum(n_dsci[s] * M["dsci"].loc[s] for s in STATES) / n_dsci[STATES].sum()}
for t in ("tenx", "HT_harmony"):
    cnt = pd.read_csv(os.path.join(OUT, t + "_meta.csv"), index_col=0)["rna_celltype"] \
            .value_counts().reindex(["OPC", "Oligodendrocyte"])
    vec[t + "_merged"] = (cnt.values[:, None] * M[t].loc[cnt.index].values).sum(0) / cnt.sum()

V = pd.DataFrame(vec, index=M["dsci"].columns)[
    ["dsci_merged", "tenx_merged", "HT_harmony_merged"]].dropna()
V.to_csv(os.path.join(OUT, "eregulon_merged_auc.csv"))
print(V.shape, "| dsciCAT mean AUC %.3f-%.3f" % (V.dsci_merged.min(), V.dsci_merged.max()))
"""))

A((CODE, r"""
# ---- 3d. the panel: two scatters, shared axes, Spearman in the title ----
from scipy.stats import spearmanr

ARMS = [("tenx_merged", "10x Multiome"), ("HT_harmony_merged", "HT-scCAT-seq (Harmony)")]
top = V.sort_values("dsci_merged", ascending=False)
top = top[~top.index.str.split(" (", regex=False).str[0].duplicated()].head(5).index.tolist()   # one per TF
star = top + [i for i in V.index if i.startswith("Mef2c (")]                        # + the story TF

fig, axes = plt.subplots(1, 2, figsize=(9, 4.4), sharex=True, sharey=True)
lim = [0, float(V.values.max()) * 1.04]
for ax, (col, name) in zip(axes, ARMS):
    x, y = V["dsci_merged"].values, V[col].values
    rho, p = spearmanr(x, y)
    ax.scatter(x, y, s=14, alpha=0.7, linewidths=0)
    ax.plot(lim, lim, "k--", linewidth=0.7, alpha=0.5)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Mean AUC, dsciCAT-seq 6-state merged")
    ax.set_title("%s    rho = %.2f, p = %.0e" % (name, rho, p))
    for i in star:
        ax.scatter([V.at[i, "dsci_merged"]], [V.at[i, col]], facecolors="none",
                   edgecolors="#b32020" if i.startswith("Mef2c") else "#151515",
                   s=90, linewidths=1.4, zorder=5)
        ax.annotate(i.split(" (")[0], (V.at[i, "dsci_merged"], V.at[i, col]),
                    fontsize=7, xytext=(4, 4), textcoords="offset points")
fig.tight_layout()
out = os.path.join(OUT, "fig5c_crosstech_scatter.pdf")
fig.savefig(out); plt.close()
print("wrote", out)
for col, name in ARMS:
    print("%-24s n=%d rho=%.3f p=%.1e" % (name, len(V), *spearmanr(V.dsci_merged, V[col])))
"""))

A((MD, r"""
**Expected**

```
10x Multiome             n=145 rho=0.811 p=4.1e-35
HT-scCAT-seq (Harmony)   n=145 rho=0.807 p=1.7e-34
```

Both arms print `rho = 0.81` as on the panel. The rings mark the five strongest
eRegulons by dsciCAT mean activity, de-duplicated to one per TF
(`Tcf12, Nfix, Rfx3, Creb5, Mbnl2`), plus `Mef2c` in red because panel B is its
network — in the published figure those labels are placed by a radial greedy search
so that nothing collides; same information, more lines of code.

### Reading the published axes back (why the caption says "its own run")

This notebook scores dsciCAT **once** and shares that x axis between the two panels,
which is the only self-consistent choice. The published figure does not: each panel
used the dsciCAT vector produced by its own script run, and those four runs
(2026-09-22/23) predate the `seed = 555` fix, so they differ by up to 0.11 AUC
(97/145 rows differ between the two arms, max |delta| 0.034, median 0.0003).
Scored against each other instead, the same y data give:

| x axis taken from | rho vs 10x | rho vs HT |
|---|---|---|
| 10x run (published 10x panel) | **0.811** | 0.763 |
| HT run (published HT panel) | 0.784 | **0.807** |
| EL14 run | 0.778 | 0.771 |
| 10x batch-1 run | 0.804 | 0.716 |

So a clean re-run lands in the 0.72-0.81 band depending on the AUCell draw — which
is exactly the instability the fixed seed removes. The conclusion (rank agreement
about 0.8, `p` far below 1e-30, no arm systematically better) is unaffected, but a
re-run is not expected to reproduce 0.811/0.807 to the third decimal. The archived
per-run vectors (`*_dsci_6state_auc.csv`, `*_merged_auc_vec.csv`, dated 2026-09-22 to
10-09) reproduce them exactly if each panel is plotted against its own dsciCAT
column — that is the sense in which the figure caption says the x axis is "taken from
that comparison's own run".

At the (cell type x eRegulon) matrix level quoted in the text rather than drawn here,
correlations are 0.83 / 0.90 / 0.94 across comparison granularities, and MOL3_6 is
the lowest state in every one of them.
"""))

# ---------------------------------------------------------------- versions --
A((MD, r"""
## 4. Software
"""))

A((CODE, r"""
import importlib.metadata as md
mods = ["numpy", "pandas", "scipy", "matplotlib", "scanpy", "anndata", "pyscenic",
        "scenicplus", "ctxcore", "plotnine", "dill", "networkx", "statsmodels"]
print("python", sys.version.split()[0])
out = []
for m in mods:
    try:
        v = md.version(m)
    except Exception:
        v = "-"
    out.append((m, v)); print("%-14s %s" % (m, v))
with open(os.path.join(OUT, "requirements.txt"), "w") as fh:
    fh.write("python>=%s\n" % sys.version.split()[0])
    for m, v in out:
        if v != "-":
            fh.write("%s==%s\n" % (m, v))
print("wrote %s/requirements.txt" % OUT)
"""))

A((MD, r"""
Versions used for the published figure (SCENIC+ environment, `pip list`):

| Package | Version | Package | Version |
|---|---|---|---|
| python | 3.9 | **scenicplus** | **1.0.1.dev4+ge4bdd9f** |
| numpy | 1.22.4 | pyscenic | 0.12.1+6.g31d51a1 |
| pandas | 1.5.1 | ctxcore | 0.2.0 |
| scipy | 1.10.1 | plotnine | 0.12.4 |
| matplotlib | 3.6.2 | dill | 0.3.6 |
| scanpy | 1.9.2 | h5py | 3.8.0 |

Panel C's export step needs R with `Seurat` + `Matrix`; the snippet handles both the
v4 (`obj[["RNA"]]$data`) and v5 (`LayerData`) slot APIs, so either major version
works. Everything else is optional — the notebook prints its own
`requirements.txt`, which is the authoritative list for whatever you ran it with.

Upstream references (not re-run here): mm39, GENCODE v45, `v10nr_clust_public`
mouse motif collection, SCENIC+ `run_scplus` defaults
(`min_target_genes = 10`, quantiles `(0.85, 0.90, 0.95)`, `top_n = (5, 10, 15)`,
BASC, `rho_threshold = 0.05`, extended kept, merged).

**Fixed seeds:** AUCell and pseudobulk `555`; GRN layout is deterministic
(`kamada_kawai_layout`).

**Deliberately absent:** panel compositing, Arial substitution, panel letters,
collision-free label placement, and the marker-gene closure bar chart that was
removed from the main figure (its numbers are in the caption).
"""))

# ------------------------------------------------------------------- emit ---
nb = {"cells": [cell(k, s) for k, s in C],
      "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python",
                                  "name": "python3"},
                   "language_info": {"name": "python", "version": "3.9"}},
      "nbformat": 4, "nbformat_minor": 5}
here = os.path.dirname(os.path.abspath(__file__)) or "."
p = os.path.join(here, "Figure_5.ipynb")
with open(p, "w") as fh:
    json.dump(nb, fh, indent=1, ensure_ascii=False)
    fh.write("\n")
nc = sum(len(s.splitlines()) for k, s in C if k == CODE)
print("wrote %s | %d cells | %d code lines" % (p, len(nb["cells"]), nc))
