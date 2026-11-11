# Figure 5 — SCENIC+ eRegulons of the oligodendrocyte lineage

`Figure_5.ipynb` reproduces the content of the three Fig 5 panels from the analysis
objects. Journal assembly (compositing, Arial substitution, panel letters,
collision-free label placement) is done by hand and is not part of this notebook.

| Panel | What it shows | Notebook section |
|---|---|---|
| A | eRegulon activity (dot size) x target expression (cell colour), 20 eRegulons x 6 lineage states | §1 |
| B | MOL3 eRegulon–target network, 44 nodes / 44 edges, *Mef2c* hub | §2 |
| C | Cross-technology reproducibility of eRegulon activity, 145 eRegulons, 10x Multiome and HT-scCAT-seq | §3 |

## Inputs (`FIG5_DATA`, default `./data`)

| File | Produced by | Used by |
|---|---|---|
| `scplus_obj_final.pkl` | SCENIC+ on the dsciCAT-seq oligodendrocyte-lineage cells — full chain in §0 of the notebook, compressed below | A, B, C |
| `10x_subset_harmony.RDS` | 10x Multiome, Harmony-integrated, 4 isolation batches | C |
| `HT_harmony.RDS` | HT-scCAT-seq, Harmony-integrated, 4 libraries | C |
| `NIHMS70326-supplement-Table_S1.xlsx` | Marques et al. 2016, 24 states × 50 up-regulated genes | B |

## Where `scplus_obj_final.pkl` comes from

**It is not a KLP table, and no KLP was ever produced here.** This SCENIC+ build
(1.0.1.dev4) has no `run_klp` step: region→gene linking is a GBM classifier inside
`run_scenicplus`, and the TF→gene edges come from the gsea-based `build_grn`. Compressed
chain (parameters, file names and dated products are in notebook §0, table of 14 steps):

1. **Seurat/Signac RDS** of the dsciCAT oligo subset (`3k_dsci_oligo_harmony.RDS` → 5,977
   lineage cells) — *not* CellRanger-ARC, *not* ArchR — exported from R to two h5ad files:
   RNA (all genes) and the **binarised peak matrix** (`atac_assay@data > 0 → 1`).
2. **cisTopic** `run_cgs_models(n_topics = [5,20,40,60,80,100,130,150,180], n_iter = 150,
   random_state = 555, alpha = 50, eta = 0.1)`; 20/40/60 compared on per-topic coherence →
   **60 topics** kept, then `binarize_topics` three ways (`otsu`, `ntop 3000`, `li`) plus
   per-state DARs (`find_marker_features`).
3. **pyCisTarget**: cistrome regions → mm10 FASTA (`bedtools getfasta`, `1kb_bg_padding.fa`),
   enrichment against `v10nr_clust_public` (`run_pycistarget(..., run_without_promoters = True,
   n_cpu = 40, annotation_version = 'v10nr_clust')`) → `CTX_topics_li_*`, `CTX_topics_top_3k_*`,
   `CTX_DARs_*` (60 `TopicN.html` each) and `menr.pkl`.
4. `create_SCENICPLUS_object(GEX_anndata = adata.raw.to_adata(), cisTopic_obj = …cisTopic60.pkl,
   menr = menr.pkl)` → `run_scenicplus(variable = ['ACC_celltype_l2'], species = 'mmusculus',
   assembly = 'mm10', tf_file = allTFs_mm.txt, upstream = downstream = [1000, 150000], n_cpu = 70)`
   → **only 10 eRegulons, 0 true repressors**, which is why the GRN is rebuilt.
5. `build_grn` → `format_egrns` → `make_rankings` (region + gene) → `score_eRegulons(auc_threshold
   = 0.05)` → pseudobulk `TF_cistrome_correlation` → std filter (`|rho| > 0.1`, > 5 targets) →
   `regulon_specificity_scores` → **`scplus_obj_final.pkl`** (8.0 GB, 181 eRegulons).

## Run

```bash
export FIG5_DATA=/path/to/objects FIG5_OUT=fig5_out
export FIG5_R=/public/tools/r-packages/R-4.5.3/Rscript   # optional, this is the default if present
jupyter nbconvert --to notebook --execute --inplace --allow-errors Figure_5.ipynb
```

Panels A and B need the SCENIC+ environment (`scenicplus`, `pyscenic`, `ctxcore`,
`plotnine`, `dill`); §3a additionally needs an **R with Seurat and Signac** — those two public
RDS objects carry a `ChromatinAssay`, and R refuses to open them without Signac. The `R` on
`PATH` here is R 4.1.3 / Seurat 4.2.1 with **no Signac and no network to install it**, so the
notebook defaults to the shared R 4.5.3 tree (Seurat 5.5.1 / SeuratObject 5.4.0 / Signac 1.17.1);
`FIG5_R` overrides it. The rest of §3 is pandas/scipy/matplotlib only. `--allow-errors` matters:
nbconvert 7 otherwise discards every cell output if one cell fails. `requirements.txt` is written
by the last cell.

## Outputs (`FIG5_OUT`)

`fig5a_eregulon_dotplot.pdf` · `fig5b_grn_{mol3,mol2,mol6,mol3_6}.pdf` ·
`fig5c_crosstech_scatter.pdf` · `eregulon_peak.csv` ·
`grn_*_cytoscape_{edges,nodes}.csv` · `dsci_6state_auc.csv` ·
`{tenx,HT_harmony}_celltype_auc.csv` · `eregulon_merged_auc.csv`

## Two parameters that decide whether the numbers reproduce

* **`seed = 555` for AUCell and pseudobulk.** `pyscenic.aucell` shuffles gene
  columns before ranking; with sparse matrices many genes tie at zero, and an
  unfixed seed moves the panel-C Spearman rho between 0.72 and 0.81.
* **One shared dsciCAT axis.** This notebook scores dsciCAT once and plots both
  panels against that vector. The published figure used each arm's own run
  (`rho = 0.811` / `0.807`), from runs that predate the seed fix; see §3 of the
  notebook for the measured difference and for why the caption says the x axis is
  "taken from that comparison's own run".

## Upstream references

**mm10** for everything in this figure - asserted, not assumed: the four Mef2c-eRegulon
regions that link to *Nrxn3* (chr12:88,707,581 / 88,952,720 / 89,284,141 / 89,376,068) only
make sense on mm10 (on mm39 *Nrxn3* sits near 56 Mb), `run_scenicplus(assembly = 'mm10')`, and
the ArchR subset that supplies the coverage tracks logs `BSgenome.Mmusculus.UCSC.mm10`.
Gene annotation for `species = 'custom'` came from the mm10 `genes.gtf`.

`v10nr_clust_public` mouse motif collection (`motifs-v10-nr.mgi-m0.00001-o0.0.tbl`,
`allTFs_mm.txt`) · SCENIC+ 1.0.1.dev4 wrapper defaults for the region-gene GBM step
(`min_target_genes = 10`, quantiles `0.85/0.90/0.95`, `top_n = 5/10/15`, BASC,
`rho_threshold = 0.05`, extended regulons kept and merged) · `kamada_kawai_layout` for the
networks · Spearman for all panel-C correlations. (GENCODE v45 / hg19-hg38 belong to the
GWAS side, Fig 6 / Fig S6 - not to this figure.)

## What changed on the published figure (2026-10-29, reviewer round 2)

| Item | Why | Script |
|---|---|---|
| A graphical key for the network was added to panel B, as **one row along the top of the panel, level with the panel letter**: `◇ TF ──▸ ◇ target gene  ─ 1.0×  ▬ 4.0×` | the panel had no key at all; the authors then asked for the explanatory sentences to live in the legend instead of the figure | `assembly/network_key.py` (`draw_key` = the multi-row form used in Figure S5B, `draw_row` = this single row), `assembly/fix_fig5_v8.py`, `assembly/fix_fig5_v9.py` |
| Everything the key does *not* say is defined in the legend: edge colour = source eRegulon, gene-node colour = subtype log2FC, "all retained links are predicted activators", and the printed width range | house rule: minimal in-figure text | manuscript legend (2026-10-30 revision round) |

The key's content was **measured from the delivered vectors**, not read off the plotting config —
`assembly/network_key.py`'s module docstring records the audit. Two config options that look like
they are in use are not: `transparency_edge_by` never took effect (no ExtGState in the output, edge
colour is plain tab10), and `shape_node_by` never took effect (every node is a diamond). The graph is
also **undirected** (`from_pandas_edgelist`, no arrows), which is why direction is keyed by node type
rather than by an arrowhead.

`assembly/fix_fig5_v8.py` / `fix_fig5_v9.py` paste the key onto the delivered page; the second round
starts from the 2026-10-21 original, which never had a key, so nothing has to be erased
(`apply_redactions` on a page holding nested Form XObjects is a hazard — see
`../00_Figure_Tools/README.md`). `assembly/finalize_fig5_s5.py` installs both figures into the
delivery block and re-exports the PNGs.
