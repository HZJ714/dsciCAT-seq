# upstream/ — how the inputs of `Figure_6.ipynb` are produced

The notebook starts from per-nucleus scDRS scores. This directory contains the chain that turns the
count matrices (CNGBdb **CNP0010263**) into those scores, in the order the scripts were run. Nothing
here is needed to *read* the notebook — it is needed to re-derive its inputs on a different machine.

Runtime of the whole chain: the export steps are minutes; the scDRS step is **444 runs**
(74 traits × 2 channels × 3 chemistries) and is the only expensive part (run it with `NP=8` on 8
cores; it took ~1 day on the author's HPC allocation, resumable — a trait already scored is skipped).

## Step 0 — where the two own-data `.h5ad` files come from

`Figure_6.ipynb` reads two h5ad files that are *not* scDRS products; they are converted from the
author's own Seurat/Signac RDS objects. Both conversions are here, and both were re-run from the RDS
and compared element-wise with the delivered files: **X (`indptr`/`indices`/`data`), `obs` and `var`
identical** in each case.

| file | source RDS | content | how |
|---|---|---|---|
| `FY/oligo/3k_dsci_oligo_subset.h5ad` | `FY/oligo/3k_dsci_oligo_subset.RDS` (794 MB, the Fig 2 lineage object) | 5,977 nuclei × 31,263 genes; `obs.cell_type` = the six subtypes (MOL6 2,574 / MOL2 1,016 / MOL3-6 1,006 / MOL3 688 / MFOL 490 / COP-NFOL 203); UMAP inherited, never recomputed | `rds_to_dump.R` → `build_h5ad.py` |
| `…/dualomics_probe/genescore_oligo_5977.h5ad` | `<online>/FY/3k_dsci_ht_10x_gs_subset.RDS` (13.2 GB, Signac multiome object) | same 5,977 nuclei in the same order × 24,333 genes; `@assays$gene_score` `data` slot | `export_genescore.R` → `build_genescore_h5ad.py` |

`D0_own_lineage_h5ad.sh` runs all four steps in order. The nucleus list (`lineage_cells_5977.txt`) is
exactly `obs_names` of the RNA object — the same 5,977 nuclei in both channels. The six subtype labels
in the chromatin file are **joined from the RNA object by barcode** (the RDS's own `cell_type` column
is the 11-class broad label — AST/EC/Ex_*/In_D1/MG/OPC/Oligo/… — and is not what the figure uses).

Slot semantics, which is the whole point of the conversion: the RNA matrix is the Seurat **`data`**
slot = `log1p(counts / colsum × 1e4)` (max 8.02, and `expm1(X)` sums to exactly 1e4 per nucleus), not
the `counts` slot — `rds_to_dump.R` writes both (`matrix.mtx`, `matrix_data.mtx`) and `build_h5ad.py`
prefers `matrix_data.mtx` and says so on stdout. The chromatin matrix is the Signac gene-score
**`data`** slot = `log1p(ArchR gene score)` (max 6.47, `expm1(X)` sums to ~4,630 per nucleus).

## ⚠ Input transform of the scDRS matrices — found by audit, fixed in code (2026-11-04)

`scdrs compute-score` is called everywhere with `--flag-raw-count False`, and that branch applies
`expm1()` to the input matrix internally (`scdrs/pp.py:613`): it therefore **requires log-transformed
input**. Broken, this promise does not raise — it prints a numpy `RuntimeWarning` into the per-batch
`compute.log` and still writes score files, which is why every downstream consistency gate kept
passing. `A6_fig6_input_transform_audit.py` walks every scDRS tree Figure 6 uses, recovers the exact
`.h5ad` its launcher named, measures the scale, and counts warning lines:

| panel / tree | matrix actually fed | max | median per-nucleus `expm1` | warn jobs/logs | verdict |
|---|---|---|---|---|---|
| A/B `scdrs_ms_out` | `3k_dsci_oligo_subset.h5ad` (Seurat `data`) | 8.02 | 1e4 | no logs retained | correct |
| C RNA `full74` | same | 8.02 | 1e4 | 0 / 25 | correct |
| **C chromatin `out_gs74`** | `genescore_oligo_5977.h5ad` (Signac `gene_score` **data** = log1p) | **6.47** | 4,630 | **0 / 25** | **correct — panel C never had this error** |
| D `dsciCAT_seq_rna` | measured RNA, logged | 8.02 | 1e4 | 0 / 25 | correct |
| D `10x_rna` | same convention | 6.73 | 1e4 | 0 / 25 | correct |
| D `HT_scCAT_seq_HT_rna` | same convention | 8.69 | 1e4 | 0 / 25 | correct |
| D `dsciCAT_seq_atac` | **ArchR `GeneScoreMatrix`, un-logged** | **1,388.6** | 3.7e34 | **31 / 31** (`overflow in expm1` → inf, then `overflow in multiply`) | **wrong** |
| D `10x_atac` | same, un-logged | 297.1 | 3.6e34 | 0 / 25 (silent) | **wrong** — `expm1(297) ≈ 1e129` |
| D `HT_scCAT_seq_HT_atac` | same, un-logged | 372.7 | 2.2e35 | **25 / 25** (`overflow in multiply`) | **wrong** |
| (abandoned) GIM RNA | ArchR `GeneIntegrationMatrix` | — | — | 0 / 26 | rejected by V1, never used |

Table = `tables/fig6_input_transform_audit_20261104.csv`. So the defect is confined to the chromatin
half of panel D; panels A/B/C are on the correct scale. Both candidate fixes converge on the same
matrix, because the raw gene scores are already 1e4-normalised: `--flag-raw-count True` would make
scDRS normalise (a no-op) and log them, and pre-applying `log1p` with the flag left at `False`
produces the identical input. We chose the second, because it keeps one rule for all six arms —
**every matrix Figure 6 scores is log-scale, and `--flag-raw-count False` never changes**.

### The fix, in code

| script | what it does |
|---|---|
| `D4b_pack_atac_log1p.py` | builds `h5ad_log1p/<tech>_atac.h5ad` from the delivered matrices (`X.data → log1p(X.data)`), asserts the `expm1` round trip (measured max abs error 4.5e-13 / 1.1e-13 / 8.5e-14), obs/var identity, and copies the cov files. The delivered un-logged files are left in place as the record of what was submitted. |
| `preflight_h5ad.py` | the gate that was missing: for a file about to be passed to `--h5ad-file`, checks `max(X) ≤ 30` and `n_inf(expm1(X)) == 0` (or, with `--raw-count`, that raw counts *are* large). On the three delivered chromatin files it prints `REJECT` (`expm1_max = inf / 1.06e129 / 7.54e161`); on the fixed ones, `OK` (`expm1_max = 1390 / 297 / 373`, 0 inf). |
| `D5_run_scdrs.sh` | now runs that preflight before launching a job and refuses the job if it fails, and post-checks each job's `compute.log` for `overflow encountered` plus batch/score.gz completeness. Also `NICE`/`IONICE` (this box has other tenants) and overridable `H5` / `SCDRS_OUT` / `LOG` so a corrected rerun can live in a parallel tree. |
| `D6_crosschem_analysis.py` | `DSCI_SCDRS_SUBDIR` / `DSCI_SCDRS_SUBDIR_ATAC` / `DSCI_TAB_SUBDIR` analyse a rerun into `tables_*` without touching the delivered tables. Only the atac arms were re-run, so `DSCI_SCDRS_SUBDIR_ATAC=scdrs_log1p_20261104` points the chromatin channel at the corrected tree while the expression channel keeps reading the delivered `scdrs/` tree. |
| `A6_fig6_scores_integrity.py` | gate before D6: every batch dir vs its `.gs` trait list — per-trait `*.score.gz` present, gzip readable, `norm_score` finite, row count == h5ad nuclei (or the single within-arm mode). Exit 2 lists dirs to delete so the D5 resume logic re-runs them. |

### Measured effect of the correction (probe, MAGMA batch 5 = 3 traits × 3 chemistries)

`V3_input_transform_probe.py`, results `V3_input_transform_probe.csv`; score files kept at
`D_crosschemistry/scdrs_V3_log1p_probe_20261104/`. Δ1 = median Oligo − median OPC.

| chemistry | trait | per-nucleus Spearman old vs new | Δ1 delivered | Δ1 with `log1p` input |
|---|---|---|---|---|
| dsciCAT-seq (own) | MS | **0.464** | **+0.262** | **−0.151** (sign flip) |
| dsciCAT-seq (own) | PBC / ReactionTime | 0.438 / 0.449 | +0.672 / −1.017 | +0.679 / −0.429 |
| 10x Multiome | MS | 0.914 | −0.430 | −0.474 |
| 10x Multiome | PBC / ReactionTime | 0.903 / 0.885 | +0.254 / −0.204 | +0.291 / −0.253 |
| HT-scCAT-seq | MS | **0.373** | **−0.297** | **+0.064** (sign flip) |
| HT-scCAT-seq | PBC / ReactionTime | 0.332 / 0.289 | +0.118 / −0.761 | +0.631 / +0.353 (flip) |

Two of the three own-arm and HT-arm contrasts change sign; 10x barely moves. That is consistent with
the magnitude of the raw matrix: the own arm contains entries above 709 (→ `inf`), HT has a long tail
out to 373, and the 10x gene scores are the flattest of the three.

What is affected is **only** the chromatin half of panel D: its 74 × 3 ranks, off-diagonal z, Δ1, the
ladder ρ_adj (delivered values 0.726 / 0.763 / 0.518) and therefore the "2 of 3 chemistries" verdict.
The rerun is **222 chromatin runs** (74 traits × 3 chemistries); the 222 expression runs and all of
panel C are untouched. Measured on this machine: HT ≈ 95 s per 3-trait batch at 0.35 GB RSS per
process, own ≈ 5 min per batch, so the full rerun is well under an hour at `NP=6`.

Status (2026-11-04): rerun **complete** into `scdrs_log1p_20261104/` (222 runs, 0 overflow warnings; integrity gate passed: 74/74 traits per arm, row counts exactly match h5ad). D6 re-analysis landed in `tables_log1p_20261104/`.

### Corrected rerun result (before → after, atac channel only)

| chemistry | ρ_adj old → new | perm P | verdict old → new | MS row rna→atac rank old → new |
|---|---|---|---|---|
| dsciCAT-seq (own) | 0.726 → **0.896** | 1.0e-04 | PASS → PASS | 7→8 → 7→37 |
| 10x Multiome | 0.763 → **0.824** | 1.0e-04 | PASS → PASS | 4→55 → 4→49 |
| HT-scCAT-seq | 0.518 → **0.630** | 1.0e-04 | fail → **PASS** | 11→47 → 11→41 |

The **"2 of 3 chemistries" verdict becomes 3 of 3**; the MS counterexample (chromatin ranks MS lower than expression) is still reproduced in all three. Per-trait turnover in the atac Δ1/ranks: own 19/74 sign flips (rank Spearman 0.73), 10x 6/74 (0.93), HT 27/74 (0.38) — the delivered atac half was materially wrong, which is exactly why the input transform had to be fixed. Full per-trait table: `D_crosschemistry/tables_log1p_20261104/fig6D_input_fix_before_after.csv`.

Author decision (2026-11-04): the delivered figure is **re-rendered with the corrected tables** and
**tentatively finalized** —
`fig6D_heatmap_pair_compact.pdf`, `fig6D_heatmap_pair.pdf`, `fig6_composed_v21.pdf`, the S6m
scatter, the S6g/S6h support panels and the Table S6/S7 ledger all draw from the corrected tree;
`A6_fig6_v21_verify.py` PASSES end-to-end (its legend expectations are derived from the corrected
MS/ladder tables, not hardcoded). `make_all_current.sh` defaults to `h5ad_log1p` + `STRICT=1` +
`TABD=tables_log1p_20261104`. **Notebook**: `Figure_6.ipynb` was regenerated from
`build_fig6_notebook.py` (defaults point at the corrected tree; the §0.1 "provisional" caveat and the
§4 "un-logged" wording are rewritten) and re-executed — its panel-D numbers match the tables
digit-for-digit (ladder 0.896/0.824/0.630, MS ranks 7→37 / 4→49 / 11→41, offdiag z 2.48/2.62/0.89),
so two independent implementations (D6 and the notebook) now agree. **Independent audit**: 37/37
checks pass (log1p bit-identical, score trees finite, Δ1/ρ_raw recomputed from the raw score trees
match the tables, RNA half bit-identical, probe cross-validation exact for 10x/HT and within 0.002
for own). The manuscript Results sentence "2 of 3 chemistries" is the author's to change (draft
replacement text in the legend file).

## Steps

| Step | Script | Input | Output | Key parameters |
|---|---|---|---|---|
| D0 | `D0_own_lineage_h5ad.sh` (runs `rds_to_dump.R`, `build_h5ad.py`, `export_genescore.R`, `build_genescore_h5ad.py`) | the two own RDS objects above | the two own `.h5ad` files | verified bit-for-bit against the delivered files; see Step 0 above |
| D1 | `D1_export_genescore_archr.R` | Fig 2 ArchR project (`DSCI_ONLINE_ROOT/FY/archr-subset-copnfol`) | `derived/gs_<tech>.mtx.gz` + `_genes.txt`, `tables/cells_<tech>.csv` | `getMatrixFromProject(useMatrix="GeneScoreMatrix")`; column names stripped of the `sample#` prefix; only OL-lineage nuclei (`Oligo`, `OPC`, `COP/NFOL`) kept; per-nucleus sum = 1e4 asserted |
| D2c | `D2c_export_cov.R` | Fig 5 meta (10x / HT) + the own multiome RDS `@meta.data` | `tables/cov_by_cell.csv` | 8,948 OL-lineage nuclei, **missing 0**; covariate column names and count copied from the C analysis so that D is not scored with a differently-specified model |
| D2e | `D2e_export_rna_seurat.R` | `FY13/FY14 3k_*_multiome.RDS` | `derived/rnaSeurat_dsciCAT_FY1{3,4}.mtx` | `GetAssayData(obj, assay="RNA", slot="data")` — **not** `LayerData()` / `$data` (see traps) |
| D2f | `D2f_pack_rna_seurat.py` | D2e product + Fig 5 `tenx_expr.mtx.gz` / `HT_harmony_expr.mtx.gz` | `derived/rna_<tech>.npz` + `_genes.txt`, `tables/rna_seurat_sources.csv`, `rna_common_genes.csv` | all three channels cut to the **gene intersection (30,559)** so "this technology measures more genes" cannot become a cross-technology confound |
| D4 | `D4_build_h5ad.py` | D1 + D2f + D2c | `h5ad/<tech>_<channel>.h5ad` + `.cov.tsv` (6 + 6) | one h5ad per technology per channel, nuclei aligned to `cells_<tech>.csv`. **Its chromatin matrices are un-logged; use D4b before scoring them (see the input-transform section).** |
| D4b | `D4b_pack_atac_log1p.py` | `h5ad/<tech>_atac.h5ad` | `h5ad_log1p/<tech>_atac.h5ad` (+ cov) | `X.data → log1p(X.data)`; asserts `expm1` round trip < 1e-8, obs/var identity, `max < 30`; then runs `preflight_h5ad.py` |
| — | `preflight_h5ad.py` | any `.h5ad` about to be scored | pass/REJECT + exit code | the scale promise of `--flag-raw-count False`, checked instead of assumed |
| — | `A6_fig6_input_transform_audit.py` | every scDRS tree + its launcher's `.h5ad` | `tables/fig6_input_transform_audit_20261104.csv` | panel-by-panel verdict; proves A/B/C were scored on correct inputs |
| D5 | `D5_run_scdrs.sh` (+ `D5b_repair_atac6.sh`) | D4 + `magma_10kb_1000.74_traits.gs.batch/` | `scdrs/<tech>_<channel>/<batch>/<trait>.score.gz` | command below, parameters copied verbatim from the own-data run so no difference between chemistries can come from the scorer |
| V1 | `V1_rna_source_scale.py` | ArchR `GeneIntegrationMatrix` vs Seurat `RNA$data` | `tables/V1_rna_source_scale.csv`, `V1_per_cell_rho.csv` | the negative result that decided the RNA channel (below) |
| D2d | `D2d_read_rna_from_arrows.py` | ArchR arrow files | `derived/rnagim_*` (demoted) | block-wise arrow reader used to *check* the GIM; rebuilds the GeneScoreMatrix and compares element-wise with D1: **max abs diff 0** |
| O1 | `O1_own_same_nucleus.py` | own 74-trait scores, both channels | `O1_own_multiome_paired/tables/T1…T7` | the own-data (panel C) cross-channel tests; Δ1 = median(MOL6) − median(COP/NFOL) |
| D6 | `D6_crosschem_analysis.py` | D5 scores + `cells_<tech>.csv` | `tables/fig6D_*` | the delivered producer of the panel-D tables; `Figure_6.ipynb` §4 reimplements the same criteria in one place and reproduces these files bit-for-bit |
| D8 | `D8_fig6C_disclosure.py` | O1 tables | `FIG6_v17/tables/fig6C_fourline_disclosure.csv`, `fig6C_sixpoint_all74.csv` | the four disclosure layers quoted in the Fig 6 legend (ladder / MS off-diagonal / six-subtype / per-nucleus) |
| A5 | `A5_fig6_v15_stats.py` | own h5ad + MS scDRS output | `FIG6_v15/tables/fig6A_stats.csv` … | the delivered producer of the panel-A statistics |

## The scDRS run (D5), verbatim

```bash
scdrs compute-score --h5ad-file <tech>_<channel>.h5ad --h5ad-species mouse \
    --gs-file <trait>.gs --gs-species human --out-folder <out>/<trait> \
    --flag-filter-data True --flag-raw-count False --n-ctrl 1000 \
    --flag-return-ctrl-raw-score False --flag-return-ctrl-norm-score True \
    --cov-file <tech>_<channel>.cov.tsv
```

The gene sets are human (MAGMA 10 kb, 1,000 genes, 74 traits) and the cells are mouse, exactly as in
the own-data analysis; `n_ctrl = 1000`, the covariate file and the runner are identical for all three
chemistries, so a difference between chemistries cannot be a difference in scoring.

## V1: why the expression channel is *measured* RNA, not ArchR's `GeneIntegrationMatrix`

`GeneIntegrationMatrix` is a CCA integration/imputation product, not measured expression. Measured on
the same 5,957 nuclei × 22,112 genes:

| evidence | Seurat `RNA$data` | `GeneIntegrationMatrix` |
|---|---|---|
| non-zero fraction | **6.0 %** | **48.3 %** |
| per-nucleus sum | median 2,300 (log-sum, naturally ragged) | **exactly 1e4** (constructively normalised) |
| per-nucleus cross-source Spearman | — | median **0.174**, while a *scrambled* null (mismatched nuclei) already gives **0.128** |
| myelination gradient | Plp1 COP 3.25 / OPC 3.50 / **Oligo 3.94** (monotonic) | Plp1 COP 45.7 / **OPC 117.2** / Oligo 109.8 (**OPC overtakes mature OL**) |
| Mbp | COP 2.99 / OPC 1.03 / Oligo 2.93 | COP 55.8 / OPC 50.0 / Oligo 54.2 (gradient flattened) |

Using the GIM as the "expression channel" would push chromatin-space information into the expression
channel and constructively inflate the cross-channel agreement — circular — and would also invert the
maturation gradient itself. The rejected GIM matrices remain in `derived/rnagim_*` and some scores in
`scdrs/ABANDONED_dsciCAT_seq_rna_GIM/`; they are evidence only and enter no figure.

## Traps that cost time here

1. **Seurat v4 vs v5 layer access.** This environment has SeuratObject 4.1.3 / Seurat 4.3.0, where
   `LayerData()` and `obj[["RNA"]]$data` do not reach `Assay@data`. `GetAssayData(assay="RNA",
   slot="data")` does. The earlier "this object has no RNA layer" conclusion was an API mistake, not
   missing data.
2. **ArchR arrow sparse blocks.** `jLengths` is the nnz per column, `jValues` are the **1-based global
   column numbers of non-empty columns** (empty columns are not stored, so the array can be shorter
   than the number of nuclei); row indices are 1-based *within* the block; block row order = the
   chromosome's genes by `FeatureDF$idx`, global row order = alphabetical by chromosome name — which
   is why it matches `getMatrixFromProject` (verified element-wise).
3. **Do not pull whole matrices into R.** `getMatrixFromProject` on 22k × 70k reached 28 GB of
   anonymous memory against a 32 GB limit and was killed (`D2b_export_rna_archr.R` is kept as the
   record of that failure); the block-wise arrow read needs 4–5 GB.
4. **`writeMM` / `mmwrite`** must write to a path, not a `gzfile` connection; and the Fig 5
   `*_expr.mtx.gz` files are *uncompressed* despite the name — check the magic bytes before opening.
5. **`pkill -f <pattern>` in a script matches the script's own command line** and kills the rest of
   it; use a pid file or a bracketed pattern.
6. **`--flag-raw-count False` is a promise about the input, not a switch about robustness.** scDRS
   `expm1()`s the matrix in that branch, so an un-logged matrix is silently accepted and then
   destroyed (`expm1(1388)` → `inf`, and `inf` in a control-cell mean propagates to that control's
   whole gene set). This is now enforced rather than remembered: `D5_run_scdrs.sh` will not launch a
   job until `preflight_h5ad.py` accepts its matrix, and it post-checks every batch log for
   `overflow encountered`. `A6_fig6_input_transform_audit.py` re-derives the whole table above from
   disk if you ever want to re-verify it.
