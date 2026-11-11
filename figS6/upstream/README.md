# upstream/ — how the inputs of `Figure_S6.ipynb` are produced

The notebook starts from delivered tables; this directory holds the scripts that created
them. Nothing here is needed to *run* the notebook — only to re-derive its inputs on a
different machine. Runtime note: the scDRS sweeps run under
`$DSCI_PY` (the analysis environment: python 3.9.23, **scdrs 1.0.2**; export it as described in
the top-level README, 《Data layout and environment variables》); the drawing
notebook does not need scdrs.

| script | produces | consumed by |
|---|---|---|
| `ms_locus_overlap.py` | `FY/oligo/scdrs_ms_out/ms_locus_overlap_table.csv` (lead-SNP ± 500 kb windows × GENCODE TSS) | Panel A, Table S6 locus block |
| `O1_own_same_nucleus.py` | `O1_own_multiome_paired/tables/` — `cell_level_table.csv` (per-nucleus RNA/ATAC MS scores, 5,962 nuclei), `T4_delta1_74traits_both_channels.csv` (per-trait Δ1/Δ1b, narrow + wide contrast), T1/T2/T3 | Panels C, F, G; Table S6 |
| `O1b_own_by_library.py` | `L1_by_library.csv` (per-library cross-channel rho, FY13/FY14) | library rungs provenance |
| `M2b_formal_authoranno.py` | `M2b_authoranno_v1/tables/M2b_per_mouse.csv` — external mouse EAE **GSE250589**, author annotation, per-mouse axis contrasts + split-half AUC (26 screened, 18 delivered; M2b rule MOL < 20 or OPC < 20 nuclei) | Panel F mouse rung; Table S7 |
| `A6_fig6_input_transform_audit.py` | audit of the **log1p-corrected** ATAC rerun (delivered 2026-11-04 by the data owners after the raw-count deflation was found): channel-wise input-transform integrity across `D_crosschemistry/scdrs_log1p_20261104/` | Panels E, F, G (ATAC side); the "log1p-corrected" data version |
| `D5_run_scdrs.sh` (+ `D5b_repair_atac6.sh`, `D5c_mem_guard.sh` in the authors' tree) | the scDRS batch sweep itself: per trait × per arm `*.score.gz` / `PASS_*.full_score.gz` trees under `D_crosschemistry/scdrs/` (RNA; ATAC rerun delivered 2026-11-04). Runs under `$DSCI_PY` (scdrs 1.0.2); interrupted batches must pass the integrity check before resuming | Panels C, F, G (raw per-nucleus scores) |
| `A6_fig6_scores_integrity.py` | completeness gate for the score tree (one `*.score.gz` per trait, gzip readable, row count == arm h5ad, finite norm_score) — protects against silently half-written batches being skipped by the resume logic | gate before any panel using the score trees |
| `D6_crosschem_analysis.py` | `D_crosschemistry/tables_log1p_20261104/fig6D_delta1_74traits_per_tech.csv` (per-chemistry, per-trait Δ1 both channels) + `fig6D_crosschem_ladder/MS_row/matched_genes` tables; arm labels per cell (`cells_*.csv`, Fig 2 `cell_type_updated` convention, joined against the delivered score tree) | Panel F dsciCAT-rung source, Panel G labels, Table S6 cross-chemistry block |
| `D3_crossconcordance_domainadj.py` | `D_crosschemistry/tables/D3_crossconcordance_domainadjust.csv` — cross-cohort/cross-channel concordance raw vs adjusted; **the O1 same-nucleus row residualises on ln n_hit (0.823), not ln n_peak** — the fact that forced the Table S6 v11 relabel | Table S6 v11 "A contrast-definition sensitivity" block |
| `D10_deconfound_policy_table.py` | `SUPPORT_v1/tables/deconfound_policy_table.csv` — 8-dataset × channel ρ/P/CI under the ln n_peak criterion (B = 10,000: O1 0.708 [0.550, 0.815], perm P 0.0001, residual confound |ρ| ≤ 0.055) | Panel F n_peak rung CI cross-check; Table S6 v11; Table Sxx criterion ledger |
| `D10_subtype_6pt.py` | `D10_subtype_6pt` six-state six-point cross-channel statistics | Panel D support values |
| `../assembly/S6v8_window_permutation.py` (lives in `assembly/`, not here - it is a panel producer and stays in one place) | `SUPPORT_v2/tables_v8/figS6v8_window_permutation.csv` + `figS6v8_window_permutation_null.csv` - lead-SNP window permutation test, `N_PERM = 10,000` draws (override with `$NPERM`), seed fixed in the script, plus the null distribution it is scored against | Panel B; notebook §0.1 row "window permutation" |

Note (2026-11-18, provenance audit S2): notebook §0.1 used to credit this table to
`upstream/S6v8_window_permutation.py`. That file does not exist; the script is
`assembly/S6v8_window_permutation.py` (md5 `77dfee4983d83b7335777e011c570527`). The label was corrected
and this row was added so the contract table covers every input the notebook reads.

The log1p-corrected score objects themselves (`scdrs_log1p_20261104/**/PASS_*.full_score.gz`,
`tables_log1p_20261104/cells_*.csv`, `fig6D_delta1_74traits_per_tech.csv`) were re-scored on
the data-owner side and delivered as analysis objects; this environment audits rather than
re-runs that step (see `A6_fig6_input_transform_audit.py`). GEO sources: GSE250589 (mouse);
the own dsciCAT-seq cohort is CNP0010263 (CNGBdb).
