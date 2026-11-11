# upstream/ — how the inputs of `Figure_7.ipynb` are produced

The notebook starts from delivered analysis tables (per-nucleus scDRS score tables and the
74-trait sweep tables). This directory contains the chain that turns the downloaded cohort
matrices into those tables, in run order. Nothing here is needed to *read* the notebook — it is
needed to re-derive its inputs on a different machine.

Runtime: the scDRS sweeps dominate. They run under the scoring kernel
`$DSCI_PY` (python 3.9.23, **scdrs 1.0.2**); the drawing
notebook (`Figure_7.ipynb`) runs under the figure kernel (python 3.8.14) and does not need scdrs.
Each cohort script scores the 74 MAGMA gene sets with
`scdrs.method.score_cell(n_ctrl=200)` (~40 s/trait on 1 core → ~50 min/cohort; A4/A7 are larger
and take a few hours). Everything else is minutes.

## Step 0 — cohort sources (GEO supplementary files)

| cohort | files | content |
|---|---|---|
| GSE277430 | `GSE277430/` (snRNA-seq of MS + Ctrl brain, 154,522 nuclei after QC) | expert annotations (`expertAnno.l1`), donor, disease |
| GSE118257 | `GSE118257_matrix.txt.gz` + `GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz` | raw counts (9 donors, 4 MS / 5 Ctrl) |
| GSE301908 | `GSE301908/` (17 samples, 233,941 OL-lineage nuclei) | `majorCluster`, `patient_info` |
| GSE231585 | mouse-free human MS cortex RNA + ATAC (separate GSEs; RNA arm here) | author-panel cell types re-annotated by `A7c_author_panel_annot.py` |
| GSE250589 | mouse EAE spinal cord paired multiome, 19 libraries (13 EAE / 6 Ctrl) | author lesion-phase annotation (Ctrl/Early/Peak/Late) |

All matrices are subset to the oligodendroglial compartment (mature OL + OPC) before scoring;
the subset dumps (`OL_csr.bin`) carry **log1p(CP10K)** for A2/A4 (Seurat `data` slot) and raw
counts for A3/A7 — the scale each cohort was actually run with is recorded in notebook §0.1 and
is frozen.

## Step 1 — per-cohort objects

```
A2_GSE277430_v1.py     # dumps OL subset (genes/barcodes/meta/csr) + per-nucleus table
A4_GSE301908_v1.py     # same for GSE301908; A4d_extract_counts.py builds the CSR dump
A7c_author_panel_annot.py  # re-annotates GSE231585 with the author panel -> derived/scored.h5ad
M1b_authoranno.py      # GSE250589: per-library delta1 of every trait (author annotation),
                       #   tables/per_lib_delta1_author_vs_self.csv
```

## Step 2 — 74-trait scDRS sweeps (one per cohort)

```
A2_GSE277430_trait_sweep74.py    # score_cell(n_ctrl=200, random_seed=0)      -> trait_sweep_74.csv
A3_GSE118257_trait_sweep74.py    # random_seed=20261029                       -> trait_sweep_74_GSE118257.csv
A4_GSE301908_trait_sweep74.py    # random_seed=0                              -> trait_sweep_74_GSE301908.csv
A7_GSE231585_delta2_full74.py    # random_seed=20261030                       -> trait_sweep_delta2_GSE231585_74.csv
```

### Step 2b - the two mirror axes (these are the tables Panel C actually reads)

Panel C's dA/dB columns do **not** come from the Step-2 sweeps. They come from the delta2 scripts
below, which re-score each cohort on the astro-immune and neuron-non-neuronal axes. (An earlier
version of this README attributed the dual-axis table to `A2_fig7_v7_dualaxis.py`; that script
*consumes* it - the producer is `A2_GSE277430_delta2_mirror.py`. The attribution was wrong until
2026-11-18, when the three producers were pulled into this bundle and their outputs vendored.)

```
A2_GSE277430_delta2_mirror.py    # n_ctrl=100, random_seed=0  -> trait_sweep_dual_axis.csv (A2)
A3_GSE118257_delta2_full74.py    # n_ctrl=100, MIN_CELL=20    -> trait_sweep_delta2_GSE118257_74.csv
A4_GSE301908_delta2_full74.py    # n_ctrl=200, CAP=12k/lineage -> trait_sweep_delta2_GSE301908_74.csv
A2_fig7_v7_dualaxis.py           # CONSUMES trait_sweep_dual_axis.csv to draw the dual-axis figure
```

The `n_ctrl` values differ between the d1 sweeps (200) and the A2/A3 mirror runs (100); each script
states the value that produced the delivered table, because that is the parameter a re-run must
reproduce, not a style choice.

Each sweep scores the same gene sets, `magma_10kb_1000.74_traits.gs` (74 MAGMA 10-kb gene sets,
1000 genes each), and writes per-trait: cohort-level Δ on the OL-maturation axis (d1) and, where
applicable, on the astro−immune (dA) and neuron−glia (dB) axes, with the descending ranks
(r1/rA/rB). GSE231585's dB axis rests on the cohort's 717 annotated neurons (1.1 % of 62,911 nuclei,
4-133 per sample); under the frozen ">= 10 nuclei per group" rule two of its 15 samples drop out of
that axis, which is why the table carries nB = 13 while n1 = nA = 14, and why Panel C marks that one
violin with `*`. Nothing is thinned or subsampled here - 717 is what the annotation yields. (An earlier
wording of this paragraph said "after thinning to 717"; that was inherited from a code comment and is
wrong: `A7_GSE231585_delta2_full74.py` applies no cap to any group, unlike A4's <= 12,000 per lineage.)

## Step 3 — derived tables the panels read

```
A6_stepC2_sizematched.py   # n_peak_raw: deduplicated ATAC peaks (public human snATAC GSE231584,
                           #   GRCh38, 121,505 peaks, TSS +/- 100 kb) covering each human gene set ->
                           #   trait_chromatin_sizeMatched.csv
A8c_confound_unified.py    # the B-panel ledger: delta1 matrix over the four cohorts, raw ranks,
                           #   stratified-by-domain-size ranks on ln(n_peak_raw), matched-window
                           #   (+/-15 %) control, confound rho -> seven4_confound_ALL_cohorts.csv
                           #   (B's face colour is rank_adj_* = the size-adjusted ranks; the raw ranks
                           #   and the matched-window control are the Fig S7 sensitivity views)
D9_fig7_assembly.py        # the A-panel ledger: unit-level delta1 for 4 human cohorts + mouse
                           #   GSE250589 (fig7A_donors_delta1.csv) and between-unit MWU P
                           #   (fig7A_between_donor_P.csv)
```

## Step 4 -- Panel D's inputs

Panel D is drawn from the drug-target audit chain that was built for Fig 7D. As of 2026-11-18 the
five producers are **inside this bundle** (`s13_two_contrasts.py`, `s15b_allgene_pop_counts.py`,
`s23_expand_new_targets.py`, `s27_target_ledger.py`, `s32_flags_v2.py`); their outputs are vendored in
small form under `../tables_v17/` so Panel D needs no 1.1 GB matrix read.

| script | produces | consumed by | note |
|---|---|---|---|
| `s13_two_contrasts.py` | `S13_localization_classes.csv`, `S13_two_contrasts.csv` | `s23`, `s27`, Panel D | two one-sided Fisher contrasts per gene x cohort (lineage vs non-lineage; mature OL vs OPC), BH within cohort and direction |
| `s15b_allgene_pop_counts.py` | `S15_det_*.csv`, `S15_n_*.csv` | `s23`, `s27`, Panel D | detected / total nuclei per gene x population over the **whole released matrix**; the only step that reads the matrices (~1.1 GB), everything downstream is table-only |
| `s23_expand_new_targets.py` | `S23_new_targets_stats.csv`, `S23_presence_in_release.csv`, `S23_gene_universe_status.csv` | `s27` | the 33-row expansion scored as atoms (its universe includes HTR2B from row 21, which the ledger's split does not) |
| `s27_target_ledger.py` | `S26_target_ledger.csv`, `S26_rowset_C12.csv` | `s32`, Panel D | the ledger the figure row set is a FUNCTION of: 83 atoms, tier, class, floor, drawn/not-drawn + reason |
| `s32_flags_v2.py` | `S29_flags_v2.csv` | Panel D | the 18 drawn rows x 3 evidence flags, author-side re-audit (openFDA label text + Europe PMC) |

`s23` and `s27` read `../data/EXPANSION_33/expansion33_rows.csv` instead of the expansion memo. The
memo (`MS_DRUG_TARGET_EXPANSION_33.md`) is an internal Chinese-language document and is not pushed;
the CSV is its ASCII extract - 33 rows, claim cell verbatim, and `fda_status_en` holding **exactly the
English strings the delivered ledger carries**. Before this existed the memo's Chinese status was
translated to English ad hoc inside a session, so the delivered `S26_target_ledger.csv` could not be
re-derived by any committed script; the CSV closes that hole and the tier rule in `s27` was rewritten
against the English vocabulary (same branch order). Re-running `s13 -> s23 -> s27 -> s32` from the
analysis-tree inputs now reproduces all five delivered tables cell-for-cell.

The donor-level re-test of Panel D's detection rule is **not** part of Figure 7: it is a robustness
check and therefore belongs to the supplement. The scripts and the tables were moved to
`figS7/handoff/donor_s1pr5_20261117/` (see the README there for the numbers and for what the S7 panel
still needs).

`Figure_7.ipynb` recomputes the panel quantities from the Step-2/3/4 outputs and draws the four panel
PDFs; it does not re-run scDRS.
