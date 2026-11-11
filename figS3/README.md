# Figure S3 — Same cells in single-omic space, cross-library cluster matching, and how deep you have to sequence to see COP/NFOL

`Figure_S3.ipynb` reproduces the measured content of Figure S3 (panels A–G) from
the files the analysis started from. Journal assembly — the swim-lane artwork of
F, the hatch/tick marks, page geometry, panel letters, compositing — is a
production step (`figures/scripts/compose_figS3*.py`,
`figures/scripts/s3_cop_stability.py`) and is deliberately **not** part of the
notebook. What is reproduced is the numbers behind those panels, and both
pictures are redrawn from them in a plain form.

| Panel | What it shows | Section | Output |
|---|---|---|---|
| A | RNA UMAP of each of the three libraries, coloured by its own cell-type call | §1, §2 | `figS3ab_umaps.pdf` (+ `.png`, top row) |
| B | ATAC UMAP of the same cells | §1, §2 | same file, bottom row |
| C–E | cluster–cluster correlation of average expression: 50k vs HT, 100k vs HT, 50k vs 100k | §3 | `figS3cde_xcorr.pdf`, `figS3_xcorr_*.csv`, `figS3_xcorr_stats.csv` |
| F | can a COP/NFOL cluster be found at all, as a function of how many lineage cells enter the clustering (9 conditions, 20 seeds each) | §4, §5.1 | `figS3f_resolvability.pdf`, `figS3_f_conditions.csv` |
| G | the seven-marker identity of the COP/NFOL / OPC / MOL cluster found in each library | §4, §5.2 | `figS3g_marker_columns.pdf`, `figS3g_columns.csv` |
| H | re-clustering the 100k COP/NFOL pool (and, as the depth control, the 50k pool) at five Leiden resolutions: does the claimed second substate appear? | §6 | `s3c_cop2/` → `FigS3c_v10_HIJK.pdf` |
| I | independent-marker check on those two subclusters: 21 frozen markers × 4 groups, BH-FDR | §6 | same |
| J | split stability: 5 whole-library runs + 100 80 % downsamples per resolution, with ARI | §6 | same |

## The chain

```
              3k_dsci_ht_10x_gs_subset.RDS     Figure 2's joint object, 70,643 nuclei
                        │
                        │  one subset() per library - written out and checked in
                        │  Figure 3 §1.0 (cells, order, RNA counts, QC metadata):
                        │    EL14_subset.RDS       batch      == 'Cortex_EL14'     4,029
                        │    3k_FY13_subset.RDS    orig.ident == 'dsciCAT-seq-50k'  17,093
                        │    3k_FY14_subset.RDS    orig.ident == 'dsciCAT-seq-100k' 30,350
                        │    HT_subset_harmony.RDS orig.ident == 'HT-scCAT-seq-HT'  14,762
                        │  each file adds the analysis run at that scope: its own
                        ▼  pca/umap/umap.atac/wnn.umap, graphs, RNA_snn_res.1, cell_type
EL14_subset.RDS   3k_FY13_subset.RDS   3k_FY14_subset.RDS        (Figure 3's fig3_light.RDS)
   ct labels,        ct labels,           ct labels,                umap.atac of all 51,472
   RNA umap,         RNA umap,            RNA umap,                 nuclei, barcodes 100%
   RNA data          RNA data             RNA data                        │
        └───────────────────────┬──────────────────────────────────────────┘
                        §1  S3_light.RDS + S3_avgexpr_{EL14,FY13,FY14}.csv
                        (cell types × the 5,597 HVGs the three libraries share)
                                │
              ┌─────────────────┴──────────────────┐
              §2  panels A/B (RNA + ATAC UMAPs)     §3  panels C/D/E (cluster-cluster r)

HT_subset_harmony.RDS + the 3 subsets ── §4.1 lineage cache (Oligo + OPC, RNA assay)
        └─ §4.2  the delivered clustering at 9 depths × 20 seeds
                 -> S3_u13_clusters.csv · S3_marker_table.csv · S3_strict4029.csv
        └─ §5.1  panel F        └─ §5.2  panel G
```

**All four input objects are selections of one joint object**, not four
datasets: Figure 3 §1.0 writes out the four `subset()` lines and checks them
against the files (same barcodes in the same order, same RNA count matrix, same
QC metadata cell for cell; `HT_subset_harmony.RDS` even agrees with the joint
object on `cell_type`, 14,762/14,762). What each file adds is the analysis run
at its own scope — the single-library or pooled-HT PCA/UMAP/WNN, graphs,
`RNA_snn_res.1` and, for the three single-library ones, the hand-written
`cell_type` table — which is what panels A–E plot and what §4.1 re-clusters, so
the files are read rather than rebuilt.

§1's average-expression matrices are **byte-equal to the cache the delivered
panels were drawn from**: `max |diff| = 0` over 22×5,597, 32×5,597 and 33×5,597
values (only the gene column order differs, which a correlation does not care
about). §3 checks this at run time when `S3_REF` points at that cache.

## Inputs

| Variable | Default | Contents |
|---|---|---|
| `FIGS3_DATA` | `<work>/FY` | the four per-library objects above (0.6 / 1.7 / 2.6 / 1.8 GB) — all four are `subset()`s of `3k_dsci_ht_10x_gs_subset.RDS` (13 GB), which is only read for `umap.atac` and only when Figure 3's cache is absent |
| `FIGS3_LIGHT` | `../fig3/fig3_out/fig3_light.RDS` | Figure 3's derived cache — it already carries `umap.atac` for all 51,472 nuclei, so §1 skips the 13 GB read |
| `FIGS3_OUT` | `figS3_out` | every panel PDF/PNG/CSV, the R sources as run, the caches |
| `S3_REF` | `<work>/figures/scripts/cache` | the type × HVG matrices the delivered C–E were drawn from (validation only) |
| `NSEED` | `20` | down-sampling draws in §4.2 (set it lower for a quick look; the published panel uses 20) |
| `FIGS3_FORCE` | unset | `1` re-runs §1 / §4.1 / §4.2 even when their outputs exist |

## Provenance: read this before comparing with the printed figure

* **The four input objects are selections of Figure 2's joint object** — see the
  chain above and Figure 3 §1.0, where the four `subset()` lines are run and
  verified. Nothing in this notebook is a separate dataset, and nothing here
  re-derives the library-level clustering.

* **A/B are the published embeddings**, not a re-clustering: RNA space is each
  library's own `umap`, ATAC space is the joint object's `umap.atac`, both
  coloured by that library's cell-type call (`cell_type`, except COP/NFOL taken
  by name from `RNA_snn_res.1` cluster 29 / 31 — the Figure 2/3 annotation).
* **C–E are computed from §1's matrices.** Type vectors are matched by Pearson
  r over the 5,597 shared HVGs; colour is the z-score along the query row; the
  support number is a 300-sample gene bootstrap. Measured here, and identical to
  the delivered annotation: COP/NFOL → **OPC** r = 0.81 (50k vs HT) and 0.80
  (100k vs HT), 300/300 bootstraps; OPC → OPC 0.95 / 0.96; 50k vs 100k gives
  COP/NFOL → **COP/NFOL** r = 0.99.
* **F and G are the delivered analysis, re-run** with its clustering parameters
  (NormalizeData → 3,000 vst HVGs → PCA 30 → 10 neighbours → SLM `algorithm = 3`
  at res 0.5/1/1.5/2, primary `RNA_snn_res.1`) and its blind COP-like rule
  (`Bmp4 ≥ 0.25 & Gpr17 ≥ 0.5 & %Gpr17 ≥ 30 & C1ql1 ≤ 0.2 & Ptprz1 ≤ 1.5 & n ≥ 10`).
  Every published number in the panel is reproduced: 144 lineage cells → **0/20**
  draws, 571 → **6/20** (14–20 cells), the strict whole-brain→4,029 thinning →
  **7/20** (12–21 cells), the 694-cell proxy → **9/20** (10–19 cells), pooled-HT
  →144 → **5/20** (16–27 cells), and the full conditions 50k → **65** cells
  (Bmp4 0.524, Gpr17 0.812), 100k → **141** (0.361, 0.877), HT pooled → **59**
  (0.623, 1.168), EL14 single library → **not found**. Panel G's twelve columns
  carry n = 65/292/205, 141/397/430, 59/28/63 and 15/36/21, with the single HT
  library's best candidate **failing the rule on C1ql1** — as printed.
* **Two bugs that changed these numbers are kept in the code, not in a
  footnote.** Cluster means must be divided column-wise
  (`sweep(ex %*% W, 2, n, "/")`); `(ex %*% W) / n` recycles element-wise and made
  every marker value wrong. And a down-sample must be materialised before
  `subset()`: `subset(o, cells = sample(...))` lets Seurat reuse the promise, so
  20 "independent" seeds drew the same cells and a "0/20" was really "0/1".
* **`ScaleData` is deliberately not called in §4.2.** §4.1 caches the RNA
  *assay*, so the lineage cells keep the whole-brain `scale.data`; `RunPCA` then
  uses the genes the subset's variable list shares with it — 646 / 783 / 391 /
  425 of 3,000 for FY13 / FY14 / HT-pooled / EL14. That is the space in which the
  published clusters were found. Re-scaling on the subset (the obvious shortcut)
  moves boundaries and flips the EL14 and HT-pooled outcomes; §4.2 prints the
  count per condition.
* **The strict and proxy versions of the 4,029-nuclei depth are two different
  experiments** and both are in `figS3_f_conditions.csv`: the proxy thins the
  lineage subset to the 694 cells that 4,029 nuclei would contribute, the strict
  one thins the whole brain to 4,029 nuclei and then takes its lineage (median
  706 cells here). The delivered F panel plots the strict one.

## Run

```bash
export FIGS3_DATA=/path/to/FY            # the four Seurat objects + the joint one
export FIGS3_LIGHT=../fig3/fig3_out/fig3_light.RDS   # optional, saves the 13 GB read
export FIGS3_OUT=figS3_out
jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=-1 Figure_S3.ipynb
```

§1, §4.1 and §4.2 are guarded by their outputs. Measured on this machine
(376 GB RAM, 96 cores), 2026-10-22: §1 ≈ 3 min (three objects, ATAC from
Figure 3's cache), §2 0.2 min, §3 ≈ 1 min, §4.1 ≈ 4 min, §4.2 **834 s** (84
re-clusterings + 20 strict draws), §5 < 1 min — **≈ 23 min cold, ≈ 3 min with
the caches in place**.

## Requirements

**Python ≥ 3.8**: `numpy`, `pandas`, `scipy`, `matplotlib`. Nothing is displayed
inline: figures are saved and reported by `show()` (the notebook stays small on
GitHub).

**R ≥ 4.1** with `Seurat` (4.x), `Signac` (loaded with it), `ggplot2`, `dplyr`,
`tidyr`, `ggrepel`, `pals`, `grid`.

Executed with R 4.2.2, Seurat 4.3.0, ggplot2 3.5.2, ggrepel 0.9.5, Python 3.8.14,
numpy 1.22.4, pandas 1.5.1, scipy 1.10.1, matplotlib 3.6.2.

**Fonts**: the delivered page is set in Arial; this machine has no Arial, so the
devices substitute the metric-compatible Liberation Sans / Helvetica.

## §6 — panels H–J: testing "the 100k COP/NFOL pool splits into two clusters" (added 2026-10-29)

Panels A–G are produced by `Figure_S3.ipynb`. **H–J are not notebook cells**: they come from a
separate, pre-registered chain in `s3c_cop2/`, because the question was raised after the notebook
was frozen and the answer had to be fixed *before* looking at the data.

### The contract comes first

`s3c_cop2/S3C_COP2_FROZEN_CRITERIA_20261029.md` was written and committed **before** any of the
numbers below existed. It fixes three independent acceptance tests (S1 independent markers,
S2 resolution stability, S3 no technical confounder) and the rule for what counts as a pass.
`s3c_cop2/S3C_COP2_RESULT_20261029.md` records the outcome against that contract, unedited.

**Outcome: S1 passes, S2 and S3 fail — so the claim is reported descriptively, not as two substates.**

| Test | Result |
|---|---|
| S1 independent markers | **pass** — 260 genes at BH-FDR < 0.05 with \|Δ\| ≥ 0.5 (146 under Bonferroni; Seurat's `p_val_adj` is Bonferroni, not BH, so both are reported). Within the frozen 21-marker panel only *Gpc5*, *Sox6*, *Mobp* qualify |
| S2 stability | **fail** — across the whole frozen grid only 41.9 % of runs recover two clusters (ARI 0.624); restricted to resolution 0.3–0.5 it is 210/210 with ARI 0.936, so the split is real *at that resolution* and resolution-dependent overall |
| S3 technical confounders | **fail** — the two clusters differ in depth: `nCount_RNA` 2,682 vs 4,292 (BH 6.0e-4, δ = −0.40) and `nFeature_RNA` 1,740 vs 2,322 (BH 2.1e-3, δ = −0.35) |

A depth-matched re-check (`s3c_depth_matched_recheck.py`) keeps 251 qualifying genes with 100 %
direction consistency (median effect 0.95×), i.e. the marker difference is not *only* depth — but the
frozen contract still fails on S2/S3, and the contract is what was committed to.
`s3c_c1ql1_sensitivity.py` shows the `C1ql1 ≤ 0.2` criterion can only cancel candidates, never create
them (it changes 1 of 24 runs and 0 of 20 draws).

### Chain

```
s3c_cop_umap_pair.R       re-cluster the 50k / 100k COP-NFOL pools, export UMAP coords + cluster ids
s3c_cop2_analysis.R       the three frozen tests (markers, stability grid, confounders)
s3c_depth_matched_recheck.py   depth-matched control  (fixes a BH-FDR bug: sorted p × unsorted rank)
s3c_c1ql1_sensitivity.py  leave-one-criterion-out sensitivity
s3c_cop_figure.py         panels H–J (matplotlib, Arial, 18.00 × 7.75 cm)
merge_figS3_v8.py         paste that row under the delivered A–G page  (offset snapped to 0.72 pt)
finalize_figS3_AJ.py      install into the delivery block, crop S3b = E–J, re-export PNGs
```

`assembly/patch_figS3_v7.py` + `verify_figS3_v7.py` are the separate A–G wording fix of the same
round: four over-asserted labels ("tracks depth, not chemistry", "never separable", "always
separable", "fails C1ql1") were replaced in the delivered PDF without touching a single other pixel
(`verify_figS3_v7.py` proves the change is confined to the four label footprints).

### Reproduction notes worth knowing before you re-run this

- **Merge offsets must land on the rasterisation lattice.** Pasting a page at an arbitrary y offset
  changes the glyph phase, and a 200/300 dpi pixel diff then reports differences inside panels that
  were not touched. `merge_figS3_v8.py` forces the offset to a multiple of 0.72 pt (= 1/100 inch),
  which is a common lattice point for both.
- **Compare against a freshly opened copy of the baseline.** Rendering the *in-memory* page object
  you just edited can differ from the file on disk (931 px of phantom difference in the Figure 4
  scatter), and comparing a product against itself always passes.
- `s3c_cop_figure.py` must call `figkit.use_matplotlib_arial()`; without it matplotlib silently falls
  back to DejaVu Sans, whose taller boxes make unrelated labels "overlap" in the ink check.
