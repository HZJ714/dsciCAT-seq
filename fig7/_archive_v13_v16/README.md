# Figure 7 — External replication of the mature-OL vs OPC scDRS contrast

`Figure_7.ipynb` reproduces the three delivered panels of **Figure 7** from the analysis objects
upstream of it: per-nucleus scDRS MS scores of four independent human MS cohorts
(GSE277430 / GSE118257 / GSE301908 / GSE231585) and of the mouse EAE spinal-cord cohort
(GSE250589), plus the 74-trait sweep tables built from the same scores. Journal assembly (the
single-column 8.0 cm sheet, panel letters, final cropping) is done outside the notebook;
`assembly/` holds those scripts and the QA gates the delivered page passes.

## Panels

| Panel | File | Size | Content |
|---|---|---|---|
| A | `fig7v13_A_donors.pdf` | 7.92 × 4.44 cm | unit-level Δ1 = median(mature OL) − median(OPC) of the MS scDRS score, per donor (A2/A3/A7) or library (A4); mouse libraries split by lesion stage Early/Peak/Late; MWU *P* per column |
| B | `fig7v16_B_rawband.pdf` | 7.84 × 3.23 cm | 74 GWAS traits × 4 human cohorts; colour = raw rank of cohort-level Δ1 (1 = strongest OL contrast), left → right = mean raw rank over the four cohorts; MS / Alzheimer / Schizophrenia marked; stratified-by-regulatory-domain-size and ±15 % matched-window sensitivity views are in Fig S7A–C |
| C | `fig7v13_C_threeaxes.pdf` | 7.90 × 5.71 cm | same cohorts 2 × 2, three cell-identity axes (OL maturation / astro−immune / neuron−glia) as robust *z* within cohort × axis; MS/SCZ marked with axis ranks; `*` = GSE231585 neuron−glia axis rests on 717 paired nuclei |

Delivered composed sheet: `fig7_composed_v16_rawband.pdf`, **7.93 × 13.98 cm** (single-column
hard limit 8.0 cm); panel letters 8 pt Arial regular; in-figure text 8 pt, 6 pt only where 8 pt
demonstrably collides.

## Pipeline

```
GEO cohorts (GSE277430 / 118257 / 301908 / 231585 / 250589)
  └── upstream/ Step 1: OL-subset dumps (log1p CP10K for A2/A4; raw for A3/A7), annotations
        └── Step 2: scDRS score_cell(n_ctrl=200, fixed seed) over magma_10kb_1000.74_traits.gs
              ├── per-nucleus MS score tables  ──────────────────┐
              ├── 74-trait sweep tables (d1 / dA / dB + ranks) ──┼──> A (via Step 3 D9 ledger)
              └── public human snATAC GSE231584 (GRCh38, 121,505 peaks): n_peak_raw ────────────> Fig S7A–C
                  stratified / matched-window sensitivity views (via Step 3 A8c ledger)
```

## Inputs

All inputs are read through environment variables; nothing is downloaded by the notebook.

| Variable | Object | Produced by |
|---|---|---|
| `FIG7_ROOT` | `routeA_out/A{2,3,4,7}_*_v1/tables/*` | `upstream/` Steps 1–3 (cohort objects, sweeps, ledgers) |
| `FIG7_ROOT` | `A6_GSE231584_atac_v1/tables/trait_chromatin_sizeMatched.csv` | `upstream/A6_stepC2_sizematched.py` (public human snATAC GSE231584, GRCh38, 121,505 deduplicated peaks, TSS ± 100 kb) |
| `FIG7_ROOT` | `M1b_authoranno_v1/tables/per_lib_delta1_author_vs_self.csv` | `upstream/M1b_authoranno.py` (GSE250589, author annotation) |
| `FIG7_PROBE` | `routeA_probe/GSE*` | GEO supplementary downloads (see `upstream/README.md` Step 0) |
| `FIG7_OUT` | panel PDFs | this notebook |
| `ARIAL_TTF` | `ARIAL.TTF` | Arial face (registered when the OS font is absent; cosmetic only) |

The scDRS score objects are large (A4 alone is 233,941 nuclei × 5.4 GB dump) and are not in the
repository; the notebook reads only the delivered tables, so an end-to-end run is ~2 min.

## Running

```bash
cd fig7
FIG7_ROOT=/path/to/routeA_out FIG7_PROBE=/path/to/routeA_probe \
  jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.timeout=-1 Figure_7.ipynb
```

Executed end-to-end on 2026-11-06 with **0 errors** (Python 3.8.14, numpy 1.22.4, pandas 1.5.1,
scipy 1.10.1, statsmodels 0.13.5, anndata 0.9.2, matplotlib 3.6.2; scdrs itself is not needed by
the notebook — the upstream sweeps ran in `$DSCI_PY` (python 3.9.23,
**scdrs 1.0.2**), and §0.2 prints both the local versions and this registration).

## Validation checklist

Every number below is printed by the notebook and matches the delivered tables bit-for-bit.

| Quantity | Notebook | Delivered |
|---|---|---|
| A: units per cohort (case/ctrl) | 13/8, 4/5, 14/3, 5/3, mouse 13/6 | `fig7A_donors_delta1.csv` |
| A: between-unit *P* | 0.500, 1.000, 0.091, 0.393; mouse overall 0.009 | `fig7A_between_donor_P.csv` |
| B: MS raw ranks (A2/A3/A4/A7) | 2 / 4 / 7 / 2 | `seven4_confound_ALL_cohorts.csv` |
| B: leftmost trait (raw #1) ranks | PLATELET 1 / 2 / 1 / 1 | idem |
| B: confound ρ (Δ1, ln n_peak) | −0.838 / −0.691 / −0.778 / −0.873 | idem (in Fig S7 legend, additional statistics — not on the face) |
| C: OL-axis ranks of MS/SCZ | identical to B raw ranks | four sweep tables (`r1`) |

## What is *not* here

- The deconfounding diagnostics (confound-ρ strip, stratified-by-regulatory-domain-size ranks,
  ±15 % matched-window controls) are supplementary (Fig S7A–C) — the main figure deliberately shows only
  the raw rank band.
- The paired multiome (ATAC-channel) probe on GSE250589 stayed negative at the targeted-trait
  level (12 traits, cross-channel Spearman ρ = −0.12 raw / 0.06 stratified by regulatory-domain
  size, n = 12) and is reported in the supplement, not in Figure 7.
