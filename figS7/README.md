# Figure S7 — rank basis, cross-channel controls and the donor-level detection test (v30 FINAL)

`Figure_S7.ipynb` reproduces the ten delivered PDFs of **Figure S7** from the tables upstream of it:
the 74-trait scDRS sweeps of four human post-mortem MS cohorts, a human snATAC cohort profiled from the
same donors (GSE231584), mouse EAE spinal-cord multiome measured in the same nuclei (GSE250589), a
target x population detection table (GSE118257), and a per-donor re-test of Fig 7D's S1PR5 detection
rule. Journal assembly (the two-column sheet, panel letters, final crop) and the QA gates run **outside**
the notebook; `assembly/` holds the panel code.

Everything the notebook reads is vendored in `data/` and `tables_v29/`, listed row by row in
`data/MANIFEST.csv` with the script that produced it. No scoring is re-run here.

## Panels

| Panel | File | Size (cm) | Content |
|---|---|---|---|
| K | `figS7_out/panels_v30/S7v30_K_keyband.pdf` | 14.13 x 0.31 | legend artwork: the seven channels in one row, each named by accession (RNA block: GSE277430 / 118257 / 301908 / 231585 / 250589; ATAC block: GSE231584 / 250589) |
| A | `..._A_ladder.pdf` | 4.38 x 4.03 | 74 gene sets x 4 human RNA cohorts ranked twice - as computed, and after the regulatory-domain-size trend is removed. Rows ordered by domain size (`n_peak_raw` desc; MS row 37, PLATELET row 56). MS ranks 2/4/7/2 and 1/1/3/1 |
| B | `..._B_tightpool.pdf` | 3.27 x 3.24 | the tight pool: MS + the 15 gene sets nearest to it on domain size, ranked within those 16; MS is rank 1 in every cohort |
| C | `..._C_human_samedonor.pdf` | 3.32 x 3.57 | RNA d1 (GSE231585) vs ATAC d1 (GSE231584, size-matched) for the 74 sets: rho 0.722 raw -> 0.137 after size residualisation (P 0.24); MS **2 -> 36** on the plotted axes |
| D | `..._D_mouse_samenucleus.pdf` | 2.91 x 3.68 | mouse EAE multiome, both channels measured in the same 14,928 nuclei (25 libraries): rho 0.248, P 0.033; MS 11 -> 36 |
| E | `..._E_ledger.pdf` | 15.17 x 3.00 | the 74 x 21 ledger: 7 channels x 3 axes, every cell a rank of 74; rows ordered by the mean of the four human cohorts' as-computed d1 rank (MS row 2) |
| RB | `..._RB_rankbar.pdf` | 4.06 x 0.56 | the rank bar shared by A, B and E |
| L2 | `..._L2_axiskey.pdf` | 8.15 x 0.31 | what d1 / d2A / d2B subtract (`d2A = AST - immune`, `d2B = neuron - non-neuron`) |
| F | `..._F_detection.pdf` | 8.35 x 3.57 | 5 Fig 7D targets x 15 populations in GSE118257 (17,799 nuclei): area = % detecting, fill = mean log2(CP10k+1), per-population n under each column |
| G | `..._G_donor.pdf` | 6.16 x 3.36 | the detection rule per donor: 27 evaluable donors, 22 meet it (8/11 GSE301908, 7/8 GSE118257, 7/8 GSE231585); y is a within-donor margin, the red dashed line is its +1.0 pp arm |

Delivered composed sheet: **17.00 x 14.28 cm** (two-column cap 17.00 cm), assembled from the ten
panels below; the composed sheet PDF is intentionally **not in git** (this repository ships single
panels + the code that draws them — composition is left to the reader/typesetter, per author policy
2026-11-25). In the working tree it sits at `figS7_out/figS7_v30.pdf`.
Layout: K above row 1; A | B | C | D on row 1, their matrices aligned to a common measured axes top; the
rank bar and the delta key in one band directly above E, left-aligned to E's **measured** matrix edge;
E full text width; F | G on row 3. Panel letters A-G are 8.0 pt Arial regular drawn by the composer into
a reserved 0.34 cm strip above each row - the panels themselves carry no letters. The face uses exactly
two tiers, 8.0 and 6.0 pt.

Legend and face-number contract: `LEGEND_S7_final.md`. It states which numbers are deliberately not on
the face (the withdrawn three-basis table, the MS ledger profile, the no-floor 6/11 donor tally, the
floor sweep, the diagnosis stratification) and maps the old S7A...S7I references onto the new letters.

## Pipeline

```
GEO cohorts  GSE277430 / 118257 / 301908 / 231585 (human RNA), GSE231584 (human snATAC),
             GSE250589 (mouse EAE spinal cord, paired multiome)
  |
  |-- upstream/*_trait_sweep74.py, human_atac_threeaxis.py, mouse_*_threeaxis.py
  |      scDRS 1.0.2 score_cell(n_ctrl=200) over magma_10kb_1000.74_traits.gs, fixed seeds
  |      -> 74-trait tables of d1 / d2A / d2B per cohort and per channel          [data/, MANIFEST]
  |
  |-- upstream/confound_all_cohorts.py   -> seven4_confound_ALL_cohorts.csv
  |      d1, ln(peak-domain size), raw ranks, ranks after residualising d1 on ln n_peak
  |      (the same file is a Fig 7 input: the two figures cannot disagree on a rank)
  |-- upstream/axis_rho_and_rank_matrix.py -> the 21-column rank matrix behind E
  |-- upstream/matched_window_table.py     -> the +-15% size window asserted by A and B
  |-- upstream/delta1_vs_peakdomain.py     -> the size-matched ATAC d1 that C plots
  |-- upstream/crosschannel_pairs_table.py -> the rho / P C is allowed to print
  |-- upstream/export_mouse_same_nucleus_delta1.py -> D's two series + rho / P / ranks / n
  |-- upstream/target_stats_GSE118257_v7.py -> F's detection and expression table
  |
  |-- handoff/donor_s1pr5_20261117/s41_donor_all_human_cohorts.py
  |      per-donor OL vs strongest-comparator margins under the >= 30-nucleus rule
  |      -> tables_v29/S41_S1PR5_per_donor.csv, S41_summary_by_cohort.csv          [panel G]
  |
  +-- assembly/make_figS7_v30_panels.py    one PDF per panel, asserted, no compose
        (Figure_S7.ipynb is generated from this file by build_figS7_notebook.py)
```

## Inputs and their provenance

`data/MANIFEST.csv` is the authoritative list (path, bytes, rows, md5, producer script, cohort, panels,
content). The short version:

| Used by | Input | Producer | Note |
|---|---|---|---|
| A, B, C, E | `data/FIG7_v1/tables/seven4_confound_ALL_cohorts.csv` | `upstream/confound_all_cohorts.py` | d1 + ln n_peak + raw and size-adjusted ranks for the four human cohorts; C plots its `A7` column |
| A, B | `data/SUPPORT_v1/tables/fig7Sb_matched_window.csv` | `upstream/matched_window_table.py` | the +-15 % window (40 traits) and rho(ln n_peak, d1), both asserted |
| E | `data/SUPPORT_v1/tables/fig7Sc_rank_matrix_4cohorts.csv` | `upstream/axis_rho_and_rank_matrix.py` | the ledger's human RNA cells; its `_1` columns are the as-computed ranks and their mean is E's row order |
| E | `S7_atac_crossaxis_human_full.csv`, `S7_atac_crossaxis_mouse_full.csv`, `S7_rna_crossaxis_mouse_full.csv` | `upstream/human_atac_threeaxis.py`, `mouse_atac_threeaxis.py` (+ composition fix), `mouse_rna_threeaxis.py` | the three non-RNA channels; mouse RNA and mouse ATAC share nuclei, labels, gene sets and seed 20261030 |
| C | `data/A6_GSE231584_atac_v1/tables/delta1_vs_peakdomain_allchannels.csv` | `upstream/delta1_vs_peakdomain.py` | `median_groupdiff`: the size-matched ATAC d1 that C plots |
| C | `data/SUPPORT_v1/tables/figS7c_crosschannel_pairs.csv` | `upstream/crosschannel_pairs_table.py` | the only rho / P the panel may print, for GSE231585 <-> GSE231584 |
| C, Table S7 | `data/A6_GSE231584_atac_v1/tables/same_donor_RNA_vs_CHROMATIN_delta1_74traits.csv` | delivered input (donor-paired ATAC chain) | the donor-paired ranks (2 -> 21). **Not plotted since v29** - see "What changed" below |
| D | `data/SUPPORT_v1/tables/S7_mouse_same_nucleus_delta1.csv` + `.json` | `upstream/export_mouse_same_nucleus_delta1.py` | per-trait d1 for both channels of the same nuclei; the JSON carries rho / P / ranks / 14,928 nuclei / 25 libraries |
| F | `data/A3_GSE118257_v1/tables/A3_GSE118257_target_stats.csv` | `upstream/target_stats_GSE118257_v7.py` | 30 targets x 15 populations, % detecting + mean log2(CP10k+1), 17,799 nuclei |
| G | `tables_v29/S41_S1PR5_per_donor.csv`, `S41_summary_by_cohort.csv` | `handoff/donor_s1pr5_20261117/s41_donor_all_human_cohorts.py` | one row per donor; the summary keeps both the >= 30-nucleus and the no-floor tallies |

Cohort objects, annotations and the scDRS gene sets behind these tables are described in
`upstream/README.md` and `handoff/donor_s1pr5_20261117/README.md`. The 74 MAGMA gene sets are the same
file Fig 6 scores, so a trait means the same thing in both figures.

## Frozen parameters

| | |
|---|---|
| scoring | scDRS **1.0.2**, `score_cell(n_ctrl=200)`, `scdrs.pp.preprocess`, seeds 0 (A2/A4), 20261029 (A3), 20261030 (A7 and both mouse channels) |
| axes | d1 = median(mature OL) - median(OPC); d2A = median(AST) - median(immune); d2B = median(neuron) - median(all non-neuronal) |
| group floor | >= 10 nuclei per group for a sweep unit to count; **>= 30 nuclei** for a population to enter a donor's comparator set in G |
| size adjustment | rank of d1 after residualising on ln(regulatory-domain size) |
| tight pool (B) | MS + 15 nearest neighbours on ln n_peak, ranked within 16 |
| detection rule (F, G) | margin >= 1.0 pp AND (OL >= 1.5x comparator OR comparator = 0) |
| statistical unit | donor where the cohort is donor-resolved, library where it is not; GSE301908's 17 released libraries are 11 donors |

## How to run

```bash
cd figS7
export DSCI_WORK_ROOT=/path/to/analysis/tree      # only for fonts/ARIAL.TTF
python3 build_figS7_notebook.py                    # regenerates Figure_S7.ipynb from assembly/
jupyter nbconvert --to notebook --execute --inplace Figure_S7.ipynb
```

The notebook writes the ten panel PDFs to `figS7_out/panels_v30/` and prints each panel's measured size
(`export()` reads the produced PDF; nominal canvas sizes are never reported). To rebuild the composed
sheet, run the composer in the analysis tree
(`figures/SUPPL_ANALYSES_20261112/FIGS7_V23_RELAYOUT_20261118/render_v30.py`), then the gate battery
(`gates_v30.py`) and the ink-equivalence check (`verify_v30_ink.py`) there. The notebook deliberately
contains no compose, no letters and no QA cells (house spec: those are assembly, not content).

## Validation checklist for this version

- [x] `Figure_S7.ipynb` executed end-to-end, **0 error cells**, English throughout, no CJK.
- [x] every pushable file of this bundle (and of `../fig7`) is free of CJK - scanned, not assumed;
      the producer scripts in `upstream/` were Englishised with an AST-skeleton guard, which also exposed
      and fixed the Chinese `kind`/`pair` labels that `delta1_vs_peakdomain.py` wrote into the D3 chain.
- [x] Notebook panel PDFs carry the same words and the same drawing geometry as the delivered ones.
- [x] `verify_v30_ink.py`: C and D re-drawn from tables reproduce the delivered v22 scatters to
      **0.03 pt / 0.02 pt** over all 74 dots; A, B, E, F, G, K, RB, L2 are ink-identical to the approved
      v28 prototype (raster pixels identical, internal geometry <= 0.02 cm).
- [x] `gates_v30.py`: font `--strict` PASS, collisions 0, text-on-data 2 (the MS callouts on their own
      dots, as in the signed sheet), ink outside page 0, ink on page border 0, run-on 0, word ledger
      167 = 160 panel words + 7 letters with the 7 letters as the only phantom, clipped ink 0 of 10,
      detection-ink audit CONSISTENT. **FAILURES: none.**
- [x] `check_palette()` asserts the colour floors at render time (adjacent 68.6, pair 40.7, semantic
      18.6 CIE76, worst of three vision models).
- [x] every face number is computed and asserted in the script, never typed: ranks 2/4/7/2 and 1/1/3/1,
      window 40, pool ranks 1/1/1/1, MS row 37/56 in A and row 2 in E, rho 0.722/0.137/0.248,
      MS 2->36 and 11->36, 17,799 nuclei in F, 27 evaluable donors and 22/27 rule met in G.
- [ ] `Table S7` written (three-basis ranks, MS ledger profile, per-donor detail with the diagnosis
      column, the floor sweep and the no-floor 6/11) - required before the manuscript cites this sheet.
- [x] face-number register: `assembly/register_figS7_faces.py` -> `figS7_out/figS7v30_face_numbers.csv`
      (38 entries; forward and reverse gates both pass, run it before any push).
- [ ] manuscript pointers S7F/S7G repointed to Table S7; Fig 7D legend gains the 8 patients + 3 controls
      disclosure; abbreviation table reconciled with the face (`NEEDS_ABBREVIATIONS.md`).

## What changed from the delivered v22, and the corrections in v29 / v30

v22 (A-I, 17.00 x 16.81 cm) is archived in `_archive_v22/` (this working tree only; just its README is versioned). The relayout withdrew three panels whose
content is a table or a row of the ledger (the three-basis rank table, the cross-channel rank profile,
the MS rank table), re-lettered the rest, and added the donor panel that used to sit in Fig 7. v28
reclaimed two pieces of dead whitespace; v29 added the two corrections below. **v30 is the pushable
version** and adds one thing: the cell-type vocabulary on the face is now the same set of symbols
Fig 7 uses (`OL/MOL - OPC`, `AST - (MG+LYM)`, `NEU - non-NEU` in the delta key, and AST / MG / LYM / NEU
as panel F's column heads for those four populations). The v29 sheet is archived in `_archive_v29/` (not versioned).

1. **C's MS annotation corrected `2 -> 21` to `2 -> 36`.** The plotted y is the size-matched ATAC series
   - the series the printed rho is computed on - and MS is rank 36 there; 21 was its rank in the
   donor-paired table, which is not plotted. Both ranks are now computed from the plotted columns.
2. **C and D drawn instead of carried over.** The prototype chain reused the delivered v22 PDFs for
   these two panels; a notebook must produce every panel, so they are re-implemented from the vendored
   tables and proven point-by-point against the signed ink.

## Directory

```
figS7/
  Figure_S7.ipynb                  generated, executed; one section per panel
  build_figS7_notebook.py          the generator (edit code in assembly/, prose here)
  assembly/make_figS7_v30_panels.py  the panel code - single source of truth
  figS7_out/panels_v30/            the ten delivered panel PDFs (+ frames JSON)
  (figS7_out/figS7_v30.pdf)        the composed sheet - working tree only, not in git
  data/                            vendored inputs + MANIFEST.csv (md5, producer, panels)
  tables_v29/                      the two donor tables (from handoff/)
  upstream/                        the scripts that produced data/ (see upstream/README.md)
  handoff/                         the donor-level pipeline this figure inherited, + its tables
  figstyle.py                      the house matplotlib style module other bundles import
  verify_text_on_data.py           the occlusion gate
  LEGEND_S7_final.md               journal legend + face-number contract
  _archive_v22/                    superseded v22 sheet, panels, composer, verifier, legend (not pushed)
```
