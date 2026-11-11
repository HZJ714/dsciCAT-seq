# Figure 7 — External replication of the mature-OL vs OPC scDRS contrast (v17 FINAL)

`Figure_7.ipynb` reproduces the four delivered panels of **Figure 7** from the analysis objects
upstream of it: per-nucleus scDRS MS scores of four independent human MS cohorts
(GSE277430 / GSE118257 / GSE301908 / GSE231585), of the mouse EAE spinal-cord cohort (GSE250589),
the 74-trait sweep tables built from the same scores, and the drug-target localisation table behind
Panel D. Journal assembly (the two-column sheet, panel letters, final crop) is done outside the
notebook; `assembly/` holds the panel code and the QA gates the composed page passes.

## Panels

| Panel | File | Size | Content |
|---|---|---|---|
| A | `fig7_out/fig7v17_A_donors.pdf` | 6.58 × 5.90 cm | unit-level Δ1 = median(mature OL) − median(OPC) of the MS scDRS score, per donor (A2/A3/A7) or library (A4); mouse libraries split by lesion stage Early/Peak/Late; Mann-Whitney *P* above each column (bare number; the "P =" prefix is in the caption) |
| B | `fig7_out/fig7v17_B_adjranks.pdf` | 3.78 × 5.85 cm | 74 GWAS gene sets × 4 human cohorts, portrait rank heatmap; colour = rank of cohort-level Δ1 **after residualising Δ1 on ln(regulatory-domain size)**, rows ordered by the mean of those four ranks; MS / AD / SCZ named. Raw ranks, the size-stratified ladder and the ±15 % matched-window control are in Fig S7 |
| C | `fig7_out/fig7v17_C_threeaxes.pdf` | 5.84 × 5.98 cm | same cohorts 2 × 2, three cell-identity axes, labelled with the manuscript's cell-type symbols — `OL−OPC`, `AST−(MG+LYM)`, `NEU−non-NEU` — as robust *z* within cohort × axis; MS/SCZ dots carry their rank on the axis they sit on; the single `*` sits **directly above GSE231585's `NEU−non-NEU` violin**; tiles are open-frame (left + bottom spines only) |
| D | `fig7_out/fig7v17_D_targets.pdf` | 17.00 × 7.90 cm | 18 audited drug targets: within-cohort log2 OR (mature OL vs OPC; OL lineage vs non-lineage) over three cohorts, drug-target flags (label / stated MOA / used in MS), and % of nuclei detecting each gene by population in GSE301908; the red box marks the only row passing the detection rule |

Composed sheet (current, **not in this repository** - composed outside the repo by the typesetting
round, per author policy 2026-11-25 this repo ships only the single panels above):
`<work>/figures/SUPPL_ANALYSES_20261112/FIG7_2COL_20261117/pdf/fig7_composed_v17d_2col.pdf`,
17.00 × 14.42 cm; the four panels that make it are `fig7_out/fig7v17_{A,B,C,D}*.pdf`. **v17d (2026-11-19, author fix round)**: all rotated axis/cohort labels re-anchored.
The v17c "half-projection right-shift" (`rot_anchor`, +tw×0.3536) made each 45°/60° label's ink centre
sit on its column but pushed the reading END 0.35–0.5 cm right of the tick — the author read panel D's
cohort names as shifted onto the next column, and panel B's 90° labels were uniformly 0.115 cm LEFT of
their ticks (no `rotation_mode="anchor"`). All labels now use the matplotlib standard: 45°/60° =
`ha="right", va="top", rotation_mode="anchor"` (upper tip ON the tick, body reads down-left);
B's 90° = `ha="right", va="center", rotation_mode="anchor"` (vertical line centred on the tick).
Measured on the composed page: B centre offset −0.115 → −0.019 cm; A/D reading-end offset → ≤0.07 cm.
Scripts `s36/s37`: shift calls removed, `rot_anchor()` kept as dead code for provenance. Gates: font
strict PASS / collide 0 / text-on-data 0 / run-on 0; composed word multiset identical to v17c (only
positions moved). v17c kept byte-identical as `fig7_composed_v17c_2col.pdf` (+ `pdf/v17c_backup_20261119/`).
v17c is v17 with the two 2026-11-18 author revisions and nothing else: panel C's
three axis labels switched to the abbreviation set (panel D of the same figure already printed AST /
MG / LYM / NEU), and panels A and C lost their top/right spines (`open_spines()`; the house style
loaded by `use_arial()` re-enables all four spines after the local rcParams, which is why the frames
were there at all). The pre-revision sheet is kept byte-identical as
`fig7_composed_v17_2col.pdf`, with `fig7_composed_v17_2col_ASTRO_LABELS_backup.pdf` as its twin.
Original delivered sheet: `figures/SUPPL_ANALYSES_20261112/FIG7_2COL_20261117/pdf/fig7_composed_v17_2col.pdf`,
**17.00 × 14.34 cm** (two-column width 17.00 cm; A | B | C on row 1, D on row 2). Panel letters are
8.0 pt Arial regular, added by the composer; the face uses exactly two tiers, 8.0 and 6.0 pt.

**The `*` in Panel C** is a data-thinness caveat, not a significance marker: GSE231585 carries 62,911
nuclei of which **717 are annotated neurons** (1.1 %; 4–133 per sample), so under the frozen "≥ 10
nuclei per group" rule two of its 15 samples drop out of the neuron−glia axis (nB = 13 while
n1 = nA = 14). The caveat is axis-specific, which is why the mark hangs over that violin rather than
on the cohort title. No other tile is marked: the other cohorts' neuron pools are 7–100× larger
(GSE277430 76,436; GSE118257 5,031; GSE301908 ≤ 12,000 per lineage by the stratified draw).

## Pipeline

```
GEO cohorts (GSE277430 / 118257 / 301908 / 231585 / 250589)
  └── upstream/ Step 1: OL-subset dumps (log1p CP10K for A2/A4; raw for A3/A7), annotations
        └── Step 2: scDRS score_cell(n_ctrl=200, fixed seed) over magma_10kb_1000.74_traits.gs
              ├── per-nucleus MS score tables  ──────────────────┐
              ├── 74-trait sweep tables (d1 / dA / dB + ranks) ──┼──> A, B, C (via Step 3 ledgers D9 / A8c)
              ├── public human snATAC GSE231584 (GRCh38, 121,505 peaks): n_peak_raw
              │        └── rank after residualising Δ1 on ln n_peak ──> B colour
              └── per-gene × per-population detection counts (S15) ──> D
```

## Inputs

All inputs are read through environment variables; nothing is downloaded by the notebook.

| Variable | Object | Produced by |
|---|---|---|
| `FIG7_DATA` | `data/A{2,3,4,7}_*_v1/tables/*.csv`, `data/FIG7_v1/`, `data/FIG7_v11/`, `data/EXPANSION_33/` | **vendored** here; producers `upstream/` Steps 1–3b + the audit chain (see `data/MANIFEST.csv`) |
| `DSCI_WORK_ROOT` | the analysis tree | needed only for the Arial TTF when drawing; the cohort tables behind A/B/C are vendored since 2026-11-18 (provenance audit G1) |
| `DSCI_WORK_ROOT` | `figures/routeA_out/A6_GSE231584_atac_v1/tables/trait_chromatin_sizeMatched.csv` | `upstream/A6_stepC2_sizematched.py` (public human snATAC GSE231584, GRCh38, 121,505 deduplicated peaks, TSS ± 100 kb) |
| `DSCI_WORK_ROOT` | `figures/routeA_out/M1b_authoranno_v1/tables/per_lib_delta1_author_vs_self.csv` | `upstream/M1b_authoranno.py` (GSE250589, author annotation) |
| `DSCI_WORK_ROOT` | `routeA_probe/GSE301908_counts/` | GEO supplementary download; read by this bundle only through the vendored `S15_*` tables (the per-donor recount lives with Fig S7) |
| `FIG7_TABLES` | `tables_v17/*.csv` | vendored: Panel D's statistics / flags / detection tables (producers: the Fig 7D audit chain, now **inside** `upstream/`, see Step 4) |
| `FIG7_OUT` | panel PDFs | this notebook writes here |

The scDRS score objects are large (A4 alone is 293,970 nuclei, a 4.4 GB dump) and are not in the
repository; the notebook reads only the delivered tables, so an end-to-end run is ~2 min. Panel D reads
the small vendored CSVs in `data/` and `tables_v17/`. Vendoring was closed on 2026-11-18: re-running
the notebook from the clone alone reproduced the four delivered panels (page size, every word, coordinates
to 0.02 cm, drawing multisets, fonts), and re-running `s13 -> s23 -> s27 -> s32` reproduces the five
delivered Panel D tables cell-for-cell.

**Inputs that come from another figure** are labelled in the notebook's §0.1 and in `upstream/README.md`:
`seven4_confound_ALL_cohorts.csv` and `trait_chromatin_sizeMatched.csv` are the *same* files Fig S7 reads
(produced by `upstream/A8c_confound_unified.py` and `A6_stepC2_sizematched.py` here, vendored into
`figS7/data/`); the 74 MAGMA gene sets are the same `magma_10kb_1000.74_traits.gs` that Fig 6 scores; the
mouse Δ1 comes from `upstream/M1b_authoranno.py`, which reuses the per-cell scores of the earlier M1
mouse run instead of re-scoring; and Panel D's tables come from the Fig 7D drug-target audit chain
(`s13 / s15b / s23 / s27 / s32`), whose outputs are vendored. The four Spearman rho values quoted for
Panel B are **recomputed inside the notebook** from the delivered Δ1 and ln n_peak columns, so this bundle
no longer reads Fig S6/S7's `D3_crossconcordance_domainadjust.csv` at all.

## Running

```bash
cd fig7
DSCI_WORK_ROOT=/path/to/work \
W_A=6.36 H_A=6.30 W_B=4.00 H_B=5.98 W_C=5.90 H_C=5.80 FIG7D_FW=17.00 \
jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 Figure_7.ipynb
```

Those `W_*/H_*` values are the delivered geometry; with no environment at all the panels still render,
just on their default canvases. `assembly/make_fig7_v17_panels.py` is the same code as a plain script
(the notebook is generated from it by `build_fig7_notebook.py`, so the two cannot drift):

```bash
python3 build_fig7_notebook.py            # regenerates Figure_7.ipynb from assembly/make_fig7_v17_panels.py
```

Executed end-to-end on 2026-11-17 with **0 errors** (18 cells: 11 markdown, 7 code); the four PDFs it
writes are text-identical to the delivered panels and agree with them to 0.01 cm (A 6.58 × 5.90,
B 3.78 × 5.85, C 5.84 × 5.90, D 17.00 × 7.90). Versions: rendering kernel python 3.8.14 / numpy 1.22.4 /
pandas 1.5.1 / scipy 1.10.1 / statsmodels 0.13.5 / anndata 0.9.2 / matplotlib 3.6.2; the sweeps that
produced every input table ran under `$DSCI_PY` (python 3.9.23, **scdrs 1.0.2**,
anndata 0.10.9, matplotlib 3.9.4) — §0.2 prints both sets.

`upstream/*.py` is English-only: the 262 Chinese comment / log lines were translated with a checked
one-for-one line map (`<work>/tools/englishize.py`, `englishize2.py`). The translation was verified
by token-stream comparison (code tokens identical, format placeholders identical in every string, only
comment and string text changed), by `ast.parse` on all 15 files, and by a CJK scan over the whole bundle
including `tables_v17/`.

## Validation checklist

Every number below is printed by the notebook and matches the delivered tables bit-for-bit.

| Quantity | Notebook | Delivered |
|---|---|---|
| A: units per cohort (case/ctrl) | 13/8, 4/5, 14/3, 5/3, mouse 13/6 | `fig7A_donors_delta1.csv` |
| A: between-unit *P* | 0.500, 1.000, 0.091, 0.393; mouse overall 0.009 | `fig7A_between_donor_P.csv` (recomputed and asserted in-cell) |
| B: MS ranks on the face (size-adjusted) | 1 / 1 / 3 / 1 | `seven4_confound_ALL_cohorts.csv` `rank_adj_*` |
| B: MS raw ranks (sensitivity, Fig S7) | 2 / 4 / 7 / 2 | idem `rank_raw_*` |
| B: confound ρ (Δ1, ln n_peak) | −0.838 / −0.691 / −0.778 / −0.873 | idem (Fig S7 legend; not on the face) |
| C: axis ranks of MS/SCZ | identical to the sweep tables' `r1/rA/rB` | four sweep tables |
| C: GSE231585 neuron pool | 717 of 62,911; nB = 13 vs n1 = nA = 14 | `trait_sweep_delta2_GSE231585_74.csv`, `logs/d2_0.log` |
| D: red box | S1PR5 only, on the pooled GSE301908 rates | `S15_det_GSE301908.csv` / `S15_n_GSE301908.csv` |
| D: "18 of 80 scored" | 83 ledger rows minus 3 non-drug atoms | `S26_target_ledger.csv` (asserted in-cell) |
| B: confound rho | −0.838 / −0.691 / −0.778 / −0.873, recomputed in-cell | `seven4_confound_ALL_cohorts.csv` (`A2..A7` vs `lnp`) |

The donor-level re-test of Panel D's detection rule is deliberately **not** part of Figure 7: a
per-donor robustness check of a pooled detection statistic belongs to the supplement. It was run on
2026-11-17 and handed over with its scripts and tables to `figS7/handoff/donor_s1pr5_20261117/`.

## What is *not* here

- The composed sheet, panel letters and the final crop (journal assembly, see the composed PDF path above).
- The QA gates (font-tier audit, collision, text-on-data, run-on): they run on the composed page.
- The raw-rank band, the size-stratified ladder, the ±15 % matched-window control and the confound-ρ
  strip: Fig S7.
- The paired multiome (ATAC-channel) probe on GSE250589, which stayed negative at the targeted-trait
  level (12 traits, cross-channel Spearman ρ = −0.12 raw / 0.06 stratified by regulatory-domain size)
  and is reported in the supplement.
- The mouse 74-trait sweep in Panel C: its neuronal layer holds 78 neurons, so the neuron−glia axis is
  not computed there; it is reported in S7.

## Directory

```
Figure_7.ipynb                 the deliverable notebook (generated; run it, do not hand-edit the code cells)
build_fig7_notebook.py         generates the notebook from assembly/make_fig7_v17_panels.py
assembly/make_fig7_v17_panels.py   panel code, one "# %% PANEL X" block per panel  (edit here)
assembly/arial_style.py        the manuscript-wide font policy (shared copy; registers Arial/Liberation)
assembly/README.md             what runs here, what stays out, and the decision history of this figure
data/                         vendored inputs for Panels A/B/C + the expansion extract, laid out like
                              routeA_out/; MANIFEST.csv gives bytes/rows/md5/producer/panel per file
tables_v17/                    vendored small inputs for Panel D (producers: the Fig 7D audit chain)
upstream/                      the scripts that produced every input table (see upstream/README.md)
fig7_out/                      the four panel PDFs written by the notebook
_archive_v13_v16/              superseded single-column scripts/outputs and the pre-v17 notebook, README
                               and builder; kept on disk for provenance, excluded from the push by the
                               repository .gitignore (`fig*/_archive*/`)
```
