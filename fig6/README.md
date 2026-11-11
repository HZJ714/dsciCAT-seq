# Figure 6 — Cross-channel (expression vs chromatin) concordance of the Oligo–OPC contrast

`Figure_6.ipynb` reproduces the delivered **Figure 6 (v36 final)** from the analysis objects
upstream of it: the scDRS per-nucleus scores of the own dsciCAT-seq multiome data (both channels
of the *same* nuclei). The delivered figure is one row — **A | B | C** — A the MS score across
the six subtypes, B the same score on the UMAP, C the 74-trait cross-channel Δ1 ladder. Per the
author instruction the notebook stops at the three panels (`fig6_out/fig6v36_{A_score,B_umap,
C_ladder}.pdf`, all at the delivered sizes); the 17.00 × 4.89 cm row
(`Figure_6_v36_ABC_row.pdf`, letters, page-width gate) is composed by the delivered driver
`assembly/fig6v36_abd_final.py` (run by `run_all.sh`), not by the notebook. The gene-level
landing panel that used to be main-figure C was pulled in v36 (no defensible positive);
`assembly/` holds the delivered render script plus the v21-era QA gates for the superseded 2 × 2.

The superseded v21 notebook and its outputs are kept in `_superseded_20261023/`.

## Panels

| Panel | File | Size | Content |
|---|---|---|---|
| A | `fig6v36_A_score.pdf` | 5.82 × 4.54 cm | scDRS MS `norm_score` per nucleus across the six Fig 2 subtypes; boxes = median/IQR, whiskers = 5–95 %, jittered dots = nuclei, MWU *P* vs COP/NFOL above |
| B | `fig6v36_B_umap.pdf` | 5.76 × 4.53 cm | the same score on the UMAP of the same nuclei (hexbin of binned mean, aspect locked to the data's 1.63 ratio) |
| C | `fig6v36_C_ladder.pdf` | 4.87 × 4.55 cm | 74 GWAS traits: Δ1 in the chromatin channel (x) vs the expression channel (y), own same-nucleus data; MS / SCZ / HEIGHT named |
| sheet | `Figure_6_v36_ABC_row.pdf` | **17.00 × 4.89 cm** | the three panels spread to exactly 17.00 cm, letters stamped, two-tier type 8.0/6.0, zero vector scaling |

## Pipeline

```
3k_dsci_oligo_subset.RDS ── rds_to_dump.R + build_h5ad.py ──┐
3k_dsci_ht_10x_gs_subset.RDS ── export_genescore.R + build_genescore_h5ad.py ─┤  (upstream/D0, verified bit-for-bit)
                                                                            ▼
own dsciCAT-seq multiome (Fig 2 annotation: 5,977 lineage nuclei)
  ├── obs.cell_type + obs.UMAP  ─┐
  ├── scDRS MS score (RNA)  ─────┴──> A, B
  ├── scDRS 74 traits, RNA channel  (full74/)   ─┐
  └── scDRS 74 traits, gene-score channel (out_gs74/) ┴──> C   (same nuclei, Δ1 per trait)

three chemistries (own dsciCAT-seq + 10x Multiome + HT-scCAT-seq, mouse OL lineage)
  ├── cells_<tech>.csv  (Fig 2 Oligo / OPC / COP-NFOL labels)
  ├── scdrs/<tech>_<channel>/**/*.score.gz   (74 traits x 2 channels x 3 techs = 444 runs)
  └── <channel>_<tech>_genes.txt  (measurable genes -> n_hit domain-size proxy)  ──> D
```

## Inputs

All inputs are read through three environment variables; nothing is downloaded by the notebook.

| Variable | Object | Produced by |
|---|---|---|
| `DSCI_FY_ROOT` | `oligo/3k_dsci_oligo_subset.h5ad` | Fig 2 lineage curation of the own dsciCAT-seq multiome, exported as a Seurat RDS; converted by `upstream/rds_to_dump.R` + `upstream/build_h5ad.py` (X = RNA `data` slot, `log1p(counts/colsum·1e4)`) |
| `DSCI_FY_ROOT` | `oligo/fig6_rerun_20261019/dualomics_probe/genescore_oligo_5977.h5ad` | Signac `gene_score` **data** slot (= `log1p(ArchR gene score)`) of `3k_dsci_ht_10x_gs_subset.RDS`, same 5,977 nuclei; converted by `upstream/export_genescore.R` + `upstream/build_genescore_h5ad.py` |

Both conversions are re-run and checked element-wise by `upstream/D0_own_lineage_h5ad.sh`; they
reproduce the delivered files exactly (`indptr`/`indices`/`data`, `obs`, `var` all identical).
| `DSCI_FY_ROOT` | `oligo/scdrs_ms_out/PASS_Multiple_sclerosis.score.gz`, `*.scdrs_group.cell_type` | `scdrs compute-score` (see `upstream/README.md`) |
| `DSCI_FY_ROOT` | `oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs` | MAGMA 10 kb / 1,000-gene gene sets for the 74 GWAS traits |
| `DSCI_FY_ROOT` | `oligo/scdrs_traits_draft/full74/`, `.../dualomics_probe/out_gs74/` | scDRS, 74 traits, expression and chromatin channel of the same nuclei |
| `DSCI_FIG_ROOT` | `routeA_out/D_crosschemistry/{tables,derived,scdrs}` | `upstream/` chain D1 → D2c/D2e/D2f → D4 → D5 |
| `DSCI_HOMOLOGENE` | `homologene.data` | Homologene (human GWAS gene sets → mouse symbols, for `n_hit`) |
| `DSCI_FONTS` | `*.ttf` | Arial (the manuscript-wide face; panel C's label anchors are tuned to its metrics) |

The 444 scDRS score files are ~6 GB and are not in the repository. The notebook reads `norm_score`
only and caches each of the six matrices as `fig6_out/scores_<job>.npz` (2–6 MB each), so the first
run spends most of its time on the gzip reads and every later run is ~2 min.

## Running

**Use the python 3.9 kernel.** The delivered panel sizes are matplotlib tight-bbox outcomes; under
matplotlib 3.6.2 (python 3.8) the three panels come out 0.05 cm short with otherwise identical ink, and
the notebook now refuses to run there (assert in the version cell). Register the kernel once, then:

```bash
mkdir -p ~/.local/share/jupyter/kernels/py39 && cat > ~/.local/share/jupyter/kernels/py39/kernel.json <<'EOF'
{"argv": ["<your analysis-environment python>", "-m", "ipykernel_launcher", "-f", "{connection_file}"],
 "display_name": "py39 (matplotlib >= 3.9)", "language": "python"}
EOF

cd fig6
DSCI_FY_ROOT=/path/to/FY DSCI_FIG_ROOT=/path/to/figures DSCI_FONTS=/path/to/fonts \
  jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.kernel_name=py39 \
  --ExecutePreprocessor.timeout=-1 Figure_6.ipynb
```

Executed end-to-end with **0 errors** on 2026-11-18 in that kernel (python 3.9.23, matplotlib 3.9.4),
reproducing the three delivered panel PDFs at 5.82 x 4.54 / 5.76 x 4.53 / 4.87 x 4.55 cm with identical
word sets and identical drawing counts. The earlier 2026-11-03 note recorded a clean run under python
3.8.14 / matplotlib 3.6.2; that run is clean but geometry-short, which is why the kernel is now gated.
The scdrs scoring upstream ran in a private conda env (python 3.9.23, **scdrs 1.0.2**); the public tree addresses that interpreter as `$DSCI_PY` —
the notebook itself does not import scdrs, and §0.2 prints the local versions plus this registration.

## Validation checklist

Every number below is printed by the notebook. The "delivered" column is the value in the tables the
published figure was drawn from, and the reproduction is bit-identical except where noted. The §D rows
were refreshed on 2026-11-18 to the corrected post-log1p chromatin tree; they had still carried the
pre-fix values (0.726 / 0.763 / 0.518, "2 of 3") after the 2026-11-04 re-run.

| § | Quantity | Notebook | Delivered |
|---|---|---|---|
| A | nuclei joined | 5,962 of 5,977 | 5,962 |
| A | n per subtype | 203 / 490 / 1,016 / 688 / 991 / 2,574 | identical |
| A | median `norm_score` per subtype | −0.4319 / −0.0502 / 0.0314 / −0.2806 / −0.0144 / 0.0658 | identical |
| A | MWU *P*, group *P*, BH *q*, Cliff δ vs COP/NFOL | see `fig6A_stats.csv` | max abs diff 0 (δ: 1e-16) |
| B | hexbin colour limit | ±0.738 | ±0.74 |
| C | same-nucleus pool / traits | 5,962 / 74 | 5,962 / 74 |
| C | Δ1 per trait, both channels, ranks, off-diagonal *z* | `fig6C_delta1_own74.csv` | max abs diff **0** |
| C | ladder ρ | 0.8191 | 0.8191 |
| C | trait absolute level ρ (control) | 0.0291, *P* = 0.8058 | 0.0291 / 0.8058 |
| C | ρ after removing trait level | 0.6790 | 0.6790 |
| C | permutation *P* of the raw ladder | 5.0e-04 (2,000 draws) | 5.0e-04 |
| C | noise ceiling √(rel_RNA·rel_ATAC) | 0.8756 | 0.8776 — **±0.002**, see note below |
| C | six-subtype profile ρ for MS | −0.1429 (*P* = 0.787) | −0.1429 |
| C | per-nucleus cross-channel ρ | 0.0790, empirical *P* = 0.0033 | 0.0790 / 0.0033 |
| C | rank labels covering another point | 0 (asserted on the rendered figure) | 0 |
| D | ρ_adj per chemistry | 0.8959 / 0.8240 / 0.6300 | identical (diff 0) |
| D | permutation *P* | 0.0001 in all three (the numerical floor at 10,000 draws) | identical |
| D | bootstrap 95 % CI | 0.806–0.944 / 0.700–0.902 / 0.451–0.768 | identical |
| D | pre-specified ladder verdict | pass / pass / pass (**3 of 3** chemistries) | identical |
| D | MS ranks (expression → chromatin, domain-adjusted) | 7→37 / 4→49 / 11→41 | identical |
| D | MS off-diagonal *z* | +2.48 / +2.62 / +0.89 | identical |
| D | MS counter-example reproduced | True in all three | identical |
| D | panel size | 10.78 × 5.39 cm | 10.78 × 5.39 cm |
| page | `assembly/compose_fig6_v21.py` | 17.86 × 11.13 cm | 17.93 × 11.18 cm (A/B tight boxes, ≤0.07 cm) |
| page | `assembly/font_gate.py` | min 5.50 pt, Arial only, 0 spans below the floor | same gates |
| page | `assembly/check_text_overlap.py` | 0 collisions | 0 collisions |
| D | `assembly/A6_fig6_D_occlusion_audit.py` | 0 rank cells covered by white ink | 0 |

**The one non-identical quantity** is the split-half noise ceiling of panel C. It is a median over 20
random draws; the delivered value came from an RNG stream that had already been consumed by the other
tests of that run, so it cannot be re-sampled bit-for-bit from a fresh notebook. The two values differ
by 0.002, against an inter-quartile spread of ±0.01 in the same statistic. Panel D's ceiling *is*
bit-identical, because that chain's helper draws from a dedicated seed.

## What the figure does and does not show

- Δ1 is a **contrast**, not an absolute score: the trait-level control (ρ = 0.029, *P* = 0.81) says the
  two channels' absolute levels are unrelated; only the lineage contrast reproduces.
- The row order of panel D is the ascending mean of the three **expression**-channel adjusted ranks.
  The chromatin channel never enters the ordering, so the smooth look of the expression rank columns is
  by construction and is not a result.
- Δ1 columns are coloured by the **within-column z** of Δ1 (coolwarm, ±1.4), because absolute Δ1
  magnitudes are never compared across technologies or species; a trait with a positive Δ1 can still
  print blue.
- **No rank digits are printed on panel D.** 74 rows occupy 4.3 cm (0.058 cm ≈ 1.64 pt per row) while
  the submission floor is 5.5 pt, so any legible digit necessarily covers neighbouring rows. The ranks
  of the three labelled traits are stated in the legend and the main text.
- The rank ramp is `cividis` reversed (darkest = rank 1). Measured as adjacent-level Lab ΔE\*76 over the
  74 quantised levels, the previous blue-grey ramp has a median of 1.09 with **0 of 73** adjacent pairs
  above the ≈2.3 just-noticeable difference; cividis reversed reaches 2.13 with 40 % above JND, is
  monotonic in lightness and colour-vision-deficiency safe. `00_Figure_Tools`-style measurement:
  `assembly/rank_contrast.py` prints the table for any candidate ramp.
- The ladder verdict below is stated for the **corrected** chromatin tree (`tables_log1p_20261104/`,
  the 2026-11-04 re-run): all three chemistries clear the pre-specified threshold
  (ρ_adj = 0.896 / 0.824 / 0.630, all ≥ 0.6 with permutation *P* = 1.0e-04), so the claim is
  **3 of 3 chemistries reproduce the ladder**. Before the input-transform fix the same statistic was
  0.726 / 0.763 / 0.518 and HT-scCAT-seq failed; the older "2 of 3" sentence in the manuscript
  draft is the author's to retire (see `upstream/README.md`, "Corrected rerun result"). The
  MS counter-example (chromatin ranks MS lower than expression) is reproduced in all three either way.
- Panel C and panel D are **not two runs of one measurement**. They share 5,957 nuclei; C scores 5,977
  nuclei of the Seurat lineage object with the six-subtype labels (Δ1 = MOL6 − COP/NFOL), D scores
  7,523 nuclei of the ArchR project with the three-class labels (Δ1 = Oligo − OPC), and 1,488 of D's
  1,490 OPC nuclei are not in C at all. Their chromatin matrices are different objects too (logged
  Signac gene score vs un-logged ArchR gene score). The two Δ1s are quoted separately and never pooled.
- **Open issue, panel-D chromatin side — found, audited, fixed in code; rerun running 2026-11-04.**
  scDRS is called with `--flag-raw-count False`, which applies `expm1()` internally and so requires
  log-transformed input. The expression matrices satisfy that, and so does panel C's chromatin matrix:
  `upstream/A6_fig6_input_transform_audit.py` walks every scDRS tree Figure 6 uses, recovers the `.h5ad`
  its launcher named and counts warning lines (`tables/fig6_input_transform_audit_20261104.csv`):
  A/B/C are clean (0 warnings in 25/25 batches, input max 8.02 / 6.47). Only the three chromatin
  matrices of panel D were 1e4-normalised but un-logged (max 1,388.6 own / 297.1 10x / 372.7 HT): the
  own arm overflowed to `inf` in 31/31 run logs, HT in 25/25 (`overflow in multiply`), 10x failed
  silently at ~1e129. `upstream/D4b_pack_atac_log1p.py` builds the corrected `log1p` matrices and
  `upstream/preflight_h5ad.py` is the gate that was missing — `D5_run_scdrs.sh` now refuses to launch
  an arm the gate rejects and post-checks every batch log. The 222 chromatin runs are being re-scored
  into `scdrs_log1p_20261104/` and analysed into `tables_log1p_20261104/`; panels A/B/C, the expression
  arms, the layout and every page gate are unaffected. The committed figure and tables still show the
  pre-rerun values until the before/after comparison is reviewed.
- Wording red line carried over from the audit: no "benchmark", "head-to-head", "outperforms" or
  "more sensitive than" — the supported statement is *reproducible across independent experiments and
  library chemistries*.

## Layout of this directory

- `Figure_6.ipynb` — the executed notebook (0 errors)
- `build_fig6_notebook.py` — the generator that produced it (edit this, re-run, notebook rewritten)
- `fig6_out/` — panel PDFs, the tables the panels are drawn from, and the `scores_*.npz` caches (gitignored)
- `upstream/` — the scripts that produce the notebook's inputs (ArchR/Seurat export → h5ad → 444 scDRS runs → tables), with `upstream/README.md`
- `assembly/` — the 2 × 2 composition and the QA gates the delivered page passes (font floor, text-collision, number traceability, heatmap occlusion audit)
- `_superseded_20261023/` — the previous Figure 6 notebook (A/B/Manhattan/74-trait strip) and its outputs

## Data availability

All inputs of Figure 6 are the authors' own libraries: dsciCAT-seq, 10x Multiome and HT-scCAT-seq
of the same mouse brain oligodendrocyte lineage. Raw FASTQ and raw count matrices are deposited at
**CNGBdb, project CNP0010263**; the upstream FASTQ → matrix workflows are in
[`../00_Preprocessing/`](../00_Preprocessing/). No third-party data set enters this figure.
