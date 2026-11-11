# Figure 4 — Chromatin–RNA dynamics along the oligodendrocyte lineage

This notebook reproduces **Figure 4** of the paper, panel by panel, from the
upstream analysis objects.

## Panels

| Panel | File | Size | Content |
|---|---|---|---|
| A | `fig4a_traj_umap.pdf` | 6.20 × 6.20 cm | scMEGA joint trajectory UMAP: 5,977 lineage cells coloured by pseudotime (0–100, blueYellow), built in R from the Seurat object (`TrajectoryPlot`) |
| B | `fig4b_dual_heatmap.pdf` | 9.80 × 6.20 cm | Dual-omics heatmap: chromatin accessibility (top) and gene expression (bottom). The x axis is the 15,607 peak-to-gene links in row-dendrogram order and the body is rotated 90°, so pseudotime runs vertically; 14 stage markers on the side |
| C | `fig4d_models.pdf` | 16.00 × 4.90 cm | Model statistics: decoupling volcano (**signed 2ΔlogL**, not a true LRT — see §3), induction timing (`93.6% of 1798`), switch intervals printed as `x/114` |
| D | `fig4c_mv_dyn.pdf` | 16.00 × 7.80 cm | MultiVelo three-gene kinetics: **Anln / Cemip2 / Tcf7l2** × (ATAC, unspliced, spliced) |

> Note: the `fig4c_*` / `fig4d_*` filename prefixes predate the 2026-10-12
> layout swap; follow the table above for the layout letters (C = models,
> D = kinetics). Updated 2026-11-05 to the 2026-10-28 final panels (Codex
> comments #43/#81/#88/#89: example genes replaced, axis renamed to signed
> 2ΔlogL, explicit denominators, p99 display clipping).

## Pipeline

```
3k_dsci_oligo_subset.RDS  ──> R (Seurat/scMEGA) ──> Panel A
scmega_obj_target.rds + scmega_p2g.csv ──> R (ComplexHeatmap) ──> Panel B
multivelo_dynamical.h5ad + results/decoupling.tsv ──> Python ──> Panels C, D
```

The notebook draws every panel from these objects; the only read-back steps
parse intermediate PDFs that the notebook itself produces (the Fig 4A
TrajectoryPlot vector points, and the 05f-equivalent raster/colorbars/labels
for Panel B).

## Inputs (upstream analysis objects)

| File | Produced by | Used in |
|---|---|---|
| `3k_dsci_oligo_subset.RDS` | cell-type annotation of paired RNA+ATAC (Seurat, 5,977 cells) | Panel A |
| `scmega_obj_target.rds` | `scMEGA_oligo_GRN.R` (chromVAR + GRN) | Panel B |
| `scmega_p2g.csv` | `SelectGenes()` (15,607 peak-to-gene links) | Panel B |
| `multivelo/results/multivelo_dynamical.h5ad` | MultiVelo `recover_dynamics_chrom` (5977 × 1809) | Panels C, D |
| `multivelo/results/decoupling.tsv` | decoupling refit (`fit_decoupling=False`; the `LRT` column is a signed 2ΔlogL) | Panel C |

Set the `DATA` variable in Section 0 of the notebook to point to a directory
containing the same layout.

## Outputs

Written to `fig4_out/`:

| File | Size |
|---|---|
| `fig4a_traj_umap.pdf` | 6.20 × 6.20 cm |
| `fig4b_dual_heatmap.pdf` | 9.80 × 6.20 cm |
| `fig4c_mv_dyn.pdf` | 16.00 × 7.80 cm |
| `fig4d_models.pdf` | 16.00 × 4.90 cm |

## Requirements

`numpy pandas scipy matplotlib pillow scikit-learn anndata pymupdf(fitz)` and
`Rscript` with `Seurat Signac scMEGA ArchR ComplexHeatmap` (panels A/B).
Package versions are printed in Section 0.2.

## Validation

Validated against the journal-submission panels at 150 dpi (ink-disagreement:
A 0.19%, B 0.005% — 2026-10-21 check of the R-native path). Panels C/D were
re-synced on 2026-11-05 to the 2026-10-28 final drawing code and reproduce the
delivered panel PDFs byte-for-byte in content (text-layer checks: signed 2ΔlogL
axis + `1434/1809 have 2ΔlogL < 0 (P = 1)` line, `93.6% of 1798`,
`x/114` denominators, `Switch intervals (114 genes)`; Anln/Cemip2/Tcf7l2 with
the old Tnr/Plp1/Cnp/Mbp/Mal labels gone — Anln and Cemip2 are haloed vector
outlines, so they are absent from the text layer by design).
