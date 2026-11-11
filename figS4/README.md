# Figure S4 — time axes, peak-to-gene links and gene-set structure behind Fig 4

This notebook reproduces the four panels of **Supplementary Figure 4** from
the upstream analysis objects.

## Panels

| Panel | File | Size | Content |
|---|---|---|---|
| A | `figS4v3a_umap.pdf` | 17.00 × 4.20 cm | MultiVelo latent time on the Fig 4A embedding (blueYellow, same ramp as Fig 4A) + cell-type composition of the latent-time axis (20 equal-width bins) |
| B | `figS4v3b_p2g.pdf` | 17.00 × 3.40 cm | The 15,607 peak-to-gene links: r (floor 0.40), peak–TSS distance (cap 250 kb), peaks per gene (S1 markers in orange) |
| C | `figS4v3c_funnel.pdf` + `figS4a_model_heatmaps.pdf` + `figS4v3c_dec.pdf` | 2.30 + 7.55 + 4.20 cm | Gene-set funnel (1,809 → 372 / 114, overlap 32, **not nested**); model-1 (1,798 genes, k-means programme order) + model-2 (11 genes) heat maps on the latent-time axis (high end labelled **MOL3/6**, not MOL6 — 2026-10-28); Nfasc / Ndrg1 direction examples |
| D | `figS4v3d_3d.pdf` | 17.00 × 8.40 cm | 3D s–u–c phase portraits for 8 genes, subtype colours, fitted trajectory (black) and switch times (magenta); non-cubic widened boxes + `tc50/ts50` definition line under the figure (2026-10-28) |

## Pipeline

```
multivelo_dynamical.h5ad ──> Panels A, C2 (model 2), D
multivelo/results/decoupling.tsv + decoupling_direction.tsv ──> Panels C1, C3
scmega_p2g.csv + scmega_p2g_S1markers.csv ──> Panel B
multivelo/plots/heatmap__model1__1798genes.pdf (wordless raster) ──> Panel C2
blueYellow ramp: paletteContinuous("blueYellow") (R) ──> Panel A
```

The notebook draws everything from these objects; the only read-back step
parses an intermediate PDF the notebook itself produces (the MultiVelo model-1
heat map raster, which is a wordless image with rows/columns documented in the
caption).

## Inputs (upstream analysis objects)

| File | Produced by | Used in |
|---|---|---|
| `multivelo_dynamical.h5ad` | MultiVelo `recover_dynamics_chrom` (5977 × 1809) | A–D |
| `multivelo/results/decoupling.tsv` + `decoupling_direction.tsv` | LRT refit | C |
| `scmega_p2g.csv` + `scmega_p2g_S1markers.csv` | `SelectGenes()` (15,607 links) | B |
| `multivelo/plots/heatmap__model1__1798genes.pdf` | MultiVelo model-1 heat map | C2 |
| blueYellow ramp | `paletteContinuous("blueYellow")` (R) | A |

Set the `DATA` variable in Section 0 of the notebook to point to a directory
containing the same layout.

## Outputs

Written to `figS4_out/`:

| File | Size |
|---|---|
| `figS4v3a_umap.pdf` | 17.00 × 4.20 cm |
| `figS4v3b_p2g.pdf` | 17.00 × 3.40 cm |
| `figS4v3c_funnel.pdf` | 17.00 × 2.30 cm |
| `figS4a_model_heatmaps.pdf` | 17.00 × 7.55 cm |
| `figS4v3c_dec.pdf` | 17.00 × 4.20 cm |
| `figS4v3d_3d.pdf` | 17.00 × 8.40 cm |

## Requirements

`numpy pandas scipy matplotlib pillow anndata pymupdf(fitz) multivelo`
(only `predict_exp` from `multivelo.dynamical_chrom_func` for the C3 curves)
plus `Rscript` with `scMEGA` / `ArchR` for the blueYellow ramp. Package
versions are printed in Section 0.2.

## Validation

Validated against the journal-submission panels at 150 dpi (ink-disagreement:
A 0.009%, B / C1 / C2 / C3 / D all 0.000 — effectively pixel-identical).

> **Tracked wording item (open).** The C3/D panel headers still print the raw
> column name `LRT`. The delivered review asked to relabel these ten strings
> `2ΔlogL` to match the main figure's `signed 2ΔlogL` axis; this requires a
> re-render of S4 (in-place PDF patching was tried upstream and is not safe
> for these overlapping small-text rows), so it is **not yet applied** — the
> figure caption carries the definition in the meantime. `FDR` values that
> underflow to 0 print as `<1e-300` (2026-10-28, applied here).
