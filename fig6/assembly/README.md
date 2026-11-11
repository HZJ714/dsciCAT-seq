# assembly/ — page composition and the gates the delivered Figure 6 passes

`Figure_6.ipynb` produces the four panels. The published page is assembled from them here, and the
same scripts are the QA gates that the delivered PDF passes. None of this re-runs an analysis: every
number on the page comes from the tables the notebook wrote into `../fig6_out/`.

## Run

```bash
cd fig6
DSCI_PYLIBS=/path/to/pymupdf DSCI_FONTS=/path/to/fonts \
  python3 assembly/compose_fig6_v21.py            # 2x2 page -> fig6_out/fig6_composed_v21.pdf
DSCI_PYLIBS=/path/to/pymupdf \
  python3 assembly/make_fig6_S6m_scatters.py      # supplementary S6m, from the same two D tables
  python3 assembly/font_gate.py fig6_out/fig6_composed_v21.pdf fig6_out/fig6S6m_crosschem_scatters.pdf
  python3 assembly/check_text_overlap.py fig6_out/fig6_composed_v21.pdf
  python3 assembly/A6_fig6_D_occlusion_audit.py fig6_out/fig6D_heatmap_pair_compact.pdf --tall 0
  python3 assembly/A6_fig6_D_rank_contrast.py
```

## Files

| File | What it does |
|---|---|
| `compose_fig6_v21.py` | splices A/B/C/D into the 2 × 2 page, stamps the panel letters, and reports the row-width/row-height mismatch it accepted. Reads `../fig6_out/*.pdf` (override with `FIG6_OUT`, `DPDF`, `OUTPDF`). |
| `make_fig6_S6m_scatters.py` | the supplementary per-trait cross-chemistry scatters, drawn from `fig6D_crosschem_ladder.csv` + `fig6D_delta1_74traits_per_tech.csv`. |
| `font_gate.py` | submission gates on the PDF text layer: every visible span ≥ 5.5 pt, Arial (or metric-compatible Liberation Sans) only, page not wider than 18 cm. Exit 1 on failure. |
| `check_text_overlap.py` | text-collision gate: reads span bboxes and reports pairs on the same line (≥ 55 % vertical overlap) that overlap horizontally by more than 0.35 pt. Vertical text is skipped, otherwise a rotated y-label is a false positive. Column widths estimated by eye always miss this, so it is run on the rendered page. |
| `A6_fig6_D_occlusion_audit.py` | measures, at 1200 dpi, how much of each rank cell is covered by white ink that should be colour, and whether any printed digit escapes its own row. It recomputes the expected cell colour from the delivered table and the colormap, so it audits the *page*, not the drawing call. `--tall` is the height of the marked rows (0 for the delivered no-digit version). |
| `A6_fig6_D_rank_contrast.py` | the sequential-ramp legibility measurement quoted in the legend: samples a ramp into 74 levels and reports the adjacent-level Lab ΔE\*76 against the ≈2.3 just-noticeable difference. Prints the table for every candidate (`cividis_r` 2.13 / 40 % above JND; the previous blue-grey ramp 1.09 / 0 %). `--pick NAME` dumps a ramp's 74 colours. |
| `A6_fig6_v21_verify.py` | the full author-side acceptance script: font floor, all-black text outside the labelled traits, cross-panel colour identity (C vs D), the "every number printed on the page traces to a delivered table" check, the caption-catchment list and the S6m font gate. It reads the working tree (tables, legend file), so it is shipped for provenance rather than as a repo-runnable test — `font_gate.py` + `check_text_overlap.py` + `A6_fig6_D_occlusion_audit.py` are the repo-runnable subset. |
| `arial_style.py` | the house matplotlib/PyMuPDF style module (Arial registration, `use_arial`, spine helpers, `patch_fitz_arial`, the renamed PDF font copies). Overlaps with `../../00_Figure_Tools/figkit.py`; kept because these scripts were written against its API. |

## Result on the notebook's own output (2026-11-03)

```
fig6_composed_v21.pdf        17.86 x 11.13 cm | 102 spans | 5.50-8.00 pt | Arial-BoldMT, ArialMT   PASS
fig6D_heatmap_pair_compact   10.78 x  5.39 cm |  23 spans | 8.00 pt      | ArialMT                 PASS
fig6S6m_crosschem_scatters   15.28 x  5.85 cm |  53 spans | 5.50-6.80 pt | Arial-BoldMT, ArialMT   PASS
text collision                    0 pairs
rank cells covered by white ink   0 rows
```

The delivered manuscript page measures 17.93 × 11.18 cm (`fig6_out/fig6_composed_v21_delivered_20261103.pdf`
is kept next to the reproduction for diffing). The 0.07 × 0.05 cm difference is entirely in the tight
bounding boxes of panels A and B; panel D is identical to the millimetre, and both pages pass the same
gates.

## Two design decisions that are recorded here because they were measured, not chosen

- **No rank digits on panel D.** 74 rows occupy 4.3 cm (0.058 cm = 1.64 pt per row) while the
  submission floor is 5.5 pt: a legible digit is necessarily taller than its row, so "digits inside the
  cells" and "no occlusion of neighbouring cells" cannot both hold at equal row height. Three variants
  were built and audited — a white underlay behind the digit (measured: 3 non-labelled rows polluted by
  5–13 % white, digit ink overflowing its row by 20 px above and 19 px below), taller marked rows with
  the rank cell left white (0 pollution, but the author rejected the look), and no digits at all
  (delivered; 0 pollution, ranks in the legend and the text).
- **`cividis` reversed for the rank ramp.** Once the digits are gone, adjacent-row colour difference is
  the only channel carrying rank; see `A6_fig6_D_rank_contrast.py` for the ΔE\*76 table that decided it.
