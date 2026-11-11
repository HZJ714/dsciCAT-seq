# Figure S1 — dsciCAT-seq enables consistent library quality across four barcode numbers and two throughput scales

`Figure_S1.ipynb` reproduces the measured content of Figure S1 from the files the
analysis started from: eight dsciCAT-seq libraries (4 barcode numbers x 2 capture
depths), each a HEK293T : NIH3T3 mixture. Journal assembly (the 18.00 x 18.40 cm
page, panel letters, the shared legend strip, the compositing) is a production
step and is deliberately **not** part of the notebook.

Executed 2026-10-22 in the environment that drew the delivered figure: every cell
ran, ~3 min in total.

## Panels

| Panel | What it shows | Section | Output | Reproducible? |
|---|---|---|---|---|
| A | RNA-channel cell counts + per-cell log10(UMI), log10(detected genes) | §1 | `S1a_RNA_UMI.png`, `S1a_RNA_genes.png` | yes, seconds |
| B | ATAC-channel cell counts + per-cell log10(unique fragments) | §2 | `S1b_ATAC_QC.png` | yes, seconds |
| C | 8 x 8 RNA profile correlation, human and mouse | §3 | `S1c_corr_{Human,Mouse}.csv`, `S1c_RNA_similarity.png` | yes, ~2 min |
| D | species call of the same nuclei in both channels, 2 x 4 facets | §4 | `S1d_joint_mix.png` | yes, seconds |

## The one step that is real computation

Panels A, B and D read per-cell tables that CellRanger / Solo / ISSAAC had already
written, and only aggregate them. Panel C is rebuilt from scratch: `matrix.mtx` +
`features.tsv` per library and genome (16 files, 3.6 GB), gene totals, eight-way
ENSG intersection (36,601 human / 32,285 mouse genes), log2 CPM, Pearson.

Two things about those matrices are easy to get wrong and are handled explicitly
in §3:

* they are stored **cells x genes**, the transpose of the CellRanger default, so
  gene totals are `sum(axis=0)` and the `features.tsv` check compares against
  `mtx.shape[1]`;
* seven libraries carry ENSG ids in `features.tsv` and one (`HG13`, Solo 0.37)
  carries symbols only. `HG13` is not part of the eight delivered libraries, so
  the symbol->ENSG path in §3 is a guard, not a workaround: it converts through
  the one-to-one map and stops if the mapped ids are not in reference order.

## Inputs (`FIGS1_DATA`, default `<raw>/HG`)

| Used by | Path (relative to the data root) |
|---|---|
| §1 | `<run>/sub/RNA/{human,mouse}_metrics.csv` (`HG9`: `RNA-sub0.3/`) |
| §2 | `<run>/sub/ATAC/{Human,Mouse}/metrics.csv` (`HG9`: `ATAC-sub0.33/`) |
| §3 | `<run>/sub/RNA/{Human,Mouse}/{matrix.mtx,features.tsv}` |
| §4 | `<run>/sub/ISSAAC_joint_metrics_sub.csv` (`HG9`, `HG11`: at the run root) |
| §3 guard | `HG13/sub/outs/Solo.out/Gene/filtered/features.tsv`, `HG11/sub/RNA/{Human,Mouse}/features.tsv` |

Raw data are deposited at CNGBdb under accession `CNP0010263` and are not in this
repository.

| condition | run | | condition | run |
|---|---|---|---|---|
| `50k-8bc` | HG9 | | `100k-8bc` | HG11 |
| `50k-16bc` | HG12 | | `100k-16bc` | HG14 |
| `50k-24bc` | HG15 | | `100k-24bc` | HG9_1 |
| `50k-48bc` | HG10_1 | | `100k-48bc` | HG12_1 |

## Run

The figures are **saved, not displayed**: each section writes a PNG to `FIGS1_OUT`
and prints the path. The numbers are the manuscript's; the colours, sizes and
fonts are the shortest form that draws the same measurement and are deliberately
not the assembled page. Check the tables the sections print, not the artwork.

```bash
python build_figS1_notebook.py            # regenerate the notebook

jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=-1 Figure_S1.ipynb
```

Outputs go to `./figS1_out/` (`FIGS1_OUT` overrides).

## Validation against the delivered figure

| Check | How | Result |
|---|---|---|
| §1 cell counts and medians | the same `*_metrics.csv` the delivered panel used | 10,738-20,884 cells (HEK293T) and 10,764-22,610 (NIH3T3) per library; median log10(UMI) 3.15-3.36, log10(detected genes) 2.88-3.13 |
| §2 cell counts and medians | as above, ATAC channel | 13,833-27,856 (HEK293T) / 15,483-33,311 (NIH3T3) nuclei; median log10(unique fragments) 4.25-4.38 |
| §3 correlation matrix | recomputed from the matrices, compared with `scripts/cache/S8_corr_{Human,Mouse}.csv` (the CSV the delivered heatmap was drawn from) | unrounded max abs delta **4.98e-4** human / **4.96e-4** mouse; the 3-decimal CSV the notebook writes is **identical on 64/64 cells** |
| §3 gene intersection | recomputed | 36,601 human / 32,285 mouse genes |
| §3 dynamic range | recomputed | human r 0.863-0.880 (median 0.871), mouse 0.871-0.886 (median 0.878) |
| §4 the eight legend lines | `value_counts` of `species_joint` | 22,387-46,355 nuclei, Human 44.2-48.3 %, Mouse 44.4-48.0 %, Mix 4.94-10.55 % — identical to the delivered panel |

## Not reproduced here

* Page geometry and compositing: `figS1_layout.py` + `compose_figS1.py`
  (18.00 x 18.40 cm, four slots, Arial 8 pt regular panel letters, 10 pt on this
  page's exception list).
* The shared legend strip and the y-label auto-alignment.
* The **12-library** variant (`FIGS1_NLIB=12` in the production scripts). That is
  an experimental layout; the delivered and reproduced figure is the eight-library
  version, and this notebook is hardcoded to it.
