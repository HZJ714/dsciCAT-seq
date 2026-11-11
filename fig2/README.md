# Figure 2 — Joint analysis of 16-barcode dsciCAT-seq with HT-scCAT-seq resolves 33 cell types and oligodendrocyte states

`Figure_2.ipynb` starts from the **SoupX-corrected RNA matrices** of the
dsciCAT-seq libraries and ends with the five panels: §1–§5 build the analysis
object, §6–§9 draw from it. Journal assembly (page geometry, panel letters,
shared colour key, compositing) is a production step and is deliberately **not**
part of the notebook.

## The chain (§1–§5)

**Only the dsciCAT-seq libraries run through this pipeline** — FY13 and FY14, the
16-barcode runs. The HT-scCAT-seq and 10x Multiome matrices join at §3 as
processed assay blocks; from §4 on every library is treated identically.

```
SoupX-corrected RNA matrix  +  peak matrix / fragment file
        │
        ├─ §1  per-library Signac ATAC object      85,255 peaks × 23,376 | 89,096 × 43,249
        ├─ §2  per-library multiome (RNA ∩ ATAC)      21,038 + 38,924 = 59,962 cells
        ├─ §3  merge across technologies              59,962 + 14,828 + 8,452 = 83,242
        ├─ §4  ArchR arrows, then its default filter   83,242 → 83,135 → 70,643 cells
        └─ §5  WNN integration (harmony)              3k_dsci_ht_10x_gs_subset.RDS

per library, §2 → §4:  10x 8,452 → 8,452 → 8,438 | HT 14,828 → 14,828 → 14,762
                       FY13 21,038 → 20,971 → 17,093 | FY14 38,924 → 38,884 → 30,350
```

§4 is the step the figure's cell count comes from. `createArrowFiles()` gives 107
of the merged cells no arrow record (67 FY13, 40 FY14 — not low-coverage cells;
their bead-level metadata reports ~1.3 × 10⁵ unique fragments each), and ArchR's
**default** filter (`minFrags = 1000`, `minTSS = 4`, nothing tuned) then removes
the 12,492 cells below `minTSS = 4`. Nothing more falls below `minFrags`: no
arrow cell is under it. The gene-score matrix ArchR writes on the survivors
(24,333 genes × 70,643 cells, every cell summing to 10,000) becomes the object's
third assay, `gene_score`.

## The step every panel depends on (§6)

The panels consume **one annotation**: the 33-type `cell_type_updated` column of
the object built in §5, in one fixed order (`dot_order`, 1 = AST … 33 = Ex_L6b)
with one fixed palette (`global_colors_full = unique(c(glasbey(32),
alphabet(26)))[1:42]`), so the UMAPs, the dot plot and the tracks all agree on
numbers and colours. §6 also writes `gs_subset_light.RDS` (~10 MB: embeddings,
metadata, the 48 marker rows), which is enough to redraw panels A–D.

| Panel | What it shows | Section | Output |
|---|---|---|---|
| A | UMAP in RNA space (Gene expression) | §7 | `fig2_umap_rna_celltype.pdf` |
| B | UMAP in ATAC space (Chromatin accessibility) | §7 | `fig2_umap_atac_celltype.pdf` |
| C | UMAP in the WNN-integrated space | §7 | `fig2_umap_wnn_celltype.pdf` |
| D | dot plot, 48 lineage markers × 33 types | §8 | `fig2_full_celltype_dotplot.pdf` |
| E | ATAC pseudo-bulk tracks at Qk, Gpr17, Plp1, Slc32a1, Slc17a7 | §9 | `fig2_track_Qk_Gpr17_Plp1_Slc32a1_Slc17a7.pdf` (5 pages) |

## Inputs

| Variable | Default | Contents |
|---|---|---|
| `FIG2_INPUT` | `<raw>` | workflow outputs: `FY/FY13`, `FY/FY14` (RNA matrices, `ATAC/out/Peak`, `ATAC/03.d2cfile/*.fragments.tsv.gz`) and `scCAT-seq-HT/Mbrain_process.rds` |
| `FIG2_DATA` | `<work>/FY` | analysis tree: `FY1x/SoupxMatrix` (the starting point), `FY1x/*_ATAC_object*.RDS`, `FY1x/*_RNA_soupx_PCA*.RDS`, the per-stage objects §2–§5 write, `archr-subset-copnfol/` (§4, §9) |
| `FIG2_OUT` | `fig2_out` | panel PDFs/PNGs, the derived cache, the per-stage hand-off files |

Each stage loads its own output when it exists, so the notebook is both a
document and a resumable pipeline.

## Provenance: read this before comparing panels

* **What is computed here.** FY13/FY14 go through §1–§5: intersection with the
  peaks, the recorded QC cuts and DoubletFinder singlet call, the merge, the
  ArchR filter, the WNN integration. The HT-scCAT-seq and 10x Multiome matrices
  enter at §3; nothing here re-derives them.
* **The SoupX correction is upstream of this notebook.** The corrected matrix
  (`$FIG2_DATA/FY1x/SoupxMatrix`, 10x format) is read as given, together with the
  peak matrix called from the same bead run.
* **33 types, not 32.** `cell_type` has 32 classes; `RNA_snn_res.1` cluster 34
  (280 cells, annotated 276 Oligo + 4 OPC) was re-annotated as **COP/NFOL** on
  Gpr17/Bmp4 co-expression evidence, giving `cell_type_updated` = 33 classes. The
  annotation is a manual step, so §5 writes the object without it and §6 applies
  the recorded labels.
* **`technique` groups the two dsciCAT-seq libraries** (FY13 → dsciCAT-seq-50k,
  FY14 → dsciCAT-seq-100k), which is why the legend counts are 47,443 / 14,762 /
  8,438.
* **Panel D's rendering parameters were recovered by measurement.** The dot plot
  saved inside the analysis notebook (`dot.scale = 16`, `8 × 8 cm`, bold row
  labels, a plot title) is not what was published: its dots are 12 mm and merge
  into a solid mass. The delivered D panel — cropped from
  `Figure S2-dot-track.pdf` — has 6 pt plain row labels, no title, a 48 × 33 grid
  at 1.86 × 2.27 mm and a largest dot of 1.37 mm; the script that drew it is not
  on this server, so §8 reproduces those measured numbers (`dot.scale = 1.5`,
  panel box 8.90 × 7.48 cm, guides above the panel).
* **Panel E uses the ArchR project the delivered figure was drawn from**
  (`archr-subset-copnfol`, ArchR 1.0.2, mm10). Its 33-cell-type grouping is
  re-derived in §9 from the Seurat object by barcode matching — the project itself
  predates the final annotation.
* **Two final-composition touches are not reproduced here.** The delivered page
  renames the `Ependyma` class to **Ependymal cells** (dot-plot row label and
  legend row 8, a production rename with zero geometry change), and recomposes
  panels A–C with a single shared legend and re-placed centroid numbers. The
  notebook reproduces the panel sources: `Ependyma` and per-panel legends.
* **Fonts.** The delivered page is set in Arial; this machine has no Arial, so the
  R device substitutes Helvetica — metrically identical, so layout is unaffected.

## Run

```bash
export FIG2_INPUT=/path/to/HZJ        # workflow outputs (RNA matrices, peaks, fragments)
export FIG2_DATA=/path/to/FY          # analysis tree (SoupX matrices, per-stage objects)
export FIG2_OUT=fig2_out
jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=-1 Figure_2.ipynb
```

Every plotting cell saves its figure and reports the path and size (nothing is
displayed inline; the notebook stays small on GitHub).

§4 (ArchR) and §5 (WNN over 70,643 cells) are the expensive stages — hours on 16
threads and tens of GB of scratch under `FIG2_DATA`; both are guarded by their
output, so a notebook run on a machine that already has
`3k_dsci_ht_10x_gs_subset.RDS` spends its time on the panels instead. Reading the
13 GB object (§5's guard and §6) costs ~13 GB of RAM and about 5 minutes; §9
tracks cost ~10 minutes and several GB more. Keep
`--ExecutePreprocessor.timeout=-1` so the kernel does not kill them.

## Requirements

**Python** (the kernel is Python 3; it only launches R and shows the results):
`numpy`.

**R ≥ 4.1**, with

* §1, §6–§9: `Seurat` (≥ 4), `Signac`, `ggplot2`, `dplyr`, `tidyr`, `pals`,
  `readr`, `rtracklayer`, `GenomeInfoDb`
* §2: `DoubletFinder` (with its `RANN`/`Matrix` back end), `glue`
* §4, §9: `ArchR` (1.0.x) + `BSgenome.Mmusculus.UCSC.mm10`
* §5: `harmony`
* `future` / `future.apply` for the parallel steps

Executed with R 4.2.2, Seurat 4.3.0, Signac 1.9.0, ArchR 1.0.2, harmony 1.2.3,
pals 1.8, and Python 3.10.

**Reference data**: the mm10 gene annotation ArchR needs, and for the rebuild
branch of §1 a mm10 `genes.gtf` (or the pre-extracted
`<work>/mm10.gene.coords.RDS`).
