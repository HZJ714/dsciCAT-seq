# Figure 3 — One identical parameter set applied to three libraries resolves the Gpr17⁺ COP/NFOL state

`Figure_3.ipynb` reproduces the measured content of Figure 3 from the files the
analysis started from. Journal assembly (page geometry, panel letters, the
hand-drawn COP/NFOL ellipse, the inset leader lines, compositing) is a production
step and is deliberately **not** part of the notebook.

| Panel | What it shows | Section | Output |
|---|---|---|---|
| A–C | WNN UMAPs of HT-scCAT-seq (4,029), dsciCAT-seq-50k (17,093) and dsciCAT-seq-100k (30,350), one shared palette, COP/NFOL ringed in B and C | §1, §2 | `fig3abc_wnn_umaps.pdf` (+ `.png`), `fig3_copnfol_hull.csv` |
| B/C insets | *Gpr17* on the re-clustered COP/NFOL cells (60 / 133), plus the 76 EL14 OPC cells as a control | §1, §2.1 | `fig3bc_gpr17_insets.pdf` (+ `.png`), `fig3_insets_summary.txt` |
| D | UMAP of the 5,977-cell oligodendrocyte-lineage subset in its six states | §3 | `fig3d_lineage_umap.pdf` (+ `.png`), `fig3d_umap_states.csv` |
| E | 17 lineage markers × six states: size = % expressing, colour = average expression | §4 | `fig3e_marker_dotplot.pdf` (+ `.png`), `fig3e_marker_stats.csv` |
| F | ATAC coverage of the five loci of E, one row per state | §5 | `fig3f_archr_tracks.pdf` (5 pages), `fig3f_archr_track_1..5.png` |
| — | validation of E and F against the vector files the delivered panel was cropped from | §6 | printed, no file |

## The step A–C and the insets depend on (§1)

### §1.0 — the three per-library objects are subsets of the joint object

Figure 2 ends in `3k_dsci_ht_10x_gs_subset.RDS` (70,643 nuclei, four libraries).
The three objects §1 reads are one `subset()` each out of it:

```r
el14 <- subset(ht, cells = Cells(ht)[which(ht$batch      == 'Cortex_EL14')])     #  4,029
fy13 <- subset(ht, cells = Cells(ht)[which(ht$orig.ident == 'dsciCAT-seq-50k')])  # 17,093
fy14 <- subset(ht, cells = Cells(ht)[which(ht$orig.ident == 'dsciCAT-seq-100k')]) # 30,350
```

`orig.ident` names the two dsciCAT-seq libraries but lumps the four HT-scCAT-seq
libraries together (`HT-scCAT-seq-HT`, 14,762), so the single-library object has
to be picked by `batch`, which still holds `Cortex_EL13 … EL16`. §1.0 runs those
three lines and checks them against the files actually used, writing
`fig3_subset_derivation.csv` and `fig3_subset_meta_equal.csv`. Measured here:
same cell sets, same order, the RNA count matrix bit-identical over 400 genes ×
every cell, the 31,263 genes identical, and all nine RNA QC columns at 1.0000
agreement. The ATAC side is identical on every peak checked, but the two
dsciCAT-seq objects carry 186,354 peaks against the joint object's 186,321, so
their `nCount_ATAC`/`nFeature_ATAC` sit higher for 58–64 % of the cells (≤ 2.1 %
per cell, 0.09 % overall) — the joint object is quantified on the shared peak
set, and no panel here uses those two columns. What a subset cannot carry is the
**library-level analysis**: each object's own PCA/UMAP/`umap.atac`/WNN
reductions, neighbour graphs, `RNA_snn_res.1` (agreement 0.16 / 0.36 / 0.52) and
the hand-written `cell_type` table (0.76 / 0.75 / 0.80). That is why §1 still
reads the three files. `FIG3_SAVE_DERIVED=1` also writes the plain subsets as
`fig3_out/{EL14,FY13,FY14}_from_joint.RDS`.

### §1 — the cache

Three per-library Seurat objects and the joint multiomic object are reduced to
one 9 MB cache, `fig3_light.RDS`: the cell-type call used by the figure
(`cell_type`, except that the COP/NFOL cluster is taken by name from
`RNA_snn_res.1` — cluster 29 in FY13, cluster 31 in FY14), the RNA / ATAC / WNN
embeddings of every cell, the raw counts of the cells the insets re-cluster, and
the shared palette. The 13 GB joint object is read **once**, and the three
libraries are matched to it by barcode — 4,029 / 17,093 / 30,350 of 4,029 /
17,093 / 30,350, i.e. 100 %, so A–C and Figure S3 A/B are literally the same
cells. §2 and §2.1 then never touch anything but the cache.

```
              3k_dsci_ht_10x_gs_subset.RDS     Figure 2's joint object, 70,643 nuclei
                        │
                        │  §1.0  subset() per library (rules above) + the library-level
                        │        run those files carry: pca / umap / umap.atac /
                        │        wnn.umap, graphs, RNA_snn_res.1, hand-written cell_type
        ┌───────────────┼───────────────┐
3k_FY13_subset.RDS 3k_FY14_subset.RDS EL14_subset.RDS   + the joint object's own
      17,093            30,350            4,029          umap.atac and wnn.umap
        └───────────────┬───────────────┘
                        §1  fig3_light.RDS
                                │
              ┌─────────────────┴──────────────────┐
              §2  panels A–C (WNN UMAP + outline)   §2.1  Gpr17 insets (subset re-cluster)
```

## Inputs

| Variable | Default | Contents |
|---|---|---|
| `FIG3_DATA` | `<work>/FY` | `3k_dsci_ht_10x_gs_subset.RDS` (13 GB) — the joint object — and the three subsets §1.0 cuts out of it, `FY13/3k_FY13_subset.RDS` (1.7 GB), `FY14/3k_FY14_subset.RDS` (2.6 GB), `HT/EL14_subset.RDS` (0.6 GB), which add the library-level analysis |
| `FIG3_OLIGO` | `<work>/FY/oligo` | `3k_dsci_oligo_subset_{umap,celltype}.csv`, `3k_dsci_oligo_subset_rna.h5ad`, and `3kdsci/archr-subset/` (ArchR project of the same cells) |
| `FIG3_OUT` | `fig3_out` | every panel PDF/PNG/CSV, the R sources as run, the derived cache |
| `FIG3_FORCE` | unset | `1` re-runs §1 even when the cache exists |
| `FIG3_SAVE_DERIVED` | unset | `1` makes §1 also write the plain subsets (≈ 9 GB) as `fig3_out/*_from_joint.RDS` |

## Provenance: read this before comparing with the printed figure

* **The three per-library objects are not separate data** — they are `subset()`s
  of Figure 2's joint object (§1.0), verified cell for cell on barcodes, order,
  RNA counts and RNA QC metadata; only the library-level analysis is theirs, and
  only that part cannot be recovered from the joint object.
* **A–C are drawn from the joint WNN embedding**, restricted to each library,
  with each library's own cell-type call — the point of the figure is that the
  three libraries were clustered and annotated with **one identical parameter
  set**, and that COP/NFOL is resolved in both dsciCAT-seq libraries (60 and 133
  cells) and never in the 4,029-nucleus HT-scCAT-seq library (0 cells).
* **The COP/NFOL marking** rings each cell of that cluster in red in B and C and
  writes the cluster's convex hull (with its cell barcodes) to
  `fig3_copnfol_hull.csv`; the delivered page uses a hand-drawn ellipse around
  the same cells.
* **The insets** re-cluster the COP/NFOL cells from their stored counts with the
  project's default parameters and one common *Gpr17* colour limit. The delivered
  figure shows the two dsciCAT insets; the notebook adds the 76 EL14 OPC cells as
  a control, and reports what it measures: *Gpr17* means 0.67 (EL14 OPC), 0.71
  (FY13), 0.93 (FY14) — the control is not *Gpr17*-negative, the difference is
  that EL14's clustering never resolves the state in the first place.
* **D–F are the same 5,977 cells** and the same six states in the same order, so
  E and F line up row by row. That lineage subset is read as published (its
  coordinates, labels and matrices share the same barcodes; §0.3 checks it); it
  is not re-derived here.
* **E, size channel is exact.** §6 measures the 102 dots of the vector file the
  delivered panel was cropped from: `radius(cm) = 0.00170 × pct + 0.0060` against
  the `pct` column of `fig3e_marker_stats.csv`, Pearson **r = 1.0000**.
* **E, colour channel is not identified.** The printed colour is a per-gene
  rescaling of an average-expression table, but not a monotone function of the one
  recomputed here (*Mbp* sits at the delivered colour minimum in MOL2 although its
  mean log-normalised RNA is the second highest of the six states). The notebook
  plots the unscaled mean and writes `avg`, `pct` and the per-gene `z` side by
  side; measured against the printed dots: Pearson 0.415 vs `avg`, 0.857 vs the
  per-gene `z`.
* **F re-runs the original computation.** The delivered F is a crop of
  `3kdsci/archr-subset/Plots/oligo_subset_celltype_track.pdf`, the output of
  `ArchR::plotBrowserTrack` on that project (ArchR 1.0.2, mm10, `groupBy =
  'cell_type'`, TSS ± 50 kb, `normMethod = 'ReadsInTSS'`). §5 calls the same
  function on the same project; §6 compares the region strings — 5 of 5 pages
  match (`chr14:46340668−46440669`, `chr14:13234779−13334780`,
  `chr9:22339205−22439206`, `chr12:88672875−88772876`, `chr5:9675351−9775352`).
  The coverage *range* printed in each page title moves by a few percent between
  runs (0–2.09 here vs 0–2.07 in the source); the tracks themselves are identical.
* **Fonts.** The delivered page is set in Arial; this machine has no Arial, so the
  devices substitute the metric-compatible Liberation Sans / Helvetica — layout
  unaffected.

## Run

```bash
export FIG3_DATA=/path/to/FY          # four Seurat objects
export FIG3_OLIGO=/path/to/FY/oligo   # lineage subset: 2 CSVs + RNA h5ad + ArchR project
export FIG3_OUT=fig3_out
PYTHONPATH=<work>/pylibs jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=-1 Figure_3.ipynb
```

`PYTHONPATH` only feeds §6 the PyMuPDF build; without it §1–§5 run normally and
§6 prints a skip message. §1 is guarded by its output: delete
`fig3_out/fig3_light.RDS` (or set `FIG3_FORCE=1`) to redo the 13 GB read.

Measured on this machine (376 GB RAM, 96 cores), 2026-10-22: §1 ≈ 10 min (the §1.0
derivation checks included; three subsets 2.6 min, the joint object 7.2 min),
§2 0.3 min, §2.1 1.0 min, §3–§4 < 1
min, §5 1 min, §6 < 1 min — **≈ 11 min** with a warm cache ≈ 4 min.

## Requirements

**Python ≥ 3.8** (the kernel; it launches R and shows the results): `numpy`,
`pandas`, `scipy`, `anndata`, `matplotlib`; §6 additionally `Pillow` + PyMuPDF.

## Orphan intermediate (not tracked)

`fig3_out/joint_meta.csv.gz` (24.5 MB) is a joint WNN metadata export from the authors'
analysis tree. **No code or notebook in this repository reads it** (verified 2026-11-11),
and it is excluded from version control (`.gitignore`). It is kept on the authors' disk
as an intermediate of the Fig 3 data-preparation step; readers who need it should export
it from the joint object (`3k_dsci_ht_10x_gs_subset.RDS`, §1.0) as part of their own
preparation, or delete it — nothing in the figure chain depends on it.

**R ≥ 4.1** with `Seurat` (4.x), `Signac` (loaded with Seurat), `ggplot2`,
`dplyr`, `tidyr`, `ggrepel`, `pals`, `scales`, `patchwork`, `grid` (§1, §2, §2.1)
and `ArchR` 1.0.x + `BSgenome.Mmusculus.UCSC.mm10` (§5).

Executed with R 4.2.2, Seurat 4.3.0, ggplot2 3.5.2, ggrepel 0.9.5, ArchR 1.0.2,
Python 3.8.14, numpy 1.22.4, pandas 1.5.1, scipy 1.10.1, anndata 0.9.2,
matplotlib 3.6.2.
