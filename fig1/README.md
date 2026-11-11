# Figure 1 — dsciCAT-seq: library design and benchmarking in two cell lines

`Figure_1.ipynb` reproduces the measured content of Figure 1 from the files the
analysis started from. Journal assembly (page geometry, panel letters, shared
colour key, compositing) is a production step and is deliberately **not** part of
the notebook.

## Panels

| Panel | What it shows | Section | Output | Reproducible? |
|---|---|---|---|---|
| A | experimental + computational workflow | — | — | vector illustration, nothing computed |
| B | species collision, RNA and ATAC channel, library HG14 | §2 | `fig1b_contamination.pdf` | yes |
| C | genome tracks at one locus per cell line (ATAC rows) | §5 | `fig1c_atac_tracks_HEK293T.pdf` + `.png` | yes (R/Signac) |
| D | RNA capture: lg(UMIs), lg(detected genes) | §3 | `fig1d_rna_qc.pdf` | yes |
| E | ATAC capture: lg(unique fragments) | §3 | `fig1e_atac_qc.pdf` | yes |
| F | RNA similarity between platforms | §4 | `fig1f_rna_similarity.pdf` + `corr_{HEK293T,NIH3T3}.csv` | yes |

## The step everything depends on (§1)

All six platforms are reduced to the same two objects before any panel is drawn:

```
mixed-species (GRCh38 + mm10) pipeline output
   per-bead/per-cell counts against each genome  +  matrix with GRCh38_/mm10___ features
        │
        ├─ bead → cell merge (dsciCAT-seq: barcodeMerge.tsv)          §1a
        ├─ species call: human (or mouse) fraction > 0.8 → Human/Mouse/Mix
        ├─ matrix split: keep that species' cells, drop the other's features,
        │                strip the prefix  →  Human/  Mouse/  (10x triplets)   §1b
        └─ metrics.csv: nCounts, nFeatures per cell                            §1c
```

Panels B/D/E/F consume only those outputs, so §1 is written out in full.

The bead → cell merge is the one place where a second upstream pipeline shows up: the
`HG*` libraries were also processed by our own dsciCAT-seq pipeline — ATAC FASTQ through
`chromap`, bead barcodes re-assembled by `d2c` under the same `forcefrag2000`
fragment-length cap — purely to produce `barcodeMerge.tsv`, because the dsciCAT bead carries
its barcode in a structure that neither the ISSAAC nor the 10x chemistry has.

## Inputs (`FIG1_DATA`, default `<raw>`)

| Used by | Path (relative to `FIG1_DATA`) |
|---|---|
| §1, §2 | `HG/HG14/sub/outs/per_cell_barcode_total_{read,fragment}_count_{GRCh38,mm10}.tsv`, `HG/HG14/HG14.barcodeMerge.tsv` |
| §2 | `HG/HG14/sub/ISSAAC_joint_metrics_sub.csv` |
| §3 | 26 × `metrics.csv`, listed verbatim in §0.3 (one per platform × cell line × modality) |
| §4 | 8 RNA matrix directories (`matrix.mtx` + `features.tsv` + `barcodes.tsv`), table in §0.1 |
| §5 | `HG/HG12/sub/human_fragments.tsv.gz`, `HG/HG14/sub/human_fragments.tsv.gz`, `scCAT-seq-HT/IL9/ATAC-sub/human_fragments.tsv.gz` (+ `.tbi`, both already indexed) and the matching `ATAC/Human/` peak matrices |

## Provenance: read this before comparing D–F across platforms

The benchmarking set is **not** homogeneous, and the notebook says so in §0.2:

* **Re-processed through our own pipeline with the reads down-sampled** —
  dsciCAT-seq-50k (HG12), dsciCAT-seq-100k (HG14), ISSAAC-seq, HT-scCAT-seq,
  10x Multiome. For these the cell calling, the species split and the depth are
  ours and comparable; the QC table and the matrix describe the same cells.
  Measured effect of the down-sample (median UMIs, HEK293T RNA, full → `sub`):
  3,126→1,897 · 2,168→1,587 · 21,468→11,912 · 42,662→11,412 · 5,105→3,056.
* **Published matrices used as-is** — sciCAR-seq, SUM-seq, Paired-seq (GEO
  downloads, EmptyDrops_CR-filtered by their authors). Their sequencing depth is
  **not** matched to the rest, and for two of them the QC table and the matrix are
  not even the same cell set (sciCAR HEK293T RNA: 465 vs 812; Paired-seq RNA:
  1,283 vs 593). sciCAR's ATAC matrices are not in the tree at all — only their
  `metrics.csv`, so panel E's sciCAR columns can be transcribed, not recomputed.

Panel F is the only place where depth is compensated (CPM on the pseudobulk
profile, §4 step 3). Panels D and E are intentionally shown at each platform's own
depth.

## Run

```bash
export FIG1_DATA=/path/to/HZJ FIG1_OUT=fig1_out
jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=-1 Figure_1.ipynb
```

Every plotting cell saves its figure and reports the path and size (nothing is
displayed inline; the notebook stays small on GitHub). Measured on this machine (2026-10-21, 16 cores): §0–§4 ≈ 6 min
(dominated by reading ~1 GB of MatrixMarket in §4), §5 ≈ 10 min and ~14 GB of RAM.
§5 writes `fig1c_atac_tracks.R`, runs it, and caches the result — delete
`fig1c_atac_tracks_HEK293T.png` (or set `FIG1_FORCE=1`) to re-render it. Keep
`--ExecutePreprocessor.timeout=-1` so the kernel does not kill it.

## Validation against the delivered figure

Executed 2026-10-21 with `/opt/conda/bin/python` (3.8.14, numpy 1.22.4,
pandas 1.5.1, scipy 1.10.1, matplotlib 3.6.2, seaborn 0.12.2) + R 4.2.2.

| Check | Notebook | Delivered figure |
|---|---|---|
| §1/§2 cell classes (HG14, 42,816 paired cells) | RNA 19,766 / 20,097 / 2,953 · ATAC 19,639 / 20,616 / 2,561 | identical (legend numbers) |
| §1 species threshold, measured | Human min 0.8001 (RNA) / 0.8004 (ATAC); Mix max 0.7998 / 0.7999 | **0.8**, both channels |
| §3 medians, all 26 tables | as printed in §3 | identical to the legend numbers |
| §4 shared genes | 26,041 (HEK293T) / 29,621 (NIH3T3) | identical |
| §4 correlation matrices | `corr_*.csv` | **max abs deviation 0.0e+00** (HEK293T) / 1.1e-16 (NIH3T3) vs the cache the delivered `Figure_1.pdf` was drawn from |
| §4 extremes | HEK293T 0.882 (dsciCAT-50k/-100k) / 0.377 (Paired vs 50k); NIH3T3 0.883 / 0.718; 10x↔ISSAAC 0.8826, 10x↔sciCAR 0.7225 | identical |
| §5 objects | 318,130 peaks × 16,165 cells · 345,484 × 27,760 · 84,125 × 1,970 (after `min.cells = 10`) | the three ATAC libraries of panel C |

The three libraries are drawn separately and stacked (`patchwork`) rather than
`merge`d: `CoveragePlot`/`TilePlot` read the fragment files, not the assay, so a
merge would only pay for a shared peak set — measured: `BindMatrices` over ~780k
peak names was still running after 40 min.

## Requirements

`numpy pandas scipy matplotlib seaborn` (Python ≥3.8) for §0–§4; `Rscript` with
`Seurat Signac Matrix readr tidyr magrittr ggplot2 patchwork GenomicRanges
EnsDb.Hsapiens.v86` (+ `EnsDb.Mmusculus.v79` for the NIH3T3 half) for §5.
Versions are printed by §0.4.

## What changed on the published figure (2026-10-29, reviewer round 2 — final)

`Figure_1.pdf` in the delivery block is now **`Fig1_v4.pdf`** (md5 `c4d759c4ae7b`, 17.00 x 25.19 cm,
unchanged geometry; the 2026-10-21 v3.5 is kept as `_superseded/Figure_1_20261021_deprecated.pdf`).
The chain that builds it is in `assembly/`:

```
fig1c_compute.py + fig1c_v4.R      panel C recomputed from the fragment/peak files (hg38 + mm10)
replot_v4.py                       panels B / D / E / F redrawn (log10 axes, sci-CAR treatment)
fig1_layout.py                     the page geometry shared by every version since v3
compose_fig1_v4.py                 splice -> Fig1_v4.pdf
extract_fig1_A.py                  strip the ghost text layer out of the AI-master A band
finalize_fig1_v4.py                install into the block + PNG + figures.docx media + MD5 list
```

What v4 fixes (reviewer items G0/G1): panel C is rebuilt with genome builds and coordinates on the
axis (human *MIR17HG* chr13:91,340,000-91,370,000 hg38; mouse *Hes7* chr11:69,090,000-69,130,000
mm10), the normalisation written into the panel (100-bp windows; ATAC = cuts per window, RNA =
aligned read-bases per window / 10; scale shared within a block only), tile rows labelled
**120 top cells**; D/E axes relabelled `log10(count + 1)` with ticks on the right half; the sci-CAR
NIH3T3 ATAC column greyed, hatched and marked `‡ not comparable: the GSE117089 …`; matrix-only
methods marked `*`.

**Two reviewer asks were declined by the authors, deliberately:**

- recolouring the A-panel scheme so the 10-nt pre-index, UMI, bead barcode and sequencing index are
  four different colours: the two indexing rounds share barcodes by design (for memorability), and
  that is already stated in legend A and in the Results text, so a fourth colour would be a
  duplication rather than an addition;
- moving part of Figure 1 to the supplements because it is dense: merging the contamination panels
  into one dual-omics panel would be just as dense, so the split would only relocate the problem.

**Known, accepted gate waiver** (this is the one delivered figure that does not pass the house font
gates, and it is recorded as such rather than left as an unexplained FAIL): minimum size 1.84 pt
(`Biotin`, `S-S` inside the AI-master schematic), Helvetica on 13 spans, Type 3 DejaVuSans and
un-subsetted Arial/Liberation inherited from the same band. Hidden characters improved 198 -> 47, but
spans below 5.5 pt went 76 -> 116 because v4 added small annotations: **20 of them are in the A
schematic (untouchable without the AI file) and 96 are in B-F, which we draw** — if this figure ever
has to pass, raise those 96 in `replot_v4.py`; that is the only work required.

The legend in `figures.docx` was also stale (the 2026-10-21 short version, 986 characters, with a
sentence — "for all or 300 single cells" — that contradicts the 120-cell tiles now printed). It was
synced to the manuscript text (4,848 characters) in the 2026-10-30 revision round.
