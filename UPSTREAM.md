# Upstream processing — FASTQ → matrices → cell-level preprocessing

The complete upstream stage of this study is documented in
**[`00_Preprocessing/README.md`](00_Preprocessing/README.md)**, which contains:

- **`00_Preprocessing/dsciCAT_workflow/`** — the four Cromwell WDL workflows
  that convert raw FASTQ into count matrices (RNA gene × cell, ATAC
  fragment/peak × cell, spliced/unspliced), one per platform / data set:

  | WDL | Data set |
  |---|---|
  | `dsciCAT-seq-V1-matrix_1.0.0.wdl` | Mouse brain (3k_dsci_oligo; FY13 + FY14, ForceFrag ≈ 3000) |
  | `ISSAAC-full-10xformat_1.0.0.wdl` | 10x Multiome and ISSAAC-seq (3 FASTQ files per library) |
  | `ISSAAC-ATAC_1.0.0.wdl` | ATAC-only variant |
  | `ISSAAC-RNA-BGI_1.0.0.wdl` | Own BGI-platform data (HG; 30 bp barcode in R2) |

- **`00_Preprocessing/RNA_soupx_singlet.ipynb`** — matrix-level RNA
  preprocessing for the mouse-brain data (SoupX ambient-RNA correction,
  QC filtering, DoubletFinder doublet removal).

- Internal runtime dependencies (private Docker images, internal HPC script
  paths, run parameters) that must be configured before re-running outside
  the original HPC environment are listed in the README.

- **Data availability**: raw FASTQ and raw count matrices of the authors'
  own (BGI-platform) dataset are deposited in CNGBdb under accession
  **CNP0010263** (https://db.cngb.org).
