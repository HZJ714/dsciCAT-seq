# 00_Preprocessing — FASTQ → count matrices → cell-level RNA preprocessing

All data in this study start from raw dual-modality libraries (scRNA-seq
cDNA + scATAC-seq). This directory contains the complete upstream stage:

1. **`scCAT_workflow/`** — four Cromwell WDL workflows that take FASTQ and
   produce count matrices (RNA gene × cell, ATAC fragment/peak × cell,
   spliced/unspliced matrices).
2. **`RNA_soupx_singlet.ipynb`** — the matrix-level RNA preprocessing for
   the mouse-brain data (SoupX ambient-RNA correction, QC filtering,
   DoubletFinder doublet removal, PCA/UMAP/clustering), applied on the
   matrices produced by `scCAT_workflow`.

```
raw FASTQ ──► scCAT_workflow (WDL) ──► count matrices ──► RNA_soupx_singlet.ipynb ──►
    (4 workflows, one per platform/data set)      (SoupX, QC, doublets)        Seurat object
```

---

## 1. scCAT_workflow — FASTQ → matrices (WDL)

### 1.1 Workflows and data sets

| WDL | Data set | Modalities | Input FASTQ per library | Barcode layout |
|---|---|---|---|---|
| `dsciCAT-seq-V1-matrix_1.0.0.wdl` | **Mouse brain** (3k_dsci_oligo; libraries FY13, FY14) | RNA + ATAC (dsciCAT-seq) | 2 files (paired-end) | 30 bp bead barcode (R1), `bc:0:29,r1:30:-1,r2:0:49` |
| `ISSAAC-full-10xformat_1.0.0.wdl` | **10x Multiome and ISSAAC-seq** samples | RNA + ATAC | **3 files** `[barcode, R1, R2]` | I2 barcode read (16 bp), UMI 10 bp |
| `ISSAAC-ATAC_1.0.0.wdl` | ATAC-only runs of the same format | ATAC | 3 files (`tenx=true`) or 2 files (`tenx=false`) | I2 read, or 30 bp R1 (dsci format) |
| `ISSAAC-RNA-BGI_1.0.0.wdl` | **Own BGI-platform data (HG)** | RNA | 2 files (paired-end) | 30 bp barcode at the start of R2 |

Two analyses (e.g. parts of Fig 1) are matrix-only; their scripts live
outside this repository and consume the matrices listed in Section 3.

### 1.2 dsciCAT-seq (mouse brain)

Workflow `iDrop_dsciCAT_1`, single species (mm10). Both mouse-brain libraries
(FY13, FY14) were run with `ForceFrag ≈ 3000` (bead-calling fragment cut-off
in `d2c merge`).

```
Data_RNA (paired FASTQ)        Data_ATAC (paired FASTQ)
        │                             │
        ▼                             ▼
cDNAanno (scStar_dsci align,    mapping (chromap --preset atac,
  Anno_dsci --intron,             bc:0:29,r1:30:-1,r2:0:49,
  PISA count)                     trim adapters, whitelist_V1_30bp)
        │                             ▼
  PISA count ──► RawMatrix       deconvolution (d2c_v1.5.3 merge,
                                      --mapq 30 --bf ≈3000 --tn5 10
        │                              --rank inflection)
  getM280UMI (bead ↔ ATAC             │
  barcode merge, tagAdd DB)           ▼
        │                        qc (TSS/fragment size) +
  count (PISA count -cb DB)      peakcount (macs2 -q 0.001 --nomodel,
        │                            C4scATAC_Cluster_AnnotationV3)
        ▼                            ▼
  FilterMatrix                 fragments.tsv.gz (bgzip+tabix)
  SpliceMatrix (--intron)      Peak / Promoter matrices
  RNAVelocityMatrix            cell & library reports
  saturation report
```

Outputs: `04.Matrix/RawMatrix`, `04.Matrix/FilterMatrix`,
`04.Matrix/SpliceMatrix`, `04.Matrix/RNAVelocityMatrix`,
`ATAC/03.d2cfile/*.fragments.tsv.gz`, `ATAC/out/Peak|Promoter`,
`05.MergeBead/*_barcodeTranslate_16.txt`.

### 1.3 ISSAAC-full-10xformat (10x Multiome / ISSAAC-seq)

Workflow `iDrop_ISSAAC_full_10x_V1`, mixed species (GRCh38-mm10,
`mixspecies = true`). Three FASTQ files per library: `[barcode, R1, R2]`.

```
[bc, R1, R2] RNA                   [bc, R1, R2] ATAC
      │                                    │
      ▼                                    ▼
prepare_cb_umi (paste bc+R2:        chromap (--preset atac,
  CB 16 bp, UMI 10 bp)                  barcode = I2 read)
      │                                    │
      ▼                                    ▼
rna_common (STAR solo:               atac_common:
  CB_UMI_Simple, EmptyDrops_CR,         convert_frag_to_reads.py + bedClip
  clip3pNbases 115)                     macs2 (--nomodel --extsize 200
      │                                    --keep-dup all --SPMR)
      ▼                                    bedtools intersect +
  Solo.out/Gene/filtered                  generate_csc_mtx.py
  (gene × cell mtx)                       STAR soloCellFiltering
                                               │
                                               ▼
                                          raw_mtx / filtered_mtx
                                          (peak × cell mtx)
  ── mixspecies: GRCh38/mm10 fragment-size & per-cell count tables ──
```

### 1.4 ISSAAC-ATAC (ATAC-only)

Workflow `iDrop_ISSAAC_ATAC_V1`: the ATAC half of 1.3 as a standalone
pipeline. `tenx = true` (3-file input) by default; `tenx = false` switches
to the dsci 2-file format (`read_format bc:0:29,r1:30:-1,r2:0:49`). Optional
down-sampling with `seqtk sample -s 42` (`down`, `down_rate`). The
alignment → peak → matrix chain is identical to 1.3.

### 1.5 ISSAAC-RNA-BGI (own data, BGI platform)

Workflow `iDrop_ISSAAC_RNA_V1`, mixed species (GRCh38-mm10). Two FASTQ files
per library; the 30 bp cell barcode is the first 30 bp of R2.

```
[FQ1, FQ2]  ──►  rna_common (STAR solo:
                   readFilesIn = FQ2 FQ1 (barcode first),
                   --soloCBstart 1 --soloCBlen 30 --soloUMIstart 31
                   --soloUMIlen 12 --soloBarcodeReadLength 0
                   --soloCellFilter EmptyDrops_CR)
                        │
                        ▼
              Solo.out/Gene/filtered (gene × cell mtx)
              + per-species (GRCh38/mm10) per-cell read counts
```

### 1.6 Running the workflows

Required: a Cromwell-compatible WDL runner and the pipeline Docker images.
Each task declares its image via `docker_url`; the reference genomes
(`refRNA`, `refATAC`), barcode whitelist and per-library FASTQ file lists
are workflow inputs.

**Internal dependencies to configure before re-running elsewhere** (they do
not affect the published outputs, which were produced with the versions
below):

| Dependency | Where | Note |
|---|---|---|
| Docker images `<registry-namespace>/*:latest` | all tasks | Private HPC registry; substitute equivalent tool images (chromap, d2c, STAR, PISA, macs2, …) |
| `<hpc-mount>/.../sta_collision_rate.R` | dsciCAT `deconvolution` task | Internal HPC script (bead-collision rate); provide your own copy |
| `/dsciRNA`, `/C4_scATAC_analysis` | dsciCAT tasks | Paths inside the image (whitelist, adapter, helper scripts) |
| `ProjectID = "<project-id>"` | dsciCAT input | Project tag of the original run |
| `ForceFrag` | dsciCAT input | `≈3000` for FY13/FY14 (mouse brain); tune per library |

Reference data: mm10 (mouse brain) and GRCh38-mm10 combined (mixed-species
samples) STAR/chromap indexes, GENCODE GTF, TSS/promoter BED, `chrom.sizes`,
and the 30 bp bead whitelist `whitelist_V1_30bp.txt`.

**Upload hygiene (2026-11-10, extended 2026-11-11):** the four WDLs here are provenance
copies of the versions that ran. Two kinds of private string are masked:

| in the pipeline as run | in this repository |
|---|---|
| `<user>_<32-hex-env-id>_private:latest` (per-user private image tag) | `<image-tag-masked-for-upload>:latest` |
| the registry namespace derived from the analysis-machine account | `<registry-namespace>` |
| `<hpc-mount>/<project-id>/<user>/shell/…` (internal helper script) | `<project-id>/<user>` placeholders |

Nothing else in the four files was changed. To actually run them you must restore the
three strings above from the pipeline's own repository (the `docker_url` lines are the
only place they appear). Site-specific path prefixes are written as the placeholder
`<hpc-mount>/` throughout: substitute your own institution's absolute prefix, and restore the
project tag where `"<project-id>"` appears.

---

## 2. Matrix-level RNA preprocessing — SoupX (mouse brain)

`RNA_soupx_singlet.ipynb` (R, Seurat/SoupX/DoubletFinder) is the post-matrix
RNA stage applied to the mouse-brain libraries. It consumes the dsciCAT-seq
outputs directly:

| Input | From | Use |
|---|---|---|
| `04.Matrix/RawMatrix` | dsciCAT `cDNAanno` | Non-cell bead counts → soup profile |
| `04.Matrix/FilterMatrix` | dsciCAT `count` | True-cell counts (`toc`) |
| `05.MergeBead/*_barcodeTranslate_16.txt` | dsciCAT `getM280UMI` | Valid cell barcodes |

Pipeline (per library, `enable_soupx = TRUE`):

1. **SoupX** — the raw matrix minus valid cell barcodes is the ambient-RNA
   (`tod`) background: `SoupChannel(tod, toc)`, clusters from a quick
   Seurat pass (`setClusters`), fixed contamination 0.2
   (`setContaminationFraction`), `adjustCounts()`; output written as a
   10x-v3 matrix (`DropletUtils::write10xCounts`).
2. **Seurat object + QC** — `percent.mt` from `^mt-`; filters
   `percent.mt < 20`, `nCount_RNA < 20000`, `200 < nFeature_RNA < 6000`;
   `NormalizeData` + VST (3000 HVGs) + PCA.
3. **DoubletFinder (10 %)** — `paramSweep_v3`/`find.pK` to pick `pK`,
   `doubletFinder_v3(pK, nExp = 0.1·nCells)`; doublets removed
   (`doublet_info == "Singlet"`).
4. **Clustering checks** — PCA elbow, UMAP, `FindClusters(resolution = 1)`;
   QC plots (`*_QC_RNA_soupx.pdf`, `*_subset_QC_RNA_soupx.pdf`,
   `*_cluster_RNA_soupx.pdf`) and the processed object
   (`*_RNA_soupx_PCA.RDS`) are written per library.

Downstream: the SoupX-corrected objects feed the integration step that
produces the shared Seurat object used by the figure notebooks
(`3k_dsci_oligo_subset.RDS` → `.h5ad`; see `UPSTREAM.md` and the
`fig4/fig6/figS4/figS6` notebooks in this repository).

---

## 3. Output → downstream analysis objects

| Pipeline output | Assembled into | Consumed by |
|---|---|---|
| dsciCAT matrices (mouse brain) → SoupX + doublet filtering | `3k_dsci_oligo_subset.RDS` → `.h5ad` | Fig 4, Fig 6, Fig S6 notebooks |
| Same matrices, via MultiVelo preprocessing (`scripts/Multivelo/00/01/04`) | `multivelo_dynamical.h5ad` | Fig S4 notebook |
| 10x / ISSAAC filtered matrices | separate analyses | Fig 1 (matrix-only scripts outside this repo) |
| BGI (HG) gene × cell matrices | separate analyses | Fig 1 (matrix-only scripts outside this repo) |

## 4. Data availability

Raw sequencing data (FASTQ) and raw count matrices for the authors' own
(BGI-platform) dataset have been deposited in the China National GeneBank
DataBase (CNGBdb, https://db.cngb.org) under accession **CNP0010263**. The
mouse-brain matrices were generated by the `dsciCAT-seq` workflow above
from the FY13/FY14 libraries.
