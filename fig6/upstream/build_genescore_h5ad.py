#!/usr/bin/env python3
"""Rebuild `genescore_oligo_5977.h5ad` — the chromatin (gene-score) channel of the own dsciCAT-seq
nuclei, i.e. the matrix that panel C's chromatin side is scored from — from the artefacts written by
`export_genescore.R`.

Inputs (all produced by export_genescore.R in the same directory):
  genescore_5977.mtx.gz    genes x nuclei sparse matrix, the Signac `gene_score` assay's `data` slot
                           = log1p(ArchR gene score normalised to 1e4 per nucleus)
  genescore_features.txt   24,333 gene names, matrix row order
  genescore_barcodes.txt   5,977 nucleus ids, matrix column order (= lineage_cells_5977.txt)
  genescore_obs.tsv.gz     52-column metadata table of those nuclei, indexed by barcode
                           (its own `cell_type` column is the 11-class *broad* label of the multiome
                           object; it is NOT what the figure uses)
  <rna h5ad>               3k_dsci_oligo_subset.h5ad — the six oligodendroglial subtype labels
                           (obs.cell_type, Fig 2 lineage curation) are joined from here by barcode

Output:
  genescore_oligo_5977.h5ad   nuclei x genes, csr, float64, obs = the 17 columns listed below

The value scale matters: scDRS is called with --flag-raw-count False, and it internally applies
expm1() to the matrix, i.e. it requires a log1p input. This file is log1p (max 6.47); the
cross-chemistry gene-score matrices in the panel-D chain are NOT log1p (max 1,388) — see README.md
§"the two own-data chromatin matrices are not the same transform".

Usage:  python build_genescore_h5ad.py <probe_dir> [out.h5ad] [rna.h5ad]
"""
import os
import sys

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import anndata as ad

OBS_COLS = ["cell_type", "nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC",
            "nCount_gene_score", "nFeature_gene_score", "TSS.enrichment", "TSS.percentile",
            "FRiP", "nucleosome_signal", "atac_celltype", "rna_celltype",
            "wsnn_res.0.5", "ATAC_snn_res.0.5", "RNA_snn_res.0.5", "batch"]

def main(probe_dir, out_path=None, rna_h5ad=None):
    gz = os.path.join(probe_dir, "genescore_5977.mtx.gz")
    if not os.path.exists(gz):
        gz = os.path.join(probe_dir, "genescore_5977.mtx")
    M = sio.mmread(gz)                        # genes x nuclei (MatrixMarket is always genes-by-cells)
    feats = [l.strip() for l in open(os.path.join(probe_dir, "genescore_features.txt")) if l.strip()]
    barc = [l.strip() for l in open(os.path.join(probe_dir, "genescore_barcodes.txt")) if l.strip()]
    obs = pd.read_csv(os.path.join(probe_dir, "genescore_obs.tsv.gz"), sep="\t", index_col=0)

    X = sp.csr_matrix(M.tocsc().T.astype(np.float64))          # -> nuclei x genes
    X.sort_indices()
    assert X.shape == (len(barc), len(feats)), (X.shape, len(barc), len(feats))
    obs = obs.reindex(barc)
    missing = [c for c in OBS_COLS if c not in obs.columns]
    assert not missing, "obs columns missing from the R export: %s" % missing
    obs = obs[OBS_COLS].copy()
    # the six subtype labels come from the RNA-side annotated object, not from the R export
    rna_h5ad = rna_h5ad or os.environ.get("DSCI_OWN_RNA_H5AD",
                                          "DSCI_FY_ROOT/3k_dsci_oligo_subset.h5ad")
    six = ad.read_h5ad(rna_h5ad, backed="r").obs["cell_type"].astype(str)
    six = six.reindex(obs.index)
    assert six.notna().all(), "some nuclei have no Fig 2 subtype label"
    obs["cell_type"] = pd.Categorical(six.values)
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(feats)))
    adata.var_names_make_unique()
    out_path = out_path or os.path.join(probe_dir, "genescore_oligo_5977.h5ad")
    adata.write_h5ad(out_path)
    print("[build] %s  X %s nnz %d  data sum %.6f  log1p scale: min %.4f max %.4f"
          % (out_path, adata.shape, X.nnz, X.data.sum(), X.data.min(), X.data.max()))
    v = np.expm1(X.data)                   # undo the log1p: the underlying ArchR gene score
    s1 = np.asarray(sp.csr_matrix((v, X.indices, X.indptr), shape=X.shape).sum(1)).ravel()
    print("[build] underlying gene score (expm1 of X): per-nucleus sum median %.1f  "
          "(cf. exactly 1e4 for the UNLOGGED cross-chemistry matrices of panel D)" % np.median(s1))
    return adata


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(*sys.argv[1:])
