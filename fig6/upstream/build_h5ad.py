#!/usr/bin/env python3
"""Assemble an AnnData .h5ad from the tables written by rds_to_dump.R.

Usage: python build_h5ad.py <dump_dir> <out.h5ad>

Figure 6 input #1: 3k_dsci_oligo_subset.h5ad (expression channel of the own dsciCAT-seq nuclei).
The matrix is taken from the Seurat *data* slot (matrix_data.mtx) when present, i.e.
log1p(counts / column sum * 1e4): scDRS is called with --flag-raw-count False and internally
applies expm1(), so a log-transformed, size-factor-normalised matrix is required. The counts slot
is exported too (matrix.mtx) but is not what the figure scores.
"""
import os
import sys

import numpy as np
import pandas as pd
import scipy.io
import anndata as ad


def main(dump_dir, out_h5ad):
    data_mtx = os.path.join(dump_dir, "matrix_data.mtx")
    if os.path.exists(data_mtx):
        m = scipy.io.mmread(data_mtx)
        print("[py] using the data slot (log1p, 1e4-normalised)")
    else:
        m = scipy.io.mmread(os.path.join(dump_dir, "matrix.mtx"))
        print("[py] WARNING: data slot absent, using the counts slot")
    genes = pd.read_csv(os.path.join(dump_dir, "genes.tsv"), header=None)[0].astype(str).values
    barcodes = pd.read_csv(os.path.join(dump_dir, "barcodes.tsv"), header=None)[0].astype(str).values

    X = m.tocsr().T                       # MatrixMarket is genes x cells -> store cells x genes
    assert X.shape[0] == len(barcodes) and X.shape[1] == len(genes), \
        "dimension mismatch %s vs %dx%d" % (X.shape, len(barcodes), len(genes))

    meta = pd.read_csv(os.path.join(dump_dir, "meta.tsv"), sep="\t", index_col=0).reindex(barcodes)
    adata = ad.AnnData(X=X, obs=meta, var=pd.DataFrame(index=genes))
    adata.var_names_make_unique()

    umap_path = os.path.join(dump_dir, "umap.tsv")
    if os.path.exists(umap_path):
        um = pd.read_csv(umap_path, sep="\t", index_col="barcode").reindex(barcodes)
        adata.obs["UMAP_1"] = um["UMAP_1"].values
        adata.obs["UMAP_2"] = um["UMAP_2"].values
        print("[py] UMAP coordinates added to obs")

    print("[py] h5ad shape:", adata.shape)
    print("[py] obs columns:", list(adata.obs.columns))
    print("[py] cell_type:", adata.obs["cell_type"].value_counts().to_dict()
          if "cell_type" in adata.obs.columns else "MISSING")
    print("[py] X integer-valued (i.e. counts, not log1p):",
          np.issubdtype(adata.X.dtype, np.integer) or bool((adata.X.data == np.floor(adata.X.data)).all()))
    adata.write_h5ad(out_h5ad)
    print("[py] written:", out_h5ad)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
