#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D4: build the scDRS inputs (h5ad + cov) for 3 chemistries x 2 channels, the conventions copied verbatim from own-C's set.

  RNA channel  X = log1p(ArchR GeneIntegrationMatrix)   cov = nCount_RNA, nFeature_RNA
  ATAC channel X = ArchR GeneScoreMatrix (linear CP10k)  cov = nCount_ATAC, nFeature_ATAC
 - the same class and parameters as the h5ad/cov that C (O1) used; neither the covariate names nor their count are re-tuned for D.
  column order = the rows of tables/cells_{tag}.csv; the two channels pair nucleus-by-nucleus (same-nucleus pairing).
"""
import os, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from scipy.io import mmread
from scipy.sparse import load_npz

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

OUT = os.path.join(os.environ.get('DSCI_FIG_ROOT', WORK + "/figures"), 'routeA_out/D_crosschemistry')
DER, TAB, H5 = OUT + "/derived", OUT + "/tables", OUT + "/h5ad"
os.makedirs(H5, exist_ok=True)
SAY = lambda fmt, *a: print("[%s] %s" % (pd.Timestamp.now().strftime("%H:%M:%S"), fmt % a if a else fmt), flush=True)
TAGS = ["dsciCAT_seq", "10x", "HT_scCAT_seq_HT"]
CH = {"rna": dict(mtx="rna_%s.npz", genes="rna_%s_genes.txt", cov=["nCount_RNA", "nFeature_RNA"]),
      "atac": dict(mtx="gs_%s.npz", genes="gs_%s_genes.txt", cov=["nCount_ATAC", "nFeature_ATAC"])}

cov = pd.read_csv(os.path.join(TAB, "cov_by_cell.csv"))
cov = cov.drop_duplicates("cell").set_index("cell")
SAY("cov_by_cell: %d rows | columns %s", len(cov), list(cov.columns))

for tag in TAGS:
    cells = pd.read_csv(os.path.join(TAB, "cells_%s.csv" % tag))
    for ch, cfg in CH.items():
        M = (load_npz(os.path.join(DER, cfg["mtx"] % tag)) if cfg["mtx"].endswith("npz")
             else mmread(os.path.join(DER, cfg["mtx"] % tag))).tocsc()      # genes x cells
        genes = [l.strip() for l in open(os.path.join(DER, cfg["genes"] % tag)) if l.strip()]
        assert M.shape[0] == len(genes), "gene-name count disagrees"
        assert M.shape[1] == len(cells), "%s/%s nuclei %d != cells table %d" % (tag, ch, M.shape[1], len(cells))
        obs = cells.copy()
        obs.index = obs.cell.astype(str)
        missing = [c for c in obs.index if c not in cov.index]
        if missing:
            SAY("  !! %d nuclei lack covariates (e.g. %s)", len(missing), missing[:3])
        for c in cfg["cov"]:
            obs[c] = pd.to_numeric(cov[c].reindex(obs.index).values, errors="coerce")
        X = sp.csc_matrix(M.T)                                            # cells x genes
        A = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=np.array(genes)))
        A.obs_names = list(obs.index); A.var_names_make_unique()
        fp = os.path.join(H5, "%s_%s.h5ad" % (tag, ch))
        A.write_h5ad(fp)
        cv = obs[[c for c in cfg["cov"]]].copy(); cv.insert(0, "index", cv.index)
        fc = os.path.join(H5, "%s_%s.cov.tsv" % (tag, ch))
        cv.to_csv(fc, sep="\t", index=False)
        v = X.data
        SAY("  %-14s %-5s | %s cells x %d genes | nnz=%d | range %.4f..%.3f | cov missing %d | %s",
            tag, ch, X.shape[0], X.shape[1], len(v), v.min(), v.max(), int(cv[cfg["cov"]].isna().sum().sum()),
            os.path.basename(fp))
        SAY("      cov: %s", fc)
SAY("ALL_FINISHED")
