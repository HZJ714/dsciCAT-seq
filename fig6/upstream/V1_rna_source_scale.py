#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V1: check whether the RNA channel's two sources are on the SAME LAYER and the SAME SCALE.

  A = the RNA used by C: X of $DSCI_WORK_ROOT/FY/oligo/3k_dsci_oligo_subset.h5ad
      (the Seurat LogNormalize data layer, log1p(CP10k), SoupX-corrected counts)
  B = the RNA used by D: the ArchR GeneIntegrationMatrix (Info/Units="NormCounts", the CCA-integrated
      linear matrix whose per-nucleus sum is exactly 1e4), exported here as rnalin_/rna_ mtx.

To answer: are the scales the same (can they be mixed directly), and are the same-source nuclei ordered the same (does switching sources move the conclusion).
Outputs routeA_out/D_crosschemistry/tables/V1_rna_source_scale.csv and the per-nucleus correlation distribution.
"""
import os, numpy as np, pandas as pd, anndata as ad
from scipy.io import mmread
from scipy.stats import spearmanr

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

OUT = os.path.join(os.environ.get('DSCI_FIG_ROOT', WORK + "/figures"), 'routeA_out/D_crosschemistry')
DER, TAB = OUT + "/derived", OUT + "/tables"
H5_R = os.path.join(os.environ.get('DSCI_FY_ROOT', WORK + "/FY"), 'oligo/3k_dsci_oligo_subset.h5ad')
SAY = lambda fmt, *a: print("[%s] %s" % (pd.Timestamp.now().strftime("%H:%M:%S"), fmt % a if a else fmt), flush=True)

C = ad.read_h5ad(H5_R)                       # 5,977 x 31,263, X = log1p(CP10k)
Xc = C.X.tocsr().astype(np.float64)
cg = [str(g) for g in C.var_names]
cc = [str(i) for i in C.obs_names]
SAY("C's RNA h5ad: %s | range %.3f..%.3f | per-nucleus non-zero-sum median %.0f", Xc.shape, Xc.data.min(), Xc.data.max(),
    np.median(np.asarray(Xc.sum(axis=1)).ravel()))

from scipy.sparse import load_npz
D = load_npz(os.path.join(DER, "rnalin_dsciCAT_seq.npz")).tocsc()   # genes x cells linear CP10k
dg = [l.strip() for l in open(os.path.join(DER, "rna_dsciCAT_seq_genes.txt")) if l.strip()]
dc = pd.read_csv(os.path.join(TAB, "cells_dsciCAT_seq.csv")).cell.tolist()
SAY("ArchR GIM (linear): %s | range %.4f..%.3f | per-nucleus sum median %.0f", D.shape, D.data.min(), D.data.max(),
    np.median(np.asarray(D.sum(axis=0)).ravel()))

gs = sorted(set(x.upper() for x in cg) & set(x.upper() for x in dg))
SAY("gene-name intersection (case-insensitive) %d | C-only %d | GIM-only %d",
    len(gs), len(set(x.upper() for x in cg) - set(x.upper() for x in dg)), len(set(x.upper() for x in dg) - set(x.upper() for x in cg)))
cs = sorted(set(cc) & set(dc))
SAY("nucleus intersection %d / C %d / GIM %d", len(cs), len(cc), len(dc))

upc = {g.upper(): i for i, g in enumerate(cg)}
upd = {g.upper(): i for i, g in enumerate(dg)}
ci = np.array([cc.index(c) for c in cs])
di = np.array([dc.index(c) for c in cs])
cix = np.array([upc[g] for g in gs]); dix = np.array([upd[g] for g in gs])

A = np.asarray(Xc[np.ix_(ci, cix)].todense(), dtype=np.float64)   # cells x genes, log1p(CP10k)
B = np.asarray(D[np.ix_(dix, di)].T.todense(), dtype=np.float64)  # cells x genes, linear CP10k
Bloc = np.log1p(B)

rows = []
for name, M in [("C=Seurat data (log1p)", A), ("GIM log1p", Bloc), ("GIM linear", B)]:
    v = M.ravel(); v = v[v != 0]
    rows.append(dict(source=name, n_cells=len(cs), n_genes=len(gs),
                     nonzero_mean=v.mean(), p50=np.percentile(v, 50), p99=np.percentile(v, 99), maxv=v.max(),
                     cellsum_median=np.median(M.sum(axis=1)), frac_nonzero=(M != 0).mean()))
tab = pd.DataFrame(rows)
SAY("scale comparison:\n" + tab.to_string(index=False))

# per-nucleus Spearman (the same nucleus: do the two sources order the same genes the same way)
def rank_sp(a, b):
    ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])
pr = np.array([rank_sp(A[i], Bloc[i]) for i in range(len(cs))])
pl = np.array([rank_sp(A[i], B[i]) for i in range(len(cs))])
# control: permute the GIM nuclei and look at the null
rng = np.random.default_rng(20261103)
pn = np.array([rank_sp(A[i], Bloc[rng.integers(0, len(cs), len(cs))]) for i in range(200)])
SAY("per-nucleus Spearman (log1p vs log1p): median %.4f | 1%%-99%% %.4f..%.4f | shuffled-null median %.4f",
    np.median(pr), np.percentile(pr, 1), np.percentile(pr, 99), np.median(pn))
SAY("per-nucleus Spearman (log1p vs linear): median %.4f (must agree; Spearman is invariant to monotone transforms)", np.median(pl))

# per gene: the cross-nucleus mean-profile correlation
gm_a = A.mean(0); gm_b = Bloc.mean(0)
SAY("gene mean-profile Pearson (log1p vs log1p) = %.4f | Spearman = %.4f",
    np.corrcoef(gm_a, gm_b)[0, 1], spearmanr(gm_a, gm_b).correlation)

extra = dict(per_cell_rho_log1p_median=float(np.median(pr)),
             per_cell_rho_log1p_q01=float(np.percentile(pr, 1)),
             per_cell_rho_log1p_q99=float(np.percentile(pr, 99)),
             perm_null_median=float(np.median(pn)),
             gene_mean_pearson=float(np.corrcoef(gm_a, gm_b)[0, 1]),
             n_gene_shared=len(gs), n_cell_shared=len(cs))
tab = pd.concat([tab, pd.DataFrame([extra])], ignore_index=True)
tab.to_csv(os.path.join(TAB, "V1_rna_source_scale.csv"), index=False)
SAY("written %s", os.path.join(TAB, "V1_rna_source_scale.csv"))
pd.DataFrame({"cell": cs, "rho_log_vs_log": pr}).to_csv(os.path.join(TAB, "V1_per_cell_rho.csv"), index=False)
