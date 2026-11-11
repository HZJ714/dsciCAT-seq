# -*- coding: utf-8 -*-
"""p0_8c_c2_stability.py - how stable is the row clustering (modules): cell-level bootstrap.

Motivation: the merged Fig6C planned to cluster the 38 rows into modules by expression profile (the desired figure style). Clustering 38 rows x 6 columns
is inherently fragile; the reproducibility rate must be measured first, else "modules" are just the shape of noise.
Method: resample cells with replacement within each subtype -> recompute the six subtype means -> row z -> ward -> cut at k ->
compare with the observed partition by ARI (adjusted Rand index, invariant to label permutation). ATAC is done the same way.
Output: figures/routeA_out/P0_8C/tables/p08c_c2_cluster_stability.csv
"""
import os
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import adjusted_rand_score as ari
import warnings; warnings.filterwarnings("ignore")

W = os.environ.get('DSCI_WORK_ROOT', '')
OUT = f"{W}/figures/routeA_out/P0_8C/tables"; os.makedirs(OUT, exist_ok=True)
CACHE = f"{OUT}/p08c_c2_celllevel38.npz"
SUBS = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
NBOOT, SEED = 500, 20261116

sub = pd.read_csv(f"{OUT}/p08c_c1_rna_atac_by_subtype.csv")
genes = list(sub.mouse.values)

if os.path.exists(CACHE):
    z = np.load(CACHE, allow_pickle=True)
    XG, AG, key, genes = z["XG"], z["AG"], z["key"], z["genes"]
    print("reading the cached cell-level 38-gene matrices", XG.shape, AG.shape)
else:
    ar = ad.read_h5ad(f"{W}/FY/oligo/3k_dsci_oligo_subset.h5ad", backed="r")
    key = ar.obs["cell_type"].astype(str).values
    av = list(map(str, ar.var_names))
    cols = np.array([av.index(g) for g in genes])
    XG = np.asarray(sp.csr_matrix(ar.X[:, cols]).todense(), dtype=np.float32)   # 5962×38
    print("RNA cell-level", XG.shape)
    aa = ad.read_h5ad(f"{W}/FY/oligo/fig6_rerun_20261019/dualomics_probe/genescore_oligo_5977.h5ad", backed="r")
    assert list(aa.obs_names) == list(ar.obs_names), "the two objects' cell orders differ; align by barcode"
    anv = list(map(str, aa.var_names))
    AG = np.zeros((XG.shape[0], len(genes)), dtype=np.float32)
    for i, g in enumerate(genes):                              # column by column, avoiding full-matrix densification (hit with 24,333 columns)
        j = anv.index(g) if g in anv else -1
        if j >= 0:
            AG[:, i] = np.asarray(aa.X[:, j].todense()).ravel()
    np.savez_compressed(CACHE, XG=XG, AG=AG, key=key, genes=np.array(genes))
    print("ATAC cell-level", AG.shape, "| columns hit", int((AG.sum(0) > 0).sum()), "/", len(genes))

rows = {s: np.flatnonzero(key == s) for s in SUBS}
print("cell counts", {s: len(v) for s, v in rows.items()})


def prof(M, idx):
    """Given cell indices, return the 38 x 6 subtype means."""
    return np.column_stack([M[idx[s]].mean(0) for s in SUBS]) if isinstance(idx, dict) else \
        np.column_stack([M[idx[SUBS.index(s)]].mean(0) for s in SUBS])


def zrow(P):
    return (P - P.mean(1, keepdims=True)) / np.maximum(P.std(1, keepdims=True), 1e-9)


def part(P, k):
    return fcluster(linkage(zrow(P), "ward"), k, "maxclust")


rng = np.random.default_rng(SEED)
obs_r, obs_a = prof(XG, rows), prof(AG, rows)
res = []
for ch, M, obs in (("RNA", XG, obs_r), ("ATAC", AG, obs_a)):
    for k in (3, 4, 5, 6):
        o = part(obs, k)
        A, per_gene = [], np.zeros((M.shape[1],))
        for b in range(NBOOT):
            idx = {s: rng.choice(rows[s], size=len(rows[s]), replace=True) for s in SUBS}
            try:
                pb = part(prof(M, idx), k)
            except Exception:
                continue
            A.append(ari(o, pb))
        A = np.array(A)
        res.append(dict(ch=ch, k=k, n_mod=np.bincount(o)[1:].tolist(),
                        ari_med=np.median(A), ari_q25=np.percentile(A, 25), ari_q75=np.percentile(A, 75),
                        frac_ari_ge70=float((A >= 0.7).mean()), frac_ari_ge50=float((A >= 0.5).mean())))
        print("%s k=%d modules%-14s ARI median %.2f [%.2f,%.2f]  ARI>=0.7 in %.0f%%  >=0.5 in %.0f%%"
              % (ch, k, str(res[-1]["n_mod"]), res[-1]["ari_med"], res[-1]["ari_q25"], res[-1]["ari_q75"],
                 100 * res[-1]["frac_ari_ge70"], 100 * res[-1]["frac_ari_ge50"]))

# joint channels (RNA+ATAC, 12 features) - the clustering convention the merged version actually uses
comb = np.column_stack([zrow(obs_r), zrow(obs_a)])
for k in (3, 4, 5):
    o = fcluster(linkage(comb, "ward"), k, "maxclust")
    A = []
    for b in range(NBOOT):
        idx = {s: rng.choice(rows[s], size=len(rows[s]), replace=True) for s in SUBS}
        c = np.column_stack([zrow(prof(XG, idx)), zrow(prof(AG, idx))])
        A.append(ari(o, fcluster(linkage(c, "ward"), k, "maxclust")))
    A = np.array(A)
    res.append(dict(ch="RNA+ATAC", k=k, n_mod=np.bincount(o)[1:].tolist(), ari_med=np.median(A),
                    ari_q25=np.percentile(A, 25), ari_q75=np.percentile(A, 75),
                    frac_ari_ge70=float((A >= 0.7).mean()), frac_ari_ge50=float((A >= 0.5).mean())))
    print("joint k=%d modules%-12s ARI median %.2f  ARI>=0.7 in %.0f%%  >=0.5 in %.0f%%"
          % (k, str(res[-1]["n_mod"]), res[-1]["ari_med"], 100 * res[-1]["frac_ari_ge70"], 100 * res[-1]["frac_ari_ge50"]))

pd.DataFrame(res).to_csv(f"{OUT}/p08c_c2_cluster_stability.csv", index=False)
print("[saved]", f"{OUT}/p08c_c2_cluster_stability.csv")
