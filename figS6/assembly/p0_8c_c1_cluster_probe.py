import os
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Feasibility probe: can C1's 38 rows take the desired figure style "columns clustered by subtype + rows clustered into modules + program annotations on the right"?
Counts / independence tests only, no figure."""
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency

TAB = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'routeA_out/P0_8C/tables')
SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
CLS = {"lipid/energy": ["APOE", "PLIN3", "EFR3B", "PHGDH", "PRR5L", "ATP1A1"],
       "myelin structure": ["CNP", "BCAS1", "PLLP"],
       "cytoskeleton/MT": ["CAMSAP2", "TUBB4A", "SCHIP1", "YPEL2", "CEP97"],
       "vesicle/membrane": ["SEC11C", "PHLDB1", "ENDOD1", "CLIC4", "TMCC3"],
       "immune/adhesion": ["LITAF", "TRIM59", "SKAP2", "AIG1", "MCAM", "IGSF8", "CD9", "ADAMTS4", "SPON1"],
       "lineage/TF": ["OLIG1", "ZNF365", "CREB5", "ZBTB38"],
       "other signalling": ["DDC", "TMEM151A", "SEZ6L2", "KCTD13", "SBF1", "RRAS2"]}
g2c = {g: c for c, gs in CLS.items() for g in gs}

S = pd.read_csv(TAB + "/p08c_c1_rna_atac_by_subtype.csv")
S["cls"] = S.gene.map(g2c)
X = S[["rna_" + s for s in SUB]].values.astype(float)
A = S[["atac_" + s for s in SUB]].values.astype(float)

# ---------- 1) does column (subtype) clustering recover the maturation axis ----------
Zc = (X - X.mean(0)) / X.std(0)
C = np.corrcoef(Zc.T)
Lc = linkage(squareform(1 - C, checks=False), "average")
ordr = [SUB[i] for i in dendrogram(Lc, no_plot=True, labels=SUB)["leaves"]]
print("column-cluster order:", " > ".join(ordr))
print("equals the maturation axis:", ordr == SUB)
print("maturation-axis adjacent correlations:", np.round([C[i, i + 1] for i in range(5)], 2), " COP/NFOL-MOL6 = %.2f" % C[0, 5])

# ---------- 2) row modules: module count and module x program independence ----------
def mods(M, k):
    Z = (M - M.mean(1, keepdims=True)) / np.maximum(M.std(1, keepdims=True), 1e-9)
    return Z, fcluster(linkage(Z, "ward"), k, "maxclust")

print()
for ch, M in (("RNA", X), ("ATAC", A)):
    for k in (3, 4, 5, 6):
        Z, m = mods(M, k)
        p = chi2_contingency(pd.crosstab(m, S.cls))[1]
        print("%s k=%d module sizes %-14s module x program chi2 P=%.3f" % (ch, k, str(np.bincount(m)[1:]), p))

# ---------- 3) look at one concrete cut (RNA k=5) and its modules ----------
Z, m = mods(X, 5)
print()
for mm in sorted(set(m)):
    r = S[S.mod == mm] if "mod" in S else S[m == mm]
    print("M%d n=%2d  z-range median %.2f  peak subtype %s" % (mm, (m == mm).sum(), np.median(np.ptp(Z[m == mm], axis=1)),
          dict(pd.Series(S.top_subtype.values[m == mm]).value_counts())))
    print("     mean z trajectory  " + "  ".join("%s %.2f" % (s, v) for s, v in zip(SUB, Z[m == mm].mean(0))))
    print("     class composition  " + str(dict(pd.Series(S.cls.values[m == mm]).value_counts())))

# ---------- 4) the "messy colours" cause: row-z amplifies low-detection noise ----------
P = S[[c for c in S.columns if c.startswith("pct")]].values
print()
print("raw log2 mean range median %.3f -> after row z %.2f (amplified %.1f x)"
      % (np.median(np.ptp(X, axis=1)), np.median(np.ptp(Z, axis=1)),
         np.median(np.ptp(Z, axis=1)) / np.median(np.ptp(X, axis=1))))
print("cells with <10%% detection %.0f%%; genes with all six subtypes >=10%% %d/38" % (100 * (P < 10).mean(), (P >= 10).all(1).sum()))
keep = (P >= 10).all(1)
Zk, mk = mods(X[keep], 4)
print("keeping only fully-detected genes n=%d: 4 modules %s, z-range median %.2f, module x program P=%.3f"
      % (keep.sum(), np.bincount(mk)[1:], np.median(np.ptp(Zk, axis=1)),
         chi2_contingency(pd.crosstab(mk, S.cls.values[keep]))[1]))

# ---------- 5) stability: do module memberships survive cell-level bootstrap resampling ----------
try:
    import anndata as ad
    print()
    print("(cell-level bootstrap stability needs the h5ad and is measured separately; this round uses the subtype-mean convention)")
except Exception:
    pass
