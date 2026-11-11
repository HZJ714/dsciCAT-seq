# -*- coding: utf-8 -*-
"""A3_GSE118257_delta2_full74.py -- mirror identity axes for GSE118257, all 74 traits, 17,799 nuclei.

Caveat carried over from the original run (and the reason this cohort's delta2 columns are read with
care in Panel C): non-OL cells are very sparse here - the smallest donor has 1 IMM nucleus, 23 NEU,
44 AST - so donor-level medians are unstable.  The numbers are used as corroboration in the supplement
and the thin-pool consequence is what Panel C marks with `*` on the cohort whose NEU pool is thinnest.

Axis definitions (this cohort's annotation vocabulary):
  d1 = median(MAT) - median(PRE)        MAT = Oligo1-6, PRE = COPs + OPCs
  dA = median(AST) - median(IMM)
  dB = median(NEU) - median(GLIA)       GLIA = every non-neuronal nucleus in this release
Frozen parameters: n_ctrl=100, random_seed=0, MIN_CELL=20 nuclei per donor per group (a donor enters
an axis only if both of its groups reach the floor), >= 4 donors else the axis is NaN.

Output: OUT/tables/trait_sweep_delta2_GSE118257_74.csv (vendored into fig7/data/ as a Panel C input).
Runtime: scDRS scoring dominates.  Run under the scoring kernel
$DSCI_PY (scdrs 1.0.2).
"""
import gzip
import os
import time

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import wilcoxon, spearmanr
import anndata as ad
import scdrs

# Root of the author's analysis tree (FY/, figures/, routeA_probe/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
MAT = WORK + "/routeA_probe/GSE118257_matrix.txt.gz"
ANN = WORK + "/routeA_probe/GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz"
GS = WORK + "/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs"
OUT = WORK + "/figures/routeA_out/A3_GSE118257_v1"
KEY = None   # None = all 74 traits
G = {"MAT": ["Oligo%d" % i for i in range(1, 7)], "PRE": ["COPs", "OPCs"],
     "AST": ["Astrocytes", "Astrocytes2"],
     "IMM": ["Microglia_Macrophages", "Macrophages", "Immune_cells"],
     "NEU": ["Neuron%d" % i for i in range(1, 6)]}
MIN_CELL = 20
t0 = time.time()


def log(*a):
    print(*a, flush=True)


with gzip.open(MAT, "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
data = pd.read_csv(MAT, sep="\t", header=None, skiprows=1, index_col=0)
data.index = data.index.astype(str)
data.columns = list(header)
ann = pd.read_csv(ANN, sep="\t")
ann["Detected"] = ann["Detected"].astype(str)
common = [b for b in header if b in set(ann["Detected"])]
ann2 = ann.set_index("Detected").loc[common].copy()
genes = data.index.tolist()
Xcsr = sparse.csr_matrix(data[common].values.astype(np.float32).T)
del data
log("matrix %d nuclei x %d genes (%.0fs)" % (Xcsr.shape[0], Xcsr.shape[1], time.time() - t0))

adata = ad.AnnData(X=Xcsr, obs=ann2, var=pd.DataFrame(index=genes))
adata.obs_names = ann2.index.tolist()
scdrs.pp.preprocess(adata)
ann2["patient"] = ann2["Sample"].str.split("_").str[0]
ct = ann2["Celltypes"].astype(str)
groups = {k: ct.isin(v).values for k, v in G.items()}
groups["GLIA"] = ~(ct.isin(G["NEU"]).values)
for k, v in groups.items():
    log("  %s: %d nuclei" % (k, v.sum()))

with open(GS) as fh:
    traits = {}
    for line in fh:
        sp = line.rstrip("\n").split("\t")
        if sp[0] == "TRAIT":
            continue
        traits[sp[0]] = [g for g in sp[1].split(",") if g in genes]


def donor_med(score, mask):
    s = pd.Series(score[mask], index=ann2["patient"][mask].values)
    s = s.groupby(level=0)
    return s.median(), s.size()


def axis_diff(score, a, b):
    x, nx = donor_med(score, groups[a]); y, ny = donor_med(score, groups[b])
    keep = [p for p in x.index.intersection(y.index)
            if nx[p] >= MIN_CELL and ny[p] >= MIN_CELL]
    if len(keep) < 4:
        return np.nan, len(keep), np.array([])
    v = x[keep].values - y[keep].values
    return float(np.median(v)), len(keep), v


keys = KEY if KEY else [k for k in traits if k != "TRAIT"]
log("n traits %d" % len(keys))
rows = []
for i, t in enumerate(keys):
    ts = time.time()
    res = scdrs.method.score_cell(adata, gene_list=traits[t], n_ctrl=100,
                                  random_seed=0, verbose=False)
    sc = res["norm_score"].values
    dA, nA, va = axis_diff(sc, "AST", "IMM")
    dB, nB, vb = axis_diff(sc, "NEU", "GLIA")
    d1, n1, v1 = axis_diff(sc, "MAT", "PRE")
    r = {"trait": t,
         "d1_OLIGO_minus_OPC": d1, "P1": wilcoxon(v1).pvalue if len(v1) > 2 else np.nan,
         "n1_donors": n1, "n1_same_sign": int((np.sign(v1) == np.sign(d1)).sum()) if len(v1) else 0,
         "dA_AST_minus_IMM": dA, "PA": wilcoxon(va).pvalue if len(va) > 2 else np.nan,
         "nA_donors": nA, "nA_same_sign": int((np.sign(va) == np.sign(dA)).sum()) if len(va) else 0,
         "dB_NEU_minus_GLIA": dB, "PB": wilcoxon(vb).pvalue if len(vb) > 2 else np.nan,
         "nB_donors": nB, "nB_same_sign": int((np.sign(vb) == np.sign(dB)).sum()) if len(vb) else 0}
    rows.append(r)
    log("[%d/%d] %s d1=%+.3f(n%d) dA=%+.3f(n%d) dB=%+.3f(n%d) (%.0fs)"
        % (i + 1, len(keys), t, d1, n1, dA, nA, dB, nB, time.time() - ts))

tab = pd.DataFrame(rows)
for c, nm in [("d1_OLIGO_minus_OPC", "r1"), ("dA_AST_minus_IMM", "rA"), ("dB_NEU_minus_GLIA", "rB")]:
    tab[nm] = tab[c].rank(ascending=False)
os.makedirs(os.path.join(OUT, "tables"), exist_ok=True)
tab.to_csv(os.path.join(OUT, "tables", "trait_sweep_delta2_GSE118257_74.csv"), index=False)
log("\ndone %.0fs" % (time.time() - t0))
log(tab[["trait", "d1_OLIGO_minus_OPC", "r1", "dA_AST_minus_IMM", "rA", "dB_NEU_minus_GLIA", "rB",
         "nA_donors", "nB_donors"]].to_string(index=False))
for a, b in [("d1_OLIGO_minus_OPC", "dB_NEU_minus_GLIA"), ("d1_OLIGO_minus_OPC", "dA_AST_minus_IMM")]:
    ok = tab[[a, b]].dropna()
    log("rho(%s,%s) = %.2f (n=%d)" % (a, b, spearmanr(ok[a], ok[b])[0], len(ok)))
