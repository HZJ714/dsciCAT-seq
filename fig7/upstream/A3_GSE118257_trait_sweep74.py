# -*- coding: utf-8 -*-
"""GSE118257, all 74 traits (reviewer gap 1): n_ctrl=200, identical to GSE277430.
Donor-level median of mature (Oligo1-6) minus precursor (COPs+OPCs), with Wilcoxon and a same-sign count.
Output: tables/trait_sweep_74_GSE118257.csv (sorted by rank)
"""
import gzip
import os
import time
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import wilcoxon, mannwhitneyu
import anndata as ad
import scdrs

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

OUT = WORK + "/figures/routeA_out/A3_GSE118257_v1"
MAT = WORK + "/routeA_probe/GSE118257_matrix.txt.gz"
ANN = WORK + "/routeA_probe/GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz"
GS = WORK + "/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs"
os.makedirs(os.path.join(OUT, "tables"), exist_ok=True)
t0 = time.time()

def log(*a):
    print(*a, flush=True)

# ---------- read the matrix and the annotations (as in A3) ----------
with gzip.open(MAT, "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
data = pd.read_csv(MAT, sep="\t", header=None, skiprows=1, index_col=0)
data.index = data.index.astype(str)
data.columns = list(header)
ann = pd.read_csv(ANN, sep="\t")
ann["Detected"] = ann["Detected"].astype(str)
common = [b for b in header if b in set(ann["Detected"])]
ann2 = ann.set_index("Detected").loc[common].copy()
X = data[common].values.astype(np.float32)
Xcsr = sparse.csr_matrix(X.T)
genes = data.index.tolist()
log("matrix %d nuclei x %d genes (%.0fs)" % (Xcsr.shape[0], Xcsr.shape[1], time.time() - t0))

adata = ad.AnnData(X=Xcsr, obs=ann2, var=pd.DataFrame(index=genes))
adata.obs_names = ann2.index.tolist()
scdrs.pp.preprocess(adata)
log("preprocess ok (%.0fs)" % (time.time() - t0))

mature = ["Oligo%d" % i for i in range(1, 7)]
prec = ["COPs", "OPCs"]
ann2["patient"] = ann2["Sample"].str.split("_").str[0]
m_mask = ann2["Celltypes"].isin(mature).values
p_mask = ann2["Celltypes"].isin(prec).values
log("mature %d / precursor %d / donors %d" % (m_mask.sum(), p_mask.sum(),
    ann2[m_mask | p_mask]["patient"].nunique()))

with open(GS) as fh:
    traits = {}
    for line in fh:
        sp = line.rstrip("\n").split("\t")
        traits[sp[0]] = [g for g in sp[1].split(",") if g in genes]

rows = []
for i, (trait, glist) in enumerate(traits.items()):
    ts = time.time()
    res = scdrs.method.score_cell(adata, gene_list=glist, n_ctrl=200,
                                  random_seed=20261029, verbose=False)
    sc = res["norm_score"].astype(float).to_numpy()
    dm = pd.Series(sc[m_mask], index=ann2["patient"][m_mask].values).groupby(level=0).median()
    dp = pd.Series(sc[p_mask], index=ann2["patient"][p_mask].values).groupby(level=0).median()
    common_d = sorted(set(dm.index) & set(dp.index))
    diff = dm[common_d].values - dp[common_d].values
    w, pw = wilcoxon(diff)
    rows.append({"trait": trait, "n_donors": len(diff),
                 "median_diff": float(np.median(diff)),
                 "P_wilcoxon": float(pw),
                 "same_sign": int((diff > 0).sum())})
    log("[%d/74] %s hit=%d diff=%+.3f P=%.2g n=%d (%.0fs)" % (
        i + 1, trait, len(glist), rows[-1]["median_diff"], pw, len(diff),
        time.time() - ts))

tab = pd.DataFrame(rows).sort_values("median_diff", ascending=False)
tab.insert(0, "rank", range(1, len(tab) + 1))
tab.to_csv(os.path.join(OUT, "tables", "trait_sweep_74_GSE118257.csv"), index=False)
ms = tab[tab["trait"] == "PASS_Multiple_sclerosis"]
log("\nMS rank = %d/74 (empirical p ~ %.3f), diff=%+.3f, P=%.2g, same sign %d/%d" % (
    ms["rank"].values[0], ms["rank"].values[0] / 74, ms["median_diff"].values[0],
    ms["P_wilcoxon"].values[0], ms["same_sign"].values[0], ms["n_donors"].values[0]))
log("done in %.0fs -> tables/trait_sweep_74_GSE118257.csv" % (time.time() - t0))
