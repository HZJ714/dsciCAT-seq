# -*- coding: utf-8 -*-
"""A2 GSE277430: the full 74-trait scDRS sweep - donor-level mature-minus-precursor difference
and the rank table. n_ctrl=200 per trait (~40 s/trait, ~50 min in total). Writes tables/trait_sweep_74.csv.
The loading and grouping logic is reused from the A2 main script.
"""
import gzip
import os
import time

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import mannwhitneyu, wilcoxon
import anndata as ad
import scdrs

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

OUT = WORK + "/figures/routeA_out/A2_GSE277430_v1"
IN = WORK + "/routeA_probe/GSE277430"
GS = WORK + "/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs"
t0 = time.time()

def log(*a):
    print(*a, flush=True)

# ---- read the OL subset (as in the A2 main script) ----
genes = pd.read_csv(os.path.join(IN, "OL_subset", "OL_genes.tsv"), header=None)[0].astype(str).tolist()
barcodes = pd.read_csv(os.path.join(IN, "OL_subset", "OL_barcodes.tsv"), header=None)[0].astype(str).tolist()
meta = pd.read_csv(os.path.join(IN, "OL_subset", "OL_meta.tsv"), sep="\t")
n_cells = len(barcodes)
n_genes = len(genes)
with open(os.path.join(IN, "OL_subset", "OL_csr.bin"), "rb") as fh:
    indptr = np.fromfile(fh, dtype=np.int32, count=n_cells + 1)
    indices = np.fromfile(fh, dtype=np.int32, count=indptr[-1])
    data = np.fromfile(fh, dtype=np.float64, count=indptr[-1])
from scipy.sparse import csc_matrix
M = csc_matrix((data, indices, indptr), shape=(n_genes, n_cells)).T.tocsr()
log("OL: %d x %d" % (n_genes, n_cells))

adata = ad.AnnData(X=M, obs=meta.set_index("barcodes").loc[barcodes],
                   var=pd.DataFrame(index=genes))
adata.obs_names = barcodes
scdrs.pp.preprocess(adata)
log("preprocess ok %.0fs" % (time.time() - t0))

meta["patient"] = meta["patient"].astype(str)
is_mat = (meta["expertAnno.l1"] == "OLIGO").values
is_pre = (meta["expertAnno.l1"] == "OPC").values

# ---- read the 74 trait gene sets ----
with open(GS) as fh:
    traits = {}
    for line in fh:
        sp = line.rstrip("\n").split("\t")
        traits[sp[0]] = [g for g in sp[1].split(",") if g in genes]

rows = []
for i, (trait, glist) in enumerate(traits.items()):
    ts = time.time()
    res = scdrs.method.score_cell(adata, gene_list=glist, n_ctrl=200,
                                  random_seed=0, verbose=False)
    score = res["norm_score"].values
    mat_med = pd.Series(score[is_mat]).groupby(meta["patient"][is_mat].values).median()
    pre_med = pd.Series(score[is_pre]).groupby(meta["patient"][is_pre].values).median()
    common = mat_med.index.intersection(pre_med.index)
    d = mat_med[common].values - pre_med[common].values
    w = wilcoxon(d) if len(d) >= 3 else None
    rows.append({"trait": trait, "n_genes_hit": len(glist),
                 "median_diff": np.median(d), "P_wilcoxon": w.pvalue if w else np.nan,
                 "n_same_sign": int((np.sign(d) == np.sign(np.median(d))).sum()), "n_donors": len(d)})
    log("[%d/74] %s hit=%d diff=%.3f P=%.2g (%.0fs)"
        % (i + 1, trait, len(glist), rows[-1]["median_diff"], rows[-1]["P_wilcoxon"],
           time.time() - ts))

tab = pd.DataFrame(rows).sort_values("median_diff", ascending=False).reset_index(drop=True)
tab.insert(0, "rank", np.arange(1, len(tab) + 1))
tab.to_csv(os.path.join(OUT, "tables", "trait_sweep_74.csv"), index=False)
log("\ndone in %.0fs -> tables/trait_sweep_74.csv" % (time.time() - t0))
print(tab[["rank", "trait", "n_genes_hit", "median_diff", "P_wilcoxon"]].head(15).to_string())
