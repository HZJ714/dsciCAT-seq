# -*- coding: utf-8 -*-
"""A4 GSE301908: the full 74-trait scDRS sweep - the trait-specificity evidence of the third cohort.
Same criterion as A2/A3: donor-level median(mature Oligo) - median(OPC), n_ctrl=200, random_seed=0.
Output tables/trait_sweep_74_GSE301908.csv (with a rank column, 74 rows, .gs header junk already dropped).
"""
import os
import time

import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.stats import wilcoxon
import anndata as ad
import scdrs

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

IN = WORK + "/routeA_probe/GSE301908"
OUT = WORK + "/figures/routeA_out/A4_GSE301908_v1"
GS = WORK + "/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs"
t0 = time.time()


def log(*a):
    print(*a, flush=True)


genes = pd.read_csv(os.path.join(IN, "OL_genes.tsv"), header=None)[0].astype(str).tolist()
barcodes = pd.read_csv(os.path.join(IN, "OL_barcodes.tsv"), header=None)[0].astype(str).tolist()
meta = pd.read_csv(os.path.join(IN, "OL_meta.tsv"), sep="\t")
n_cells, n_genes = len(barcodes), len(genes)
with open(os.path.join(IN, "OL_csr.bin"), "rb") as fh:
    indptr = np.fromfile(fh, dtype=np.int32, count=n_cells + 1)
    indices = np.fromfile(fh, dtype=np.int32, count=indptr[-1])
    data = np.fromfile(fh, dtype=np.float64, count=indptr[-1])
M = csc_matrix((data, indices, indptr), shape=(n_genes, n_cells)).T.tocsr()
log("OL: %d genes x %d nuclei (nnz=%d)" % (n_genes, n_cells, len(data)))

ct = meta["majorCluster"].astype(str)
is_mat = ct.str.lower().str.contains("oligo").values
is_pre = ct.str.lower().str.contains("opc").values
log("mature nuclei %d, precursor nuclei %d" % (is_mat.sum(), is_pre.sum()))

raw = meta["orig.ident"].astype(str)
meta["patient"] = raw.str.extract(r"^(ms\d+)")[0].fillna(raw)
log("%d donors" % meta["patient"].nunique())

adata = ad.AnnData(X=M, obs=meta.copy(), var=pd.DataFrame(index=genes))
adata.obs_names = barcodes
scdrs.pp.preprocess(adata)
log("preprocess ok %.0fs" % (time.time() - t0))

traits = {}
with open(GS) as fh:
    for line in fh:
        sp = line.rstrip("\n").split("\t")
        if sp[0] == "TRAIT":
            continue
        traits[sp[0]] = [g for g in sp[1].split(",") if g in genes]
log("number of traits = %d" % len(traits))

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
                 "n_same_sign": int((np.sign(d) == np.sign(np.median(d))).sum()),
                 "n_donors": len(d)})
    log("[%d/%d] %s hit=%d diff=%.3f P=%.2g same sign %d/%d (%.0fs)"
        % (i + 1, len(traits), trait, len(glist), rows[-1]["median_diff"],
           rows[-1]["P_wilcoxon"], rows[-1]["n_same_sign"], len(d), time.time() - ts))
    pd.DataFrame(rows).sort_values("median_diff", ascending=False).reset_index(
        drop=True).assign(rank=lambda x: np.arange(1, len(x) + 1)).to_csv(
        os.path.join(OUT, "tables", "trait_sweep_74_GSE301908.csv"), index=False)

tab = pd.DataFrame(rows).sort_values("median_diff", ascending=False).reset_index(drop=True)
tab.insert(0, "rank", np.arange(1, len(tab) + 1))
tab.to_csv(os.path.join(OUT, "tables", "trait_sweep_74_GSE301908.csv"), index=False)
log("\ndone in %.0fs -> tables/trait_sweep_74_GSE301908.csv" % (time.time() - t0))
log(tab[["rank", "trait", "median_diff", "P_wilcoxon"]].head(12).to_string())
