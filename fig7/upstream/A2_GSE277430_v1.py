# -*- coding: utf-8 -*-
"""A2 GSE277430, the human main analysis v1: scDRS MS scoring of the OL lineage, then donor-level statistics.
Inputs: OL_matrix.mtx.gz + OL_genes.tsv + OL_barcodes.tsv + OL_meta.tsv (exported by A2_GSE277430_extract.R)
Design (written before running; criteria fixed in the HANDOVER A2 note):
  primary: mature (OLIGO) vs precursor (OPC), donor-level paired Wilcoxon over 21 donors
  A2: MS_Demye / MS_Mye vs CTRL, donor-level MWU + Cliff delta (unpaired; facility / PMI / age described)
  within MS: Demye vs Mye (7 vs 6 donors)
  specificity: 5 unrelated traits + random-set percentile (same panel as GSE118257)
"""
import gzip
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.stats import mannwhitneyu, wilcoxon
import anndata as ad

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

OUT = WORK + "/figures/routeA_out/A2_GSE277430_v1"
IN = WORK + "/routeA_probe/GSE277430"
RDS_PATH = WORK + "/routeA_probe/GSE277430/GSE277430_normCount_GEOsubmission.rds.gz"
GS = WORK + "/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs"
os.makedirs(os.path.join(OUT, "tables"), exist_ok=True)
os.makedirs(os.path.join(OUT, "pdf"), exist_ok=True)
t0 = time.time()

def log(*a):
    print(*a, flush=True)

# ---- read the OL subset (CSC binary exported by R: indptr int32 / indices int32 / data float64) ----
n_genes = len(pd.read_csv(os.path.join(IN, "OL_subset", "OL_genes.tsv"), header=None)[0].astype(str).tolist())
genes = pd.read_csv(os.path.join(IN, "OL_subset", "OL_genes.tsv"), header=None)[0].astype(str).tolist()
barcodes = pd.read_csv(os.path.join(IN, "OL_subset", "OL_barcodes.tsv"), header=None)[0].astype(str).tolist()
meta = pd.read_csv(os.path.join(IN, "OL_subset", "OL_meta.tsv"), sep="\t")
n_cells = len(barcodes)
binf = os.path.join(IN, "OL_subset", "OL_csr.bin")
with open(binf, "rb") as fh:
    indptr = np.fromfile(fh, dtype=np.int32, count=n_cells + 1)
    indices = np.fromfile(fh, dtype=np.int32, count=indptr[-1])
    data = np.fromfile(fh, dtype=np.float64, count=indptr[-1])
from scipy.sparse import csc_matrix
M = csc_matrix((data, indices, indptr), shape=(n_genes, n_cells)).T.tocsr()  # cells x genes
log("OL subset: %d genes x %d nuclei (nnz=%d)" % (n_genes, n_cells, len(data)))
log("meta %d rows; barcodes aligned: %s" % (len(meta), meta["barcodes"].iloc[0] == barcodes[0]))

import scdrs
adata = ad.AnnData(X=M, obs=meta.set_index("barcodes").loc[barcodes],
                   var=pd.DataFrame(index=genes))
adata.obs_names = barcodes
# normCount data (not raw counts), so preprocess only computes gene/cell statistics and matches the control sets - it does not normalise
scdrs.pp.preprocess(adata)
log("preprocess ok %.0fs" % (time.time() - t0))

# ---- grouping ----
OL_mature = meta["expertAnno.l1"] == "OLIGO"
OL_prec = meta["expertAnno.l1"] == "OPC"
meta["patient"] = meta["patient"].astype(str)
meta["group"] = np.where(meta["disease"] == "MS",
                         np.where(meta["pathology_class"] == "CX_Demye", "MS_Demye",
                                  np.where(meta["pathology_class"] == "CX_Mye", "MS_Mye", "MS_other")),
                         np.where(meta["pathology_class"] == "WM_Ctrl", "CTRL_WM", "CTRL_CX"))
log("nuclei per group:")
log(meta.groupby("group").size().to_string())

with open(GS) as fh:
    gs = {}
    for line in fh:
        line = line.rstrip("\n")
        if line:
            n, g = line.split("\t", 1)
            gs[n] = g.split(",")
ms = gs["PASS_Multiple_sclerosis"]
log("MS set hits %d / %d" % (len([g for g in ms if g in set(genes)]), len(ms)))

def score(gl, n_ctrl=200, seed=20261029):
    gl = [g for g in gl if g in set(genes)]
    r = scdrs.method.score_cell(adata, gl, n_ctrl=n_ctrl, random_seed=seed)
    return r["norm_score"].to_numpy()

s_ms = score(ms)
meta["MS_norm_score"] = s_ms

def donor_median(col):
    return meta.groupby("patient")[col].median()

def cliffs(a, b):
    aa = np.asarray(a)[:, None]; bb = np.asarray(b)[None, :]
    return ((aa > bb).sum() - (aa < bb).sum()) / (aa.size * bb.size)

# primary criterion: mature vs precursor, donor-level paired
mat_pat = meta[OL_mature].groupby("patient")["MS_norm_score"].median()
pre_pat = meta[OL_prec].groupby("patient")["MS_norm_score"].median()
common = sorted(set(mat_pat.index) & set(pre_pat.index))
d = mat_pat[common].values - pre_pat[common].values
w, p_main = wilcoxon(d)
log("\n[primary criterion] mature (OLIGO) vs precursor (OPC), donor-level paired: n=%d donors, P=%g, median difference %+.3f, all same sign=%s"
    % (len(common), p_main, np.median(d), bool((np.sign(d) == np.sign(np.median(d))).all())))
cd_main = cliffs(s_ms[OL_mature.values], s_ms[OL_prec.values])
log("  nucleus-level Cliff delta=%.3f (mature %d vs precursor %d)" % (cd_main, OL_mature.sum(), OL_prec.sum()))

# A2: donor-level Demye/Mye vs CTRL (unpaired MWU + Cliff delta)
res_groups = []
for grp in ["MS_Demye", "MS_Mye"]:
    for ctl in ["CTRL_CX", "CTRL_WM"]:
        g1 = meta[meta["group"] == grp].groupby("patient")["MS_norm_score"].median()
        g0 = meta[meta["group"] == ctl].groupby("patient")["MS_norm_score"].median()
        u, p = mannwhitneyu(g1, g0, alternative="two-sided")
        cd = cliffs(g1.values, g0.values)
        res_groups.append({"test": "%s_vs_%s" % (grp, ctl), "n1": len(g1), "n0": len(g0),
                           "median1": round(g1.median(), 3), "median0": round(g0.median(), 3),
                           "P": p, "cliffs_d": round(cd, 3)})
        log("  %s vs %s: n=%d/%d  P=%.3g  Cliff delta=%.3f  (median %.3f vs %.3f)"
            % (grp, ctl, len(g1), len(g0), p, cd, g1.median(), g0.median()))
# within MS: Demye vs Mye
g1 = meta[meta["group"] == "MS_Demye"].groupby("patient")["MS_norm_score"].median()
g0 = meta[meta["group"] == "MS_Mye"].groupby("patient")["MS_norm_score"].median()
u, p = mannwhitneyu(g1, g0, alternative="two-sided")
res_groups.append({"test": "MS_Demye_vs_MS_Mye", "n1": len(g1), "n0": len(g0),
                   "median1": round(g1.median(), 3), "median0": round(g0.median(), 3),
                   "P": p, "cliffs_d": round(cliffs(g1.values, g0.values), 3)})
log("  MS_Demye vs MS_Mye: n=%d/%d  P=%.3g" % (len(g1), len(g0), p))

# specificity: 5 unrelated traits (same 74-set file) + myelin markers + random sets
comp_traits = ["PASS_Alzheimers_Jansen2019", "PASS_Schizophrenia_Pardinas2018",
               "PASS_Rheumatoid_Arthritis", "PASS_Type_2_Diabetes", "PASS_Intelligence_SavageJansen2018"]
myelin = ["MBP","PLP1","MOG","MOBP","MAG","CNP","CLDN11","MAL","OPALIN","BCAS1",
          "MYRF","SOX10","OLIG1","OLIG2","NKX2-2","TMPRSS5","FA2H","CD9","GSN","QKI",
          "ST18","ERMN","PLLP","GJC2","FOLH1","CSPG4","PDGFRA","GPR17","TCF7L2"]
sp = []
for nm, gl in [("MS", ms), ("Myelin_markers", myelin)] + [(t, gs[t]) for t in comp_traits]:
    s = score(gl)
    meta["_sc"] = s
    mp_ = meta[OL_mature].groupby("patient")["_sc"].median()
    pp_ = meta[OL_prec].groupby("patient")["_sc"].median()
    dd = (mp_ - pp_).dropna()
    wp_, pv = wilcoxon(dd)
    sp.append({"trait": nm, "n_genes": len([g for g in gl if g in set(genes)]),
               "cliffs_d": round(cliffs(s[OL_mature.values], s[OL_prec.values]), 3),
               "P_wilcoxon": round(pv, 6), "n_patients": int(dd.notna().sum()),
               "median_diff": round(dd.median(), 3),
               "all_same_sign": bool((np.sign(dd) == np.sign(dd.median())).all())})
    log("  specificity: %-26s delta=%+.3f P=%g median difference %+.3f" % (nm, sp[-1]["cliffs_d"], pv, dd.median()))

# random-set null (size = the number of MS hits)
rng = np.random.default_rng(11)
geneset = np.array(genes)
gl_sz = len([g for g in ms if g in set(genes)])
rnd = []
for i in range(50):
    gl = list(rng.choice(geneset, gl_sz, replace=False))
    s = score(gl, n_ctrl=100, seed=i)
    meta["_sc"] = s
    mp_ = meta[OL_mature].groupby("patient")["_sc"].median()
    pp_ = meta[OL_prec].groupby("patient")["_sc"].median()
    rnd.append((mp_ - pp_).median())
rnd = np.array(rnd)
ms_diff = sp[0]["median_diff"]
log("\nrandom sets (n=50): median %.3f p90 %.3f max %.3f -> MS=%.3f, percentile %.1f%%"
    % (np.median(rnd), np.percentile(rnd, 90), rnd.max(), ms_diff, 100 * (rnd < ms_diff).mean()))

# ---- save ----
meta.to_csv(os.path.join(OUT, "tables", "OL_meta_with_MS_score.csv"))
pd.DataFrame({"patient": common, "mature_minus_precursor": d}).to_csv(
    os.path.join(OUT, "tables", "patient_mature_minus_precursor.csv"), index=False)
pd.DataFrame(res_groups).to_csv(os.path.join(OUT, "tables", "group_comparisons.csv"), index=False)
pd.DataFrame(sp).to_csv(os.path.join(OUT, "tables", "specificity_controls.csv"), index=False)
pd.DataFrame({"random_set_median_diff": rnd}).to_csv(
    os.path.join(OUT, "tables", "random_sets_diffs.csv"), index=False)
with open(os.path.join(OUT, "provenance.md"), "w") as fh:
    fh.write("""# A2 GSE277430, human main analysis v1

- time: %s
- rds: %s (909,468,360 B) + metadata 9,306,745 B
- OL subset: %d genes x %d nuclei (expertAnno.l1 in OLIGO/OPC; R extraction script A2_GSE277430_extract.R)
- gene set: PASS_Multiple_sclerosis (%d/%d genes hit), scDRS 1.0.2 n_ctrl=200
- primary criterion: mature (OLIGO) vs precursor (OPC), donor-level paired Wilcoxon P=%.3g, median difference %+.3f, n=%d donors
- A2 groups: unpaired donor-level MWU (Demye/Mye vs CTRL_CX/CTRL_WM, Demye vs Mye)
- specificity: 5 unrelated traits + myelin markers + 50 random sets (percentile %.1f%%)
- seeds: 20261029 / 11
""" % (time.strftime("%Y-%m-%d %H:%M:%S"), RDS_PATH, len(genes), len(barcodes),
       len([g for g in ms if g in set(genes)]), len(ms), p_main, np.median(d), len(common),
       100 * (rnd < ms_diff).mean()))
log("saved in %.0f s" % (time.time() - t0))
