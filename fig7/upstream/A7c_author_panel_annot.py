# -*- coding: utf-8 -*-
"""A7c: re-call GSE231585 (tier B) with the published author marker panel and the authors own 8-class
semantics, then compare d1 with the old call.

Author facts (paper plus github.com/AnneHartebrodt/ms-omics, both checked): the snRNA labels were
assigned by hand per cluster over 8 major classes using marker genes; the snATAC labels were
transferred from the RNA labels with FindTransferAnchors; the authors never released per-cell labels
(GEO holds raw h5 only and the three Zenodo 17733765 tar.gz files are Cell Ranger output). So this
step re-calls our own clusters with the authors panel and the authors class semantics; it does not
restore the authors labels.

Two traps, both found by measurement and fixed here:
  1) d1 must use the scDRS norm_score (the delivered basis), not the raw_score column of
     scored.h5ad (recorded in the A7b docstring). MS is re-scored here as a single trait and the
     per-nucleus scores are cached as nuc_MS_normscore.csv, so later changes of grouping take seconds.
  2) calling rule: the author panel markers are z-scored BETWEEN clusters, a class score is the mean
     of its member z, then argmax plus a margin. An absolute threshold cannot be set high: 13 of the
     19 clusters here are OL, so the OL markers own between-cluster mean is inflated by themselves
     and the real OL clusters only reach z = 0.17-1.39. The main rule R1 is therefore argmax with
     z > 0 and margin >= 0.3; the stricter R2 (z >= 0.3) is run as a sensitivity analysis (result: no
     sign flip in d1, and the pooled four-cohort statistic is completely unchanged).
"""
import os, gc, time, resource
import numpy as np
import pandas as pd
import h5py
import scipy.sparse as sp
from scipy import stats

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

BASE = WORK + "/figures/routeA_out/A7_GSE231585_rna_v1"
TAB, MD, DERV = BASE + "/tables", BASE + "/md", BASE + "/derived"
os.makedirs(MD, exist_ok=True)
H5 = DERV + "/scored.h5ad"
GS = WORK + "/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs"
TRAIT = "PASS_Multiple_sclerosis"
N_CTRL, SEED = 200, 20261030
CACHE = TAB + "/nuc_MS_normscore.csv"

PANEL = {                                            # the authors 8 classes with the panel as published
    "Oligodendrocytes":  ["PLP1", "MAG", "MOG", "OPALIN", "MBP"],
    "OPCs":              ["PDGFRA", "BCAN"],
    "Astrocytes":        ["GFAP", "AQP4", "ALDH1L1", "GJA1"],
    "Immune cells":      ["CD68", "ITGAM", "AIF1", "TMEM119", "P2RY12"],
    "Neurons":           ["RBFOX3"],
    "B-cells":           ["IGHG1", "MS4A1", "CD79A", "PAX5", "CD19"],
    "Endothelial cells": ["ERG", "VWF", "EMCN", "CLDN5", "EGFL7"],
    "Pericytes":         ["PDGFRB", "NOTCH3", "TBX18"],
}
RULES = [("R1_argmax", 0.0, 0.30), ("R2_strict", 0.30, 0.30)]     # (name, z floor, floor on the gap to the runner-up)
OLD2NEW = {"OLIGO": "Oligodendrocytes", "OPC": "OPCs", "AST": "Astrocytes", "MG": "Immune cells",
           "NEU": "Neurons", "END": "Endothelial cells", "LYM": "B-cells", "Unclear": "Unclear"}
t0 = time.time()
log_lines = []
def log(*a):
    s = "[%6ds] %s" % (time.time() - t0, " ".join(str(x) for x in a))
    print(s, flush=True); log_lines.append(s)
def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6

log("A7c: author-panel re-call + d1 comparison (delivered norm_score basis, two-rule sensitivity)")

# ---------- [1] obs + the 30 marker columns (read in chunks, never transposed) ----------
log("[1] reading the h5ad: obs + the 30 markers of the author panel (in chunks)...")
with h5py.File(H5, "r") as f:
    vn = np.asarray([x.decode() if isinstance(x, bytes) else x for x in f["var/_index"][:]])
    want = [g for gl in PANEL.values() for g in gl]
    pos = {g: int(np.where(vn == g)[0][0]) for g in want if g in set(vn)}
    gi = sorted(pos.values()); gidx = {g: gi.index(p) for g, p in pos.items()}
    indptr, data, indices = f["X/indptr"][:], f["X/data"][:], f["X/indices"][:]
    from anndata.experimental import read_elem
    obsd = {}
    for k in ["sample", "barcode", "cl", "celltype"]:
        v = read_elem(f["obs"][k]); v = v.cat if hasattr(v, "cat") else v
        obsd[k] = np.asarray(v.categories[v.codes]) if hasattr(v, "categories") else np.asarray(v)
    n = len(indptr) - 1
    M = np.zeros((n, len(gi)), dtype=np.float32)
    for i0 in range(0, n, 2048):
        i1 = min(i0 + 2048, n)
        blk = sp.csr_matrix((data[indptr[i0]:indptr[i1]], indices[indptr[i0]:indptr[i1]],
                             indptr[i0:i1 + 1] - indptr[i0]), shape=(i1 - i0, len(vn)))
        d = blk.toarray(); M[i0:i1] = d[:, gi]; del d, blk
    log("  panel hits %d/%d, missing %s; peak RSS %.2f GB" % (len(pos), len(set(want)), sorted(set(want) - set(pos)), rss()))
del data, indices, indptr; gc.collect()
obs = pd.DataFrame({"samp": obsd["sample"], "barcode": obsd["barcode"], "cl": obsd["cl"].astype(int),
                    "celltype_old": obsd["celltype"]})
cl = obs["cl"].values
samp, donr = obs["samp"].values, np.array([s.split("_")[0] for s in obs["samp"].values])

# ---------- [2] between-cluster z -> class score ----------
log("[2] z-scoring the author panel markers between clusters; class score = mean of the member z...")
ncl = int(cl.max()) + 1
ME = np.zeros((ncl, len(gi)))
order = np.argsort(cl, kind="stable"); cs = cl[order]
bnd = np.searchsorted(cs, np.arange(ncl + 1))
for k in range(ncl):
    idx = order[bnd[k]:bnd[k + 1]]
    if len(idx):
        ME[k] = M[idx].mean(0)
Zc = (ME - ME.mean(0, keepdims=True)) / (ME.std(0, keepdims=True) + 1e-9)
cls_names = list(PANEL.keys())
CS = pd.DataFrame([[float(np.mean([Zc[k, gidx[g]] for g in gl if g in gidx])) for gl in PANEL.values()]
                   for k in range(ncl)], index=range(ncl), columns=cls_names)

def build_table(z_min, marg):
    rows = []
    for k in range(ncl):
        sel = cl == k
        if not sel.sum():
            continue
        s = CS.loc[k].sort_values(ascending=False)
        call = s.index[0] if (s.iloc[0] > z_min and (len(s) < 2 or s.iloc[0] - s.iloc[1] >= marg)) else "Unclear"
        rows.append({"cl": int(k), "n": int(sel.sum()), "old_call": obs.loc[sel, "celltype_old"].iloc[0],
                     "new_call": call, "class_best": s.index[0], "z_best": round(float(s.iloc[0]), 3),
                     "runner_up": "%s(%.3f)" % (s.index[1], s.iloc[1]),
                     "margin": round(float(s.iloc[0] - s.iloc[1]), 3),
                     "same_as_old": "same" if call == OLD2NEW.get(obs.loc[sel, "celltype_old"].iloc[0], "?") else "DIFF",
                     **{c: round(float(v), 3) for c, v in s.items()}})
    return pd.DataFrame(rows).sort_values("cl").reset_index(drop=True)

# ---------- [3] MS norm_score (delivered basis) ----------
log("[3] scDRS %s per-nucleus norm_score (n_ctrl=%d, seed=%d)..." % (TRAIT, N_CTRL, SEED))
CACHE_OK = os.path.exists(CACHE) and len(pd.read_csv(CACHE)) == n
if CACHE_OK:
    sc = pd.read_csv(CACHE)["norm_score"].values
    log("  using the cache %s (%d rows)" % (CACHE, len(sc)))
else:
    if os.path.exists(CACHE):
        log("  the cache has the wrong number of rows; deleting and re-scoring..."); os.remove(CACHE)
    import anndata as ad, scdrs
    adata = ad.read_h5ad(H5); scdrs.pp.preprocess(adata)
    sets = {}
    for l in open(GS):
        f2 = l.rstrip("\n").split("\t")
        if len(f2) < 2 or f2[0] == "TRAIT": continue
        sets[f2[0]] = f2[1].split(",")
    use = [g for g in sets[TRAIT] if g in set(adata.var_names)]
    log("  gene-set hits %d genes, %d nuclei..." % (len(use), adata.n_obs))
    res = scdrs.method.score_cell(adata, use, n_ctrl=N_CTRL, random_seed=SEED)
    col = "norm_score" if "norm_score" in res.columns else res.columns[0]
    sc = res[col].to_numpy()
    pd.DataFrame({"samp": obs["samp"], "barcode": obs["barcode"], "cl": cl,
                  "celltype_old": obs["celltype_old"], "norm_score": sc}).to_csv(CACHE, index=False)
    del adata, res; gc.collect()
    log("  done in %.0fs; peak RSS %.2f GB" % (time.time() - t0, rss()))
s_ = sc.astype(float)
obs["norm_score_MS"] = s_

def delta1(mat, pre, keys, minn=10):
    out = {}
    for k in pd.unique(keys):
        m = keys == k
        a, b = s_[m & mat], s_[m & pre]
        if a.size >= minn and b.size >= minn:
            out[k] = (float(np.median(a) - np.median(b)), int(a.size), int(b.size))
    return out

mat_o = (obs["celltype_old"] == "OLIGO").values
pre_o = (obs["celltype_old"] == "OPC").values
So, Do = delta1(mat_o, pre_o, samp), delta1(mat_o, pre_o, donr)
W = pd.read_csv(TAB + "/sample_delta1_by_trait.csv", index_col=0).loc[TRAIT]
dev = max(abs(So[s][0] - float(W[s])) for s in So if s in W.index)
ok = "reproduced (the difference is only float rounding in the cached CSV)" if dev < 1e-5 else "FAIL: there is a discrepancy, investigate"
log("  * the old grouping reproduces the delivered table: max deviation %.2e (%d/14 samples) -> %s" % (dev, len(So), ok))

B = WORK + "/figures/routeA_out"
v3 = []
for f_ in ["A3_GSE118257_v1/tables/patient_pair_scores.csv", "A2_GSE277430_v1/tables/patient_pair_scores.csv",
           "A4_GSE301908_v1/tables/patient_pair_scores.csv"]:
    v3 += list(pd.read_csv(B + "/" + f_)["diff"].values)
log("[4] baseline with the old call: A7 %d donors all %s | four cohorts n=%d, same sign=%d, P=%.3g"
    % (len(Do), "positive" if all(x[0] > 0 for x in Do.values()) else "mixed",
       len(v3) + len(Do), sum(1 for x in v3 + [y[0] for y in Do.values()] if x > 0),
       stats.binomtest(sum(1 for x in v3 + [y[0] for y in Do.values()] if x > 0), len(v3) + len(Do), 0.5,
                       alternative="greater").pvalue))

# ---------- [5] call with each rule and compute d1 ----------
summary = []
for name, z_min, marg in RULES:
    log("\n===== rule %s (z>%.2f and margin>=%.2f) =====" % (name, z_min, marg))
    T = build_table(z_min, marg)
    T.to_csv(TAB + "/cluster_annotation_authorpanel_%s.csv" % name, index=False)
    if name == RULES[0][0]:
        T.to_csv(TAB + "/cluster_annotation_authorpanel.csv", index=False)
    log(T[["cl", "n", "old_call", "new_call", "class_best", "z_best", "runner_up", "margin", "same_as_old"]].to_string(index=False))
    nd = T[T.same_as_old == "DIFF"]
    log("  differs from the old call: %d clusters / %d nuclei -> %s" % (len(nd), nd.n.sum(), nd[["cl", "old_call", "new_call"]].to_dict("records")))
    lab = pd.Series(T.set_index("cl")["new_call"]).reindex(cl).values
    obs["celltype_authorpanel_" + name] = lab
    mat_n, pre_n = (lab == "Oligodendrocytes"), (lab == "OPCs")
    log("  nuclei: Oligodendrocytes=%d OPCs=%d | old OLIGO=%d OPC=%d" % (mat_n.sum(), pre_n.sum(), mat_o.sum(), pre_o.sum()))
    Sn, Dn = delta1(mat_n, pre_n, samp), delta1(mat_n, pre_n, donr)
    tb = pd.DataFrame({"samp": sorted(So)})
    tb["d1_old"] = [So[s][0] for s in tb.samp]; tb["n_mat_old"] = [So[s][1] for s in tb.samp]; tb["n_pre_old"] = [So[s][2] for s in tb.samp]
    tb["d1_" + name] = [Sn.get(s, (np.nan, 0, 0))[0] for s in tb.samp]
    tb["n_mat_" + name] = [Sn.get(s, (0, 0, 0))[1] for s in tb.samp]
    tb["n_pre_" + name] = [Sn.get(s, (0, 0, 0))[2] for s in tb.samp]
    tb["delta"] = tb["d1_" + name] - tb.d1_old
    tb["flip"] = np.where(np.sign(tb.d1_old) != np.sign(tb["d1_" + name]), "FLIP", "")
    tb.round(4).to_csv(TAB + "/delta1_MS_authorpanel_samples_%s.csv" % name, index=False)
    log(tb.round(3).to_string(index=False))
    dd = pd.DataFrame({"donor": sorted(Do)})
    dd["d1_old"] = [Do[d][0] for d in dd.donor]
    dd["d1_" + name] = [Dn.get(d, (np.nan,))[0] for d in dd.donor]
    dd["delta"] = dd["d1_" + name] - dd.d1_old
    dd["flip"] = np.where(np.sign(dd.d1_old) != np.sign(dd["d1_" + name]), "FLIP", "")
    dd.round(4).to_csv(TAB + "/delta1_MS_authorpanel_donors_%s.csv" % name, index=False)
    log("  donor-level d1:\n" + dd.round(3).to_string(index=False))
    vv = np.array(v3 + [x[0] for x in Dn.values()])
    nn, kk = len(vv), int((vv > 0).sum())
    p = stats.binomtest(kk, nn, 0.5, alternative="greater").pvalue
    log("  four cohorts pooled at donor level: n=%d, same sign=%d, opposite=%d, one-sided binomial P=%.3g (delivered value 49/48, P=8.88e-14)" % (nn, kk, nn - kk, p))
    summary.append({"rule": name, "n_OL_nuc": int(mat_n.sum()), "n_OPC_nuc": int(pre_n.sum()),
                    "clusters_diff": len(nd), "donor_flip": int((dd.flip == "FLIP").sum()),
                    "max_abs_delta_donor": float(np.nanmax(np.abs(dd.delta.values))),
                    "cross_cohort_n": nn, "cross_cohort_same_sign": kk, "cross_cohort_P": p})
obs.to_csv(TAB + "/nuc_annotation_authorpanel.csv", index=False)
SM = pd.DataFrame(summary)
SM.to_csv(TAB + "/authorpanel_summary.csv", index=False)
log("\n[5] rule comparison summary:\n" + SM.to_string(index=False, formatters={"cross_cohort_P": lambda x: "%.3g" % x}))
open(MD + "/A7c_authorpanel.md", "w").write("\n".join(log_lines) + "\n")
log("[6] done; peak RSS %.2f GB" % rss())
