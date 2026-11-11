# -*- coding: utf-8 -*-
"""M1 with the author annotation (item 1 of the outstanding list in the review receipt):
the per-cell scDRS scores are the 6 obs columns M1 already computed (nothing is re-scored); only the
grouping labels come from cellType in the authors annotation.csv.gz: mature = MOL1+MOL2+MOL56, precursor = OPC+COP;
the disease / phase labels also come from the authors orig.ident_merge. Library -> author mouse is mapped by barcode containment.
Criteria identical to M1: >= 10 nuclei on both sides of a library, d1 = median(mature) - median(precursor), binomial same-sign test, EAE vs Ctrl MWU.
Output: $DSCI_WORK_ROOT/figures/routeA_out/M1b_authoranno_v1/tables/
"""
import os, io, time, gzip, json
import numpy as np
import pandas as pd
import h5py
from scipy.stats import binomtest, mannwhitneyu

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

ANNO = WORK + "/routeA_probe/mouse_EAE/annotation.csv.gz"
CACHE = WORK + "/routeA_probe/mouse_EAE/anno_min.pkl"
H5 = WORK + "/figures/routeA_out/M1_mouse_EAE_rna_v1/derived/annotated.h5ad"
OUT = WORK + "/figures/routeA_out/M1b_authoranno_v1"
os.makedirs(OUT + "/tables", exist_ok=True)
t0 = time.time()
TRAITS = ["PASS_Multiple_sclerosis", "PASS_Schizophrenia_Pardinas2018", "PASS_Intelligence_SavageJansen2018",
          "PASS_SWB", "UKB_460K.blood_PLATELET_COUNT", "UKB_460K.body_HEIGHTz"]
MAT = ["MOL1", "MOL2", "MOL56"]
PRE = ["OPC", "COP"]


def log(*a):
    print(*a, flush=True)


# ---------- 1) author annotation (cached persistently) ----------
if os.path.exists(CACHE):
    A = pd.read_pickle(CACHE)
    log("author annotation read from cache: %d rows" % len(A))
else:
    A = pd.read_csv(ANNO, usecols=lambda c: c in ["Unnamed: 0", "orig.ident", "orig.ident_merge", "model",
                                                  "cellType", "subcellType", "batch", "Sex"])
    A["cell"] = A["Unnamed: 0"].astype(str)
    A["bc16"] = A["cell"].str.slice(0, 16)
    A["suffix"] = A["cell"].str.extract(r"-(\d+)$")[0]
    A.to_pickle(CACHE)
    log("author annotation parsed and cached: %d rows, %.0fs" % (len(A), time.time() - t0))
log("cellType distribution: %s" % A["cellType"].value_counts().head(9).to_dict())
log("orig.ident_merge: %s" % A["orig.ident_merge"].value_counts().to_dict())

# ---------- 2) per-cell scores of M1 + our own labels (obs columns only) ----------
with h5py.File(H5, "r") as f:
    def col(name):
        ds = f["obs"][name]
        if isinstance(ds, h5py.Group):        # categorical encoding: codes + categories
            codes = ds["codes"][:]
            cats = [c.decode() if isinstance(c, bytes) else c for c in ds["categories"][:]]
            return np.array([cats[c] if c >= 0 else None for c in codes], dtype=object)
        arr = ds[:]
        if arr.dtype.kind in "uf":
            return arr
        return np.array([x.decode() if isinstance(x, bytes) else x for x in arr], dtype=object)
    barcode = col("barcode"); lib = col("lib"); self_ct = col("celltype")
    m1_lib_phase = {}
    scols = {}
    for tr in TRAITS:
        k = "scDRS_" + tr
        if k in f["obs"]:
            scols[tr] = f["obs"][k][:]
        else:
            alt = [x for x in f["obs"] if tr in x]
            scols[tr] = f["obs"][alt[0]][:] if alt else None
    n = len(barcode)
log("M1 cells %d, score columns %d/%d, read in %.0fs" % (n, sum(v is not None for v in scols.values()), len(TRAITS), time.time() - t0))

# ---------- 3) library -> author mouse mapping (bc16 containment) ----------
M = pd.DataFrame({"barcode": barcode, "lib": lib, "self_ct": self_ct})
M["bc16"] = M["barcode"].astype(str).str.slice(0, 16)
amap = {m: set(g["bc16"]) for m, g in A.groupby("orig.ident")}
rows_map = []
for L, g in M.groupby("lib"):
    bs = set(g["bc16"])
    sc = {m: len(bs & s) / max(len(bs), 1) for m, s in amap.items()}
    best = max(sc, key=sc.get)
    rows_map.append({"lib": L, "n_cells": len(g), "author_mouse": best, "containment": sc[best],
                     "second": sorted(sc.values(), reverse=True)[1] if len(sc) > 1 else 0.0})
MAP = pd.DataFrame(rows_map)
MAP.to_csv(OUT + "/tables/lib_to_author_mouse.csv", index=False)
log(MAP.round(3).to_string(index=False))
good = MAP[(MAP.containment >= 0.7) & (MAP.second < 0.1)]  # unique mapping; the ~0.9 ceiling from the authors dropping nuclei in QC does not change the assignment
log("usable libraries: %d/%d (containment >= 0.7 and runner-up < 0.1); unique mice %d; duplicated assignments %s" % (len(good), len(MAP), good.author_mouse.nunique(), good.author_mouse.duplicated().sum()))

# ---------- 4) attach the author labels exactly through (lib -> mouse, bc16) ----------
A_idx = A.set_index(["orig.ident", "bc16"])[["cellType", "subcellType"]].groupby(level=[0, 1]).first()
key_lib2mouse = dict(zip(MAP.lib, MAP.author_mouse))
M["mouse"] = M["lib"].map(key_lib2mouse)
jidx = pd.MultiIndex.from_arrays([M["mouse"].astype(str), M["bc16"].astype(str)])
res = A_idx.reindex(jidx)
M["author_ct"] = res["cellType"].values
matched = M["author_ct"].notna().values
log("per-nucleus annotation matched: %d/%d (%.1f%%)" % (matched.sum(), n, 100 * matched.mean()))
ph = A.groupby("orig.ident")["orig.ident_merge"].first().to_dict()
M["author_phase"] = M["mouse"].map(ph)
M.to_csv(OUT + "/tables/cells_annotated.csv", index=False)

# ---------- 5) d1 per library (author grouping vs our own, on the same scores) ----------
lib_ok = set(good.lib)
per_lib, comp = [], []
for tr in TRAITS:
    s = scols[tr]
    if s is None:
        continue
    for kind, ctcol in [("author", "author_ct"), ("self", "self_ct")]:
        if kind == "author":
            ismat = M[ctcol].isin(MAT).values
            ispre = M[ctcol].isin(PRE).values
        else:
            ismat = (M[ctcol] == "OLIGO").values
            ispre = (M[ctcol] == "OPC").values
        vals = []
        for L in sorted(lib_ok):
            m = (M["lib"] == L).values
            a, b = s[m & ismat], s[m & ispre]
            if len(a) >= 10 and len(b) >= 10:
                vals.append((L, M.loc[m, "author_phase"].iloc[0], float(np.median(a) - np.median(b)), len(a), len(b)))
        for L, p, d, na, nb in vals:
            per_lib.append({"trait": tr, "anno": kind, "lib": L, "phase": p, "delta1": d, "n_mat": na, "n_pre": nb})
        arr = np.array([v[2] for v in vals])
        np_ = int((arr > 0).sum())
        comp.append({"trait": tr, "anno": kind, "n_lib": len(vals), "median_diff": float(np.median(arr)),
                     "n_pos": np_, "P_sign": binomtest(np_, len(vals), 0.5).pvalue if len(vals) else np.nan})
P = pd.DataFrame(per_lib)
P.to_csv(OUT + "/tables/per_lib_delta1_author_vs_self.csv", index=False)
C = pd.DataFrame(comp)
C.to_csv(OUT + "/tables/trait_summary_author_vs_self.csv", index=False)
log("\n=== 6 traits: author annotation vs our annotation (same per-cell scores)")
print(C.round(5).to_string(index=False))

# ---------- 6) disease axis (under the author labels) ----------
ms = P[(P.trait == "PASS_Multiple_sclerosis") & (P.anno == "author")]
eae = ms[ms.phase.isin(["Early", "Late", "Peak"])].delta1.values
ctl = ms[ms.phase == "Ctrl"].delta1.values
out = {"n_eae_libs": len(eae), "n_ctrl_libs": len(ctl), "median_eae": float(np.median(eae)),
       "median_ctrl": float(np.median(ctl)), "delta": float(np.median(eae) - np.median(ctl)),
       "MWU_P": mannwhitneyu(eae, ctl).pvalue if len(eae) and len(ctl) else np.nan,
       "n_pos_all": int((ms.delta1 > 0).sum()), "n_lib_all": len(ms)}
for tr in TRAITS:
    t = P[(P.trait == tr) & (P.anno == "author")]
    e = t[t.phase.isin(["Early", "Late", "Peak"])].delta1.values
    c = t[t.phase == "Ctrl"].delta1.values
    if len(e) and len(c):
        out["%s__delta" % tr] = float(np.median(e) - np.median(c))
        out["%s__P" % tr] = mannwhitneyu(e, c).pvalue
pd.DataFrame([out]).to_csv(OUT + "/tables/disease_axis_author.csv", index=False)
log("\n=== disease axis (author-annotation MS d1): EAE %d libraries vs Ctrl %d libraries, d=%+.3f P=%.4f; same sign %d/%d libraries"
    % (out["n_eae_libs"], out["n_ctrl_libs"], out["delta"], out["MWU_P"], out["n_pos_all"], out["n_lib_all"]))
log("all done in %.0fs" % (time.time() - t0))
