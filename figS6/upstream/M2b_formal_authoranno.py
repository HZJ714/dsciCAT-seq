# -*- coding: utf-8 -*-
"""M2b formal version: mouse GSE250589 same-nucleus multiome, all annotations switched to the AUTHORS' annotation.csv.gz.
Design (receipt §6 + the §54.3 pilot conclusion):
  1) cell-state groups = the authors' RNA annotation: MOL (MOL1+MOL2+MOL56) vs OPC (OPC+COP), per mouse (within-mouse split-half)
  2) ATAC-channel test: within-mouse split-half differential-peak panel (training half selects peaks -> testing half scores)
     positive-control bar AUC >= 0.6 or permutation P < 0.01 (QC failing -> M2 enters no figure)
  3) disease axis: per-mouse median ATAC axis (MOL-OPC differential-peak accessibility) -> donor-level MWU of EAE (Early+Peak+Late) vs Ctrl
  4) same-nucleus RNA x ATAC cross-channel: per-mouse median RNA axis (the authors' annotation on the same nuclei, OL/OPC marker-gene panel) vs the ATAC axis, Spearman + permutation
  5) also outputs the authors'-annotation RNA disease axis (the M1 control)
Exact peak-name matching across libraries is ~0 (Cell Ranger peak-boundary differences) -> the panel must be defined within each mouse; exact cross-library matching is impossible.
"""
import os, sys, re, time, glob
import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csc_matrix
from scipy.stats import mannwhitneyu, spearmanr

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

D = WORK + "/routeA_probe/mouse_EAE"
OUT = WORK + "/figures/routeA_out/M2b_authoranno_v1"
os.makedirs(OUT + "/tables", exist_ok=True)
os.makedirs(OUT + "/panels", exist_ok=True)
t0 = time.time()


def log(*a):
    print(*a, flush=True)


A = pd.read_pickle("/tmp/anno_meta.pkl")
A["bc16"] = A["bc16"].astype(str)
files = sorted(glob.glob(D + "/*_gex.h5"))
OL_MARK = ["Mbp", "Plp1", "Mag", "Mog", "Cldn11", "Aspa", "Galc", "Ugt8", "Opalin",
           "Myrf", "Ermn", "Cnp", "Slc44a1", "Fa2h", "St18"]
OP_MARK = ["Pdgfra", "Cspg4", "Abca1", "Ptprz1", "Tcf4", "Cspg5", "Six1", "Hes1",
           "Sox10", "Olig1", "Olig2", "Egfr", "Dcx", "Meox1", "Tubb2b"]

# ---- 1) suffix -> gex file (the 100%-containment rule) ----
sfx_n = {s: (A["suffix"] == s).sum() for s in A["suffix"].unique()}
suffix2file = {}
for f in files:
    with h5py.File(f, "r") as h:
        g0 = h["matrix"] if "matrix" in h else h[list(h.keys())[0]]
        bc16 = set(re.sub(r"-\d+$", "", b.decode() if isinstance(b, bytes) else b)
                   for b in g0["barcodes"][:])
    ov = {s: len(bc16 & set(A.loc[A["suffix"] == s, "bc16"])) for s in sfx_n}
    s = max(ov, key=ov.get)
    if ov[s] == sfx_n[s] and s not in suffix2file:
        suffix2file[s] = f
    elif ov[s] == sfx_n[s]:
        log("warning: suffix %s has multiple files 100%% overlapping" % s)
log("suffix->gex mapping: %d/26 unique" % len(suffix2file))

rows = []
for mouse, sub in A.groupby("orig.ident"):
    suffix = sub["suffix"].iloc[0]
    if suffix not in suffix2file:
        log("skipping %s: no gex file" % mouse)
        continue
    f = suffix2file[suffix]
    want = set(sub["bc16"])
    with h5py.File(f, "r") as h:
        g0 = h["matrix"] if "matrix" in h else h[list(h.keys())[0]]
        fg = g0["features"]
        if isinstance(fg, h5py.Group):
            tps = [x.decode() for x in fg["feature_type"][:]]
            nms = [x.decode() for x in fg["name"][:]]
        else:
            tps = [x[2].decode() for x in fg[:]]
            nms = [x[1].decode() for x in fg[:]]
        pk = np.array([t == "Peaks" for t in tps])
        ge = np.array([t == "Gene Expression" for t in tps])
        pk_names = np.array(nms)[pk]
        ge_names = np.array(nms)[ge]
        bcs = [re.sub(r"-\d+$", "", b.decode() if isinstance(b, bytes) else b)
               for b in g0["barcodes"][:]]
        col = [i for i, b in enumerate(bcs) if b in want]
        n_feat = len(nms)
        M = csc_matrix((g0["data"][:].astype(np.float32), g0["indices"][:], g0["indptr"][:]),
                       shape=(n_feat, len(bcs)))
        Mpk = M[pk, :][:, col].tocsr()          # peaks x nuclei
        Mge = M[ge, :][:, col].tocsr()          # genes x nuclei
    P = Mpk.T.tocsr()                            # nuclei x peaks
    G = Mge.T.tocsr()                            # nuclei x genes
    del M, Mpk, Mge
    bc2row = {bcs[i]: k for k, i in enumerate(col)}
    rows_idx = [bc2row[b] for b in sub["bc16"]]
    sub = sub.iloc[rows_idx].reset_index(drop=True)
    P = P[rows_idx]
    G = G[rows_idx]
    depP = np.asarray(P.sum(1)).ravel()
    depG = np.asarray(G.sum(1)).ravel()
    X = (P / np.maximum(depP, 1)[:, None] * 1e4).tocsr()      # log1p CPM peaks
    X = X.log1p()
    R = (G / np.maximum(depG, 1)[:, None] * 1e4).tocsr()      # log1p CPM genes
    R = R.log1p()
    # the RNA axis (the same nuclei)
    gi_ol = [np.where(ge_names == g)[0] for g in OL_MARK if g in ge_names]
    gi_op = [np.where(ge_names == g)[0] for g in OP_MARK if g in ge_names]
    gi_ol = [i[0] for i in gi_ol if len(i)]
    gi_op = [i[0] for i in gi_op if len(i)]
    if not gi_ol or not gi_op:
        log("skipping %s: marker genes missing" % mouse); continue
    rna_axis = np.asarray(R[:, gi_ol].mean(1)).ravel() - np.asarray(R[:, gi_op].mean(1)).ravel()
    grp = np.where(sub["cellType"].isin(["MOL1", "MOL2", "MOL56"]))[0]
    gop = np.where(sub["cellType"].isin(["OPC", "COP"]))[0]
    if len(grp) < 20 or len(gop) < 20:
        log("skipping %s: groups too small MOL%d OPC%d" % (mouse, len(grp), len(gop))); continue
    rng = np.random.default_rng(0)
    g1 = rng.permutation(grp)[: len(grp) // 2]
    o1 = rng.permutation(gop)[: len(gop) // 2]
    g2 = np.setdiff1d(grp, g1)
    o2 = np.setdiff1d(gop, o1)
    # the training half selects peaks
    tr = np.concatenate([g1, o1]); tl = np.concatenate([np.ones(len(g1)), np.zeros(len(o1))])
    Xt = X[tr]
    m1 = np.asarray(Xt[tl == 1].mean(0)).ravel()
    m0 = np.asarray(Xt[tl == 0].mean(0)).ravel()
    v1 = np.asarray(Xt[tl == 1].power(2).mean(0)).ravel() - m1 ** 2
    v0 = np.asarray(Xt[tl == 0].power(2).mean(0)).ravel() - m0 ** 2
    n1, n0 = (tl == 1).sum(), (tl == 0).sum()
    z = (m1 - m0) / np.sqrt(v1 / max(n1, 1) + v0 / max(n0, 1) + 1e-12)
    z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
    topM = np.argsort(-z)[:500]
    topO = np.argsort(z)[:500]
    # the testing half scores AUC
    wM = np.zeros(X.shape[1]); wM[topM] = 1.0 / 500
    wO = np.zeros(X.shape[1]); wO[topO] = 1.0 / 500
    score_all = np.asarray(X.dot(wM)).ravel() - np.asarray(X.dot(wO)).ravel()
    # within-lineage contrast convention (same as M1/Fig 7): mature-group median minus precursor-group median; the all-nucleus median is composition reference only
    atac_contrast = np.median(score_all[grp]) - np.median(score_all[gop])
    rna_contrast = np.median(rna_axis[grp]) - np.median(rna_axis[gop])
    atac_allcell = np.median(score_all)
    rna_allcell = np.median(rna_axis)
    a, b = score_all[g2], score_all[o2]
    u = mannwhitneyu(a, b, alternative="two-sided")
    auc = u.statistic / (len(a) * len(b))
    # permutation P (subsampling 500+500, 100 draws)
    sa = rng.choice(a, size=min(len(a), 500), replace=False)
    sb = rng.choice(b, size=min(len(b), 500), replace=False)
    obs = (np.median(sa) - np.median(sb))
    cnt_perm = 0
    for _ in range(100):
        s = np.concatenate([sa, sb]); rng.shuffle(s)
        cnt_perm += (np.median(s[:len(sa)]) - np.median(s[len(sa):])) >= obs
    p_perm = (cnt_perm + 1) / 101
    qc = (auc >= 0.6) or (p_perm < 0.01)
    rows.append(dict(mouse=mouse, group=sub["orig.ident_merge"].iloc[0],
                     n_MOL=len(grp), n_OPC=len(gop),
                     auc_test=auc, p_perm=p_perm, qc_pass=qc,
                     atac_axis_median=atac_allcell, rna_axis_median=rna_allcell,
                     atac_axis_contrast=atac_contrast, rna_axis_contrast=rna_contrast,
                     n_panel=1000))
    with open(OUT + "/panels/%s_topMOL.txt" % mouse, "w") as fo:
        fo.write("\n".join(pk_names[topM]))
    with open(OUT + "/panels/%s_topOPC.txt" % mouse, "w") as fo:
        fo.write("\n".join(pk_names[topO]))
    log("%-14s MOL%5d OPC%4d AUC=%.3f Pperm=%.3f QC=%s" %
        (mouse, len(grp), len(gop), auc, p_perm, qc))
    del P, X, G, R

df = pd.DataFrame(rows)
df.to_csv(OUT + "/tables/M2b_per_mouse.csv", index=False)
log("\nQC passed: %d/%d mice" % (df["qc_pass"].sum(), len(df)))

# ---- the disease axis (donor level, contrast convention) ----
for ch, col in [("ATAC", "atac_axis_contrast"), ("RNA", "rna_axis_contrast")]:
    d = df[df["qc_pass"]].copy() if ch == "ATAC" else df.copy()
    eae = d[d["group"] != "Ctrl"][col].values
    ctrl = d[d["group"] == "Ctrl"][col].values
    u = mannwhitneyu(eae, ctrl, alternative="two-sided")
    log("[%s axis] EAE n=%d vs Ctrl n=%d: Δ=%.4f (EAE-Ctrl), U P=%.3g" %
        (ch, len(eae), len(ctrl), np.median(eae) - np.median(ctrl), u.pvalue))

# ---- same-nucleus cross-channel (contrast convention) ----
d = df[df["qc_pass"]]
rho, pp = spearmanr(d["atac_axis_contrast"], d["rna_axis_contrast"])
rng = np.random.default_rng(1)
null = []
for _ in range(10000):
    y = rng.permutation(d["rna_axis_contrast"].values)
    null.append(spearmanr(d["atac_axis_contrast"].values, y).statistic)
p_perm = (np.abs(np.array(null)) >= abs(rho)).mean()
log("\nsame-nucleus cross-channel (per-mouse contrast): rho=%.3f perm P=%.4f (n=%d)" % (rho, p_perm, len(d)))
pd.DataFrame({"channel": ["ATAC", "RNA"],
              "delta_EAE_minus_Ctrl": [np.median(df[df["qc_pass"] & (df.group != "Ctrl")]["atac_axis_contrast"]) -
                                       np.median(df[df["qc_pass"] & (df.group == "Ctrl")]["atac_axis_contrast"]),
                                       np.median(df[df.group != "Ctrl"]["rna_axis_contrast"]) -
                                       np.median(df[df.group == "Ctrl"]["rna_axis_contrast"])]},
             index=[0, 1]).to_csv(OUT + "/tables/M2b_disease.csv", index=False)
pd.DataFrame({"n_mice": len(d), "spearman_rho": rho, "perm_P": p_perm}, index=[0]).to_csv(
    OUT + "/tables/M2b_crosschannel.csv", index=False)
log("done in %.0fs -> %s" % (time.time() - t0, OUT))
