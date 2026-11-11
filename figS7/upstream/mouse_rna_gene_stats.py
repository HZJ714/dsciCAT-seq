#!/usr/bin/env python3
"""mouse_rna_gene_stats.py — whole-dataset GENE_STATS for the mouse RNA channel (scdrs.method.score_cell, n_ctrl = 200), 28 libraries.

M7_gene_stats.py - mouse RNA channel full-data GENE_STATS (replicates the M5 formula, 28 libraries)

Background: M7_rna_mouse_threeaxis needs the GENE_STATS used by scdrs.method.score_cell(n_ctrl=200).
  M6 (the ATAC arm) reuses M5's S7_atac_fullstats_GENE_STATS.csv (ATAC gene activity, 28 libraries).
  M1 (the delivered RNA pipeline) used scdrs.pp.preprocess full-data statistics. M7 must use full-data
  statistics on the RNA channel too (28 libraries x all nuclei of the Gene Expression matrix), or it

Basis (structurally identical to M5):
  RNA rows = Cell Ranger feature_type=="Gene Expression";
  X = CP10K + log1p (each nucleus rescaled to a count sum of 1e4, then log1p); ct = expm1(X);
  moments accumulated per gene: sum(X), sum(X^2), sum(ct), sum(ct^2), n (all nuclei, 28 libraries);
  GS columns = mean / var / var_tech / ct_mean / ct_var / ct_var_tech / mean_var
    (var_tech from a loess(span=0.3, degree=2) fit of log10(ct_var)~log10(ct_mean), mean_var from
     pd.qcut 20x20 bins, as in scdrs.pp.compute_stats - the verified M5 block is copied directly).
  Gene space = the M5 UNION (index of S7_atac_fullstats_GENE_STATS.csv, all official RNA symbols),
  identical to M6, so the ATAC and RNA arms are compared on one gene set.

Output: figures/routeA_out/SUPPORT_v1/tables/S7_rna_fullstats_GENE_STATS.csv
      + a P2 check against Late5 all-nuclei scdrs.pp.preprocess, column by column."""
import os, re, sys, time, glob, json
import numpy as np, pandas as pd
import h5py
from scipy import sparse as sp
import anndata as ad
import scdrs
from skmisc.loess import loess

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

t0 = time.time()
def log(*a): print("[%4.0fs]" % (time.time() - t0), *a, flush=True)

MOUSE_EAE = WORK + "/routeA_probe/mouse_EAE"
SUP = WORK + "/figures/routeA_out/SUPPORT_v1/tables"
OUT = SUP + "/S7_rna_fullstats_GENE_STATS.csv"
GENE_STATS_ATAC = SUP + "/S7_atac_fullstats_GENE_STATS.csv"
VERIFY_LIB = "Late5" if os.path.exists(MOUSE_EAE + "/Late5_gex.h5") else None

# ---------- gene space: reuse the M5 UNION (same set as the ATAC arm, so both arms stay comparable) ----------
GA = pd.read_csv(GENE_STATS_ATAC, index_col=0)
UNION = sorted(GA.index)
log("reusing the M5 UNION: %d genes (same gene set for ATAC/RNA)" % len(UNION))

# ---------- P1: streaming moment accumulation over all nuclei of the 28 RNA libraries ----------
files = sorted(glob.glob(MOUSE_EAE + "/*_gex.h5"))
log("28 library H5: %d found, starting the RNA streaming accumulation (union gene space %d)" % (len(files), len(UNION)))
MOM = {g: np.zeros(5, dtype=np.float64) for g in UNION}
lib_n = {}
for f in files:
    with h5py.File(f, "r") as h:
        g = h["matrix"]
        nf, nb = int(g["shape"][0]), int(g["shape"][1])
        gt = np.array([x.decode() if isinstance(x, bytes) else x for x in g["features"]["feature_type"][:]])
        nm = np.array([x.decode() if isinstance(x, bytes) else x for x in g["features"]["name"][:]])
        full = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=(nf, nb))
        gi = np.where(gt == "Gene Expression")[0]
        X = full[gi].T.tocsr().astype(np.float32)             # nuclei x genes (all nuclei)
        names_g = nm[gi]
    # row indices of the UNION official symbols in the RNA matrix (a max over duplicates is not needed:
    want = np.isin(names_g, UNION)
    idx = np.where(want)[0]
    syms = names_g[idx]
    Xg = X[:, idx]                                            # nuclei x UNION intersect present
    s = Xg.sum(1).astype(np.float32)
    from scipy.sparse import diags as _diags
    Gn = _diags((1e4 / np.maximum(np.asarray(s), 1)).ravel()).dot(Xg)               # CP10K (sparse diagonal scaling, avoids np.matrix overflow)
    Gn = np.log1p(Gn.toarray()).astype(np.float32)            # dense, for the moment accumulation
    ct = np.expm1(Gn)
    pos = {s_: i for i, s_ in enumerate(syms)}
    un = np.array([i for i, g_ in enumerate(UNION) if g_ in pos], dtype=np.int64)
    col = np.array([pos[UNION[u]] for u in un], dtype=np.int64)
    if len(un):
        MOMv = np.stack([Gn[:, col].sum(0), (Gn[:, col] ** 2).sum(0),
                         ct[:, col].sum(0), (ct[:, col] ** 2).sum(0),
                         np.full(len(col), Gn.shape[0])], axis=1)
        for u, mv in zip(un, MOMv):
            MOM[UNION[u]] += mv
    lib_n[os.path.basename(f)] = Gn.shape[0]
    log("  %-24s nuclei %6d | union genes %5d | moments accumulated (running nuclei %d)" % (
        os.path.basename(f), Gn.shape[0], len(col), sum(lib_n.values())))
    del Xg, Gn, ct
N_TOT = sum(lib_n.values())
log("P1 done: RNA moments over %d libraries and %d nuclei" % (len(lib_n), N_TOT))

# ---------- build GENE_STATS (exact replica of scdrs.pp.compute_stats / M5) ----------
GS = pd.DataFrame(index=UNION, columns=["mean", "var", "var_tech", "ct_mean", "ct_var", "ct_var_tech", "mean_var"])
MM = np.array([MOM[g] for g in UNION])
n, sm, sq, ctm, ctq = MM[:, 4], MM[:, 0], MM[:, 1], MM[:, 2], MM[:, 3]
GS["mean"] = sm / n
GS["var"] = sq / n - GS["mean"] ** 2
GS["ct_mean"] = ctm / n
GS["ct_var"] = ctq / n - GS["ct_mean"] ** 2
not_const = GS["ct_var"].values > 0
estimat_var = np.zeros(len(UNION), dtype=np.float64)
y = np.log10(GS["ct_var"].values[not_const]); x = np.log10(GS["ct_mean"].values[not_const])
model = loess(x, y, span=0.3, degree=2); model.fit()
estimat_var[not_const] = model.outputs.fitted_values
GS["ct_var_tech"] = 10 ** estimat_var
GS["var_tech"] = GS["var"] * GS["ct_var_tech"] / GS["ct_var"]
GS.loc[GS["var_tech"].isna(), "var_tech"] = 0
n_bin_max = int(np.floor(np.sqrt(len(UNION) / 10)))
nb = min(20, n_bin_max)
v_mean_bin = pd.qcut(GS["mean"], nb, labels=False, duplicates="drop")
GS["mean_var"] = ""
for bin_ in sorted(set(v_mean_bin)):
    ind = v_mean_bin == bin_
    v_var_bin = pd.qcut(GS.loc[ind, "var"], nb, labels=False, duplicates="drop")
    GS.loc[ind, "mean_var"] = ["%d.%d" % (bin_, x) for x in v_var_bin]
log("GENE_STATS built (%d genes, %d x %d bins), non-empty mean_var %d" % (
    len(GS), nb, nb, (GS.mean_var != "").sum()))
GS.to_csv(OUT)
log("wrote %s" % OUT)

# ---------- P2: column-by-column check against Late5 all-nuclei scdrs.pp.preprocess ----------
if VERIFY_LIB:
    f = sorted(glob.glob(MOUSE_EAE + "/%s_gex.h5" % VERIFY_LIB))
    if f:
        f = f[0]
        with h5py.File(f, "r") as h:
            g = h["matrix"]
            gt = np.array([x.decode() if isinstance(x, bytes) else x for x in g["features"]["feature_type"][:]])
            nm = np.array([x.decode() if isinstance(x, bytes) else x for x in g["features"]["name"][:]])
            full = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]),
                                 shape=(int(g["shape"][0]), int(g["shape"][1])))
            gi = np.where(gt == "Gene Expression")[0]
            X = full[gi].T.tocsr().astype(np.float32)
            names_g = nm[gi]
        want = np.isin(names_g, UNION)
        idx = np.where(want)[0]
        Xg = X[:, idx]
        s = Xg.sum(1).astype(np.float32)
        from scipy.sparse import diags as _diags
        Gn = np.log1p(_diags((1e4 / np.maximum(np.asarray(s), 1)).ravel()).dot(Xg).toarray()).astype(np.float32)
        vv = names_g[idx]
        a = ad.AnnData(X=Gn)
        a.var_names = pd.Index(vv); a.obs_names = pd.Index(["c%d" % i for i in range(Gn.shape[0])])
        a.var_names_make_unique()
        scdrs.pp.preprocess(a)
        R = a.uns["SCDRS_PARAM"]["GENE_STATS"]
        C = GS.loc[list(vv)]
        same_mean = float(np.max(np.abs(R["mean"].values - C["mean"].values)))
        same_var = float(np.max(np.abs(R["var"].values - C["var"].values)))
        same_ct = float(np.max(np.abs(R["ct_mean"].values - C["ct_mean"].values)))
        log("P2 check (%s all nuclei): max|delta mean|=%.2e  max|delta var|=%.2e  max|delta ct_mean|=%.2e  %s" % (
            VERIFY_LIB, same_mean, same_var, same_ct,
            "OK" if max(same_mean, same_var, same_ct) < 1e-5 else "FAIL"))
