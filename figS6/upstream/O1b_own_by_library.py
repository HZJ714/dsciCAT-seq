#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""O1b: split O1 BY LIBRARY (author correction: our own data is not one library but two, 50k=FY13 and 100k=FY14)

Per $DSCI_WORK_ROOT/FY/oligo/fig6_rerun_20261019/ATAC_GWAS_dualomics_feasibility_20261022.md:326
  library = dsciCAT-seq-50k (FY13) / dsciCAT-seq-100k (FY14), with library == donor confounded 1:1.
The cell_id prefix is the library name (FY13_/FY14_), so stratification is direct.

Three things:
 L1 per-library 74-trait Δ1 (both channels) + cross-channel rho + split-half reliability -> can the concordance hold in each library
 L2 per-library per-nucleus cross-channel rho (the T1/T2 convention) -> is the weak coupling present in both
 D1 redundancy check: how correlated is the current Fig 6D trait vector (MOL6 selectivity, row z) with O1's RNA Δ1
"""
import os, glob, time
import numpy as np, pandas as pd
from scipy.stats import spearmanr

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

t0 = time.time()
def log(*a): print("[%6ds]" % (time.time() - t0), *a, flush=True)
D_RNA  = WORK + "/FY/oligo/scdrs_traits_draft/full74"
D_ATAC = WORK + "/FY/oligo/fig6_rerun_20261019/dualomics_probe/out_gs74"
H5_R   = WORK + "/FY/oligo/3k_dsci_oligo_subset.h5ad"
H5_A   = WORK + "/FY/oligo/fig6_rerun_20261019/dualomics_probe/genescore_oligo_5977.h5ad"
TAB    = WORK + "/figures/routeA_out/O1_own_multiome_paired/tables"
MAT, COP, MS, SCZ = "MOL6", "COP/NFOL", "PASS_Multiple_sclerosis", "PASS_Schizophrenia_Pardinas2018"
rng = np.random.default_rng(7)

def load_scores(root):
    out = {}
    for f in glob.glob(f"{root}/*/*.score.gz"):
        tr = os.path.basename(f).replace(".score.gz", "")
        if tr in out: continue
        s = pd.read_csv(f, sep="\t", index_col=0)["norm_score"]; s.name = tr; out[tr] = s
    return pd.DataFrame(out)

Sr, Sa = load_scores(D_RNA), load_scores(D_ATAC)
import anndata as ad
obsR = ad.read_h5ad(H5_R, backed="r").obs
cells = Sr.index.intersection(Sa.index)
Sr, Sa = Sr.loc[cells], Sa.loc[cells].astype(float)
traits = [t for t in Sr.columns if t in Sa.columns]
Sr, Sa = Sr[traits].astype(float), Sa[traits]
lib = pd.Series([str(c).split("_")[0] for c in cells], index=cells, name="lib")
ct  = obsR.loc[cells, "cell_type"].astype(str)
log("library distribution: " + lib.value_counts().to_dict().__str__())
log("library x subtype:\n" + pd.crosstab(lib.values, ct.values).to_string())

rows, nuc = [], []
per_lib = {}
for L in sorted(lib.unique()):
    m = (lib == L).values
    sub_ct = ct[m]
    ok = all((sub_ct == k).sum() >= 20 for k in (MAT, COP))
    Xr, Xa = Sr.values[m], Sa.values[m]
    d1r = np.median(Xr[sub_ct == MAT], 0) - np.median(Xr[sub_ct == COP], 0)
    d1a = np.median(Xa[sub_ct == MAT], 0) - np.median(Xa[sub_ct == COP], 0)
    rho = spearmanr(d1r, d1a).correlation
    null = [spearmanr(d1r, d1a[rng.permutation(len(d1a))]).correlation for _ in range(2000)]
    P = (sum(abs(v) >= abs(rho) for v in null) + 1) / (len(null) + 1)
    T = pd.DataFrame({"trait": traits, "d1_rna": d1r, "d1_atac": d1a})
    rk = lambda c, tr: int(T[c].rank(ascending=False)[T.trait == tr].iloc[0])
    import statsmodels.api as sm
    fitL = sm.OLS(np.asarray(d1r).ravel(), sm.add_constant(np.asarray(d1a).ravel())).fit()
    residL = np.asarray(d1r - np.asarray(fitL.fittedvalues).ravel()).ravel()
    rkS = pd.Series(residL, index=traits).rank(ascending=False)
    off_ms = float(np.ravel(rkS.loc[MS])[0])
    off_z = float(((residL - residL.mean()) / residL.std())[traits.index(MS)])
    # split-half reliability (within library)
    rel = {}
    for ch, X in [("rna", Xr), ("atac", Xa)]:
        v = []
        for _ in range(10):
            pa = rng.permutation(np.where(sub_ct == MAT)[0]); pb = rng.permutation(np.where(sub_ct == COP)[0])
            a1, a2 = pa[:len(pa)//2], pa[len(pa)//2:]; b1, b2 = pb[:len(pb)//2], pb[len(pb)//2:]
            d1 = np.median(X[a1], 0) - np.median(X[b1], 0); d2 = np.median(X[a2], 0) - np.median(X[b2], 0)
            v.append(spearmanr(d1, d2).correlation)
        rel[ch] = float(np.median(v))
    # per-nucleus cross-channel rho (the library's 74-trait vectors) + the MS/SCZ single traits
    rr = np.array([spearmanr(Xr[i], Xa[i]).correlation for i in range(0, m.sum(), max(1, m.sum() // 2500))])
    i_ms, i_sc = traits.index(MS), traits.index(SCZ)
    for nm, v in [("MS", Xr[:, i_ms]), ("SCZ", Xr[:, i_sc])]:
        nuc.append({"lib": L, "item": "per_nucleus_rho_" + nm, "value": spearmanr(v, Xa[:, i_sc if nm == "SCZ" else i_ms]).correlation})
    per_lib[L] = pd.DataFrame({"trait": traits, "d1_rna_" + L: d1r, "d1_atac_" + L: d1a,
                               "rank_rna_" + L: pd.Series(d1r, index=traits).rank(ascending=False).astype(int).values,
                               "rank_atac_" + L: pd.Series(d1a, index=traits).rank(ascending=False).astype(int).values})
    rows.append({"lib": L, "n_cells": int(m.sum()), "n_MOL6": int((sub_ct == MAT).sum()), "n_COP": int((sub_ct == COP).sum()),
                 "usable_gap": ok, "rho_xchan": rho, "P_perm": P, "null975": float(np.quantile(np.abs(null), .975)),
                 "d1_rna_MS": d1r[i_ms], "rank_rna_MS": rk("d1_rna", MS), "d1_atac_MS": d1a[i_ms], "rank_atac_MS": rk("d1_atac", MS),
                 "d1_rna_SCZ": d1r[i_sc], "rank_rna_SCZ": rk("d1_rna", SCZ), "d1_atac_SCZ": d1a[i_sc], "rank_atac_SCZ": rk("d1_atac", SCZ),
                 "MS_offdiag_rank": off_ms, "MS_offdiag_z": off_z,
                 "split_half_rna": rel["rna"], "split_half_atac": rel["atac"],
                 "median_per_nucleus_rho_all74": float(np.median(rr)), "frac_pos_per_nucleus": float((rr > 0).mean())})
R = pd.DataFrame(rows); R.to_csv(TAB + "/L1_by_library.csv", index=False)
pd.DataFrame(nuc).to_csv(TAB + "/L1b_by_library_per_nucleus.csv", index=False)
log("\nby library:\n" + R.round(3).to_string(index=False))
log("\nper-nucleus single-trait rho:\n" + pd.DataFrame(nuc).round(3).to_string(index=False))

# ---------- D1 redundancy check: the current Fig 6D vectors vs O1's RNA Δ1 ----------
d1r_all = np.median(Sr.values[ct.values == MAT], 0) - np.median(Sr.values[ct.values == COP], 0)
d1a_all = np.median(Sa.values[ct.values == MAT], 0) - np.median(Sa.values[ct.values == COP], 0)
C6D = WORK + "/FY/oligo/scdrs_traits_draft/full74/arialfix_trait_ranking.csv"
log("\n6D source table %s columns: %s" % (os.path.basename(C6D), list(pd.read_csv(C6D).columns)))
F = pd.read_csv(C6D)
F.index = F["Unnamed: 0"].astype(str)
F = F[F.index.isin(traits)]
d = pd.DataFrame({"d1_rna": pd.Series(d1r_all, index=traits), "d1_atac": pd.Series(d1a_all, index=traits)}).join(F, how="inner")
out = []
for c in d.columns:
    if c in ("d1_rna", "d1_atac", "Unnamed: 0") or not pd.api.types.is_numeric_dtype(d[c]): continue
    r = spearmanr(d["d1_rna"], d[c]).correlation
    out.append({"6D_column": c, "spearman_vs_O1_RNA_delta1": r})
D = pd.DataFrame(out).sort_values("spearman_vs_O1_RNA_delta1", key=lambda s: -s.abs())
D.to_csv(TAB + "/D1_redundancy_fig6D_vs_O1.csv", index=False)
log("\nD1 redundancy check (current Fig 6D columns vs O1's RNA Δ1, n=%d traits):\n%s" % (len(d), D.round(3).to_string(index=False)))

# merge the per-library Δ1 columns back into T4 (the figure script reads T4's d1_rna_FY13 / rank_atac_FY14 etc. directly)
T4f = TAB + "/T4_delta1_74traits_both_channels.csv"
T4 = pd.read_csv(T4f)
for L, PL in per_lib.items():
    T4 = T4.merge(PL, on="trait", how="left")
T4.to_csv(T4f, index=False)
log("per-library Δ1 columns merged back into T4: %s" % [c for c in T4.columns if "FY" in c])

# ---------- T7 subtype-level (cell-state identity) cross-channel concordance ----------
ORDER = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
def zmat(M):   # z-score the six subtypes within each trait (removes each trait's own level)
    return M.sub(M.mean(1), axis=0).div(M.std(1) + 1e-12, axis=0)
t7 = []
for L, m in [("pooled", np.ones(len(cells), bool)), ("FY13", (lib == "FY13").values), ("FY14", (lib == "FY14").values)]:
    ctv = ct.values[m]
    Mr = pd.DataFrame({k: Sr.values[m][ctv == k].mean(0) for k in ORDER}, index=traits)   # traits x subtypes
    Ma = pd.DataFrame({k: Sa.values[m][ctv == k].mean(0) for k in ORDER}, index=traits)
    zr, za = zmat(Mr), zmat(Ma)
    pool = spearmanr(zr.values.ravel(), za.values.ravel()).correlation
    per = [spearmanr(zr.loc[t].values, za.loc[t].values).correlation for t in zr.index]
    axis_r, axis_a = Mr.mean(0), Ma.mean(0)          # each subtype's mean z across the 74 traits (the maturation ladder)
    t7.append({"lib": L, "rho_pooled_6x74": pool, "median_rho_per_trait_6pts": float(np.median(per)),
               "frac_pos_per_trait": float(np.mean([v > 0 for v in per])),
               "rho_subtype_axis_6pts": spearmanr(axis_r.values, axis_a.values).correlation,
               "P_subtype_axis": spearmanr(axis_r.values, axis_a.values).pvalue,
               **{("RNA_" + k): float(axis_r[k]) for k in ORDER}, **{("ATAC_" + k): float(axis_a[k]) for k in ORDER}})
T7 = pd.DataFrame(t7); T7.to_csv(TAB + "/T7_subtype_level_concordance.csv", index=False)
log("\nT7 subtype-level cross-channel concordance:\n" + T7.round(3).to_string(index=False))

log("done")
