#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""O1 own dsciCAT-seq same-nucleus paired RNA x ATAC cross-channel test (Route A's "own-data M2").

Background: Fig 6's primary RNA evidence comes from our own dsciCAT-seq mouse OL multiome (5,977 nuclei,
RNA and ATAC gene scores from the SAME NUCLEI, obs_names identical). The 10-23 S6H/I only did
"each channel's own shape along the maturation ladder"; it never tested same-NUCLEUS-level cross-channel concordance, nor brought
Route A's measured "gene-set regulatory-domain size" confounding back onto our own data. This script adds those three things.

Four tests (all on the same 5,977 nuclei only):
 T1 same-nucleus cross-channel concordance: per nucleus, the Spearman of its RNA norm_score vector vs its ATAC norm_score vector
    across the 74 traits (raw / double-centred), null = permuting the ATAC vectors across nuclei (breaks pairing but keeps
    each channel's marginal distribution).
 T2 single-trait (MS) per-nucleus correlation: raw / within-subtype / after depth-residualization (nCount_RNA, nFeature_RNA,
    nCount_ATAC, nFeature_ATAC spline/linear residuals).
 T3 state-level concordance: the six subtypes' mean profiles (74-dimensional) cross-channel Spearman; plus Δ1 along the maturation ladder.
 T4 74-trait Δ1 cross-channel concordance: Δ1 = median(MOL6) - median(COP/NFOL) (the main convention),
    the sensitivity convention = (MOL6+MOL2) - (MFOL+COP/NFOL). Raw ranks + ranks after deconfounding on log(n_hit)
    (n_hit = the trait's gene set actually measurable in that channel = the within-channel proxy for regulatory-domain size).

Memory: only reads the 74 x 5,977 score table; peak < 2 GB.
"""
import os, glob, gzip, time
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, binomtest
import statsmodels.api as sm

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

t0 = time.time()
def log(*a): print("[%6ds]" % (time.time() - t0), *a, flush=True)

D_RNA  = os.path.join(os.environ.get('DSCI_FY_ROOT', WORK + "/FY"), 'oligo/scdrs_traits_draft/full74')
D_ATAC = os.path.join(os.environ.get('DSCI_FY_ROOT', WORK + "/FY"), 'oligo/fig6_rerun_20261019/dualomics_probe/out_gs74')
GS     = os.path.join(os.environ.get('DSCI_FY_ROOT', WORK + "/FY"), 'oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs')
H5_R   = os.path.join(os.environ.get('DSCI_FY_ROOT', WORK + "/FY"), 'oligo/3k_dsci_oligo_subset.h5ad')
H5_A   = os.path.join(os.environ.get('DSCI_FY_ROOT', WORK + "/FY"), 'oligo/fig6_rerun_20261019/dualomics_probe/genescore_oligo_5977.h5ad')
OUT    = os.path.join(os.environ.get('DSCI_FIG_ROOT', WORK + "/figures"), 'routeA_out/O1_own_multiome_paired')
TAB, DER = OUT + "/tables", OUT + "/derived"
os.makedirs(TAB, exist_ok=True); os.makedirs(DER, exist_ok=True)
MAT, COP = "MOL6", "COP/NFOL"
DEP_R = ["nCount_RNA", "nFeature_RNA"]
DEP_A = ["nCount_ATAC", "nFeature_ATAC"]
DEP = DEP_R + DEP_A
rng = np.random.default_rng(20261030)

def load_scores(root):
    """{trait: pd.Series(norm_score, index=cell)}; same-named files are unique within a batch directory"""
    out = {}
    for f in glob.glob(f"{root}/*/*.score.gz"):
        tr = os.path.basename(f).replace(".score.gz", "")
        if tr in out:
            continue
        s = pd.read_csv(f, sep="\t", index_col=0)["norm_score"]
        s.name = tr
        out[tr] = s
    return pd.DataFrame(out)

log("reading RNA per-nucleus scores ..."); Sr = load_scores(D_RNA)
log("reading ATAC per-nucleus scores ..."); Sa = load_scores(D_ATAC)
log("RNA %s | ATAC %s" % (Sr.shape, Sa.shape))

import anndata as ad
obsR = ad.read_h5ad(H5_R, backed="r").obs
obsA = ad.read_h5ad(H5_A, backed="r").obs
cells = Sr.index.intersection(Sa.index)
log("same-nucleus intersection (both channels have scores) %d" % len(cells))
Sr = Sr.loc[cells]; Sa = Sa.loc[cells]
meta = obsR.loc[cells, ["cell_type"] + DEP_R].join(obsA.loc[cells, DEP_A])
meta["cell_type"] = meta["cell_type"].astype(str)
for c in DEP:
    meta[c] = pd.to_numeric(meta[c], errors="coerce").astype(float)
meta["log_depR"] = np.log1p(meta.nCount_RNA) ; meta["log_depA"] = np.log1p(meta.nCount_ATAC)
log("subtype nucleus counts:\n" + meta.cell_type.value_counts().to_string())

traits = [t for t in Sr.columns if t in Sa.columns]
log("shared traits %d" % len(traits))
Sr, Sa = Sr[traits], Sa[traits]

# ---------- the gene-set measurable count n_hit (= the within-channel proxy for regulatory-domain size) ----------
sets = {}
for l in open(GS):
    f = l.rstrip("\n").split("\t")
    if len(f) < 2 or f[0] == "TRAIT":
        continue
    sets[f[0]] = [x.strip().upper() for x in f[1].split(",") if x.strip()]
vr = set(pd.read_hdf if False else ad.read_h5ad(H5_R, backed="r").var_names.str.upper())
va = set(ad.read_h5ad(H5_A, backed="r").var_names.str.upper())
n_hit = {}
for t in traits:
    g = sets.get(t, [])
    n_hit[t] = {"n_gs": len(g), "n_hit_rna": sum(x in vr for x in g), "n_hit_atac": sum(x in va for x in g)}
log("n_hit median RNA %.0f / ATAC %.0f" % (np.median([v["n_hit_rna"] for v in n_hit.values()]),
                                        np.median([v["n_hit_atac"] for v in n_hit.values()])))

# =============== T1 same-nucleus cross-channel concordance ===============
Xr, Xa = Sr[traits].values, Sa[traits].values
def per_cell_rho(A, B):
    """per nucleus: Spearman of its A vector vs its B vector across the 74 traits"""
    ra = np.argsort(np.argsort(A, axis=1), axis=1); rb = np.argsort(np.argsort(B, axis=1), axis=1)
    a = ra - ra.mean(1, keepdims=True); b = rb - rb.mean(1, keepdims=True)
    num = (a * b).sum(1); den = np.sqrt((a ** 2).sum(1) * (b ** 2).sum(1))
    return num / np.maximum(den, 1e-12)
rho_raw = per_cell_rho(Xr, Xa)
# double centring (removes "this nucleus runs high" / "this trait runs high")
Zr = Xr - Xr.mean(1, keepdims=True); Zr = Zr - Zr.mean(0, keepdims=True)
Za = Xa - Xa.mean(1, keepdims=True); Za = Za - Za.mean(0, keepdims=True)
rho_dc = per_cell_rho(Zr, Za)
# the null: ATAC vectors permuted across nuclei (keeps both channels' marginals, breaks only the pairing)
n_null = 300
null_raw = np.array([per_cell_rho(Xr, Xa[rng.permutation(len(Xa))]).mean() for _ in range(n_null)])
null_dc = np.array([per_cell_rho(Zr, Za[rng.permutation(len(Za))]).mean() for _ in range(n_null)])
T1 = pd.DataFrame([
    {"test": "T1 per-cell cross-channel rho (raw)", "median_rho": np.median(rho_raw),
     "mean_rho": rho_raw.mean(), "frac_pos": (rho_raw > 0).mean(),
     "null_mean": null_raw.mean(), "null_sd": null_raw.std(),
     "Z_vs_null": (rho_raw.mean() - null_raw.mean()) / null_raw.std(),
     "P_emp": (np.sum(null_raw >= rho_raw.mean()) + 1) / (n_null + 1), "n_cell": len(cells)},
    {"test": "T1 per-cell cross-channel rho (double-centered)", "median_rho": np.median(rho_dc),
     "mean_rho": rho_dc.mean(), "frac_pos": (rho_dc > 0).mean(),
     "null_mean": null_dc.mean(), "null_sd": null_dc.std(),
     "Z_vs_null": (rho_dc.mean() - null_dc.mean()) / null_dc.std(),
     "P_emp": (np.sum(null_dc >= rho_dc.mean()) + 1) / (n_null + 1), "n_cell": len(cells)},
])
log("\nT1:\n" + T1.round(4).to_string(index=False))
meta["xchan_rho_raw"] = rho_raw; meta["xchan_rho_dc"] = rho_dc
T1.to_csv(TAB + "/T1_same_nucleus_crosschannel_rho.csv", index=False)

# =============== T2 MS per-nucleus correlation ===============
def partial(v1, v2, cov):
    C = sm.add_constant(cov)
    r1 = sm.OLS(v1, C, missing="drop").fit().resid
    r2 = sm.OLS(v2, C, missing="drop").fit().resid
    m = ~(np.isnan(r1) | np.isnan(r2))
    return spearmanr(r1[m], r2[m])
rows2 = []
for tr in ["PASS_Multiple_sclerosis", "PASS_Schizophrenia_Pardinas2018", "UKB_460K.body_HEIGHTz"]:
    if tr not in traits:
        continue
    v1, v2 = Sr[tr].values, Sa[tr].values
    r_all = spearmanr(v1, v2)
    r_dc_, _ = None, None
    rd = partial(v1, v2, meta[["log_depR", "log_depA"]].values)
    sub = []
    for c, d in meta.groupby("cell_type"):
        m = (meta.cell_type == c).values
        if m.sum() >= 50:
            r = spearmanr(v1[m], v2[m])
            sub.append("%s %.3f" % (c, r.correlation))
    rows2.append({"trait": tr, "n_cell": int(m.sum() and len(cells)),
                  "rho_raw": r_all.correlation, "P_raw": r_all.pvalue,
                  "rho_depth_adjusted": rd.correlation, "P_depth_adjusted": rd.pvalue,
                  "within_subtype": "; ".join(sub)})
    # kept for the figure
    meta["_".join([tr, "rna"])] = v1; meta["_".join([tr, "atac"])] = v2
T2 = pd.DataFrame(rows2)
log("\nT2:\n" + T2[["trait", "rho_raw", "P_raw", "rho_depth_adjusted", "P_depth_adjusted"]].round(4).to_string(index=False))
T2.to_csv(TAB + "/T2_MS_per_nucleus_rna_vs_atac.csv", index=False)

# =============== T3 state-level profile concordance ===============
sub_order = [MAT, "MOL2", "MOL3/6", "MOL3", "MFOL", COP]
mr_ = {c: (meta.cell_type == c).values for c in sub_order}
prof_r = np.array([Xr[mr_[c]].mean(0) for c in sub_order])
prof_a = np.array([Xa[mr_[c]].mean(0) for c in sub_order])
r3 = spearmanr(prof_r.ravel(), prof_a.ravel())
r3t = [spearmanr(prof_r[:, j], prof_a[:, j]) for j in range(len(traits))]     # per trait: the 6-subtype profiles
T3 = pd.DataFrame([{"test": "T3a subtype x trait full matrix (%d x %d)" % (len(sub_order), len(traits)),
                    "rho": r3.correlation, "P": r3.pvalue,
                    "note": "includes shared inter-trait variation; not independent units"},
                   {"test": "T3b within-trait 6-subtype profile cross-channel rho",
                    "rho": np.median([x.correlation for x in r3t]),
                    "P": np.median([x.pvalue for x in r3t]),
                    "frac_pos": float(np.mean([x.correlation > 0 for x in r3t])),
                    "note": "6 points/trait, direction-finding only"}])
log("\nT3:\n" + T3.round(4).to_string(index=False))
T3.to_csv(TAB + "/T3_state_profile_concordance.csv", index=False)

# =============== T4 74-trait Δ1 cross-channel (the main analysis) ===============
def delta1(X, mask_mat, mask_cop, w=0):
    if w == 0:
        return np.median(X[mask_mat], 0) - np.median(X[mask_cop], 0)
    return X[mask_mat].mean(0) - X[mask_cop].mean(0)

d1r = delta1(Xr, mr_[MAT], mr_[COP]); d1a = delta1(Xa, mr_[MAT], mr_[COP])
g1r = delta1(Xr, np.isin(meta.cell_type, [MAT, "MOL2"]), np.isin(meta.cell_type, ["MFOL", COP]), w=1)
g1a = delta1(Xa, np.isin(meta.cell_type, [MAT, "MOL2"]), np.isin(meta.cell_type, ["MFOL", COP]), w=1)
T = pd.DataFrame({"trait": traits, "d1_rna": d1r, "d1_atac": d1a, "d1b_rna": g1r, "d1b_atac": g1a,
                  "n_gs": [n_hit[t]["n_gs"] for t in traits],
                  "n_hit_rna": [n_hit[t]["n_hit_rna"] for t in traits],
                  "n_hit_atac": [n_hit[t]["n_hit_atac"] for t in traits],
                  "n_cell_mat": mr_[MAT].sum(), "n_cell_cop": mr_[COP].sum()})

# check MS against the delivered S6H numbers (same nuclei, same convention: they must agree)
i = T.index[T.trait == "PASS_Multiple_sclerosis"][0]
log("\nMS Δ1: RNA %+.4f | ATAC %+.4f | sensitivity RNA %+.4f / ATAC %+.4f" %
    (T.d1_rna[i], T.d1_atac[i], T.d1b_rna[i], T.d1b_atac[i]))

def adjust(col, by):
    x = np.log1p(T[by].values.astype(float)); y = T[col].values
    res = sm.OLS(y, sm.add_constant(x)).fit()
    return res.resid, res            # deconfounding = regression residuals (previously wrongly written as fitted values; fixed)

T["rank_rna_raw"] = T.d1_rna.rank(ascending=False).astype(int)
T["rank_atac_raw"] = T.d1_atac.rank(ascending=False).astype(int)
for ch, by in [("rna", "n_hit_rna"), ("atac", "n_hit_atac")]:
    adj, fit = adjust("d1_" + ch, by)
    T["d1_" + ch + "_adj"] = adj
    T["rank_" + ch + "_adj"] = pd.Series(adj, index=T.index).rank(ascending=False).astype(int)
    log("deconfounded d1_%s ~ log(n_hit_%s): slope %+.4f, R2 %.3f, P %.3g" %
        (ch, ch, fit.params[1], fit.rsquared, fit.pvalues[1]))
# ---- convention C: residualize each nucleus's scores on depth first, then Δ1 (removes the "mature nuclei have different ATAC/RNA detection" artefact) ----
def depth_resid(X, cov):
    C = sm.add_constant(cov)
    out = np.empty_like(X)
    for j in range(X.shape[1]):
        out[:, j] = sm.OLS(X[:, j], C).fit().resid
    return out
Xd_r = depth_resid(Xr, meta[["log_depR"]].values)
Xd_a = depth_resid(Xa, meta[["log_depA"]].values)
T["d1_rna_depthres"] = delta1(Xd_r, mr_[MAT], mr_[COP])
T["d1_atac_depthres"] = delta1(Xd_a, mr_[MAT], mr_[COP])
T["rank_rna_depthres"] = T.d1_rna_depthres.rank(ascending=False).astype(int)
T["rank_atac_depthres"] = T.d1_atac_depthres.rank(ascending=False).astype(int)
for ch in ["rna", "atac"]:
    adj, fit = adjust("d1_" + ch + "_depthres", "n_hit_" + ch)
    T["d1_" + ch + "_dr_adj"] = adj
    T["rank_" + ch + "_dr_adj"] = pd.Series(adj, index=T.index).rank(ascending=False).astype(int)
VAR = [("raw", "d1_rna", "d1_atac"),
       ("depth-residualised", "d1_rna_depthres", "d1_atac_depthres"),
       ("depth-resid + n_hit-adjusted", "d1_rna_dr_adj", "d1_atac_dr_adj"),
       ("n_hit-adjusted only", "d1_rna_adj", "d1_atac_adj"),
       ("sens: (MOL6+MOL2)-(MFOL+COP)", "d1b_rna", "d1b_atac")]
rr = []
for nm, a, b in VAR:
    rho, pv = spearmanr(T[a], T[b])
    # the permutation null: shuffle the ATAC trait order wholesale (breaks the "same trait" correspondence)
    nullv = []
    for _ in range(2000):
        o = rng.permutation(len(T))
        nullv.append(spearmanr(T[a].values, T[b].values[o]).correlation)
    nullv = np.array(nullv)
    rr.append({"variant": nm, "rho_cross_channel": rho, "P": pv,
               "null_mean": nullv.mean(), "null_975": np.quantile(nullv, .975),
               "P_perm": (np.sum(np.abs(nullv) >= abs(rho)) + 1) / 2001,
               "MS_rank_rna": int(T[a and a and "rank_" + a.split("_")[1] + ("_raw" if a in ("d1_rna", "d1_atac") else "")] if False else 0) if False else 0})
RHO = pd.DataFrame(rr).drop(columns=["MS_rank_rna"])
RHO.to_csv(TAB + "/T4b_cross_channel_rho_variants.csv", index=False)
log("\nT4b cross-channel Δ1 concordance (5 conventions + permutation null):\n" + RHO.round(4).to_string(index=False))
T.to_csv(TAB + "/T4_delta1_74traits_both_channels.csv", index=False)
T["rho_raw_all"] = spearmanr(T.d1_rna, T.d1_atac).correlation
T["rho_adj_all"] = spearmanr(T.d1_rna_adj, T.d1_atac_adj).correlation
T.to_csv(TAB + "/T4_delta1_74traits_both_channels.csv", index=False)

r_raw = spearmanr(T.d1_rna, T.d1_atac); r_adj = spearmanr(T.d1_rna_adj, T.d1_atac_adj)
r_b = spearmanr(T.d1b_rna, T.d1b_atac)
log("\nT4 74-trait cross-channel Δ1 concordance: raw rho %+.3f (P %.3g) -> after removing n_hit confounding rho %+.3f (P %.3g)"
    % (r_raw.correlation, r_raw.pvalue, r_adj.correlation, r_adj.pvalue))
log("sensitivity convention (MOL6+MOL2 - MFOL+COP/NFOL) rho %+.3f (P %.3g)" % (r_b.correlation, r_b.pvalue))
log("MS rank: RNA raw %d/%d (deconfounded %d) | ATAC raw %d/%d (deconfounded %d)"
    % (T.rank_rna_raw[i], len(T), T.rank_rna_adj[i], T.rank_atac_raw[i], len(T), T.rank_atac_adj[i]))

# per-channel sign test (trait level): the number of traits with Δ1 > 0
sym = pd.DataFrame([{"channel": "RNA", "n_pos_d1gt0": int((T.d1_rna > 0).sum()), "median_d1": float(np.median(T.d1_rna)),
                     "P_sign": binomtest(int((T.d1_rna > 0).sum()), len(T), .5).pvalue,
                     "spread": float(np.std(T.d1_rna))},
                    {"channel": "ATAC", "n_pos_d1gt0": int((T.d1_atac > 0).sum()), "median_d1": float(np.median(T.d1_atac)),
                     "P_sign": binomtest(int((T.d1_atac > 0).sum()), len(T), .5).pvalue,
                     "spread": float(np.std(T.d1_atac))}])
log("\nwithin-channel Δ1 distribution (trait level):\n" + sym.round(4).to_string(index=False))
sym["abs_ratio_RNA_over_ATAC"] = sym.spread.iloc[0] / sym.spread.iloc[1]
sym.to_csv(TAB + "/T4_channel_dispersion.csv", index=False)


# =============== T5 split-half reliability (the noise ceiling) + independent-nucleus-subset cross-channel concordance ===============
# split the nuclei randomly into A/B halves (MOL6 and COP/NFOL each stratified), each half computes the 74-trait Δ1 once.
# the same-channel two-half correlation = the reliability ceiling of that channel's Δ1 estimate; if RNA 0.9 / ATAC 0.5, a cross-channel 0.82 is already near the ceiling.
idx_mat = np.where(mr_[MAT])[0]; idx_cop = np.where(mr_[COP])[0]
rep = {}
for it in range(20):
    pa = rng.permutation(len(idx_mat)); ha1, ha2 = idx_mat[pa[:len(pa) // 2]], idx_mat[pa[len(pa) // 2:]]
    pb = rng.permutation(len(idx_cop)); kb1, kb2 = idx_cop[pb[:len(pb) // 2]], idx_cop[pb[len(pb) // 2:]]
    for nm, (m1, m2) in {"half1": ((ha1, kb1), None), "half2": ((ha2, kb2), None)}.items():
        pass
    d1r1 = np.median(Xr[ha1], 0) - np.median(Xr[kb1], 0); d1r2 = np.median(Xr[ha2], 0) - np.median(Xr[kb2], 0)
    d1a1 = np.median(Xa[ha1], 0) - np.median(Xa[kb1], 0); d1a2 = np.median(Xa[ha2], 0) - np.median(Xa[kb2], 0)
    rep.setdefault("rna_within_channel", []).append(spearmanr(d1r1, d1r2).correlation)
    rep.setdefault("atac_within_channel", []).append(spearmanr(d1a1, d1a2).correlation)
    rep.setdefault("xchan_half1_vs_half2", []).append(spearmanr(d1r1, d1a2).correlation)
    rep.setdefault("xchan_same_half", []).append(spearmanr(d1r1, d1a1).correlation)
T5 = pd.DataFrame([{"metric": k, "median_split_half_rho": float(np.median(v)),
                    "q25": float(np.quantile(v, .25)), "q75": float(np.quantile(v, .75)), "n_draw": len(v)}
                   for k, v in rep.items()])
T5.to_csv(TAB + "/T5_split_half_reliability.csv", index=False)
log("\nT5 split-half reliability (20 draws, ~1,288 MOL6 / ~101 COP nuclei each):\n" + T5.round(3).to_string(index=False))


# =============== T6 "how far off the diagonal": does MS break out only in the expression channel ===============
# fit d1_rna ~ d1_atac with rank regression (outlier-robust), look at each trait's residual rank; then give the split-half ceiling
fit_r = sm.RLM(T.d1_rna.values, sm.add_constant(T.d1_atac.values), M=sm.robust.norms.HuberT()).fit()
res = T.d1_rna.values - fit_r.fittedvalues
T["offdiag_resid"] = res
T["offdiag_rank"] = pd.Series(res, index=T.index).rank(ascending=False).astype(int)
T["offdiag_z"] = (res - res.mean()) / res.std()
i_ms = int(np.where(T.trait.values == "PASS_Multiple_sclerosis")[0][0])
ceil_ = float(np.sqrt(T5.median_split_half_rho[0] * T5.median_split_half_rho[1]))
corr_att = float(float(T.rho_raw_all.iloc[0]) / ceil_)
T.to_csv(TAB + "/T4_delta1_74traits_both_channels.csv", index=False)
summ = pd.DataFrame([{"item": "cross-channel rho (raw)", "value": float(T.rho_raw_all.iloc[0])},
                     {"item": "noise ceiling sqrt(relRNA*relATAC)", "value": ceil_},
                     {"item": "attenuation-corrected rho", "value": corr_att},
                     {"item": "independent-nucleus-halves rho (median)", "value": float(T5[T5.metric == "xchan_half1_vs_half2"].median_split_half_rho.iloc[0])},
                     {"item": "MS off-diagonal residual z", "value": float(T.offdiag_z[i_ms])},
                     {"item": "MS off-diagonal rank of 74", "value": float(T.offdiag_rank[i_ms])}])
summ.to_csv(TAB + "/T6_cross_channel_summary.csv", index=False)
log("\nT6 summary:\n" + summ.round(3).to_string(index=False))
log("MS off-diagonal residual z=%+.2f (rank %d of 74); mirror control SCZ residual z=%+.2f (rank %d)"
    % (T.offdiag_z[i_ms], T.offdiag_rank[i_ms],
       T.offdiag_z[int(np.where(T.trait.values == "PASS_Schizophrenia_Pardinas2018")[0][0])],
       T.offdiag_rank[int(np.where(T.trait.values == "PASS_Schizophrenia_Pardinas2018")[0][0])]))

meta.reset_index().rename(columns={"index": "cell_id"}).to_csv(TAB + "/cell_level_table.csv", index=False)
np.save(DER + "/prof_rna.npy", prof_r); np.save(DER + "/prof_atac.npy", prof_a)
pd.Series(sub_order).to_csv(DER + "/sub_order.csv", index=False)
log("written to %s" % TAB)
