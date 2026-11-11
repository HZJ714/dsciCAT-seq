# -*- coding: utf-8 -*-
"""Fig 6 v15 numbers to ledger: every statistic printed on the face goes into tables, and the drawing scripts only read tables (review receipt §5-1 / §2 / §13)
Output $DSCI_WORK_ROOT/figures/routeA_out/FIG6_v15/tables/
  fig6A_stats.csv      per-subtype n / mean / median / whisker (5-95) / MWU P / group P / BH q (5-contrast family) / Cliff's d / scDRS internal FDR passes
  fig6C_stats.csv      level vs contrast: rho, P, trait-level-removed rho, noise ceiling, attenuation-corrected rho, MS off-diag z/rank, each channel's top5
  fig6D_stats.csv      disease two-arm Δ / MWU P / Cliff's d / Hedges g; cross-channel rho / 95% CI (Fisher z and bootstrap) / perm P (n_perm, seed) / AUC interval
  fig6_all_numbers.csv flat key-value (for bit-by-bit checking)
"""
import os, sys

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
FY = os.environ.get("DSCI_FY_ROOT", WORK + "/FY")
CONFIG = {"H5AD": os.path.join(FY, "oligo", "3k_dsci_oligo_subset.h5ad"),
          "DATA": os.path.join(FY, "oligo", "scdrs_ms_out")}
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, rankdata
from config import CONFIG
import scanpy as sc   # CONFIG is set from DSCI_FY_ROOT above; no external config module

FIGR = os.environ.get("DSCI_FIG_ROOT", WORK + "/figures")
OUTD = os.path.join(FIGR, "routeA_out/FIG6_v15/tables")
os.makedirs(OUTD, exist_ok=True)
T4D = os.path.join(FIGR, "routeA_out/O1_own_multiome_paired/tables")
M2B = os.path.join(FIGR, "routeA_out/M2b_authoranno_v1/tables")
SEED = 20261102
NPERM = 10000
flat = {}

# ---------------- A ----------------
adata = sc.read_h5ad(CONFIG["H5AD"])
score = pd.read_csv(os.path.join(CONFIG["DATA"], "PASS_Multiple_sclerosis.score.gz"), sep="\t", index_col=0)
df = adata.obs[["cell_type"]].join(adata.obs[["UMAP_1", "UMAP_2"]]).join(score[["norm_score"]], how="inner").dropna()
order = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
st = pd.read_csv(os.path.join(CONFIG["DATA"], "MS_scDRS_stats.csv"))
grp = pd.read_csv(os.path.join(CONFIG["DATA"], "PASS_Multiple_sclerosis.scdrs_group.cell_type"),
                  sep=None, engine="python", index_col="group")


def cliffs(x, y):
    gt = sum((a > b) + 0.5 * (a == b) for a in x for b in y)
    return 2.0 * gt / (len(x) * len(y)) - 1.0


def hedges_g(x, y):
    nx, ny = len(x), len(y)
    sp = np.sqrt(((nx - 1) * np.var(x, ddof=1) + (ny - 1) * np.var(y, ddof=1)) / (nx + ny - 2))
    d = (np.mean(x) - np.mean(y)) / sp
    J = 1 - 3 / (4 * (nx + ny) - 9)
    return J * d


ref = df.loc[df.cell_type == "COP/NFOL", "norm_score"].values
rows = []
for ct in order:
    v = df.loc[df.cell_type == ct, "norm_score"].values
    s = st[st.cell_type == ct].iloc[0]
    p = float(s["MWU_p_vs_lowest_median"]) if ct != "COP/NFOL" else np.nan
    rows.append(dict(cell_type=ct, n=len(v), mean=v.mean(), median=np.median(v),
                     sem=float(s["sem"]), whisker_lo5=np.percentile(v, 5), whisker_hi95=np.percentile(v, 95),
                     minv=v.min(), maxv=v.max(),
                     mwu_p_vs_COPNFOL=p,
                     group_P=float(grp.loc[ct, "assoc_mcp"]),
                     cliffs_delta=(0.0 if ct == "COP/NFOL" else cliffs(v, ref)),
                     scdrs_n_fdr_005=float(grp.loc[ct, "n_fdr_0.05"]),
                     scdrs_n_fdr_01=float(grp.loc[ct, "n_fdr_0.1"]),
                     scdrs_n_fdr_02=float(grp.loc[ct, "n_fdr_0.2"])))
A = pd.DataFrame(rows)
# BH q: family = the 5 mature-subtype contrasts (reference COP/NFOL), applied once to MWU P and once to group P
def bh(p):
    """BH q: rank ascending -> multiply by m/r -> cumulative min taken from the largest -> put back (the family counts only non-NaN members)"""
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    idx = np.where(~np.isnan(p))[0]
    m = len(idx)
    if m == 0:
        return q
    o = idx[np.argsort(p[idx])]
    qs = p[o] * m / np.arange(1, m + 1)
    qs = np.minimum.accumulate(qs[::-1])[::-1]
    q[o] = np.minimum(qs, 1.0)
    return q


A["mwu_q_BH5"] = bh(A.mwu_p_vs_COPNFOL.values)  # the MWU column's reference row is NaN by design
# family = the 5 mature-subtype contrasts (the COP/NFOL reference does not enter the correction)
_gp = A.group_P.where(A.cell_type != "COP/NFOL")
A["group_q_BH5"] = bh(_gp.values)
A.to_csv(OUTD + "/fig6A_stats.csv", index=False)
flat["A_n_total_joined"] = len(df)
flat["A_n_h5ad"] = len(adata)
flat["A_mwu_q_max"] = np.nanmax(A.mwu_q_BH5.values)
flat["A_group_q_min"] = np.nanmin(A.group_q_BH5.values)
flat["A_scdrs_fdr_pass_any"] = int(A[["scdrs_n_fdr_005", "scdrs_n_fdr_01", "scdrs_n_fdr_02"]].sum().sum())

# ---------------- C ----------------
T = pd.read_csv(T4D + "/T4_delta1_74traits_both_channels.csv")
T4b = pd.read_csv(T4D + "/T4b_cross_channel_rho_variants.csv")
T4c = pd.read_csv(T4D + "/T4c_level_effect_control.csv")
T6 = pd.read_csv(T4D + "/T6_cross_channel_summary.csv").set_index("item")["value"]
lvl = pd.read_csv(T4D + "/T4_channel_dispersion.csv") if os.path.exists(T4D + "/T4_channel_dispersion.csv") else None
ms = T[T.trait.str.contains("Multiple_sclerosis", regex=False)].iloc[0]
C = dict(
    rho_trait_level=float(T4c.rho_trait_means_cross_channel.iloc[0]),
    P_trait_level=float(T4c.P_trait_means.iloc[0]),
    rho_delta1_raw=float(T4c.rho_delta1_raw.iloc[0]),
    rho_delta1_level_removed=float(T4c.rho_delta1_after_removing_trait_level.iloc[0]),
    perm_P_delta1=float(T4b.loc[T4b.variant == "raw", "P_perm"].iloc[0]),
    noise_ceiling=float(T6["noise ceiling sqrt(relRNA*relATAC)"]),
    attenuation_corrected_rho=float(T6["attenuation-corrected rho"]),
    independent_halves_rho=float(T6["independent-nucleus-halves rho (median)"]),
    MS_offdiag_z=float(T6["MS off-diagonal residual z"]),
    MS_offdiag_rank=float(T6["MS off-diagonal rank of 74"]),
    MS_rank_rna=int(ms.rank_rna_raw), MS_rank_atac=int(ms.rank_atac_raw),
    MS_d1_rna=float(ms.d1_rna), MS_d1_atac=float(ms.d1_atac),
    n_traits=len(T),
    n_cell_mat=int(ms.n_cell_mat), n_cell_cop=int(ms.n_cell_cop),
    top5_rna="; ".join(T.sort_values("d1_rna", ascending=False).head(5).trait.str.replace("PASS_", "", regex=False)),
    top5_atac="; ".join(T.sort_values("d1_atac", ascending=False).head(5).trait.str.replace("PASS_", "", regex=False)),
    rho_sens_twopole=float(T4b.loc[T4b.variant.str.startswith("sens"), "rho_cross_channel"].iloc[0]),
)
pd.DataFrame([C]).to_csv(OUTD + "/fig6C_stats.csv", index=False)
flat.update({"C_" + k: v for k, v in C.items() if not isinstance(v, str)})

# C-i plotting data: per-trait "absolute level" (each channel averaged over all same nuclei) - the same convention as T4c
import glob


def load_score(root):
    o = {}
    for f in glob.glob(root + "/*/*.score.gz"):
        tr_ = os.path.basename(f).replace(".score.gz", "")
        if tr_ not in o:
            o[tr_] = pd.read_csv(f, sep="\t", index_col=0)["norm_score"]
    return pd.DataFrame(o)


CL = pd.read_csv(T4D + "/cell_level_table.csv")
Sr = load_score(os.path.join(FY, "oligo/scdrs_traits_draft/full74")).reindex(CL.cell_id)
Sa = load_score(os.path.join(FY, "oligo/fig6_rerun_20261019/dualomics_probe/out_gs74")).reindex(CL.cell_id)
trl = list(T.trait)
LV = pd.DataFrame({"trait": trl,
                   "level_rna": Sr.reindex(columns=trl).mean().values,
                   "level_atac": Sa.reindex(columns=trl).mean().values})
LV = LV.merge(T[["trait", "d1_rna", "d1_atac", "rank_rna_raw", "rank_atac_raw", "offdiag_z", "offdiag_rank"]], on="trait")
rho_lv = spearmanr(LV.level_rna, LV.level_atac)
LV.to_csv(OUTD + "/fig6C_trait_levels.csv", index=False)
C["rho_trait_level_reproduced"] = float(rho_lv.correlation)
C["P_trait_level_reproduced"] = float(rho_lv.pvalue)
C["n_nuclei_for_level"] = int(len(CL))
pd.DataFrame([C]).to_csv(OUTD + "/fig6C_stats.csv", index=False)
flat["C_rho_trait_level_reproduced"] = float(rho_lv.correlation)
flat["C_n_nuclei_for_level"] = len(CL)

# ---------------- D ----------------
d = pd.read_csv(M2B + "/M2b_per_mouse.csv"); dq = d[d.qc_pass].copy()
e = dq[dq.group != "Ctrl"]; cc = dq[dq.group == "Ctrl"]
rows = []
for ch, col in [("RNA", "rna_axis_contrast"), ("ATAC", "atac_axis_contrast")]:
    u = mannwhitneyu(e[col], cc[col], alternative="two-sided")
    rows.append(dict(channel=ch, median_EAE=float(np.median(e[col])), median_Ctrl=float(np.median(cc[col])),
                     delta_EAE_minus_Ctrl=float(np.median(e[col]) - np.median(cc[col])),
                     mwu_P=float(u.pvalue), cliffs_delta=float(cliffs(e[col].values, cc[col].values)),
                     hedges_g=float(hedges_g(e[col].values, cc[col].values)),
                     n_EAE=len(e), n_Ctrl=len(cc)))
D = pd.DataFrame(rows)
# cross-channel rho + CI + permutation (with n_perm/seed; the figure reads only here)
rho, _ = spearmanr(dq.atac_axis_contrast, dq.rna_axis_contrast)
n = len(dq)
z = np.arctanh(rho); se = 1 / np.sqrt(n - 3)
ci_lo, ci_hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
rng = np.random.default_rng(SEED)
bs = []
xy = np.c_[dq.atac_axis_contrast.values, dq.rna_axis_contrast.values]
for _ in range(NPERM):
    s = rng.integers(0, n, n)
    if len(np.unique(xy[s, 0])) > 2:
        bs.append(spearmanr(xy[s, 0], xy[s, 1]).correlation)
bs = np.array(bs, float)
boot_lo, boot_hi = np.nanpercentile(bs, 2.5), np.nanpercentile(bs, 97.5)
null = np.array([spearmanr(xy[rng.permutation(n), 0], xy[:, 1]).correlation for _ in range(NPERM)])
perm_p = float((np.sum(np.abs(null) >= abs(rho)) + 1) / (len(null) + 1))
X = pd.read_csv(M2B + "/M2b_crosschannel.csv").iloc[0]
extra = dict(rho=rho, n_mice=n, perm_P_new=perm_p, n_perm=NPERM, seed=SEED,
             perm_P_delivered=float(X.perm_P), rho_delivered=float(X.spearman_rho),
             ci_fisher_lo=ci_lo, ci_fisher_hi=ci_hi, ci_boot_lo=boot_lo, ci_boot_hi=boot_hi,
             auc_min=float(dq.auc_test.min()), auc_max=float(dq.auc_test.max()),
             auc_perm_P_mode=float(dq.p_perm.mode()), n_pass=int(dq.qc_pass.sum()), n_total=int(len(d)),
             n_dropped=int(len(d) - len(dq)), phase_trend_rho=None)
# staging trend (CE8: descriptive only, not a trend)
ph = {"Ctrl": 0, "Early": 1, "Peak": 2, "Late": 3}
dq["phase_ord"] = dq.group.map(ph)
tr = {("phase_trend_rho_" + c): float(spearmanr(dq.phase_ord, dq[c]).correlation) for c in ["rna_axis_contrast", "atac_axis_contrast"]}
tr.update({("phase_trend_P_" + c): float(spearmanr(dq.phase_ord, dq[c]).pvalue) for c in ["rna_axis_contrast", "atac_axis_contrast"]})
pd.DataFrame([{**extra, **tr}]).to_csv(OUTD + "/fig6D_crosschannel.csv", index=False)
D.to_csv(OUTD + "/fig6D_stats.csv", index=False)
dq.to_csv(OUTD + "/fig6D_per_mouse.csv", index=False)
flat.update({"D_" + k: v for k, v in extra.items() if not isinstance(v, str)})
flat.update({"D_" + k: v for k, v in tr.items()})
for _, r in D.iterrows():
    flat["D_%s_delta" % r.channel] = r.delta_EAE_minus_Ctrl
    flat["D_%s_mwuP" % r.channel] = r.mwu_P
    flat["D_%s_cliff" % r.channel] = r.cliffs_delta
    flat["D_%s_hedges" % r.channel] = r.hedges_g

pd.DataFrame(sorted(flat.items()), columns=["key", "value"]).to_csv(OUTD + "/fig6_all_numbers.csv", index=False)
print("A total nuclei=%d (h5ad %d) ; MWU q max=%.3f ; group q min=%.3f ; scDRS internal FDR passes=%d"
      % (len(df), len(adata), flat["A_mwu_q_max"], flat["A_group_q_min"], flat["A_scdrs_fdr_pass_any"]))
print("C level rho=%.3f (P=%.2f) | contrast rho=%.3f | level-removed %.3f | ceiling %.3f | MS offdiag z=%.2f rank %d"
      % (C["rho_trait_level"], C["P_trait_level"], C["rho_delta1_raw"], C["rho_delta1_level_removed"],
         C["noise_ceiling"], C["MS_offdiag_z"], C["MS_offdiag_rank"]))
print(D.to_string(index=False))
print("D cross-channel rho=%.4f  95%%CI(Fisher) %.2f-%.2f  CI(boot) %.2f-%.2f  perm P=%.4f (n_perm=%d, seed=%d)  [delivered table %.3f]"
      % (rho, ci_lo, ci_hi, boot_lo, boot_hi, perm_p, NPERM, SEED, X.perm_P))
print("D AUC %.3f-%.3f, %d/%d pass the bar; staging trend RNA rho=%.2f P=%.2f / ATAC rho=%.2f P=%.2f"
      % (extra["auc_min"], extra["auc_max"], extra["n_pass"], extra["n_total"],
         tr["phase_trend_rho_rna_axis_contrast"], tr["phase_trend_P_rna_axis_contrast"],
         tr["phase_trend_rho_atac_axis_contrast"], tr["phase_trend_P_atac_axis_contrast"]))
print("table dir:", OUTD)
