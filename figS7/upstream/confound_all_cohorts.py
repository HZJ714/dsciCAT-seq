#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""confound_all_cohorts.py — regulatory-domain-size confound audit of the 74-trait d1 specificity across the four human cohorts (raw vs size-stratified vs size-matched controls).

A8: audit of the **peak-set size confound** in the 74-trait d1 specificity across the four human
cohorts (cross-channel, with n_peak-matched controls). d1 correlates with the peak count of the
+/-100 kb chromatin domain (A2 rho=0.84, A3 0.69, A4 ~0.7, chromatin 0.71; all P<1e-11), so claims
like "MS ranks 2/74" must be shown not to proxy gene-set regulatory-domain size. Method:
  1) per cohort regress median_diff on log(n_peak) (or + n_genes_hit) and re-rank the residuals;
  2) matched control: read MS's position only inside the trait subset whose n_peak falls within
     +/-15% of MS (model-free);  3) write tables/seven4_crosschannel_confound_table.csv."""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
import os
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

B = WORK + "/figures/routeA_out"
T = pd.read_csv(B + "/A6_GSE231584_atac_v1/tables/trait_chromatin_sizeMatched.csv")[["trait", "n_peak_raw", "median_groupdiff"]]
a2 = pd.read_csv(B + "/A2_GSE277430_v1/tables/trait_sweep_74.csv")[["trait", "median_diff", "n_genes_hit"]].rename(columns={"median_diff": "A2"})
a3 = pd.read_csv(B + "/A3_GSE118257_v1/tables/trait_sweep_74_GSE118257.csv")[["trait", "median_diff"]].rename(columns={"median_diff": "A3"})
a4 = pd.read_csv(B + "/A4_GSE301908_v1/tables/trait_sweep_74_GSE301908.csv")[["trait", "median_diff"]].rename(columns={"median_diff": "A4"})
# A7 is **donor-level** (8 donors; several lesion domains of one donor are not independent, so the sample-level 14/14 P would be optimistic - Fig 7 receipt section 29)
a7 = pd.read_csv(B + "/A7_GSE231585_rna_v1/tables/trait_sweep_74_GSE231585_donorLevel.csv")[["trait", "median_diff_donor"]].rename(columns={"median_diff_donor": "A7"})
m = T.merge(a2, on="trait").merge(a3, on="trait").merge(a4, on="trait", how="left").merge(a7, on="trait", how="left")
m["lnp"] = np.log(m.n_peak_raw)
m = m.rename(columns={"median_groupdiff": "chr"})
CH = ["A2", "A3", "A4", "A7", "chr"]
MS = "PASS_Multiple_sclerosis"
ims = m.index[m.trait == MS][0]


def rank_desc(x):
    return (-pd.Series(x)).rank(method="min").astype(int).values


def resid(y, cols, df):
    Xd = sm.add_constant(df[cols].astype(float))
    return sm.OLS(pd.Series(y, dtype=float), Xd).fit().resid.values


print("n=%d traits; MS n_peak=%d (cohort median %d), genes hit %d"
      % (len(m), m.n_peak_raw[ims], int(m.n_peak_raw.median()), m.n_genes_hit[ims]))
print("\nrho(log n_peak, d1 value):", {c: round(spearmanr(m[c].dropna(), m.loc[m[c].notna(), "lnp"]).statistic, 3) for c in CH})
print("rho(log n_peak, d1 rank):", {c: round(spearmanr(m[c].dropna().rank(ascending=False), m.loc[m[c].notna(), "lnp"].rank(ascending=False)).statistic, 3) for c in CH})
rows = []
for c in CH:
    d = m[m[c].notna()].copy()
    raw = rank_desc(d[c])
    r_lnp = rank_desc(resid(d[c].values, ["lnp"], d))
    r_both = rank_desc(resid(d[c].values, ["lnp", "n_genes_hit"], d))
    idx = d.index.get_loc(ims)
    # matched control (model-free)
    win = m[(m.n_peak_raw >= m.n_peak_raw[ims] * .85) & (m.n_peak_raw <= m.n_peak_raw[ims] * 1.15) & m[c].notna()]
    rows.append({"channel": c, "n_trait": len(d),
                 "MS_rank_raw": int(raw[idx]), "MS_rank_resid_lnp": int(r_lnp[idx]),
                 "MS_rank_resid_lnp+hit": int(r_both[idx]),
                 "MS_rank_in_matched_window": int(rank_desc(win[c].values)[win.index.get_loc(ims)]),
                 "n_in_window": len(win),
                 "rho_rank_vs_lnp": round(spearmanr(raw, d.lnp.values).statistic, 3),
                 "rho_resid_vs_lnp": round(spearmanr(r_lnp, d.lnp.values).statistic, 3)})
R = pd.DataFrame(rows)
print("\n=== MS position per channel (raw rank vs rank after removing the domain-size confound) ===")
print(R.to_string(index=False))

for c in CH:
    m["rank_raw_" + c] = np.nan
    m["rank_adj_" + c] = np.nan
d = m[m.A2.notna()]
for c in CH:
    dd = m[m[c].notna()]
    m.loc[dd.index, "rank_raw_" + c] = rank_desc(dd[c].values)
    m.loc[dd.index, "rank_adj_" + c] = rank_desc(resid(dd[c].values, ["lnp"], dd))
m.sort_values("rank_adj_A2").to_csv(B + "/A6_GSE231584_atac_v1/tables/seven4_crosschannel_confound_table.csv", index=False)
m.sort_values("rank_adj_A2").to_csv(B + "/FIG7_v1/tables/seven4_confound_ALL_cohorts.csv", index=False)
print("\n=== top 15 of the A2 channel after de-confounding ===")
print(m.sort_values("rank_adj_A2").head(15)[["trait", "n_peak_raw", "A2", "rank_raw_A2", "rank_adj_A2", "rank_adj_A3", "rank_adj_chr"]].round(3).to_string(index=False))
print("\n=== MS row ===")
print(m[m.trait == MS][["trait", "n_peak_raw", "A2", "A3", "A4", "A7", "chr"] + [c for c in m.columns if c.startswith("rank_raw_")] + [c for c in m.columns if c.startswith("rank_adj_")]].round(3).to_string(index=False))
print("\n=== raw / adjusted rank of MS per channel + genes hit ===")
print(R.to_string(index=False))
print("\ncross-channel rank correlation (after de-confounding):")
for x, y in [("A2", "A3"), ("A2", "chr"), ("A3", "chr"), ("A2", "A4"), ("A2", "A7"), ("A3", "A4"), ("A3", "A7"), ("A4", "A7"), ("A4", "chr"), ("A7", "chr")]:
    dd = m[m[x].notna() & m[y].notna()]
    print("  %s-%s raw %.3f | both residualised %.3f P=%.2g (n=%d)"
          % (x, y, spearmanr(dd["rank_raw_" + x], dd["rank_raw_" + y]).statistic,
             spearmanr(dd["rank_adj_" + x], dd["rank_adj_" + y]).statistic,
             spearmanr(dd["rank_adj_" + x], dd["rank_adj_" + y]).pvalue, len(dd)))
