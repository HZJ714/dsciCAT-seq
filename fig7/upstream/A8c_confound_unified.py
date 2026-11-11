#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A8c: regulatory-domain-size confounding check for the 74-trait d1, unified-basis version
(author decision 11-05: every number on the main figure must come from one self-consistent basis).

Background: main-figure B (rank / rho strip) and C (three axes) used to read two different d1 -
B from seven4 (A4 scored on the OL-lineage object, A7 at donor level) and C from the 11-04
whole-lineage three-axis tables - so one cohort could show two d1 side by side (A4: 0.963/#5 vs
0.707/#7). This script moves every derived statistic of B onto the tables C already uses: the
method is untouched, only the inputs are swapped, and seven4 is then rewritten isomorphically.

Inputs (exactly what main-figure C uses):
  A2: trait_sweep_dual_axis.csv            d1_OLmat           (whole lineage, 154,522 nuclei)
  A3: trait_sweep_delta2_GSE118257_74.csv  d1_OLIGO_minus_OPC (whole-lineage draw, 9 donors)
  A4: trait_sweep_delta2_GSE301908_74.csv  d1                 (whole-lineage draw, <= 12k per lineage, 11 donors)
  A7: trait_sweep_delta2_GSE231585_74.csv  d1                 (whole lineage, 62,911 nuclei, 15 samples)
  chr: A6 trait_chromatin_sizeMatched.csv  median_groupdiff   (A6 unchanged)

Method (word for word as A8b): rank_desc, OLS residualisation of d1 on ln(n_peak_raw), and the
MS +-15% matched window.

Output: FIG7_v1/tables/seven4_confound_ALL_cohorts.csv (written in place; the previous table is
archived first as *_20261030_OLobject.csv).
"""
import shutil
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
OLD = B + "/FIG7_v1/tables/seven4_confound_ALL_cohorts.csv"
ARCHIVE = B + "/FIG7_v1/tables/seven4_confound_ALL_cohorts_20261030_OLobject.csv"
shutil.copy(OLD, ARCHIVE)
print("archiving the old table ->", ARCHIVE)

T = pd.read_csv(B + "/A6_GSE231584_atac_v1/tables/trait_chromatin_sizeMatched.csv")[["trait", "n_peak_raw", "median_groupdiff"]]
a2 = pd.read_csv(B + "/A2_GSE277430_v1/tables/trait_sweep_dual_axis.csv")[["trait", "d1_OLmat"]].rename(columns={"d1_OLmat": "A2"})
a3 = pd.read_csv(B + "/A3_GSE118257_v1/tables/trait_sweep_delta2_GSE118257_74.csv")[["trait", "d1_OLIGO_minus_OPC"]].rename(columns={"d1_OLIGO_minus_OPC": "A3"})
a4 = pd.read_csv(B + "/A4_GSE301908_v1/tables/trait_sweep_delta2_GSE301908_74.csv")[["trait", "d1", "n_genes_hit"]].rename(columns={"d1": "A4"})
a7 = pd.read_csv(B + "/A7_GSE231585_rna_v1/tables/trait_sweep_delta2_GSE231585_74.csv")[["trait", "d1"]].rename(columns={"d1": "A7"})
m = T.rename(columns={"median_groupdiff": "chr"}).merge(a2, on="trait").merge(a3, on="trait").merge(a4, on="trait").merge(a7, on="trait")
m["lnp"] = np.log(m.n_peak_raw)
CH = ["A2", "A3", "A4", "A7", "chr"]
MS = "PASS_Multiple_sclerosis"
ims = m.index[m.trait == MS][0]


def rank_desc(x):
    return (-pd.Series(x)).rank(method="min").astype(int).values


def resid(y, cols, df):
    Xd = sm.add_constant(df[cols].astype(float))
    return sm.OLS(pd.Series(y, dtype=float), Xd).fit().resid.values


print("n=%d traits; the MS n_peak=%d (overall median %d)"
      % (len(m), m.n_peak_raw[ims], int(m.n_peak_raw.median())))
print("rho(log n_peak, d1 value): ", {c: round(spearmanr(m[c].dropna(), m.loc[m[c].notna(), "lnp"]).statistic, 3) for c in CH})
rows = []
for c in CH:
    d = m[m[c].notna()].copy()
    raw = rank_desc(d[c])
    r_lnp = rank_desc(resid(d[c].values, ["lnp"], d))
    idx = d.index.get_loc(ims)
    win = m[(m.n_peak_raw >= m.n_peak_raw[ims] * .85) & (m.n_peak_raw <= m.n_peak_raw[ims] * 1.15) & m[c].notna()]
    rows.append({"channel": c, "n_trait": len(d),
                 "MS_rank_raw": int(raw[idx]), "MS_rank_resid_lnp": int(r_lnp[idx]),
                 "MS_rank_in_matched_window": int(rank_desc(win[c].values)[win.index.get_loc(ims)]),
                 "n_in_window": len(win),
                 "rho_rank_vs_lnp": round(spearmanr(raw, d.lnp.values).statistic, 3),
                 "rho_resid_vs_lnp": round(spearmanr(r_lnp, d.lnp.values).statistic, 3)})
print("\n=== MS position in each channel (unified basis) ===")
print(pd.DataFrame(rows).to_string(index=False))

for c in CH:
    m["rank_raw_" + c] = np.nan
    m["rank_adj_" + c] = np.nan
for c in CH:
    dd = m[m[c].notna()]
    m.loc[dd.index, "rank_raw_" + c] = rank_desc(dd[c].values)
    m.loc[dd.index, "rank_adj_" + c] = rank_desc(resid(dd[c].values, ["lnp"], dd))
out = m.sort_values("rank_adj_A2")[["trait", "n_peak_raw", "chr", "A2", "n_genes_hit", "A3", "A4", "A7", "lnp"]
                                   + ["rank_raw_%s" % c for c in CH] + ["rank_adj_%s" % c for c in CH]]
out.to_csv(OLD, index=False)
print("\n=== MS row (unified basis) ===")
print(out[out.trait == MS].round(4).to_string(index=False))
print("\n=== cross-channel rank correlation (after removing the size confound, unified basis) ===")
for x, y in [("A2", "A3"), ("A2", "chr"), ("A3", "chr"), ("A2", "A4"), ("A2", "A7"), ("A3", "A4"), ("A3", "A7"), ("A4", "A7"), ("A4", "chr"), ("A7", "chr")]:
    dd = out[out[x].notna() & out[y].notna()]
    print("  %s-%s raw %.3f | both sides residualised %.3f P=%.2g (n=%d)"
          % (x, y, spearmanr(dd["rank_raw_" + x], dd["rank_raw_" + y]).statistic,
             spearmanr(dd["rank_adj_" + x], dd["rank_adj_" + y]).statistic,
             spearmanr(dd["rank_adj_" + x], dd["rank_adj_" + y]).pvalue, len(dd)))
