# -*- coding: utf-8 -*-
"""matched_window_table.py — model-free matched n_peak window (no regression): MS rank inside the window, per cohort.

(A8_confound_matched_rank.py printed this window but never wrote it to disk). This script applies
**the same rule** - rank by descending d1 inside the +/-15% n_peak window around MS - to the
delivered seven4_confound_ALL_cohorts.csv, and carries D3 rho(ln n_peak, d1) in, so the face has a table."""
import sys, os

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
import os
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
sys.path.insert(0, WORK + "/scripts/routeA")
from support_common import *
import numpy as np, pandas as pd

S4 = read("FIG7_v1/tables/seven4_confound_ALL_cohorts.csv")
D3 = read("D_crosschemistry/tables/D3_crossconcordance_domainadjust.csv")
MS = "PASS_Multiple_sclerosis"
NAME = {"A2": "GSE277430", "A3": "GSE118257", "A4": "GSE301908", "A7": "GSE231585"}


def rho_from_d3(tag):
    r = D3[D3["pair"].astype(str).str.contains(tag, regex=False)]
    return float(r.iloc[0, list(D3.columns).index("rho_raw")]) if len(r) else np.nan


ms_np = float(S4.loc[S4.trait == MS, "n_peak_raw"].iloc[0])
lo, hi = ms_np * 0.85, ms_np * 1.15
rows = []
for c, nm in NAME.items():
    d = S4[S4[c].notna()]
    win = d[(d.n_peak_raw >= lo) & (d.n_peak_raw <= hi)]
    rk = (-win[c].astype(float)).rank(method="min").astype(int)
    rows.append(dict(cohort=c, study=nm, n_traits=len(d), MS_n_peak=ms_np,
                     window_lo=lo, window_hi=hi, n_in_window=len(win),
                     MS_rank_in_window=int(rk[win.index[win.trait == MS][0]]),
                     rho_lnp_vs_d1=rho_from_d3("%s: Spearman" % nm),
                     rank_raw=int(S4.loc[S4.trait == MS, "rank_raw_" + c].iloc[0]),
                     rank_adj=int(S4.loc[S4.trait == MS, "rank_adj_" + c].iloc[0]),
                     source="FIG7_v1/tables/seven4_confound_ALL_cohorts.csv + D3 confound-check rows"))
out = pd.DataFrame(rows)
out.to_csv(os.path.join(TABDIR, "fig7Sb_matched_window.csv"), index=False)
say(out.to_string(index=False))
