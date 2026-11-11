#!/usr/bin/env python3
"""export_mouse_same_nucleus_delta1.py — compact panel-E input: per-trait d1 for both channels of the same nuclei.

export_mouse_same_nucleus_delta1.py - generator of the panel-D input table (compact form).

Aggregates the per-nucleus scDRS scores (one *.scores.csv.gz per library, with a cellType column)
into "one row per trait, d1 for both channels", which the drawing script reads directly:
  delta1 = median(mature-OL nuclei) - median(precursor-OL nuclei)
Difference from S7_rna_crossaxis_mouse_full.csv / S7_atac_crossaxis_mouse_full.csv: this table
pairs **the same nuclei across the two channels** (same nuclei, same labels, same seed) and is what
panel D scatters; those two are three-axis tables (d1/d2A/d2B, all 74 traits x cell types).

usage:  S7_SCORES=/path/to/m3_probe S7_DATA=../data python export_mouse_same_nucleus_delta1.py"""
import os, glob, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

MOL = {"MOL1", "MOL2", "MOL56"}
OPC_ = {"OPC", "COP"}
MS = "PASS_Multiple_sclerosis"
SC = os.environ.get("S7_SCORES", WORK + "/routeA_probe/m3_probe")
OUT = os.environ.get("S7_DATA", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "SUPPORT_v1", "tables"))


def d1_series(S):
    a = S.loc[S.cellType.isin(MOL).values]
    b = S.loc[S.cellType.isin(OPC_).values]
    traits = [c for c in S.columns if c not in ("cellType", "model")]
    out = pd.Series(np.nan, index=traits)
    for tr in traits:
        x, y = a[tr].dropna(), b[tr].dropna()
        if min(len(x), len(y)) >= 10:
            out[tr] = x.median() - y.median()
    return out


def load(pattern):
    df = pd.concat([pd.read_csv(f, index_col=0)
                    for f in sorted(glob.glob(os.path.join(SC, pattern)))])
    return df[~df.index.duplicated()]


Sm6 = load("m6_scores/*.scores.csv.gz")     # ATAC channel
Sm7 = load("m7_scores/*.scores.csv.gz")     # RNA channel
assert set(Sm6.index) == set(Sm7.index), "M6/M7 nucleus sets differ (should be the same 14,928)"
E1 = d1_series(Sm6)                          # ATAC d1
E2 = d1_series(Sm7)                          # RNA d1
tab = pd.DataFrame({"trait": E1.index, "delta1_atac": E1.values, "delta1_rna": E2.values})
tab.to_csv(os.path.join(OUT, "S7_mouse_same_nucleus_delta1.csv"), index=False)
rho = spearmanr(E1, E2, nan_policy="omit")
info = {"n_nuclei": int(len(Sm6)), "n_libs": int(Sm6.index.str.split("|").str[0].nunique()),
        "n_traits": int(E1.notna().sum()),
        "spearman_rho": float(rho.statistic), "spearman_p": float(rho.pvalue),
        "ms_rank_rna": int((E2 > E2[MS]).sum()) + 1, "ms_rank_atac": int((E1 > E1[MS]).sum()) + 1,
        "seed_note": "same nuclei, same cellType labels, same gene sets in both channels"}
json.dump(info, open(os.path.join(OUT, "S7_mouse_same_nucleus_delta1.json"), "w"), indent=1)
print(json.dumps(info, indent=1))
