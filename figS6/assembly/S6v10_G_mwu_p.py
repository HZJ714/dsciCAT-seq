#!/usr/bin/env python
"""S6 v10 panel G: two-sided Mann-Whitney U test of OPC vs Oligo per-nucleus norm_score.

Author 2026-11-05 instruction: the ±Δ1 in each G cell's lower-right corner becomes a P value; the caption states the method.
Δ1 itself is not deleted; it retreats to Table S6 (which already carries the MS rows of delta1_74traits_per_tech).

The convention matches make_S6_v10_panels.py's ms_scores() exactly:
  RNA  = the norm_score column of scdrs/<tech>/*/PASS_Multiple_sclerosis.full_score.gz
  ATAC = the same column of scdrs_log1p_20261104/<tech>/*/... (the log1p-corrected lineage, 2026-11-04)
  subtype labels = the label column of tables_log1p_20261104/cells_<tech>.csv, intersection only
Test: scipy.stats.mannwhitneyu(two-sided, default method='auto'), no permutation, no seed (deterministic).
Output tables_v8/figS6v10_G_opc_oligo_mwu.csv
"""
import glob
import os
import pandas as pd
from scipy.stats import mannwhitneyu

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

R = WORK + "/figures/routeA_out"
V8 = R + "/../routeA_out/SUPPORT_v2/tables_v8"
V8 = WORK + "/figures/routeA_out/SUPPORT_v2/tables_v8"
CELLF = {"dsciCAT-seq": "cells_dsciCAT_seq.csv", "10x Multiome": "cells_10x.csv",
         "HT-scCAT-seq": "cells_HT_scCAT_seq_HT.csv"}
DIRF = {"dsciCAT-seq": ("dsciCAT_seq_rna", "dsciCAT_seq_atac"),
        "10x Multiome": ("10x_rna", "10x_atac"),
        "HT-scCAT-seq": ("HT_scCAT_seq_HT_rna", "HT_scCAT_seq_HT_atac")}


def scores(tech, chan):
    rn, an = DIRF[tech]
    base = (R + "/D_crosschemistry/scdrs/" + rn if chan == "rna"
            else R + "/D_crosschemistry/scdrs_log1p_20261104/" + an)
    fs = sorted(glob.glob(base + "/*/PASS_Multiple_sclerosis.full_score.gz"))
    df = pd.concat([pd.read_csv(f, sep="\t", index_col=0, low_memory=False) for f in fs])
    lb = pd.read_csv(R + "/D_crosschemistry/tables_log1p_20261104/" + CELLF[tech]).set_index("cell")["label"]
    df = df[df.index.isin(lb.index)].copy()
    df["label"] = lb.loc[df.index].values
    return df


def fig_p(p):
    """The on-face P format (goes into the ledger as p_fig; the panel prints this exact string so the reverse gate matches it character by character):
    P < 1e-3 -> two-significant-digit scientific notation; 1e-3 <= P < 0.01 -> four decimals; otherwise two decimals."""
    if p < 1e-3:
        return "%.1e" % p
    if p < 0.01:
        return "%.4f" % p
    return "%.2f" % p


rows = []
for tech in CELLF:
    for chan in ["rna", "atac"]:
        df = scores(tech, chan)
        op = df.loc[df.label == "OPC", "norm_score"].values
        mo = df.loc[df.label == "Oligo", "norm_score"].values
        u, p = mannwhitneyu(mo, op, alternative="two-sided")
        rows.append(dict(tech=tech, channel=chan, n_opc=len(op), n_oligo=len(mo),
                         median_opc=round(float(pd.Series(op).median()), 4),
                         median_oligo=round(float(pd.Series(mo).median()), 4),
                         median_diff=round(float(pd.Series(mo).median() - pd.Series(op).median()), 4),
                         mwu_U=float(u), p_two_sided=float(p), p_fig=fig_p(p),
                         p_label="P = " + fig_p(p),
                         method="two-sided Mann-Whitney U on per-cell scDRS norm_score"))
        print("%-14s %-5s  OPC n=%5d med=%+.3f | Oligo n=%5d med=%+.3f | diff=%+.3f  U=%.3g  P=%.3g"
              % (tech, chan, len(op), pd.Series(op).median(), len(mo), pd.Series(mo).median(),
                 pd.Series(mo).median() - pd.Series(op).median(), u, p))
out = V8 + "/figS6v10_G_opc_oligo_mwu.csv"
pd.DataFrame(rows).to_csv(out, index=False)
print("[written]", out)
