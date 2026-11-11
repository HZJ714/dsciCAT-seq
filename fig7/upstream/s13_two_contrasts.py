#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s13_two_contrasts.py -- the two contrasts Fig7D should be built on.

Reference Fig7A/C define one axis: scDRS delta1 = median(mature OL) - median(OPC), and identity axes
OL-maturation vs astro-immune vs neuron-glia. Fig7D therefore asks, for each curated gene, the SAME
two questions in detection-frequency terms, within each cohort:

  C1 lineage vs rest : OL lineage (OL + iOL + OPC)  vs  non-lineage (AST MG LYM VASC NEU)
  C2 within lineage  : mature OL vs OPC                     (the delta1 analogue)

Both are one-sided Fisher on detected/not-detected nuclei counts, BH across the 65 genes within
cohort and within direction, so an enrichment and a depletion never share a q-value ladder.
A2 (GSE277430) is an OL-subset release: it has no non-lineage cells, so C1 is undefined there and
only C2 is drawn. A3's ImOlGs/COPs are folded into iOL exactly as in s10b.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# Root of the author's analysis tree (figures/, routeA_out/, routeA_probe/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
B = WORK + "/figures/SUPPL_ANALYSES_20261112"
LIN, NON = ["OL", "iOL", "OPC"], ["AST", "MG", "LYM", "VASC", "NEU"]
FULL = ["GSE118257", "GSE301908", "GSE231585"]
MAP = {
 "GSE118257": {**{f"Oligo{i}": "OL" for i in range(1, 7)}, "ImOlGs": "iOL", "COPs": "iOL", "OPCs": "OPC",
               "Astro": "AST", "Micro/Mac": "MG", "Immune": "LYM", "Endo": "VASC", "Pericyte/SMC": "VASC", "Neuron": "NEU"},
 "GSE301908": {"Oligo": "OL", "OPC": "OPC", "Astro": "AST", "Micro": "MG", "Lym": "LYM", "Endo": "VASC", "Neuron": "NEU"},
 "GSE231585": {"OLIGO": "OL", "OPC": "OPC", "AST": "AST", "MG": "MG", "LYM": "LYM", "END": "VASC", "NEU": "NEU", "Unclear": "QC"},
 "GSE277430": {"OLIGO": "OL", "OPC": "OPC"}}


def cell(a, n):
    """detected / not-detected counts, clipped to be non-negative."""
    a = int(round(a)); n = int(round(n))
    return max(a, 0), max(n - a, 0)


def stat(dl, nl, dn, nn):
    """log2 odds ratio of detection in group1 vs group2 plus both one-sided p-values."""
    if nl <= 0 or nn <= 0:
        return np.nan, np.nan, np.nan
    d1, u1 = cell(dl, nl); d2, u2 = cell(dn, nn)
    if d1 + u1 + d2 + u2 == 0:
        return np.nan, np.nan, np.nan
    lor = np.log2(((d1 + .5) / max(u1, .5)) / ((d2 + .5) / max(u2, .5)))
    p_hi = fisher_exact([[d1, u1], [d2, u2]], alternative="greater")[1]
    p_lo = fisher_exact([[d1, u1], [d2, u2]], alternative="less")[1]
    return lor, p_hi, p_lo


d = pd.read_csv(B + "/tables/S10_all_cells.csv")
d["det"] = d.n * d.pct / 100.0
d["cp"] = [MAP[r.dataset].get(r.group, "?") for r in d.itertuples()]
rows = []
for ds, g in d.groupby("dataset"):
    for gene, gg in g.groupby("gene"):
        def ag(pops):
            x = gg[gg["cp"].isin(pops)]
            return x["n"].sum(), x["det"].sum()
        NL, DL = ag(LIN); NN, DN = ag(NON); NO, DO = ag(["OL"]); NP_, DP = ag(["OPC"])
        l1, h1, o1 = stat(DL, NL, DN, NN)
        l2, h2, o2 = stat(DO, NO, DP, NP_)
        rows.append([ds, gene, NL, round(100 * DL / NL, 3) if NL else np.nan, NN,
                     round(100 * DN / NN, 3) if NN else np.nan, NO, round(100 * DO / NO, 3) if NO else np.nan,
                     NP_, round(100 * DP / NP_, 3) if NP_ else np.nan, l1, h1, o1, l2, h2, o2])
r = pd.DataFrame(rows, columns=["ds", "gene", "N_lin", "pct_lin", "N_non", "pct_non", "N_OL", "pct_OL", "N_OPC", "pct_OPC",
                                "lor_lineage", "p_lin_hi", "p_lin_lo", "lor_OLvsOPC", "p_OL_hi", "p_OL_lo"])
for ds in r.ds.unique():
    m = r.ds == ds
    for pc, qc in [("p_lin_hi", "q_lin_hi"), ("p_lin_lo", "q_lin_lo"), ("p_OL_hi", "q_OL_hi"), ("p_OL_lo", "q_OL_lo")]:
        r.loc[m, qc] = multipletests(np.nan_to_num(r.loc[m, pc].values, nan=1.0), method="fdr_bh")[1]
# p_hi = "lineage/OL odds greater than comparator" -> enrichment; p_lo = depletion
r["LIN_ENR"] = (r.q_lin_hi < 0.05) & (r.lor_lineage > 1)
r["NON_ENR"] = (r.q_lin_lo < 0.05) & (r.lor_lineage < -1)
r["OL_ENR"] = (r.q_OL_hi < 0.05) & (r.lor_OLvsOPC > 1)
r["OPC_ENR"] = (r.q_OL_lo < 0.05) & (r.lor_OLvsOPC < -1)
r.to_csv(B + "/tables/S13_two_contrasts.csv", index=False)

s = r[r.ds.isin(FULL)]
ev = pd.read_csv(B + "/tables/S10_evidence_rows.csv").set_index("gene")
gr = pd.read_csv(B + "/tables/S12_target_curation_grades.csv").set_index("gene")
S7I = {"S1PR1", "S1PR5", "S1PR3", "S1PR4", "ITGA4", "CD52", "NFE2L2"}


def grade(g):
    v = gr.grade.get(g, None)
    if v in ("A", "B", "C", "D"):
        return "A" if (v == "D" and g in S7I) else v
    return "A" if g in S7I else "n/a"


summ = []
for g in s.gene.unique():
    x = s[s.gene == g]
    summ.append([g, int(x.LIN_ENR.sum()), int(x.NON_ENR.sum()), int(x.OL_ENR.sum()), int(x.OPC_ENR.sum()),
                 round(np.nanmax(x.lor_lineage), 2) if x.lor_lineage.notna().any() else np.nan,
                 round(np.nanmax(x.lor_OLvsOPC), 2) if x.lor_OLvsOPC.notna().any() else np.nan,
                 round(np.nanmin(x.lor_OLvsOPC), 2) if x.lor_OLvsOPC.notna().any() else np.nan,
                 ev.peak_pct.get(g, np.nan), ev.rv1_rank.get(g, np.nan), grade(g)])
S = pd.DataFrame(summ, columns=["gene", "coh_lineage_enr", "coh_nonlineage_enr", "coh_OLgtOPC", "coh_OPCgtOL",
                                "max_lor_lineage", "max_lor_OLgtOPC", "min_lor_OLvsOPC", "peak_pct", "rv1_rank", "grade"])
S.to_csv(B + "/tables/S13_gene_summary.csv", index=False)
pd.set_option("display.width", 250)
print("== OL-lineage-enriched in >=2/3 full-atlas cohorts ==")
print(S[S.coh_lineage_enr >= 2].sort_values("max_lor_lineage", ascending=False).to_string(index=False))
print("\n== mature-OL > OPC in >=2/3 ==")
print(S[S.coh_OLgtOPC >= 2].sort_values("max_lor_OLgtOPC", ascending=False).to_string(index=False))
print("\n== OPC > mature-OL in >=2/3 ==")
print(S[S.coh_OPCgtOL >= 2].sort_values("min_lor_OLvsOPC").to_string(index=False))
print("\n== neither contrast (flat / peripherally restricted) among A/B-grade targets ==")
print(S[(S.coh_lineage_enr == 0) & (S.coh_OLgtOPC == 0) & (S.coh_OPCgtOL == 0) & (S.grade.isin(["A", "B"]))]
      .sort_values("peak_pct", ascending=False).to_string(index=False))
