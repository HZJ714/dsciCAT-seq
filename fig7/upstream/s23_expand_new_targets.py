#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s23_expand_new_targets.py -- score the 33-row parallel-env target list with the SAME two
contrasts as s13, using ONLY local tables (no network): tables/S15_det_<ds>.csv and
S15_n_<ds>.csv hold detected / total nuclei counts for every gene of every released matrix.

C1 = OL lineage (OL+iOL+OPC) vs non-lineage (AST MG LYM VASC NEU)
C2 = mature OL vs OPC
one-sided Fisher in both directions; BH within cohort and within direction, across the union of
the previously curated 65 genes and the 38 atomic genes implied by the 33 new rows.
Outputs: tables/S23_new_targets_stats.csv, tables/S23_gene_universe_status.csv
"""
import os
import numpy as np, pandas as pd, re, os
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# Root of the author's analysis tree (figures/, routeA_out/, routeA_probe/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
B = WORK + "/figures/SUPPL_ANALYSES_20261112"; T = B + "/tables"
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
EXP_CSV = os.environ.get("FIG7_EXPANSION_CSV", os.path.join(REPO, "data", "EXPANSION_33", "expansion33_rows.csv"))
DS = ["GSE118257", "GSE301908", "GSE231585", "GSE277430"]
LAB = dict(zip(DS, ["A3", "A4", "A7", "A2"]))
LIN, NON = ["OL", "iOL", "OPC"], ["AST", "MG", "LYM", "VASC", "NEU"]

# --- gene lists -----------------------------------------------------------
OLD = pd.read_csv(T + "/S12_target_curation_audit.csv")
oldg = set(OLD.gene.astype(str))
NEWROWS = {}
# Vendored ASCII extract of the author's expansion memo (33 rows).  The memo itself is an internal
# Chinese-language document and is not pushed.  claim_cell is kept verbatim and THIS script's own
# atom expansion is applied to it: the universe here is not the ledger's (row 21 also yields HTR2B),
# so the two scripts must not share one pre-expanded atom list.
X = pd.read_csv(EXP_CSV)
assert len(X) == 33, len(X)
for r in X.itertuples():
    # explicit expansion of the receipt's shorthand cells; regex on "IMPDH1/2" loses IMPDH2
    EXP = {5:["IMPDH1","IMPDH2"], 14:["CD80","CD86"], 21:["SLC6A4","SLC6A2","HTR2B"],
           24:["SCN1A","SCN8A"], 28:["ESR1","ESR2"]}
    syms = EXP.get(int(r.exp_no), [x for x in re.split(r"[;,/] ", r.claim_cell) if re.match(r"^[A-Z][A-Z0-9-]{1,15}$", x)])
    if int(r.exp_no) == 19: syms = ["TNFSF13"]
    if int(r.exp_no) == 30: syms = ["IL7R"]
    NEWROWS[int(r.exp_no)] = dict(claim_cell=r.claim_cell, syms=syms, drug=r.drug, fda=r.fda_status_en,
                                  year_j=r.year_journal, title=r.title, doi=r.doi, point="")
newatoms = sorted({s for v in NEWROWS.values() for s in v["syms"]})
UNIVERSE = sorted(oldg | set(newatoms))
print("curated 65 | 33 rows -> %d atomic new symbols | union universe = %d" % (len(newatoms), len(UNIVERSE)))

# --- per cohort, per gene --------------------------------------------------
def stat(dl, nl, dn, nn):
    if nl <= 0 or nn <= 0: return (np.nan,) * 3
    d1, u1 = max(int(round(dl)), 0), max(int(round(nl - dl)), 0)
    d2, u2 = max(int(round(dn)), 0), max(int(round(nn - dn)), 0)
    if d1 + u1 + d2 + u2 == 0: return (np.nan,) * 3
    lor = np.log2(((d1 + .5) / max(u1, .5)) / ((d2 + .5) / max(u2, .5)))
    tab = [[d1, u1], [d2, u2]]
    return lor, fisher_exact(tab, alternative="greater")[1], fisher_exact(tab, alternative="less")[1]

FULL3 = ["GSE118257", "GSE301908", "GSE231585"]
rows, status = [], []
for ds in DS:
    det = pd.read_csv(f"{T}/S15_det_{ds}.csv", index_col=0).fillna(0)
    n = pd.read_csv(f"{T}/S15_n_{ds}.csv", index_col=0)["n"]
    universe_all = set(det.index)
    for g in UNIVERSE:
        present = g in universe_all
        status.append([ds, g, present])
        if not present: continue
        r = det.loc[g]
        if isinstance(r, pd.DataFrame): r = r.max()
        NL = sum(n[p] for p in LIN if p in n.index and p in r.index)
        DL = sum(r[p] for p in LIN if p in n.index and p in r.index)
        NN = sum(n[p] for p in NON if p in n.index and p in r.index)
        DN = sum(r[p] for p in NON if p in n.index and p in r.index)
        NO, DO = (n.get("OL", 0), r.get("OL", 0)) if "OL" in r.index else (0, 0)
        NP, DP = (n.get("OPC", 0), r.get("OPC", 0)) if "OPC" in r.index else (0, 0)
        l1, h1, o1 = stat(DL, NL, DN, NN) if ds in FULL3 else (np.nan,) * 3
        l2, h2, o2 = stat(DO, NO, DP, NP)
        rows.append([ds, g, NL, round(100 * DL / NL, 3) if NL else np.nan, NN,
                     round(100 * DN / NN, 3) if NN else np.nan, NO,
                     round(100 * DO / NO, 3) if NO else np.nan, NP,
                     round(100 * DP / NP, 3) if NP else np.nan, l1, h1, o1, l2, h2, o2])
R = pd.DataFrame(rows, columns=["ds", "gene", "N_lin", "pct_lin", "N_non", "pct_non", "N_OL", "pct_OL",
                               "N_OPC", "pct_OPC", "lor_lineage", "p_lin_hi", "p_lin_lo",
                               "lor_OLvsOPC", "p_OL_hi", "p_OL_lo"])
for ds in R.ds.unique():
    m = R.ds == ds
    for pc, qc in [("p_lin_hi", "q_lin_hi"), ("p_lin_lo", "q_lin_lo"), ("p_OL_hi", "q_OL_hi"), ("p_OL_lo", "q_OL_lo")]:
        R.loc[m, qc] = multipletests(np.nan_to_num(R.loc[m, pc].values, nan=1.0), method="fdr_bh")[1]
R["LIN_ENR"] = (R.q_lin_hi < .05) & (R.lor_lineage > 1)
R["NON_ENR"] = (R.q_lin_lo < .05) & (R.lor_lineage < -1)
R["OL_ENR"] = (R.q_OL_hi < .05) & (R.lor_OLvsOPC > 1)
R["OPC_ENR"] = (R.q_OL_lo < .05) & (R.lor_OLvsOPC < -1)
R.to_csv(f"{T}/S23_new_targets_stats.csv", index=False)
ST = pd.DataFrame(status, columns=["ds", "gene", "present_in_release"])
ST.to_csv(f"{T}/S23_presence_in_release.csv", index=False)
print("\nabsent from a released matrix:", ST[~ST.present_in_release].groupby("ds").gene.apply(list).to_dict())

S = pd.read_csv(T + "/S13_localization_classes.csv").set_index("gene")
summ = []
for g in R.gene.unique():
    x = R[(R.gene == g) & R.ds.isin(FULL3)]
    y = R[(R.gene == g) & (R.ds == "GSE277430")]
    n1, n2, n3, n4 = int(x.LIN_ENR.sum()), int(x.NON_ENR.sum()), int(x.OL_ENR.sum()), int(x.OPC_ENR.sum())
    # class numbering kept identical to tables/S13_localization_classes.csv
    if n1 >= 2 and n3 >= 2: cl = "1 lineage + mature-OL"
    elif n1 >= 2 and n4 >= 2: cl = "2 lineage + precursor/OPC"
    elif n1 >= 2: cl = "3 lineage, no maturation contrast"
    elif n2 >= 2: cl = "4 non-lineage restricted"
    elif n4 >= 2: cl = "5 no contrast (ubiquitous or undetectable)" if False else "6 OPC>OL but not lineage-enriched"
    else: cl = "5 no contrast (ubiquitous or undetectable)"
    pp = R[R.gene == g].pct_OL.max()
    summ.append([g, "both" if g in oldg and g in newatoms else ("old65" if g in oldg else "new33"),
                 n1, n2, n3, n4, cl,
                 round(np.nanmax(x.lor_lineage), 2) if x.lor_lineage.notna().any() else np.nan,
                 round(np.nanmax(x.lor_OLvsOPC), 2) if x.lor_OLvsOPC.notna().any() else np.nan,
                 round(pp, 2) if pd.notna(pp) else np.nan,
                 S["class"].get(g, "")])
U = pd.DataFrame(summ, columns=["gene", "list", "coh_lineage_enr", "coh_nonlin_enr", "coh_OLgtOPC",
                                "coh_OPCgtOL", "class_new", "max_lor_lineage", "max_lor_OLvsOPC",
                                "max_pct_OL", "class_old65"])
U.to_csv(f"{T}/S23_gene_universe_status.csv", index=False)
pd.set_option("display.width", 250)
nnew = U[U.list == "new33"]
print("\n== the %d atomic NEW genes, classified on the same axes ==" % len(nnew))
print(nnew.sort_values(["class_new", "max_pct_OL"], ascending=[True, False]).to_string(index=False))
print("\nnew-gene class distribution:\n", nnew.class_new.value_counts().to_string())
print("\nclass 1/2/3 (OL-side signal) among new:", sorted(nnew[nnew.class_new.str.startswith(("1 ", "2 ", "3 "))].gene))
print("\nuniverse total=%d | informative-class total=%d" %
      (len(U), int(U.class_new.str.startswith(("1 ", "2 ", "3 ")).sum())))
