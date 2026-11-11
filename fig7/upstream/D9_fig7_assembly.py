#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D9: the assembly tables for Fig 7 - the delivered numbers that live in different cohorts, collected
into the two tables the figure reads.

fig7A_donors_delta1.csv     donor / library-level d1 of the 5 studies (the dots), the condition and
                            the between-donor P printed above each column of panel A:
                            A2 GSE277430 (donor 13/8) | A3 GSE118257 (donor 4/5) |
                            A4 GSE301908 (library 14/3) | A7 GSE231585 (donor 4/4, our annotation) |
                            M GSE250589 (19 libraries, author annotation, EAE/Ctrl)
fig7E_disease_block.csv     three layers of disease evidence (between-donor / within-lesion / mouse
                            EAE) plus the composition-confounding caveats. Everything is read from
                            delivered tables, nothing is recomputed here; A7 always uses
                            trait_sweep_74_..._donorLevel.csv together with
                            delta1_MS_authorpanel_donors.csv (handover section 6.4 item 2: the old
                            raw table is retired).

Unit disclosure: the A4 main analysis has donor=11 but its case-control contrast uses library=17
(14/3); GSE277430 donor 21 (13/8); GSE118257 donor 9 (4/5); A7 donor 8 (4/4); GSE250589 19
libraries / 18 mice.
"""
import os, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

R = WORK + "/figures/routeA_out"
OUT = R + "/FIG7_v11/tables"
os.makedirs(OUT, exist_ok=True)
SAY = lambda fmt, *a: print("[%s] %s" % (pd.Timestamp.now().strftime("%H:%M:%S"), fmt % a if a else fmt), flush=True)
MS = "PASS_Multiple_sclerosis"
rows = []

# ---- A2 GSE277430: donor-level d1 (median MS norm_score difference OLIGO - OPC); condition taken from disease in meta ----
a2 = pd.read_csv(R + "/A2_GSE277430_v1/tables/OL_meta_with_MS_score.csv",
                 usecols=["expertAnno.l1", "patient", "disease", "MS_norm_score"])
dm = a2[a2["expertAnno.l1"] == "OLIGO"].groupby("patient")["MS_norm_score"].median()
dp = a2[a2["expertAnno.l1"] == "OPC"].groupby("patient")["MS_norm_score"].median()
dis = a2.groupby("patient")["disease"].first()
cm = dm.index.intersection(dp.index)
v2 = (dm[cm] - dp[cm])
for p, val in v2.items():
    rows.append(dict(study="A2 GSE277430", unit="donor", id=str(p), cond=str(dis[p]), d1=float(val)))
# ---- A3 GSE118257 ----
a3 = pd.read_csv(R + "/A3_GSE118257_v1/tables/patient_pair_scores.csv")
for _, r in a3.iterrows():
    rows.append(dict(study="A3 GSE118257", unit="donor", id=str(r.patient), cond=str(r.Condition), d1=float(r["diff"])))
# ---- A4 GSE301908 (library level, patient_info pure per library) ----
a4 = pd.read_csv(R + "/A4_GSE301908_v1/tables/OL_meta_with_MS_score.csv")
a4["patient"] = a4["orig.ident"].astype(str)
dm4 = a4[a4.majorCluster == "Oligo"].groupby("patient")["MS_norm_score"].median()
dp4 = a4[a4.majorCluster == "OPC"].groupby("patient")["MS_norm_score"].median()
st4 = a4.groupby("patient")["patient_info"].first()
cm4 = dm4.index.intersection(dp4.index)
v4 = dm4[cm4] - dp4[cm4]
for p, val in v4.items():
    rows.append(dict(study="A4 GSE301908", unit="library", id=str(p), cond={"patient": "MS", "ctrl": "Ctrl"}[str(st4[p])], d1=float(val)))
# ---- A7 GSE231585 (donor-level d1 under the new author annotation; the old raw table is retired) ----
a7 = pd.read_csv(R + "/A7_GSE231585_rna_v1/tables/delta1_MS_authorpanel_donors.csv")
for _, r in a7.iterrows():
    rows.append(dict(study="A7 GSE231585", unit="donor", id=str(r.donor), cond="MS" if str(r.donor).startswith("MS") else "Ctrl",
                     d1=float(r.d1_new)))
# ---- M GSE250589 (author annotation, library level; phase = Ctrl/Early/Peak/Late) ----
mm = pd.read_csv(R + "/M1b_authoranno_v1/tables/per_lib_delta1_author_vs_self.csv")
mm = mm[(mm.trait == MS) & (mm.anno == "author")]
for _, r in mm.iterrows():
    rows.append(dict(study="M GSE250589", unit="library", id=str(r.lib), cond="Ctrl" if r.phase == "Ctrl" else "EAE", d1=float(r.delta1)))

A = pd.DataFrame(rows)
A["cond"] = A.cond.replace({"MS_Mye": "MS", "MS_Demye": "MS", "CTRL": "Ctrl", "CONTROL": "Ctrl"})
A["cond"] = A.cond.str.strip().str.replace("^MS$", "MS", regex=True)
A.to_csv(OUT + "/fig7A_donors_delta1.csv", index=False)
SAY("fig7A_donors_delta1.csv: %d rows | %s", len(A),
    A.groupby(["study", "cond"]).size().to_dict())

# between-donor P for every column (computed in this same table; the face only reads it)
pr = []
for st, g in A.groupby("study"):
    x, y = g[g.cond.isin(["MS", "EAE"])].d1.values, g[g.cond == "Ctrl"].d1.values
    u = mannwhitneyu(x, y, alternative="two-sided")
    pr.append(dict(study=st, unit=g.unit.iloc[0], n_case=len(x), n_ctrl=len(y),
                   median_case=float(np.median(x)), median_ctrl=float(np.median(y)),
                   delta=float(np.median(x) - np.median(y)), P=float(u.pvalue),
                   frac_pos_case=float((x > 0).mean()), frac_pos_ctrl=float((y > 0).mean())))
P = pd.DataFrame(pr).sort_values("study")
P.to_csv(OUT + "/fig7A_between_donor_P.csv", index=False)
SAY("\n" + P.round(4).to_string(index=False))

# ---- three layers of the disease block ----
cc = pd.read_csv(R + "/A5_case_control/tables/case_control_delta1.csv")
gr = pd.read_csv(R + "/A2_GSE277430_v1/tables/group_comparisons.csv")
dz = pd.read_csv(R + "/M1b_authoranno_v1/tables/disease_axis_author.csv").iloc[0]
sw = pd.read_csv(R + "/A7_GSE231585_rna_v1/tables/trait_sweep_74_GSE231585_donorLevel.csv")
# 11-05 unified basis: the sweep rows now use the sample-level three-axis tables that B/C use too (the old donor-level table is kept as a reference)
a7new = pd.read_csv(R + "/A7_GSE231585_rna_v1/tables/trait_sweep_delta2_GSE231585_74.csv")
a7ms = a7new[a7new.trait == MS].iloc[0]
comp = pd.read_csv(R + "/D_crosschemistry/tables/fig6D_composition.csv")
E = []
for _, r in cc.iterrows():
    E.append(dict(layer="1 between-donor MS vs Ctrl", study=r.cohort, unit=r.unit, n="%s/%s" % (r.n_MS, r.n_Ctrl),
                  value=r.median_MS_minus_Ctrl, P=r.P_ranksum, effect=r.rank_biserial,
                  source="A5_case_control/tables/case_control_delta1.csv",
                  note="all negative: nothing measurable at the diagnosis level"))
for _, r in gr.iterrows():
    E.append(dict(layer="2 within-lesion (same donor)", study="GSE277430", unit="lesion", n="%s/%s" % (r.n1, r.n0),
                  value=float(r.median1) - float(r.median0), P=r.P, effect=r.cliffs_d,
                  source="A2_GSE277430_v1/tables/group_comparisons.csv", note=r.test))
E.append(dict(layer="3 mouse EAE (author annotation)", study="GSE250589", unit="library", n="%s/%s" % (dz.n_eae_libs, dz.n_ctrl_libs),
              value=float(dz.delta), P=float(dz.MWU_P), effect=np.nan,
              source="M1b_authoranno_v1/tables/disease_axis_author.csv",
              note="EAE below Ctrl (printed as the table says; criteria unchanged); MS d1 positive in %d/%d libraries" % (dz.n_pos_all, dz.n_lib_all)))
E.append(dict(layer="3b mouse specificity", study="GSE250589", unit="library", n="6 traits",
              value=float(dz.PASS_Schizophrenia_Pardinas2018__delta), P=float(dz.PASS_Schizophrenia_Pardinas2018__P),
              effect=np.nan, source="as above", note="the disease effect of SCZ (control comparison, expected none)"))
E.append(dict(layer="4 composition confound", study="GSE250589", unit="mouse", n="18",
              value=float(comp[comp.quantity.str.startswith("composition")].delta.iloc[0]),
              P=float(comp[comp.quantity.str.startswith("composition")].mwu_P.iloc[0]), effect=np.nan,
              source="D_crosschemistry/tables/fig6D_composition.csv",
              note="n_MOL/n_OPC itself differs between EAE and Ctrl; adjusted RNA-axis P=%.3f (see the same table)"
                   % float(comp[comp["quantity"].str.contains("after adjusting") & (comp.axis == "RNA")].mwu_P.iloc[0])))
E.append(dict(layer="5 A7 sample-level sweep", study="GSE231585", unit="sample (donor×lesion domain)", n="%d" % int(a7ms.n1),
              value=float(a7ms.d1), P=float(a7ms.P1),
              effect=np.nan, source="A7_GSE231585_rna_v1/tables/trait_sweep_delta2_GSE231585_74.csv",
              note="sample level (donor x lesion domain, %d samples pass the >= 10 nuclei per group gate), median over the 74-trait scan %+.4f - not a case-vs-control contrast; the same tables as B/C"
                   % (int(a7ms.n1), float(a7ms.d1))))
pd.DataFrame(E).to_csv(OUT + "/fig7E_disease_block.csv", index=False)
SAY("\nfig7E_disease_block.csv written with %d rows" % len(E))
SAY(pd.DataFrame(E)[["layer", "study", "n", "value", "P"]].round(4).to_string(index=False))
