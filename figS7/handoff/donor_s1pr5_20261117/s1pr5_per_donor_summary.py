#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s1pr5_per_donor_summary.py -- corrected donor-level statistics for Fig 7D (reads s40 tables).

Fixes one bug in s40: the evaluable-donor test leaked a loop variable, so every donor was judged on the
nuclei counts of ms8 and all 11 came out "evaluable".  Here the floor is applied per donor properly, and
two additional views are added:
  * an exact one-sided Fisher test per donor (OL vs each non-lineage group, same library, so depth cancels
    inside the donor) -- this is what a small comparator pool (9 nuclei in ms10 LYM) should be allowed to say;
  * a floor sweep (30 / 50 / 100 nuclei per donor x population) so the reader can see which statement is
    robust to the floor and which is not.
"""
import numpy as np
import pandas as pd
from scipy.stats import binomtest, fisher_exact

import os
T = os.environ.get("S7_DONOR_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "tables"))
NL = ["AST", "MG", "LYM", "VASC", "NEU"]
GENES = ["S1PR5", "IFNAR2", "S1PR2", "BTK", "CD40LG", "CD52", "DCK", "ITGA4", "MS4A1", "S1PR1",
         "S1PR3", "S1PR4", "DHODH", "IFNAR1", "MC2R", "NFE2L2", "NR3C1", "TOP2A", "PLP1", "PDGFRA"]
L = pd.read_csv(T + "/S40_det_donor_pop.csv")
R = L.pivot(index=["gene", "pop"], columns="donor", values="rate")
Dt = L.pivot(index=["gene", "pop"], columns="donor", values="det")
Nn = L.pivot(index=["gene", "pop"], columns="donor", values="n")
DON = list(R.columns)

# nuclei per donor x population is gene-independent; take it from any row
nDP = Nn.loc["S1PR5"].loc[NL]
print("comparator pool sizes (min across the 5 non-lineage groups, per donor):")
print(nDP.min(axis=0).astype(int).to_string(), "\n")

rows = []
for g in GENES:
    ol = R.loc[g].loc["OL"]; nl = R.loc[g].loc[NL]
    mx = nl.max(); who = nl.idxmax()
    for floor in (30, 50, 100):
        ev = [d for d in DON if nDP[d].min() >= floor]
        ok = [d for d in ev if ol[d] > mx[d] + 1.0 and ol[d] > 1.5 * mx[d]]
        hi = [d for d in ev if ol[d] > mx[d]]
        # per-donor exact test: OL beats every comparator at p<0.05 (one-sided, greater)
        sig = []
        for d in ev:
            good = True
            for p in NL:
                a, b = int(Dt.loc[(g, "OL"), d]), int(Dt.loc[(g, p), d])
                na, nb = int(Nn.loc[(g, "OL"), d]), int(Nn.loc[(g, p), d])
                if b / nb >= a / na:
                    good = False; break
                pv = fisher_exact([[a, na - a], [b, nb - b]], alternative="greater")[1]
                if pv >= 0.05:
                    good = False; break
            if good:
                sig.append(d)
        n = len(ev)
        rows.append(dict(gene=g, floor=floor, n_evaluable=n, OL_highest=len(hi),
                         p_highest=binomtest(len(hi), n, .5, alternative="greater").pvalue if n else np.nan,
                         rule_pass=len(ok), p_rule=binomtest(len(ok), n, .5, alternative="greater").pvalue if n else np.nan,
                         fisher_all5=len(sig), med_fold=round(float(np.median([ol[d] / max(mx[d], 1e-9) for d in ev])), 2)))
S = pd.DataFrame(rows)
S.to_csv(T + "/S40b_summary_floors.csv", index=False)
for fl in (30, 50, 100):
    print("\n=== floor %d nuclei per donor x population ===" % fl)
    print(S[S.floor == fl].sort_values(["rule_pass", "OL_highest"], ascending=False)
          .head(8).to_string(index=False))

d5 = R.loc["S1PR5"]
det = pd.DataFrame({"OL_pct": d5.loc["OL"].round(2), "binding_group": d5.loc[NL].idxmax(),
                    "binding_pct": d5.loc[NL].max().round(2),
                    "n_binding": [int(Nn.loc[("S1PR5", d5.loc[NL].idxmax()[d]), d]) for d in DON],
                    "fold": (d5.loc["OL"] / d5.loc[NL].max()).round(2)})
det["rule_1pp_1p5x"] = (det.OL_pct > det.binding_pct + 1.0) & (det.OL_pct > 1.5 * det.binding_pct)
det.to_csv(T + "/S40b_S1PR5_per_donor.csv")
print("\nS1PR5 per donor:\n", det.to_string())
