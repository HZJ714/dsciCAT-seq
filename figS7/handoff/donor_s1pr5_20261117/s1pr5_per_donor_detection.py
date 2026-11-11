#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s1pr5_per_donor_detection.py -- Fig 7D's detection block re-tested at the DONOR unit (author 2026-11-17).

Why: the "% nuclei (log)" block of Fig 7D pools 293,970 nuclei of GSE301908 and draws the red box on
S1PR5 with a pooled rule (OL beats every non-lineage group by >=1.0 pp AND >=1.5x).  A pooled rate
inherits the house rule that the statistical unit is the donor, not the nucleus: the same rule has to
be re-applied inside each donor and counted.  Two further facts make this non-cosmetic:
  * the 17 released libraries come from only 11 donors (library_group_map.csv; ms1/ms4/ms6/ms12/ms15/
    ms16 contribute two libraries each), and three libraries GEO labels "Ctrl" (ms6r1, ms12r1, ms16r1)
    are second libraries of donors that also appear as MS -- so "Ctrl" here is a second site of the same
    subject, not a different subject;
  * population sizes per donor swing from 5 nuclei (Ctrl1 Lym) to 29,531 (MS9 Oligo), so a per-donor
    rate needs a floor before it can be compared at all.

Outputs (tables only, no figure): S40_det_donor_pop.csv, S40_rate_donor_pop.csv, S40_rule_per_donor.csv,
S40_summary.csv.  Read-only on all inputs.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import binomtest

_ROOT = os.environ.get("DSCI_WORK_ROOT", os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
D = os.environ.get("FIG7_PROBE", os.path.join(_ROOT, "routeA_probe", "GSE301908_counts")) + os.sep
OUT = os.environ.get("S7_DONOR_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "tables"))
os.makedirs(OUT, exist_ok=True)
MAP = {"Oligo": "OL", "OPC": "OPC", "Astro": "AST", "Micro": "MG", "Lym": "LYM",
       "Endo": "VASC", "Neuron": "NEU"}
POPS = ["OL", "OPC", "AST", "MG", "LYM", "VASC", "NEU"]
NONLIN = ["AST", "MG", "LYM", "VASC", "NEU"]          # the 7D rule's comparator set (lineage groups
MIN_POP = 30                                          # excluded from "non-lineage" by construction)
TARGETS = ["S1PR5", "IFNAR2", "S1PR2", "BTK", "CD40LG", "CD52", "DCK", "ITGA4", "MS4A1", "S1PR1",
           "S1PR3", "S1PR4", "DHODH", "IFNAR1", "MC2R", "NFE2L2", "NR3C1", "TOP2A"]
MARKERS = ["PLP1", "PDGFRA", "GAPDH", "ACTB"]         # lineage / housekeeping controls
GENES = TARGETS + MARKERS

meta = pd.read_csv(D + "full_meta.tsv", sep="\t", dtype=str)
genes = [l.rstrip("\n") for l in open(D + "full_genes.tsv")]
assert len(meta) == 293970, len(meta)
donor = meta.grp.str.extract(r"^(ms\d+)")[0].values          # msNN prefix = subject
pop = meta.majorCluster.map(MAP).values
donor = pd.Series(donor).fillna(meta.grp).values             # anything unprefixed keeps its library id
ud, up = sorted(set(donor)), POPS
di = pd.Series({v: i for i, v in enumerate(ud)}).reindex(donor).to_numpy(np.int16)
pi = pd.Series({v: i for i, v in enumerate(up)}).reindex(pop).to_numpy(np.float32)
pi = np.where(np.isnan(pi), -1, pi).astype(np.int16)
gslot = np.full(len(genes), -1, np.int16)
for s, g in enumerate(GENES):
    gslot[genes.index(g)] = s
nD, nP, nG = len(ud), len(up), len(GENES)
tot = np.zeros((nD, nP), np.int64)
np.add.at(tot, (di, np.where(pi >= 0, pi, 0)), np.where(pi >= 0, 1, 0))

ip = np.fromfile(D + "full_csr.bin", dtype=np.int32, count=len(meta) + 1)
nnz = int(ip[-1])
det = np.zeros((nG, nP, nD), np.int64)
CH = 40_000_000
for s in range(0, nnz, CH):
    e = min(s + CH, nnz)
    ix = np.memmap(D + "full_csr.bin", dtype=np.int32, mode="r", offset=(len(meta) + 1) * 4,
                   shape=(nnz,))[s:e]
    va = np.memmap(D + "full_csr.bin", dtype=np.float32, mode="r",
                   offset=(len(meta) + 1) * 4 + nnz * 4, shape=(nnz,))[s:e]
    pos = np.nonzero(np.asarray(va) > 0)[0]
    cell = np.searchsorted(ip, pos + s, side="right") - 1
    g = gslot[np.asarray(ix)[pos]]
    keep = (g >= 0) & (pi[cell] >= 0)
    flat = g[keep] * (nP * nD) + pi[cell][keep] * nD + di[cell][keep]
    det += np.bincount(flat, minlength=det.size).reshape(det.shape)
    print("  %d/%d nonzero (%.1f GB read)" % (e, nnz, (e * 8 + (len(meta) + 1) * 4) / 1e9), flush=True)

rows = []
for gi_, g in enumerate(GENES):
    for p_, p in enumerate(up):
        for d_, d in enumerate(ud):
            rows.append(dict(gene=g, pop=p, donor=d, det=det[gi_, p_, d_], n=tot[d_, p_],
                             rate=100.0 * det[gi_, p_, d_] / max(tot[d_, p_], 1)))
L = pd.DataFrame(rows)
L.to_csv(f"{OUT}/S40_det_donor_pop.csv", index=False)
R = L.pivot_table(index=["gene", "pop"], columns="donor", values="rate")
R.to_csv(f"{OUT}/S40_rate_donor_pop.csv")

# donors where every comparator of the rule is evaluable (>= MIN_POP nuclei in that donor x population)
ok_donor = {d: all(tot[d_, up.index(p)] >= MIN_POP for p in NONLIN) and tot[d_, up.index("OL")] >= MIN_POP
            for d, dd in zip(ud, range(nD))}
print("\ndonor nuclei (OL / evaluable):", {d: int(tot[ud.index(d), 0]) for d in ud})
print("evaluable donors (all 5 non-lineage groups >= %d nuclei): %d/%d -> %s"
      % (MIN_POP, sum(ok_donor.values()), len(ud), [d for d in ud if ok_donor[d]]))

out = []
for g in GENES:
    rr = R.loc[g]
    ol = rr.loc["OL"]
    nl = rr.loc[[p for p in NONLIN if p in rr.index]]
    mx = nl.max()
    d_ok = [d for d in ud if ok_donor[d] and d in ol.index and not np.isnan(ol[d])]
    hi = [d for d in d_ok if ol[d] > mx[d]]                                   # OL highest of the 5
    # exactly the rule that draws the red box in 7D (s36_render_C18.py:161), applied inside each donor
    rule = [d for d in d_ok if ol[d] > mx[d] + 1.0 and ol[d] > 1.5 * mx[d] and ol[d] > 0]
    pp1 = [d for d in d_ok if ol[d] > mx[d] + 1.0]
    fc = [d for d in d_ok if ol[d] > 1.5 * mx[d]]
    both = rule                                                               # the 7D red-box rule
    n = len(d_ok)
    bt = binomtest(len(hi), n, 0.5, alternative="greater") if n else None
    bt2 = binomtest(len(both), n, 0.5, alternative="greater") if n else None
    out.append(dict(gene=g, n_donors=n, OL_highest=len(hi), pass_1pp=len(pp1), pass_1p5x=len(fc),
                    pass_rule=len(both), p_highest=(bt.pvalue if bt else np.nan),
                    p_rule=(bt2.pvalue if bt2 else np.nan),
                    med_OL=round(float(np.median([ol[d] for d in d_ok])), 2),
                    med_max_nonlin=round(float(np.median([mx[d] for d in d_ok])), 2),
                    med_fold=round(float(np.median([ol[d] / max(mx[d], 1e-9) for d in d_ok])), 2)))
S = pd.DataFrame(out).sort_values("pass_rule", ascending=False)
S.to_csv(f"{OUT}/S40_summary.csv", index=False)
print("\n=== donor-level (GSE301908, %d evaluable donors, floor %d nuclei/donor x population) ==="
      % (sum(ok_donor.values()), MIN_POP))
print(S.to_string(index=False))

# per-donor detail for S1PR5, and the same rule pooled over nuclei (what 7D actually shows)
d5 = R.loc["S1PR5"]
pd.DataFrame({"OL": d5.loc["OL"], "max_nonlin": d5.loc[NONLIN].max(),
              "fold": d5.loc["OL"] / d5.loc[NONLIN].max().clip(lower=1e-9)}).round(2).to_csv(
    f"{OUT}/S40_rule_per_donor.csv")
print("\nS1PR5 per donor (%% detected):\n", d5.round(2).to_string())

pooled = pd.read_csv(os.path.join(OUT, "S15_det_GSE301908.csv"),
                     index_col=0)
np_ = pd.read_csv(os.path.join(OUT, "S15_n_GSE301908.csv"),
                  index_col=0)["n"]
PR = pooled.div(np_.reindex(pooled.columns).values.astype(float), axis=1) * 100
pl = [g for g in GENES if g in PR.index]
pool_rule = {g: bool(PR.loc[g, "OL"] > np.nanmax(PR.loc[g, NONLIN].values) + 1.0 and
                     PR.loc[g, "OL"] > 1.5 * np.nanmax(PR.loc[g, NONLIN].values)) for g in pl}
print("\npooled (as drawn in 7D) red-box rule:", {k: v for k, v in pool_rule.items() if v})
S.assign(pooled_rule=[pool_rule.get(g, False) for g in S.gene]).to_csv(f"{OUT}/S40_summary.csv",
                                                                       index=False)
print("\nwrote 4 tables to", OUT)
