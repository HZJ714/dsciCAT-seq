#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A6-step C2: remove the regulatory-domain-size confound of step C by downsampling EVERY peak set
to KPEAK peaks and averaging over RREP draws.

Motivation, measured in step C: the rank correlation of median_groupdiff with peak-set size is
-0.715 (more peaks -> more mean regression -> a more diluted group difference), which inflates the
cross-channel 74-trait rank correlation - rho = 0.76 falls to 0.46 (A2) and 0.19 (A3) once
log(n_peak) is controlled.

Method: the trait peak sets and the null-distribution peak sets are both randomly reduced to the
same size, the per-nucleus score is (trait mean - null mean) / SD(null), and step C then runs its
unchanged criteria (sign test of the group difference over the 6 blinded libraries, and per
library Spearman(score, maturation axis)).
"""
import os
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr, binomtest

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

OUTD = WORK + "/figures/routeA_out/A6_GSE231584_atac_v1"
DERV = os.path.join(OUTD, "derived")
TAB = os.path.join(OUTD, "tables")
GS = WORK + "/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs"
KPEAK = 6000
RREP = 12
NNULL = 40
t0 = time.time()


def log(*a):
    print("[%6ds]" % (time.time() - t0), *a, flush=True)


R = np.load(DERV + "/peak_resid_f16.npy", mmap_mode="r")     # peaks x nuclei (float16)
ann = pd.read_pickle(DERV + "/annot5.pkl")
PS = pd.read_csv(TAB + "/panel_scores_peak.csv", index_col=0)
n = R.shape[1]
pos = np.arange(n)


def set_score(pidx):
    out = np.zeros(n)
    for s in range(0, len(pidx), 20000):
        blk = pidx[s:s + 20000]
        out += np.sum(np.asarray(R[blk]), axis=0, dtype=np.float64)
    return out / len(pidx)


# ---------- peak-set mapping (as in step C) ----------
import h5py
peaks = pd.read_pickle(DERV + "/peaks.pkl")
gene = pd.read_pickle(DERV + "/gene_v2.pkl")
porder = peaks.sort_values(["chrom", "center"]).reset_index(drop=False)
gsorted = gene.sort_values(["chrom", "start"]).reset_index(drop=True)
WIN = 100_000


def gene_peaks(symbols):
    selg = gsorted[gsorted.symbol.isin(set(symbols))]
    out = []
    for ch, gp in selg.groupby("chrom"):
        sub = porder[porder.chrom == ch]
        pc = sub.center.values
        for s, e in zip(gp.start.values - WIN, gp.end.values + WIN):
            lo, hi = np.searchsorted(pc, s), np.searchsorted(pc, e)
            if hi > lo:
                out.append(sub["index"].values[lo:hi])
    return np.unique(np.concatenate(out)) if out else np.array([], dtype=int)


sets = {}
for l in open(GS):
    f = l.rstrip("\n").split("\t")
    if len(f) < 2 or f[0] == "TRAIT":
        continue
    sets[f[0]] = [x for x in f[1].split(",") if x]
allg = gsorted.symbol.values
pk_all = {tr: gene_peaks(gl) for tr, gl in sets.items()}
small = [t for t, p in pk_all.items() if len(p) < KPEAK]
log("%d traits; %d of them have fewer than %d peaks and need topping up: %s" % (len(sets), KPEAK, len(small), small[:8]))
nullgl = [np.random.default_rng(7000 + j).choice(allg, size=1000, replace=False) for j in range(NNULL)]
pk_null = [gene_peaks(g) for g in nullgl]
log("null sets: median number of peaks %.0f" % np.median([len(p) for p in pk_null]))

# ---------- scores after equalising the peak-set size ----------
rng = np.random.default_rng(20261030)


def score_sets(pidx_list, reps=RREP):
    """down/up-sample every peak set to KPEAK, reps times, and return the per-nucleus mean score matrix"""
    acc = np.zeros(n)
    for _ in range(reps):
        if len(pidx_list) >= KPEAK:
            pick = rng.choice(pidx_list, size=KPEAK, replace=False)
        else:                                    # too few peaks -> top up to KPEAK with replacement (same dilution everywhere)
            pick = rng.choice(pidx_list, size=KPEAK, replace=True)
        acc += set_score(np.sort(pick))
    return acc / reps


NULL = np.vstack([score_sets(p, reps=4) for p in pk_null])
nmu, nsd = NULL.mean(0), NULL.std(0) + 1e-9
log("null sets (same %d peaks): median per-nucleus mean %.3f | median SD %.3f" % (KPEAK, np.median(nmu), np.median(nsd)))

axis = PS.axis.values if "axis" in PS else (PS.MAT - PS.PRE).values
lab = PS.lab.values
lib = PS.lib.values
summary = pd.DataFrame({"pos": pos, "lib": [str(x) for x in lib], "lab": lab, "axis": axis})
log("OL nuclei %d (mature %d / precursor %d) in %d libraries" % ((lab != "nonOL").sum(), (lab == "mature").sum(),
     (lab == "precursor").sum(), summary.lib.nunique()))

rows, per_lib = [], {}
for i, tr in enumerate(sorted(sets), 1):
    pk = pk_all[tr]
    if len(pk) < 300:
        continue
    sc = (score_sets(pk) - nmu) / nsd
    dif, rho, diflab = [], [], []
    for L, d in summary.groupby("lib"):
        m = (d.lab != "nonOL").values
        if m.sum() < 30:
            continue
        a = sc[d.pos.values[(d.lab == "mature").values]]
        b = sc[d.pos.values[(d.lab == "precursor").values]]
        if len(a) >= 10 and len(b) >= 10:
            dif.append(np.median(a) - np.median(b)); diflab.append(str(L))
        rho.append(spearmanr(sc[d.pos.values[m]], d.axis.values[m]).statistic)
    if len(dif) < 4:
        continue
    dif = np.array(dif)
    npos = int((dif > 0).sum())
    rows.append({"trait": tr, "n_peak_raw": len(pk), "n_lib": len(dif), "pos_libs": npos,
                 "P_sign": binomtest(npos, len(dif), 0.5).pvalue, "median_groupdiff": float(np.median(dif)),
                 "rho_median": float(np.median(rho)), "rho_pos": int(np.sum(np.array(rho) > 0))})
    per_lib[tr] = {"lib" + L: v for L, v in zip(diflab, dif)}
    if i % 10 == 0:
        log("  [%d/%d] %s group difference %+.3f, positive %d/%d" % (i, len(sets), tr, np.median(dif), npos, len(dif)))
T = pd.DataFrame(rows)
T["rank"] = (-T.median_groupdiff).rank(method="min").astype(int)
T["rank_rho"] = (-T.rho_median).rank(method="min").astype(int)
T.to_csv(TAB + "/trait_chromatin_sizeMatched.csv", index=False)
pd.DataFrame(per_lib).T.to_csv(TAB + "/trait_groupdiff_by_library_matched.csv")
log("\nMS:\n" + T[T.trait == "PASS_Multiple_sclerosis"].round(4).to_string(index=False))
log("residual size confound: rho(group difference, original peak count) = %.3f P=%.3g" % spearmanr(T.median_groupdiff, T.n_peak_raw)[:2])
log("\ntop 10:\n" + T.sort_values("rank").head(10)[["trait", "median_groupdiff", "pos_libs", "P_sign", "n_peak_raw"]].round(3).to_string(index=False))

# ---------- align with the RNA channel ----------
a2 = pd.read_csv(WORK + "/figures/routeA_out/A2_GSE277430_v1/tables/trait_sweep_dual_axis.csv")[["trait", "r1"]].rename(columns={"r1": "A2_d1"})
a3 = pd.read_csv(WORK + "/figures/routeA_out/A3_GSE118257_v1/tables/trait_sweep_delta2_GSE118257_74.csv")[["trait", "r1"]].rename(columns={"r1": "A3_d1"})
a4f = WORK + "/figures/routeA_out/A4_GSE301908_v1/tables/trait_sweep_74_GSE301908.csv"
m = T.merge(a2, on="trait").merge(a3, on="trait")
m["r_chr"] = m["rank"]; m["r_A2"] = m.A2_d1; m["r_A3"] = m.A3_d1
log("\n=== cross-channel 74-trait rank agreement (peak-set sizes matched) n=%d ===" % len(m))
for x, y in [("r_chr", "r_A2"), ("r_chr", "r_A3"), ("r_A2", "r_A3")]:
    r, p = spearmanr(m[x], m[y])[:2]
    log("  rho(%s,%s) = %.3f  P=%.2g" % (x, y, r, p))
if os.path.exists(a4f):
    a4 = pd.read_csv(a4f)
    rc = "rank" if "rank" in a4.columns else a4.columns[0]
    m2 = m.merge(a4[["trait", rc]].rename(columns={rc: "r_A4"}), on="trait", how="inner")
    r, pv = spearmanr(m2.r_chr, m2.r_A4)[:2]
    log("  rho(chr,A4) = %.3f  P=%.2g (n=%d)" % (r, pv, len(m2)))
open(TAB + "/summary_stepC2.md", "w").write(
    "\n".join(["# A6 step C2: chromatin channel, 74 traits after matching the peak-set size",
               "", "Every peak set was reduced or topped up to %d peaks and averaged over %d draws (the null distribution at the same size; %d random gene sets)." % (KPEAK, RREP, NNULL),
               "OL nuclei = %d (mature %d / precursor %d); unit = %d blinded libraries." % ((lab != "nonOL").sum(), (lab == "mature").sum(), (lab == "precursor").sum(), summary.lib.nunique()),
               "", "rho(group difference, original peak count) = %.3f" % spearmanr(T.median_groupdiff, T.n_peak_raw).statistic,
               "", "MS:", T[T.trait == "PASS_Multiple_sclerosis"].round(4).to_string(index=False),
               "", "cross-channel rank correlations: chr-A2 %.3f / chr-A3 %.3f / A2-A3 %.3f" % (
                   spearmanr(m.r_chr, m.r_A2).statistic, spearmanr(m.r_chr, m.r_A3).statistic, spearmanr(m.r_A2, m.r_A3).statistic),
               "", T.sort_values("rank").round(4).to_string(index=False)]) + "\n")
log("wrote trait_chromatin_sizeMatched.csv / summary_stepC2.md")
