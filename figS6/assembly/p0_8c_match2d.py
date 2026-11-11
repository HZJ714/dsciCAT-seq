# -*- coding: utf-8 -*-
"""p0_8c_match2d.py - upgrade of the top-1000 landing test's controls: two-dimensional matching (detection breadth x mean expression)

scDRS's own control sets match on "mean expression + detection breadth"; a one-dimensional match (detection rate only) may be insufficient.
Here the background controls upgrade to 5x5 two-dimensional strata (detection breadth deciles -> quintiles x the channel's mean-signal quintiles),
one control drawn per stratum, 10,000 permutations, testing the set-mean Δ. Positive controls = the 8 myelin/OL anchor genes.
"""
import os
import numpy as np, pandas as pd
TAB = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'routeA_out/P0_8C/tables')
M = pd.read_csv(TAB + "/p08c_dual_landing_genomewide.csv")
L = pd.read_csv(TAB + "/p08c_genelevel_all6_pool_allgenes.csv")
ms875 = set(L.gene.astype(str)); ms223 = set(L[~L.near_zero].gene.astype(str))
cis38 = set(pd.read_csv(TAB + "/p08c_c1_rna_atac_by_subtype.csv").mouse.astype(str))
PLAT = os.path.join(os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')), 'oligo/scdrs_ms_out/genes_of_interest.txt')
anch = ["Mbp", "Plp1", "Mog", "Mag", "Cnp", "Mobp", "Sox10", "Olig2"]
SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
P331 = set(x.strip().upper() for x in open(PLAT) if x.strip())
M["gu"] = M.gene.astype(str).str.upper()          # platform/window lists are uppercase; M.gene keeps the human symbols as-is
P293 = set(M.gu) & (P331 - set(str(x).upper() for x in cis38))
CONTRAST = os.environ.get("M2_CONTRAST", "pool")
OUTF = os.environ.get("M2_OUT", "p08c_dual_setlevel_matched2d.csv")
NPERM, SEED = 10000, 20261116
out = []
for tag in ("rna", "atac"):
    d = M[tag + "_d"].values
    if CONTRAST == "mol6":      # the main-figure convention: Δ = MOL6 − COP/NFOL (author 2026-11-17)
        d = M[tag + "_" + SUB[5]].values - M[tag + "_" + SUB[0]].values
    pct = M[tag + "_pct"].values
    mean = M[[tag + "_" + s for s in SUB]].mean(1).values
    pool = (pct >= 10) & np.isfinite(d) & np.isfinite(mean)
    isbg = pool & ~M.gene.isin(ms875).values
    # two-dimensional strata: detection breadth quintiles x mean-signal quintiles (quantiles fixed on the background pool only, so the risk sets cannot steer the strata)
    qb = np.quantile(pct[isbg], [.2, .4, .6, .8]); qm = np.quantile(mean[isbg], [.2, .4, .6, .8])
    cellb = np.digitize(pct, qb) * 6 + np.digitize(mean, qm)
    cells = {c: np.flatnonzero(isbg & (cellb == c)) for c in np.unique(cellb[isbg])}
    # same-panel control: the background pool becomes "the 293 platform genes outside the cis windows", with the SAME strata and permutation as the main test.
    # it used to live in p0_8c_p2_38vs293.py (4x4 strata, RNA P=0.0155) and p0_8c_p2_robust.py (4x4, P=0.0318),
    # two stratum schemes gave two P values - the face/caption can carry only one, so it merged into this file, sharing the stratum convention with the main test (SKILL: single source).
    isbg293 = pool & M.gu.isin(P293).values
    cells293 = {c: np.flatnonzero(isbg293 & (cellb == c)) for c in np.unique(cellb[isbg293])}
    print("\n=== %s: background pool %d, non-empty 2D strata %d ===" % (tag.upper(), isbg.sum(), len(cells)), flush=True)
    for nm, gs in (("myelin anchors (positive control)", anch), ("top1000 measurable 223", ms223), ("cis38", cis38),
                   ("cis38 vs platform 293", cis38)):
        gi = np.flatnonzero(pool & M.gene.isin(gs).values)
        cb = cellb[gi]
        pool_cells = cells293 if nm.endswith("vs platform 293") else cells
        empty = np.array([c not in pool_cells or len(pool_cells[c]) == 0 for c in cb])
        gi2, cb2 = gi[~empty], cb[~empty]
        if gi2.size < 5:
            print("%-14s too few in stratum %d" % (nm, gi2.size)); continue
        obs = d[gi2].mean()
        rng = np.random.default_rng(SEED)
        null = np.empty(NPERM)
        for k in range(NPERM):
            pick = np.array([pool_cells[c][rng.integers(len(pool_cells[c]))] for c in cb2])
            null[k] = d[pick].mean()
        p_hi = (np.sum(null >= obs) + 1) / (NPERM + 1); p_lo = (np.sum(null <= obs) + 1) / (NPERM + 1)
        z = (obs - null.mean()) / null.std()
        print("%-14s n=%4d (%d not stratified) Δ%+.4f matched null %+.4f±%.4f z=%+.2f P_hi=%.4f P_lo=%.4f"
              % (nm, gi2.size, empty.sum(), obs, null.mean(), null.std(), z, p_hi, p_lo), flush=True)
        out.append(dict(channel=tag, set=nm, n=gi2.size, dropped=int(empty.sum()), mean_d=obs,
                        null_mean=null.mean(), null_sd=null.std(), z=z, p_greater=p_hi, p_less=p_lo))
        loo_p, loo_gene = np.nan, ""
        if os.environ.get("M2_LOO", "0") == "1" and nm in ("cis38", "cis38 vs platform 293"):
            # leave-one-out: the WORST P over dropping any single gene. Guards against "one gene carrying the result".
            for k in range(gi2.size):
                m = np.ones(gi2.size, bool); m[k] = False
                o = d[gi2[m]].mean(); c = cb2[m]; r2 = np.random.default_rng(SEED + k)
                nn = np.array([d[[pool_cells[q][r2.integers(len(pool_cells[q]))] for q in c]].mean()
                               for _ in range(NPERM)])
                pp = (np.sum(nn >= o) + 1) / (NPERM + 1)
                if pp > loo_p or not loo_gene:
                    loo_p, loo_gene = pp, str(M.gene.iloc[gi2[k]])
            print("%-14s LOO worst P=%.4f (dropping %s)" % (nm, loo_p, loo_gene), flush=True)
        out[-1].update(loo_worst_p=loo_p, loo_drop=loo_gene)
pd.DataFrame(out).to_csv(TAB + "/" + OUTF, index=False)
print("\n[written] %s (Δ convention = %s)" % (OUTF, CONTRAST))
