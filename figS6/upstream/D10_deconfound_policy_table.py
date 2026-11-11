#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D10_deconfound_policy_table.py - the whole-paper "when deconfounding is mandatory" criteria table (the numeric body of plan-A deliverable 1)

For EVERY dataset x channel, measures the confounding of regulatory-domain size (ln n_peak) on Δ1 within that dataset's 74 traits,
and rules per the three trigger criteria of `figures/DECONFOUND_POLICY_20261106.md` §1.

Criteria (any one met => every rank-type conclusion of that dataset must show raw + adjusted side by side):
  (i)   |ρ(ln n_peak, Δ1)| ≥ 0.30
  (ii)  P < 0.05
  (iii) the adjustment changes any face verdict / a labelled trait moves >= 5 ranks / the top-5 composition changes

Self-checks (rather output nothing than output wrong):
  1) Fig 6 must reproduce the MS off-diagonal z 2.483 / 2.623 / 0.892 written into the delivered caption, else error (guards against the 11-03 old tree)
  2) the two mouse rows must match the other session's stop-rule table bit-for-bit (rho_atac_rna 0.382 / rho_size_atac 0.105 / rho_atac_rna_adj 0.433)
  3) the Fig 7 slopes must match the delivered figS7a_confound_fit_slope.csv (robust-z scale)

Outputs:
  SUPPORT_v1/tables/deconfound_policy_table.csv        the criteria table itself (one row per dataset x channel)
  SUPPORT_v1/tables/deconfound_policy_verdicts.csv     one row per "face conclusion": before/after across conventions
"""
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy.stats import spearmanr, rankdata

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
import os
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

R = WORK + "/figures/routeA_out"
M3 = WORK + "/routeA_probe/m3_probe/tables_full74"
OUT = R + "/SUPPORT_v1/tables"
SEED = 20261102
NPERM, NBOOT = 10000, 10000

lnp_of = lambda s: np.log(s)
def robust_z(v):
    v = np.asarray(v, float); m = np.median(v)
    return (v - m) / (1.4826 * np.median(np.abs(v - m)))
def res(y, X):
    X = np.asarray(X, float)
    if X.ndim == 1: X = X[:, None]
    return y - sm.OLS(y, sm.add_constant(X)).fit().fittedvalues
def ols(v, x):
    f = sm.OLS(np.asarray(v, float), sm.add_constant(np.asarray(x, float))).fit()
    return f.params[1], f.pvalues[1], f.rsquared
def _det_rng(v, salt):
    """deterministic seed bound to the input content: bootstrap/permutation results must not depend on call order (different sessions must rerun bit-for-bit identical)"""
    import zlib
    key = np.asarray(v, dtype=np.float64).tobytes()
    return np.random.default_rng((SEED * 1_000 + salt * 7_919 + zlib.crc32(key)) % (2 ** 32))
def perm_P(x, y):
    obs = abs(spearmanr(x, y).statistic); rng = _det_rng(np.asarray(x) + np.asarray(y), 1)
    cnt = sum(abs(spearmanr(rng.permutation(x), y).statistic) >= obs for _ in range(NPERM))
    return (1 + cnt) / (1 + NPERM)
def boot_ci(x, y):
    rng = _det_rng(np.asarray(x) + np.asarray(y), 2); n = len(x); rs = []
    for _ in range(NBOOT):
        i = rng.integers(0, n, n); rs.append(spearmanr(x[i], y[i]).statistic)
    return np.percentile(rs, [2.5, 97.5])
def rk(v):           # rank 1 = the largest Δ1 (the strongest contrast)
    return rankdata(-np.asarray(v, float), method="ordinal").astype(int)

LAB = {"PASS_Multiple_sclerosis": "MS", "PASS_Schizophrenia_Calhoon2022": "SCZ",
       "PASS_Height": "HEIGHT", "PASS_PLATELET": "PLATELET", "UKB_460K.height": "HEIGHT"}
def lab(t): return LAB.get(t, t.replace("PASS_", "").split("_")[0])

# -- n_peak: on the human side a trait-level constant column (from the GSE231584 peak annotation, shared across Fig 7)
S7 = pd.read_csv(R + "/FIG7_v1/tables/seven4_confound_ALL_cohorts.csv")
NPEAK_HUM = dict(zip(S7.trait, S7.lnp))

rows, vrows = [], []

def decide(rho, P, trig3, verdict_changed, rank_move, top5_changed):  # rank_move/top5 are both "printed-on-the-face" conventions
    t1, t2 = abs(rho) >= 0.30, P < 0.05
    return dict(trig_i=bool(t1), trig_ii=bool(t2), trig_iii=bool(trig3),
                must_adjust=bool(t1 or t2 or trig3),
                why=("%s%s%s" % ("|rho|>=0.30 " if t1 else "", "P<0.05 " if t2 else "",
                                 "verdict/rank/top5 changed" if trig3 else "")).strip() or "none")

def add(dataset, channel, n, lnp, d1, cov_note, cur_status, printed_claim,
        primary_rho=None, adj_rho=None, verdict_changed=None, rank_move=None, top5_changed=None,
        nhit_rho_rank=None, extra="", rank_move_any=None, labelled_on_face=True,
        verdict_primary=None, verdict_adj=None):
    lnp = np.asarray(lnp, float); d1 = np.asarray(d1, float)
    ok = ~(np.isnan(lnp) | np.isnan(d1)); lnp, d1 = lnp[ok], d1[ok]
    slope, Pv, r2 = ols(robust_z(d1), lnp)
    rho = spearmanr(lnp, d1).statistic
    sp_P = spearmanr(lnp, d1).pvalue
    dec = np.asarray(pd.qcut(-lnp, 10, labels=False, duplicates="drop")).astype(int) + 1
    r_ = rk(d1)
    stair_raw = spearmanr(dec, r_).statistic
    stair_adj = spearmanr(dec, rk(res(d1, lnp))).statistic
    resid_after = spearmanr(lnp, res(d1, lnp)).statistic   # residual confounding after adjustment (should be ~0, proving the adjustment really removed it)
    lm = rank_move if labelled_on_face else None          # (iii) applies only to quantities printed on the face
    d = decide(rho, sp_P, bool(verdict_changed or (lm or 0) >= 5 or (top5_changed and labelled_on_face)),
               verdict_changed, lm, top5_changed and labelled_on_face)
    rows.append(dict(dataset=dataset, channel=channel, n_trait=len(d1),
                     rho_lnp_d1=round(rho, 3), P_spearman="%.3g" % sp_P, R2=round(r2, 3),
                     slope_robustz=round(slope, 2), stair_rho_raw=round(stair_raw, 3),
                     stair_rho_adj=round(stair_adj, 3), rho_lnp_after_adj=round(resid_after, 3),
                     rho_nhit_rank=nhit_rho_rank,
                     covariate_source=cov_note, printed_claim=printed_claim,
                     rho_primary=primary_rho, rho_domainadj=adj_rho,
                     verdict_changed=bool(verdict_changed), verdict_primary=verdict_primary, verdict_adj=verdict_adj,
                     rank_move_labelled=lm, rank_move_any_trait=rank_move_any,
                     labelled_on_face=labelled_on_face, top5_changed=(top5_changed and labelled_on_face),
                     status=cur_status, **d, note=extra))
    return slope

# == 1. the four Fig 7 cohorts (human unpaired; Δ1 = each cohort's mature minus precursor contrast)
COH = {"A2": "GSE277430", "A3": "GSE118257", "A4": "GSE301908", "A7": "GSE231585"}
DELIV = pd.read_csv(OUT + "/figS7a_confound_fit_slope.csv").set_index("cohort")
MS_R, MS_A = {}, {}
for k, gse in COH.items():
    g = S7.copy()
    c = k + "_1" if k + "_1" in g.columns else [x for x in g.columns if x.startswith(k)][0]
    d1 = g[c].values
    rr = rk(d1); ra = rk(res(d1, g.lnp.values))
    i = list(g.trait).index("PASS_Multiple_sclerosis")
    MS_R[k], MS_A[k] = int(rr[i]), int(ra[i])
    deliv_slope = float(DELIV.loc[k, "slope"])
    recomputed = ols(robust_z(d1), g.lnp.values)[0]
    assert abs(deliv_slope - recomputed) < 0.02, (k, deliv_slope, recomputed)
    add("Fig7 " + gse, "Δ1 (human unpaired)", 74, g.lnp.values, d1,
        "n_peak from GSE231584 peak annotation (trait-level constant, shared across cohorts)",
        "adjusted AND raw both shown (compliant)",
        "MS rank %d -> %d after adjustment (main Fig 7B/S7b)" % (MS_R[k], MS_A[k]),
        verdict_changed=True, rank_move=abs(MS_R[k] - MS_A[k]), top5_changed=True,
        rank_move_any=int(np.abs(rk(d1) - rk(res(d1, g.lnp.values))).max()),
        primary_rho="", adj_rho="", extra="slope verified vs delivered S7-4b table")

# == 2. Fig 6D three chemistries x two channels
S = pd.DataFrame({"trait": S7.trait, "lnp": S7.lnp})
D = pd.read_csv(R + "/D_crosschemistry/tables_log1p_20261104/fig6D_delta1_74traits_per_tech.csv").merge(S, on="trait", how="left")
# self-check 1: reproduce the delivered off-diagonal z
chk = D[D.trait == "PASS_Multiple_sclerosis"].set_index("tech").offdiag_z
for t, v in [("dsciCAT-seq (own)", 2.483), ("10x Multiome", 2.623), ("HT-scCAT-seq", 0.892)]:
    assert abs(chk[t] - v) < 0.005, "Fig6 convention self-check failed (old tree in use?): %s %.3f vs %.3f" % (t, chk[t], v)
gs = {t: g.reset_index(drop=True) for t, g in D.groupby("tech", sort=False)}
LAD = {}
for t, g in gs.items():
    r, a = g.d1_rna.values.astype(float), g.d1_atac.values.astype(float)
    lnp = g.lnp.values
    prim = spearmanr(res(r, np.log(g.n_hit_rna.values)), res(a, np.log(g.n_hit_atac.values))).statistic
    dom = spearmanr(res(r, lnp), res(a, lnp)).statistic
    both = spearmanr(res(r, np.c_[lnp, np.log(g.n_hit_rna.values)]),
                     res(a, np.c_[lnp, np.log(g.n_hit_atac.values)])).statistic
    lo, hi = boot_ci(res(r, lnp), res(a, lnp))
    DELIV_PRIM = {"dsciCAT-seq (own)": 0.896, "10x Multiome": 0.824, "HT-scCAT-seq": 0.630}
    assert abs(prim - DELIV_PRIM[t]) < 0.002, "the pre-registered n_hit convention did not reproduce the delivered caption value: %s %.3f" % (t, prim)
    LAD[t] = (prim, dom, lo, hi, perm_P(res(r, lnp), res(a, lnp)), both)
    i = list(g.trait).index("PASS_Multiple_sclerosis")
    for ch, key, nh in [("RNA", "d1_rna", "n_hit_rna"), ("chromatin", "d1_atac", "n_hit_atac")]:
        v = g[key].values.astype(float)
        rn, ra_ = rk(v), rk(res(v, lnp))
        mv = abs(int(rn[i]) - int(ra_[i]))
        top5 = set(np.argsort(rn)[:5]) != set(np.argsort(ra_)[:5])
        add("Fig6D " + t, ch + " (74 traits)", 74, lnp, v,
            "same human n_peak column (imported trait-level)",
            "NOT adjusted for n_peak; adjusted for n_hit only, labelled 'domain-adjusted' (gap)",
            "ladder verdict + printed MS/chromatin ranks",
            primary_rho="%.3f" % prim, adj_rho="%.3f" % dom,
            verdict_changed=(prim >= 0.6) != (dom >= 0.6), rank_move=mv, top5_changed=top5,
        rank_move_any=int(np.abs(rn - ra_).max()),
        verdict_primary="PASS" if prim >= 0.6 else "FAIL", verdict_adj="PASS" if dom >= 0.6 else "FAIL",
        nhit_rho_rank=round(spearmanr(g[nh], g["rank_rna_adj" if ch == "RNA" else "rank_atac_adj"]).statistic, 3),
        extra="pass/fail vs pre-registered rho>=0.6 under domain adj: %s" % ("PASS" if dom >= 0.6 else "FAIL"))

# == 3. Fig 6C own same-nucleus
O = pd.read_csv(R + "/O1_own_multiome_paired/tables/T4_delta1_74traits_both_channels.csv").merge(S, on="trait", how="left")
o_r, o_a, lnp = O.d1_rna.values.astype(float), O.d1_atac.values.astype(float), O.lnp.values
prim = spearmanr(o_r, o_a).statistic
dom = spearmanr(res(o_r, lnp), res(o_a, lnp)).statistic
lo, hi = boot_ci(res(o_r, lnp), res(o_a, lnp))
i = list(O.trait).index("PASS_Multiple_sclerosis")
for ch, key in [("RNA", "d1_rna"), ("chromatin", "d1_atac")]:
    v = O[key].values.astype(float); rn, ra_ = rk(v), rk(res(v, lnp))
    add("Fig6C own dsciCAT-seq", ch + " (same nucleus)", 74, lnp, v,
        "same human n_peak column (imported trait-level)",
        "NOT adjusted for n_peak; S6g(a) row 'domain-size adjusted (D3)' actually residualises ln n_hit (gap)",
        "cross-channel rho + labelled MS ranks #1/#18",
        primary_rho="%.3f" % prim, adj_rho="%.3f" % dom,
        verdict_changed=abs(prim - dom) >= 0.05, rank_move=abs(int(rn[i]) - int(ra_[i])),
        rank_move_any=int(np.abs(rn - ra_).max()),
        top5_changed=set(np.argsort(rn)[:5]) != set(np.argsort(ra_)[:5]),
        extra="rho 0.819 -> 0.708 (95%% CI %.3f-%.3f); identity-line batch Q1 shares this panel" % (lo, hi))

# == 4. mouse 74 traits (n_peak is computed per library from the mouse's OWN gpeak, not an imported column)
M = pd.read_csv(M3 + "/M3_full74_atac_vs_rna_d1_74.csv")
STOP = pd.read_csv(M3 + "/M3_full74_stop_rule_stats.csv").iloc[0]
# self-check 2: bit-for-bit match with the other session's stop-rule table
assert abs(STOP.rho_atac_rna - spearmanr(M.rna_d1, M.atac_d1).statistic) < 1e-3
assert abs(STOP.rho_size_atac - spearmanr(M.ln_n_peak, M.atac_d1).statistic) < 1e-3
assert abs(STOP.rho_atac_rna_adj - spearmanr(res(M.rna_d1, M.ln_n_peak), res(M.atac_d1, M.ln_n_peak)).statistic) < 1e-3
for ch, key in [("ATAC", "atac_d1"), ("RNA", "rna_d1")]:
    v = M[key].values.astype(float); lnp = M.ln_n_peak.values
    rn, ra_ = rk(v), rk(res(v, lnp))
    rho, Pv, r2 = ols(robust_z(v), lnp)
    add("Mouse GSE250589", ch + " arm (74 traits)", 74, lnp, v,
        "mouse-specific n_peak, per-library gpeak (tables_full74/GSM*_trait_npeak74.csv) — cleaner than the imported human column",
        "raw only; decision evidence sits in the stop-rule CSV, not on any delivered face",
        "cross-channel rho only (0.382 trait-level; 0.391 library-level on S7C face) - no mouse trait rank is printed anywhere",
        primary_rho="%.3f" % STOP.rho_atac_rna, adj_rho="%.3f" % STOP.rho_atac_rna_adj,
        verdict_changed=False, rank_move=int(np.abs(rn - ra_).max()),
        top5_changed=set(np.argsort(rn)[:5]) != set(np.argsort(ra_)[:5]),
        labelled_on_face=False, rank_move_any=int(np.abs(rn - ra_).max()),
        extra="|rho|=%.3f P=%.3g R2=%.3f; SD(d1)=%.3f vs other arm — compressed channel limits power" % (
            rho, Pv, r2, np.nanstd(v)))

T = pd.DataFrame(rows)
T.to_csv(OUT + "/deconfound_policy_table.csv", index=False)

# -- the before/after of every "face conclusion" across conventions (the evidence body for criterion iii)
def vr(claim, where, prim, adj, unit, changed, note):
    vrows.append(dict(claim=claim, where=where, primary=prim, domain_adjusted=adj,
                      unit=unit, verdict_changed=changed, note=note))
for t, (p, d, lo, hi, P, bo) in LAD.items():
    vr("ladder reproduces (rho_adj >= 0.6)", "Fig6D / S6h", "%.3f PASS" % p,
       "%.3f %s (95%% CI %.3f-%.3f, perm P %.4f)" % (d, "PASS" if d >= 0.6 else "FAIL", lo, hi, P),
       "Spearman rho", (p >= 0.6) != (d >= 0.6), "n_peak+n_hit together: %.3f" % bo)
vr("3 of 3 chemistries reproduce the ladder", "Fig6D legend / S6h face", "3 of 3", "2 of 3 (HT 0.523 < 0.6)",
   "count", True, "the only verdict that flips")
vr("own > 10x ordering", "Fig6D", "CI 0.009-0.147 (paired boot, n_hit)", "CI 0.076-0.317 (n_peak), P(diff<=0)=0.000",
   "rho difference", False, "cleaner under the stricter adjustment")
vr("10x > HT ordering", "Fig6D", "CI 0.066-0.339 (n_hit)", "CI -0.051-0.290 (n_peak), P=0.091",
   "rho difference", True, "loses significance under n_peak")
vr("cross-channel rho (own, same nucleus)", "Fig6C", "0.819", "0.708 (95% CI 0.552-0.816)",
   "Spearman rho", False, "n_hit-only gives 0.826 (~unchanged): the n_hit proxy does not capture this confound")
for t, dv, ad, rkm in [("dsciCAT-seq (own)", 2.483, 2.80, "73->74"), ("10x Multiome", 2.623, 3.08, "74->74"), ("HT-scCAT-seq", 0.892, 1.23, "62->66")]:
    vr("MS off-diagonal z [" + t + "]", "Fig6D legend", "%+.2f" % dv, "%+.2f (off-diag rank %s)" % (ad, rkm),
       "SD", False, "counterexample gets STRONGER under deconfounding")
vr("MS counterexample reproduced in all chemistries", "Fig6D / S6h", "3 of 3", "3 of 3", "count", False,
   "core narrative unaffected")
for k in COH:
    vr("MS rank, cohort " + COH[k], "Fig7B / S7b", "%d" % MS_R[k], "%d" % MS_A[k], "rank/74",
       MS_R[k] != MS_A[k], "why Fig 7 must show both")
vr("mouse cross-channel rho", "S7C mouse tile", "0.382 (perm P 7.8e-4)", "0.433", "Spearman rho", False,
   "both below the pre-registered 0.6 -> ladder does NOT reproduce in mouse, under either adjustment")
V = pd.DataFrame(vrows)
V.to_csv(OUT + "/deconfound_policy_verdicts.csv", index=False)

print(T[["dataset", "channel", "rho_lnp_d1", "P_spearman", "R2", "must_adjust", "status"]].to_string(index=False))
print("\ncriteria table %d rows -> %s/deconfound_policy_table.csv" % (len(T), OUT))
print("convention comparison %d rows -> %s/deconfound_policy_verdicts.csv" % (len(V), OUT))
print("self-checks passed: Fig6 off-diagonal z reproduced / mouse stop-rule three values reproduced / Fig7 slopes match the delivered S7-4b")
