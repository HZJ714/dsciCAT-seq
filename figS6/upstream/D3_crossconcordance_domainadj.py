# -*- coding: utf-8 -*-
"""D3: cross-channel / cross-cohort Δ1 concordance - the formal ledger table, raw vs adjusted for gene-set regulatory-domain size (ln n_peak).
Why it is mandatory:
  - Fig 7 receipt §5-C / §7 requires "raw -> adjusted for domain size: does it survive?" as a checkable table (10 human pairs + 1 own same-nucleus row);
  - Fig 6 receipt §10.3 and Fig 7 receipt §10.5 retract the claim "74-trait ranking is highly concordant across channels",
    the own same-nucleus O1 row is the only positive; it may go on a figure only after recomputation with the exact adjustment used for the retraction.
Conventions (identical to receipt §23):
  - human: OLS-residualize Δ1 on ln(n_peak_raw), both sides residualized separately, then Spearman;
  - O1 same-nucleus: both sides residualized on ln(n_hit_atac) (n_gs is the constant 1000 across the 74 traits, so it cannot confound);
  - permutation: random permutation of one adjusted vector 10,000 times (seed=20261102), P = (1+#{|rho_perm|>=|rho|})/(1+n_perm);
  - bootstrap: 10,000 trait resamples (same seed), 2.5/97.5 percentiles.
Output D_crosschemistry/tables/D3_crossconcordance_domainadjust.csv
"""
import numpy as np, pandas as pd, os
from scipy.stats import spearmanr, norm
from statsmodels.regression.linear_model import OLS
from statsmodels.tools import add_constant as AC

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
import os
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

BASE = WORK + "/figures/routeA_out"
OUT = os.path.join(BASE, "D_crosschemistry", "tables")
os.makedirs(OUT, exist_ok=True)
SEED, NPERM, NBOOT = 20261102, 10000, 10000

s4 = pd.read_csv(f"{BASE}/FIG7_v1/tables/seven4_confound_ALL_cohorts.csv")
pk = pd.read_csv(f"{BASE}/A6_GSE231584_atac_v1/tables/delta1_vs_peakdomain_allchannels.csv")
t4 = pd.read_csv(f"{BASE}/O1_own_multiome_paired/tables/T4_delta1_74traits_both_channels.csv")

hum = s4.merge(pk[["trait", "median_groupdiff"]], on="trait", how="left")
hum["A6"] = hum["median_groupdiff"]
assert hum[["A2", "A3", "A4", "A7", "A6", "n_peak_raw"]].notna().all().all(), "the human table has gaps"
hum["lnp"] = np.log(hum["n_peak_raw"])

COH = ["A2", "A3", "A4", "A7"]
LBL = {"A2": "GSE277430", "A3": "GSE118257", "A4": "GSE301908", "A7": "GSE231585", "A6": "GSE231584(snATAC)"}


def resid(y, x):
    return np.asarray(OLS(np.asarray(y, float), AC(np.asarray(x, float))).fit().resid)


def one_pair(x, y, adj_x, name, kind, ms_rank=None):
    n = len(x)
    r_raw, p_raw = spearmanr(x, y)
    rx = resid(x, adj_x) if adj_x is not None else np.asarray(x, float)
    ry = resid(y, adj_x) if adj_x is not None else np.asarray(y, float)
    r_adj, p_adj = spearmanr(rx, ry)
    rng = np.random.default_rng(SEED)
    obs = abs(r_adj)
    cnt = 0
    for _ in range(NPERM):
        rr, _ = spearmanr(rx, rng.permutation(ry))
        if abs(rr) >= obs - 1e-15:
            cnt += 1
    perm_P = (1 + cnt) / (1 + NPERM)
    rng2 = np.random.default_rng(SEED + 1)
    bs = []
    for _ in range(NBOOT):
        i = rng2.integers(0, n, n)
        bs.append(spearmanr(rx[i], ry[i])[0])
    lo, hi = np.percentile(bs, [2.5, 97.5])
    z = 0.5 * np.log((1 + r_raw) / (1 - r_raw))
    se = 1 / np.sqrt(n - 3)
    return dict(pair=name, kind=kind, n=n,
                rho_raw=r_raw, P_raw=p_raw,
                rho_adj=r_adj, P_adj=p_adj,
                perm_P_adj=perm_P, n_perm=NPERM, seed=SEED,
                ci_lo_adj=lo, ci_hi_adj=hi,
                ci_lo_raw_fisher=z - 1.96 * se, ci_hi_raw_fisher=z + 1.96 * se,
                rho_retained=r_adj / r_raw if r_raw else np.nan,
                MS_rank_raw=ms_rank.get("raw") if ms_rank else "",
                MS_rank_adj=ms_rank.get("adj") if ms_rank else "")


rows = []
# --- the 6 RNA-to-RNA pairs (between cohorts) ---
for i in range(len(COH)):
    for j in range(i + 1, len(COH)):
        a, b = COH[i], COH[j]
        sub = hum[[a, b, "lnp"]].dropna()
        msr = {"raw": int(s4.loc[s4.trait.str.contains("Multiple_sclerosis"), f"rank_raw_{a}"].iloc[0]),
               "adj": int(s4.loc[s4.trait.str.contains("Multiple_sclerosis"), f"rank_adj_{a}"].iloc[0])}
        rows.append(one_pair(sub[a], sub[b], sub["lnp"], f"{LBL[a]} ↔ {LBL[b]}", "RNA↔RNA", msr))
# --- the 4 RNA-to-chromatin pairs (A6, per dataset) ---
for a in COH:
    sub = hum[[a, "A6", "lnp"]].dropna()
    rows.append(one_pair(sub[a], sub["A6"], sub["lnp"], f"{LBL[a]} <-> {LBL['A6']}", "RNA<->chromatin (separate libraries)"))
# --- own same-nucleus RNA<->ATAC ---
o = t4[["d1_rna", "d1_atac", "n_hit_rna", "n_hit_atac"]].dropna()
ms = t4[t4.trait.str.contains("Multiple_sclerosis")].iloc[0]
rows.append(one_pair(o["d1_rna"], o["d1_atac"], np.log(o["n_hit_atac"]),
                     "O1 own dsciCAT-seq same-nucleus RNA<->gene score", "RNA<->chromatin (same-nucleus paired)",
                     {"raw": int(ms["rank_rna_raw"]), "adj": int(ms["rank_atac_raw"])}))
# robustness: both sides residualized on ln n_hit_rna instead; and no adjustment (equivalent to raw)
rows.append(one_pair(o["d1_rna"], o["d1_atac"], np.log(o["n_hit_rna"]),
                     "O1 same-nucleus (adjustment variable switched to ln n_hit_rna)", "sensitivity"))
# the domain-size vs Δ1 correlation (each channel's own diagnostic row)
for a in COH + ["A6"]:
    r, p = spearmanr(hum["lnp"], hum[a])
    rows.append(dict(pair=f"{LBL[a]}: Spearman(ln n_peak, Δ₁)", kind="confound-check", n=len(hum),
                     rho_raw=r, P_raw=p, rho_adj=np.nan, P_adj=np.nan, perm_P_adj=np.nan,
                     n_perm=NPERM, seed=SEED))
r, p = spearmanr(np.log(o["n_hit_atac"]), o["d1_atac"])
rows.append(dict(pair="O1 same-nucleus: Spearman(ln n_hit_atac, Δ1_atac)", kind="confound-check", n=len(o),
                 rho_raw=r, P_raw=p, rho_adj=np.nan, P_adj=np.nan, perm_P_adj=np.nan,
                 n_perm=NPERM, seed=SEED))

d = pd.DataFrame(rows)
d.to_csv(f"{OUT}/D3_crossconcordance_domainadjust.csv", index=False)
pd.set_option("display.width", 250, "display.max_columns", 30)
print(d[["pair", "kind", "n", "rho_raw", "rho_adj", "P_adj", "perm_P_adj", "ci_lo_adj", "ci_hi_adj", "rho_retained"]].to_string(index=False))
print("\nwritten:", f"{OUT}/D3_crossconcordance_domainadjust.csv", d.shape)
