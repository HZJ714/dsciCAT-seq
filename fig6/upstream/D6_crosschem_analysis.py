#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D6: the formal statistics of Fig 6 D - the same Δ1 cross-channel test rerun on each of the three library chemistries.

Criteria (syntax identical to C / D3, nothing new):
  Δ1(trait, tech, channel) = median(norm_score | Oligo) - median(| OPC)   <- Fig 2's cell_type_updated
  rho_raw = Spearman(Δ1_rna, Δ1_atac) over the 74 traits
  rho_adj = Spearman after OLS-residualizing both sides on ln(n_hit) (n_hit = genes of that trait's set measurable in that channel)
            - the same function and seed as D3_crossconcordance_domainadjust.py
  perm_P_adj: one-side random permutation 10,000 times after adjustment, seed=20261102, P=(1+#{|rho|>=obs})/(1+n)
  CI: bootstrap by trait 10,000 times, 2.5/97.5 percentiles
  noise ceiling: the geometric mean of within-channel split-half reliabilities (20 splits, each half of Oligo/OPC)
  off-diagonal: the residual z of a Huber robust regression Δ1_rna ~ Δ1_atac and its rank within the 74 (same method as O1/C)
  effect size: Cliff's δ of Oligo vs OPC per tech per channel (the MS row + the 74-trait median)
          - across chemistries ONLY effect sizes and ranks are compared, never P values (review receipt §10.10)
Pre-registered criteria (fixed before the run, handover sheet §3):
  ladder reproduces    = rho_adj >= 0.6 and perm_P_adj < 0.05 (all three must pass for "reproduced")
  MS counter-example reproduces = rank_atac_adj > rank_rna_adj and offdiag_z > 0

V2 bridge (handover sheet §2, mandatory): recompute the dsciCAT arm with C's original criteria (MOL6 vs COP/NFOL, labels from
C's h5ad obs["cell_type"], intersection of both batches' nuclei) and compare the ladder ρ with the delivered 0.819;
within 0.819 ± 0.05 => swapping the matrix does not move the conclusion; report as the same family.
"""
import os, glob, json, shutil, numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.regression.linear_model import OLS
from statsmodels.tools import add_constant as AC
import statsmodels.api as sm
import anndata as ad

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

OUT = os.path.join(os.environ.get('DSCI_FIG_ROOT', WORK + "/figures"), 'routeA_out/D_crosschemistry')
# DSCI_SCDRS_SUBDIR / DSCI_TAB_SUBDIR let a corrected rerun (h5ad_log1p -> scdrs_log1p_*) be analysed
# into a parallel table directory without overwriting the delivered numbers. defaults = delivered tree.
SCDRS_SUBDIR = os.environ.get('DSCI_SCDRS_SUBDIR', 'scdrs')
# 2026-11-04 corrected re-run: only the ATAC arm needs recomputation from h5ad_log1p -> scdrs_log1p_*; the RNA arm stays in the delivered tree.
# so the ATAC arm can point at the re-run tree alone while the RNA arm stays under DSCI_SCDRS_SUBDIR.
SCDRS_SUBDIR_ATAC = os.environ.get('DSCI_SCDRS_SUBDIR_ATAC', SCDRS_SUBDIR)
TAB, DER = OUT + "/" + os.environ.get('DSCI_TAB_SUBDIR', 'tables'), OUT + "/derived"
os.makedirs(TAB, exist_ok=True)
if TAB != OUT + "/tables":
    for _f in os.listdir(OUT + "/tables"):
        if _f.startswith("cells_") and not os.path.exists(os.path.join(TAB, _f)):
            shutil.copy(os.path.join(OUT + "/tables", _f), os.path.join(TAB, _f))
GSF = os.path.join(os.environ.get('DSCI_FY_ROOT', WORK + "/FY"), 'oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs')
H5C = os.path.join(os.environ.get('DSCI_FY_ROOT', WORK + "/FY"), 'oligo/3k_dsci_oligo_subset.h5ad')
T6C = os.path.join(os.environ.get('DSCI_FIG_ROOT', WORK + "/figures"), 'routeA_out/O1_own_multiome_paired/tables/T6_cross_channel_summary.csv')
SEED, NPERM, NBOOT, NHALF = 20261102, 10000, 10000, 20
HOM = os.environ.get("DSCI_HOMOLOGENE", WORK + "/routeA_probe/homologene.data")
_hg = pd.read_csv(HOM, sep="\t", header=None, names=["hid", "tax", "gid", "symbol", "gi", "acc"])
_hh = _hg[_hg.tax == 9606].set_index("hid").symbol.rename("human")
_mm = _hg[_hg.tax == 10090].set_index("hid").symbol.rename("mouse")
_pair = pd.concat([_hh.groupby(level=0).first(), _mm.groupby(level=0).first()], axis=1).dropna()
_pair["human"] = _pair["human"].str.upper()
H2M = _pair.groupby("human").mouse.apply(list).to_dict()
MAT, OPC, COP = "Oligo", "OPC", "COP/NFOL"
TECHS = ["dsciCAT_seq", "10x", "HT_scCAT_seq_HT"]
LBL = {"dsciCAT_seq": "dsciCAT-seq (own)", "10x": "10x Multiome", "HT_scCAT_seq_HT": "HT-scCAT-seq"}
SAY = lambda fmt, *a: print("[%s] %s" % (pd.Timestamp.now().strftime("%H:%M:%S"), fmt % a if a else fmt), flush=True)


def load_scores(job):
    sub = SCDRS_SUBDIR_ATAC if job.endswith("_atac") else SCDRS_SUBDIR
    out = {}
    for f in glob.glob(os.path.join(OUT, sub, job, "**", "*.score.gz"), recursive=True):
        if f.endswith(".full_score.gz"):
            continue
        tr = os.path.basename(f).replace(".score.gz", "")
        if tr in out:
            continue
        out[tr] = pd.read_csv(f, sep="\t", index_col=0)["norm_score"]
    return pd.DataFrame(out)


def resid(y, x):
    return np.asarray(OLS(np.asarray(y, float), AC(np.asarray(x, float))).fit().resid)


def cliffs(x, y):
    x, y = np.asarray(x), np.asarray(y)
    if len(x) == 0 or len(y) == 0:
        return np.nan
    return float((x[:, None] > y[None, :]).mean() - (x[:, None] < y[None, :]).mean())


def ladder(x, y, adj_x, adj_y, rng):
    n = len(x)
    r_raw, p_raw = spearmanr(x, y)
    rx, ry = resid(x, adj_x), resid(y, adj_y)
    r_adj, p_adj = spearmanr(rx, ry)
    rr = np.random.default_rng(SEED)
    obs = abs(r_adj)
    cnt = sum(1 for _ in range(NPERM) if abs(spearmanr(rx, rr.permutation(ry))[0]) >= obs - 1e-15)
    perm_P = (1 + cnt) / (1 + NPERM)
    rb = np.random.default_rng(SEED + 1)
    bs = [spearmanr(rx[i], ry[i])[0] for i in (rb.integers(0, n, n) for _ in range(NBOOT))]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    z = 0.5 * np.log((1 + r_raw) / (1 - r_raw)); se = 1 / np.sqrt(n - 3)
    return dict(rho_raw=r_raw, P_raw=p_raw, rho_adj=r_adj, P_adj=p_adj, perm_P_adj=perm_P,
                ci_lo_adj=lo, ci_hi_adj=hi, ci_lo_raw_fisher=float(np.tanh(z - 1.96 * se)), ci_hi_raw_fisher=float(np.tanh(z + 1.96 * se)),
                rho_retained=r_adj / r_raw if r_raw else np.nan, n_perm=NPERM, seed=SEED)


def split_half_ceiling(Xr, Xa, m_mat, m_opc):
    rng = np.random.default_rng(SEED + 7)
    im, io = np.where(m_mat)[0], np.where(m_opc)[0]
    vals = {"rna": [], "atac": []}
    for _ in range(NHALF):
        for key, X in [("rna", Xr), ("atac", Xa)]:
            pm, po = rng.permutation(im), rng.permutation(io)
            h1m, h2m = pm[:len(pm) // 2], pm[len(pm) // 2:]
            h1o, h2o = po[:len(po) // 2], po[len(po) // 2:]
            d1 = np.median(X[h1m], 0) - np.median(X[h1o], 0)
            d2 = np.median(X[h2m], 0) - np.median(X[h2o], 0)
            vals[key].append(spearmanr(d1, d2).correlation)
    rel = {k: float(np.median(v)) for k, v in vals.items()}
    return float(np.sqrt(max(rel["rna"], 1e-6) * max(rel["atac"], 1e-6))), rel


def n_hit_per_trait(genes):
    gs = {}
    for l in open(GSF):
        f = l.rstrip("\n").split("\t")
        if len(f) < 2 or f[0] == "TRAIT":
            continue
        gs[f[0]] = [x.strip().upper() for x in f[1].split(",") if x.strip()]
    have = set(g.upper() for g in genes)
    rows = []
    for t, g in gs.items():
        n_mat = sum(any((m.upper() in have) for m in H2M.get(gg, [])) for gg in g)
        rows.append({"trait": t, "n_gs": len(g), "n_hit": n_mat})
    return pd.DataFrame(rows).set_index("trait")


lad, msr, mg, D1, extra = [], [], [], [], {}
for tag in TECHS:
    cells = pd.read_csv(os.path.join(TAB, "cells_%s.csv" % tag)).set_index("cell")
    Sr, Sa = load_scores(tag + "_rna"), load_scores(tag + "_atac")
    if len(Sr) == 0 or len(Sa) == 0:
        SAY("!! %s scores not on disk yet (rna %d / atac %d) - skipping", tag, len(Sr), len(Sa)); continue
    common = [c for c in cells.index if c in Sr.index and c in Sa.index]
    traits = sorted(set(Sr.columns) & set(Sa.columns))
    Sr, Sa = Sr.loc[common, traits], Sa.loc[common, traits]
    lab = cells.loc[common, "cell_type_updated"] if "cell_type_updated" in cells else cells.loc[common, "label"]
    m_mat, m_opc = (lab == MAT).values, (lab == OPC).values
    Xr, Xa = Sr.values, Sa.values
    SAY("%s | nuclei %d (Oligo %d / OPC %d / COP %d) | traits %d", LBL[tag], len(common), m_mat.sum(), m_opc.sum(),
        (lab == COP).sum(), len(traits))
    d1r = np.median(Xr[m_mat], 0) - np.median(Xr[m_opc], 0)
    d1a = np.median(Xa[m_mat], 0) - np.median(Xa[m_opc], 0)
    gr = [l.strip() for l in open(os.path.join(DER, "rna_%s_genes.txt" % tag)) if l.strip()]
    ga = [l.strip() for l in open(os.path.join(DER, "gs_%s_genes.txt" % tag)) if l.strip()]
    nr, na = n_hit_per_trait(gr).reindex(traits), n_hit_per_trait(ga).reindex(traits)
    L = ladder(d1r, d1a, np.log(nr.n_hit.values), np.log(na.n_hit.values), None)
    ceil_, rel = split_half_ceiling(Xr, Xa, m_mat, m_opc)
    ok = bool(L["rho_adj"] >= 0.6 and L["perm_P_adj"] < 0.05)
    lad.append(dict(tech=LBL[tag], n_oligo=int(m_mat.sum()), n_opc=int(m_opc.sum()), n_cop=int((lab == COP).sum()),
                    n_traits=len(traits), **L, noise_ceiling=ceil_, rel_rna=rel["rna"], rel_atac=rel["atac"],
                    ladder_reproduced=ok))
    T = pd.DataFrame({"trait": traits, "d1_rna": d1r, "d1_atac": d1a,
                      "n_hit_rna": nr.n_hit.values, "n_hit_atac": na.n_hit.values})
    rr, ra_ = resid(T.d1_rna, np.log(T.n_hit_rna)), resid(T.d1_atac, np.log(T.n_hit_atac))
    T["rank_rna_raw"] = pd.Series(T.d1_rna).rank(ascending=False).astype(int).values
    T["rank_atac_raw"] = pd.Series(T.d1_atac).rank(ascending=False).astype(int).values
    T["rank_rna_adj"] = pd.Series(rr).rank(ascending=False).astype(int).values
    T["rank_atac_adj"] = pd.Series(ra_).rank(ascending=False).astype(int).values
    fit = sm.RLM(T.d1_rna.values, AC(T.d1_atac.values), M=sm.robust.norms.HuberT()).fit()
    off = T.d1_rna.values - fit.fittedvalues
    T["offdiag_resid"] = off; T["offdiag_z"] = (off - off.mean()) / off.std()
    T["offdiag_rank"] = pd.Series(off).rank(ascending=False).astype(int).values
    # effect sizes (only these are compared, never P)
    T["cliff_delta_rna"] = [cliffs(Xr[m_mat, j], Xr[m_opc, j]) for j in range(len(traits))]
    T["cliff_delta_atac"] = [cliffs(Xa[m_mat, j], Xa[m_opc, j]) for j in range(len(traits))]
    T.insert(1, "tech", LBL[tag]); T["n_nuclei"] = len(common)
    i = int(np.where(T.trait == "PASS_Multiple_sclerosis")[0][0])
    msr.append(dict(tech=LBL[tag], trait="PASS_Multiple_sclerosis", d1_rna=T.d1_rna[i], d1_atac=T.d1_atac[i],
                    cliff_rna=T.cliff_delta_rna[i], cliff_atac=T.cliff_delta_atac[i],
                    rank_rna_raw=int(T.rank_rna_raw[i]), rank_rna_adj=int(T.rank_rna_adj[i]),
                    rank_atac_raw=int(T.rank_atac_raw[i]), rank_atac_adj=int(T.rank_atac_adj[i]),
                    offdiag_z=T.offdiag_z[i], offdiag_rank=int(T.offdiag_rank[i]),
                    ms_counterexample_reproduced=bool(T.rank_atac_adj[i] > T.rank_rna_adj[i] and T.offdiag_z[i] > 0),
                    direction=("chromatin ranks MS lower (counter-example reproduced)" if T.rank_atac_adj[i] > T.rank_rna_adj[i] else
                               "chromatin does not rank MS lower (counter-example not reproduced)")))
    D1.append(T)
    mg.append(pd.DataFrame({"trait": traits, "tech": LBL[tag], "n_genes_rna": len(gr), "n_genes_atac": len(ga),
                            "n_hit_rna": nr.n_hit.values, "n_hit_atac": na.n_hit.values}))
    extra[LBL[tag]] = dict(n_cells=len(common), n_traits=len(traits))

if not D1:
    SAY("no scores on disk yet - run this after D5 finishes"); raise SystemExit(0)
DL = pd.DataFrame(lad); DM = pd.DataFrame(msr); DT = pd.concat(D1, ignore_index=True); MG = pd.concat(mg, ignore_index=True)
gr = {t: set(x.strip() for x in open(os.path.join(DER, "rna_%s_genes.txt" % t)) if x.strip()) for t in TECHS}
ga = {t: set(x.strip() for x in open(os.path.join(DER, "gs_%s_genes.txt" % t)) if x.strip()) for t in TECHS}
MG["n_genes_all3_rna"] = len(set.intersection(*gr.values()))
MG["n_genes_all3_atac"] = len(set.intersection(*ga.values()))
SAY("\n=== ladder ===\n" + DL.round(4).to_string(index=False))
SAY("\n=== the MS row ===\n" + DM.round(4).to_string(index=False))
DL.to_csv(os.path.join(TAB, "fig6D_crosschem_ladder.csv"), index=False)
DM.to_csv(os.path.join(TAB, "fig6D_MS_row.csv"), index=False)
DT.to_csv(os.path.join(TAB, "fig6D_delta1_74traits_per_tech.csv"), index=False)
MG.to_csv(os.path.join(TAB, "fig6D_matched_genes.csv"), index=False)
json.dump(extra, open(os.path.join(TAB, "fig6D_inputs.json"), "w"), ensure_ascii=False, indent=1)

# ---------------- V2 bridge: recompute the dsciCAT ladder with C's original criteria ----------------
try:
    Sr = load_scores("dsciCAT_seq_rna"); Sa = load_scores("dsciCAT_seq_atac")
    Oc = ad.read_h5ad(H5C, backed="r").obs
    both = [c for c in Sr.index if c in Oc.index] + [c for c in Sa.index if c in Oc.index and c not in Sr.index]
    both = sorted(set(Sr.index) & set(Sa.index) & set(Oc.index))
    traits = sorted(set(Sr.columns) & set(Sa.columns))
    labc = Oc.loc[both, "cell_type"].astype(str).values
    m6, cop = labc == "MOL6", labc == "COP/NFOL"
    Xr, Xa = Sr.loc[both, traits].values, Sa.loc[both, traits].values
    d1r = np.median(Xr[m6], 0) - np.median(Xr[cop], 0)
    d1a = np.median(Xa[m6], 0) - np.median(Xa[cop], 0)
    rho, p = spearmanr(d1r, d1a)
    ref = float(pd.read_csv(T6C).set_index("item").loc["cross-channel rho (raw)", "value"])
    grc = [str(g) for g in ad.read_h5ad(H5C, backed="r").var_names]
    nr_c, na_c = n_hit_per_trait(grc).reindex(traits), n_hit_per_trait([l.strip() for l in open(os.path.join(DER, "gs_dsciCAT_seq_genes.txt")) if l.strip()]).reindex(traits)
    L = ladder(d1r, d1a, np.log(nr_c.n_hit.values), np.log(na_c.n_hit.values), None)
    V2 = pd.DataFrame([dict(item="C delivered ladder ρ (3k_dsci_oligo_subset.h5ad, 5,962 nuclei, MOL6 vs COP/NFOL)", value=ref, n=5962),
                       dict(item="V2 dsciCAT recomputed on D's RNA matrix (raw)", value=float(rho), n=len(both)),
                       dict(item="V2 same (after domain-size adjustment)", value=float(L["rho_adj"]), n=len(both)),
                       dict(item="difference (V2 raw − C)", value=float(rho - ref), n=np.nan),
                       dict(item="criterion |difference| <= 0.05 ?", value=np.nan, n=np.nan)])
    V2["pass_0.05"] = ["", "", "", "YES" if abs(rho - ref) <= 0.05 else "NO", ""]
    V2.to_csv(os.path.join(TAB, "V2_bridge_ladder_dsciCAT.csv"), index=False)
    SAY("\n=== V2 bridge ===\n" + V2.round(4).to_string(index=False))
except Exception as e:
    SAY("V2 cannot be computed yet (scores incomplete?): %s", e)
SAY("written -> %s", TAB)
