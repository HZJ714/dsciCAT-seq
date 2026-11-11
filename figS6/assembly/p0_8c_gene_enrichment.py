#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
P0-8c  Fig6C pre-computation: one-sided Monte Carlo test of MS-enriched genes (scDRS gene-level fit)

Spec: FIG7_DRUG_HEATMAP_RECEIPT_20261112.md §8.3 P0-8c v3
  mechanism (the derivation of PMID 38092915 Fig 5d's six genes, self-contained):
    per-gene Pearson(gene expression, norm_score) within the enriched subtype's cells;
    null = the same gene's correlation with the same 1000 ctrl_norm_score columns;
    one-sided MC P = (#{ctrl_corr >= obs} + 1) / 1001        <- add-one, mandatory
    z = (obs - mean(ctrl)) / sd(ctrl) -> the theoretical one-sided P; BH FDR across genes.
   Discipline:
    * gene-level enrichment is a DECOMPOSITION/FIT of the cell-level signal, not an independent replication -> "driver" is banned.
    * the enrichment list is computed on the RNA side only (scDRS is a transcriptomic tool; ATAC has no norm_score);
      the ATAC heatmap = the same list's gene activity, a "two-channel confirmation", not "two-omics independent enrichment".
    * the orthologue table uses only scDRS's built-in load_homolog_mapping (skill §15.6); hand-rolled upper() is banned.
  Two conventions (both run; used by S6C):
    convention 1 = the enriched-subtype pool (Fig6A group_q<0.05 subtypes)
    convention 2 = the all-six-subtype pool
"""
import os, sys, json
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad

WORK = os.environ.get('DSCI_WORK_ROOT', '')
OUT = WORK + "/figures/routeA_out/P0_8C"
os.makedirs(OUT + "/tables", exist_ok=True); os.makedirs(OUT + "/panels", exist_ok=True)

RNA_H5 = WORK + "/FY/oligo/3k_dsci_oligo_subset.h5ad"
ATAC_H5 = WORK + "/FY/oligo/fig6_rerun_20261019/dualomics_probe/genescore_oligo_5977.h5ad"
RROOT = WORK + "/FY/oligo/scdrs_traits_draft/full74"
TRAIT = "PASS_Multiple_sclerosis"
STATS = WORK + "/GitHub/fig6/fig6_out/fig6A_stats.csv"
LOGO = WORK + "/figures/routeA_out/GENELEVEL_v1/logo_delta1_contributions.csv"
SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
ANCHOR = ["Mbp", "Plp1", "Mog", "Cnp", "Mobp", "Mag", "Olig2", "Sox10"]   # mouse symbols
PCT_MIN = 10.0          # near-zero exclusion: mean positive rate across the six subtypes < 10% (applied before any top-N cut)
NCTRL_EXPECT = 1000

sys.path.insert(0, os.path.expanduser('"~/.local/lib/python3.9/site-packages'))
from scdrs.util import load_homolog_mapping

# ---------------------------------------------------------------- cell convention (same source as 6A)
def score_frame():
    import glob
    fs = sorted(glob.glob(RROOT + "/*/" + TRAIT + ".full_score.gz"))
    if not fs:
        raise RuntimeError("cannot find " + TRAIT + ".full_score.gz")
    return pd.read_csv(fs[0], sep="\t", index_col=0), fs[0]

S, SRC = score_frame()
ctrl_cols = [c for c in S.columns if c.startswith("ctrl_norm_score_")]
assert len(ctrl_cols) == NCTRL_EXPECT, "control score columns %d != %d" % (len(ctrl_cols), NCTRL_EXPECT)
print("[1] scDRS score table %s | nuclei %d | norm_score + %d ctrl_norm_score columns"
      % (os.path.basename(SRC), len(S), len(ctrl_cols)))

a = ad.read_h5ad(RNA_H5)
key = (pd.Series(a.obs["cell_type"].astype(str)).str.upper()
       .str.replace("/", "", regex=False).str.replace("-", "", regex=False)).values
cells = set(S.index)
mask = a.obs_names.isin(cells)
grp = {s: np.where(mask & (key == s.replace("/", "")))[0] for s in SUB}
cnt = "/".join(str(len(grp[s])) for s in SUB)
print("[2] cell convention %s (fig6A_stats.csv: 203/490/1016/688/991/2574)" % cnt)
assert cnt == "203/490/1016/688/991/2574", "cell convention differs from Fig6A's batch!"

# ---------------------------------------------------------------- enriched-subtype delimitation (convention 1)
st = pd.read_csv(STATS)
enrich_q = st[st.group_q_BH5 < 0.05].cell_type.tolist()
enrich_w = st[st.mwu_q_BH5 < 0.05].cell_type.tolist()
print("[3] enriched subtypes: group_q<0.05 -> %s | mwu_q<0.05 -> %s" % (enrich_q, enrich_w))
if not enrich_q:
    raise RuntimeError("convention 1 finds no enriched subtype")
POOL = {"enriched_pool": [s for s in SUB if s in enrich_q], "all6_pool": list(SUB)}

# ---------------------------------------------------------------- gene list (the MS 1000, built-in orthologue table)
H = load_homolog_mapping("hsapiens", "mmusculus")
L = pd.read_csv(LOGO, usecols=["gene"])
pairs, unmapped = [], []
for g in L.gene:
    mo = H.get(g, None)
    if mo is None:
        unmapped.append(g); continue
    pairs.append((g, mo))
pos = {str(v): i for i, v in enumerate(a.var_names)}
used = [(g, mo) for g, mo in pairs if mo in pos]
print("[4] MS gene set %d -> mapped by the built-in table %d (no mouse symbol %d) -> present in the RNA matrix %d"
      % (len(L), len(pairs), len(unmapped), len(used)))

def densify(rows, cols_idx):
    X = a.X[np.ix_(rows, cols_idx)]
    return np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)

# ---------------------------------------------------------------- per-subtype statistics (RNA + ATAC)
atac = ad.read_h5ad(ATAC_H5)
apos = {str(v): i for i, v in enumerate(atac.var_names)}
amask = atac.obs_names.isin(cells)
akey = (pd.Series(atac.obs["cell_type"].astype(str)).str.upper()
        .str.replace("/", "", regex=False).str.replace("-", "", regex=False)).values
agrp = {s: np.where(amask & (akey == s.replace("/", "")))[0] for s in SUB}
print("[5] ATAC gene-activity matrix %d genes; %d columns hit in the same nuclei" % (atac.shape[1], len(used) - sum(1 for _, mo in used if mo not in apos)))

# the anchor genes also need per-subtype statistics (outside the MS gene set; they enter the heatmap, not the MC test)
anch_in = [(mo, mo) for mo in ANCHOR if mo in pos]
used_all = used + [x for x in anch_in if x[1] not in {m for _, m in used}]
print("[5b] statistics list = enrichment candidates %d + anchor fill %d = %d" % (len(used), len(used_all) - len(used), len(used_all)))
rows = []
for s in SUB:
    r = grp[s]
    E = densify(r, [pos[mo] for _, mo in used_all])
    for k, (hu, mo) in enumerate(used_all):
        rec = dict(gene_human=hu, gene=mo, subtype=s,
                   rna_mean=float(E[:, k].mean()),
                   rna_pct=float((E[:, k] > 0).mean() * 100), n=int(len(r)))
        if mo in apos:
            j = apos[mo]
            col = atac.X[agrp[s], j]
            col = np.asarray(col.todense()).ravel() if sp.issparse(atac.X) else np.asarray(col).ravel()
            rec["atac_mean"] = float(col.mean()); rec["atac_pct"] = float((col > 0).mean() * 100)
        else:
            rec["atac_mean"] = np.nan; rec["atac_pct"] = np.nan
        rows.append(rec)
T = pd.DataFrame(rows)
T.to_csv(OUT + "/tables/p08c_subtype_stats_long.csv", index=False)
print("[6] per-subtype statistics written, %d rows" % len(T))

# ---------------------------------------------------------------- MC test (both conventions)
def mc_test(pool_subtypes):
    r = np.concatenate([grp[s] for s in pool_subtypes])
    E = densify(r, [pos[mo] for _, mo in used])
    y = S.loc[a.obs_names[r], "norm_score"].values
    C = S.loc[a.obs_names[r], ctrl_cols].values
    Ec = E - E.mean(0); Cc = C - C.mean(0); yc = y - y.mean()
    sdE = np.sqrt((Ec ** 2).sum(0)); sdC = np.sqrt((Cc ** 2).sum(0)); sdy = np.sqrt((yc ** 2).sum())
    obs = (Ec * yc[:, None]).sum(0) / (sdE * sdy + 1e-12)                    # (n_genes,)
    null = (Ec.T @ Cc) / (sdE[:, None] * sdC[None, :] + 1e-12)               # (n_genes, 1000)
    mu, sd = null.mean(1), null.std(1, ddof=1)
    mc_p = (np.sum(null >= obs[:, None], axis=1) + 1) / (null.shape[1] + 1)  # add-one
    z = (obs - mu) / (sd + 1e-12)
    theo = 1.0 - __import__("scipy.stats", fromlist=["norm"]).norm.cdf(z)
    return obs, mc_p, z, theo

def bh(p):
    o = np.argsort(p); n = len(p); q = np.empty(n)
    qs = p[o] * n / np.arange(1, n + 1)
    q[o] = np.minimum.accumulate(qs[::-1])[::-1]      # cumulative min taken from the largest P backwards (the original omitted the reversal)
    return np.clip(q, 0, 1)

results = {}
for name, pool in POOL.items():
    obs, mc_p, z, theo = mc_test(pool)
    G = pd.DataFrame(dict(gene_human=[g for g, _ in used], gene=[m for _, m in used],
                          r=obs, mc_p=mc_p, mc_z=z, theo_p=theo))
    mean_pct = T[T.gene.isin(G.gene)].groupby("gene").rna_pct.mean()
    G["mean_pct"] = G.gene.map(mean_pct)
    G["near_zero"] = G.mean_pct < PCT_MIN
    keep = G[~G.near_zero].copy()
    keep["fdr_bh"] = bh(keep.mc_p.values)
    keep = keep.sort_values("mc_p")
    results[name] = dict(kept=keep, dropped=int(G.near_zero.sum()), all=G)
    print("[%s] pool=%s | near-zero excluded %d | kept %d | FDR<0.05 %d | MC P<0.05 %d"
          % (name, "+".join(pool), G.near_zero.sum(), len(keep),
             int((keep.fdr_bh < 0.05).sum()), int((keep.mc_p < 0.05).sum())))

# ---------------------------------------------------------------- the anchor block (anti-circularity: C ∩ E = empty)
sig = results["enriched_pool"]["kept"][results["enriched_pool"]["kept"].fdr_bh < 0.05]
anchor_rows = []
for mo in ANCHOR:
    hu = [k for k, v in H.items() if v == mo]
    anchor_rows.append(dict(gene=mo, gene_human=(hu[0].upper() if hu else mo.upper()),
                            in_rna=mo in pos, in_ms_set=mo in set(L.gene)))
AR = pd.DataFrame(anchor_rows)
AR["in_enriched_list"] = AR.gene.isin(set(sig.gene))
print("[7] anchor block 8 rows | intersection with the enrichment list = %d (must be 0; anti-circularity)" % int(AR.in_enriched_list.sum()))
AR.to_csv(OUT + "/tables/p08c_anchor_block.csv", index=False)

# ---------------------------------------------------------------- delivery
for name in POOL:
    results[name]["kept"].to_csv(OUT + "/tables/p08c_genelevel_%s.csv" % name, index=False)
    results[name]["all"].to_csv(OUT + "/tables/p08c_genelevel_%s_allgenes.csv" % name, index=False)
inter = set(results["enriched_pool"]["kept"].head(20).gene) & set(results["all6_pool"]["kept"].head(20).gene)
print("[8] top-20 intersection of the two conventions %d/20 (used by S6C)" % len(inter))
K = results["enriched_pool"]["kept"]
print("[9] significance-rate check: MC P<0.05 %d/%d = %.1f%% | FDR<0.05 %d (%.1f%%) | r quantiles %s"
      % (int((K.mc_p < 0.05).sum()), len(K), 100 * (K.mc_p < 0.05).mean(),
         int((K.fdr_bh < 0.05).sum()), 100 * (K.fdr_bh < 0.05).mean(),
         np.round(np.percentile(K.r, [5, 25, 50, 75, 95]), 3)))
json.dump(dict(spec="P0-8c v3", source_score_file=os.path.basename(SRC), n_ctrl=NCTRL_EXPECT,
               add_one=True, pct_min=PCT_MIN, cells=cnt,
               enriched_pool=POOL["enriched_pool"], enriched_by_group_q=enrich_q,
               enriched_by_mwu_q=enrich_w, ms_set_size=int(len(L)),
               homolog_mapped=len(pairs), rna_hit=len(used),
               near_zero_dropped={k: results[k]["dropped"] for k in POOL},
               fdr_lt05={k: int((results[k]["kept"].fdr_bh < 0.05).sum()) for k in POOL},
               anchor_intersection=int(AR.in_enriched_list.sum()),
               top20_two_pool_overlap=len(inter),
               note=("gene-level enrichment is a decomposition/concordance of the cell-level "
                     "signal, not an independent replication; never say driver")),
          open(OUT + "/tables/p08c_meta.json", "w"), indent=1)
print("\nconvention-1 enriched genes (FDR<0.05):")
print(sig[["gene_human", "gene", "r", "mc_p", "mc_z", "fdr_bh", "mean_pct"]].head(20)
      .round(4).to_string(index=False))
