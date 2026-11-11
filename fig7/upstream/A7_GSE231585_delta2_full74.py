# -*- coding: utf-8 -*-
"""GSE231585 (human snRNA, FANS oligodendrocyte enrichment): three axes x all 74 traits.

The object and the protocol are exactly those of the delivered A7b d1 sweep:
  derived/scored.h5ad (62,911 nuclei x 36,591 genes; celltype from the author panel annotation)
  scdrs.method.score_cell(n_ctrl=200, random_seed=20261030)
  unit = sample (donor x lesion domain, 15 of them); a group needs >= 10 nuclei to count in a sample
so d1 reconciles directly with trait_sweep_74_GSE231585.csv and all three axes share one scoring
run, which keeps their internal scale identical.

The three axes are defined exactly as in A2:
  d1  = median(OLIGO) - median(OPC)
  d2A = median(AST)   - median(MG + LYM)
  d2B = median(NEU)   - median(all non-neuronal)   # only 717 annotated neurons, thin; see the n columns

Usage:  python A7_GSE231585_delta2_full74.py --shard 0 --nshards 2
        python A7_GSE231585_delta2_full74.py --merge
Output: tables/a7_delta2_shard{N}.csv  ->  tables/trait_sweep_delta2_GSE231585_74.csv
"""
import argparse
import os
import time

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

OUTD = WORK + "/figures/routeA_out/A7_GSE231585_rna_v1"
TAB = os.path.join(OUTD, "tables")
GS = WORK + "/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs"
N_CTRL, SEED = 200, 20261030
MIN_CELL = 10
t0 = time.time()


def log(*a):
    print("[%5ds]" % (time.time() - t0), *a, flush=True)


ap = argparse.ArgumentParser()
ap.add_argument("--shard", type=int, default=0)
ap.add_argument("--nshards", type=int, default=1)
ap.add_argument("--merge", action="store_true")
A = ap.parse_args()

sets = {}
for l in open(GS):
    f = l.rstrip("\n").split("\t")
    if len(f) < 2 or f[0] == "TRAIT":
        continue
    sets[f[0]] = f[1].split(",")
traits = sorted(sets)

if A.merge:
    d = pd.concat([pd.read_csv(os.path.join(TAB, "a7_delta2_shard%d.csv" % i))
                   for i in range(A.nshards)], ignore_index=True)
    d = d[d.trait != "TRAIT"].sort_values("median_diff_d1", ascending=False).reset_index(drop=True)
    d["r1"] = d.d1.rank(ascending=False, method="min").astype(int)
    d["rA"] = d.dA.rank(ascending=False, method="min").astype(int)
    d["rB"] = d.dB.rank(ascending=False, method="min").astype(int)
    d.to_csv(os.path.join(TAB, "trait_sweep_delta2_GSE231585_74.csv"), index=False)
    log("merged %d rows -> tables/trait_sweep_delta2_GSE231585_74.csv" % len(d))
    for pat, tag in [("Multiple_sclerosis", "MS"), ("Schizophrenia", "SCZ")]:
        r = d[d.trait.str.contains(pat)]
        if len(r):
            r = r.iloc[0]
            log("%s: Δ1 %+.3f #%d (n=%d) | Δ2A %+.3f #%d (n=%d) | Δ2B %+.3f #%d (n=%d)"
                % (tag, r.d1, r.r1, r.n1, r.dA, r.rA, r.nA, r.dB, r.rB, r.nB))
    raise SystemExit(0)

import anndata as ad
import scdrs

mine = [t for i, t in enumerate(traits) if i % A.nshards == A.shard]
log("shard %d/%d handles %d traits" % (A.shard, A.nshards, len(mine)))

adata = ad.read_h5ad(os.path.join(OUTD, "derived", "scored.h5ad"))
log("loaded %s | %s" % (adata.shape, adata.obs.celltype.value_counts().to_dict()))
scdrs.pp.preprocess(adata)
log("preprocess done")

genes = set(adata.var_names)
ct = adata.obs["celltype"].astype(str)
samp = adata.obs["sample"].astype(str).values
G = {
    "MAT": (ct == "OLIGO").values,
    "PRE": (ct == "OPC").values,
    "AST": (ct == "AST").values,
    "IMM": ct.isin(["MG", "LYM"]).values,
    "NEU": (ct == "NEU").values,
    "REST": (ct != "NEU").values,
}
for k, v in G.items():
    log("  %-4s %6d nuclei" % (k, v.sum()))
samples = sorted(set(samp))


def axis_vals(s, a, b):
    """Within one sample: median(group a) - median(group b); each group needs >= MIN_CELL nuclei."""
    out = []
    for sm in samples:
        m = samp == sm
        x, y = s[m & a], s[m & b]
        if len(x) >= MIN_CELL and len(y) >= MIN_CELL:
            out.append(float(np.median(x) - np.median(y)))
    return np.array(out)


rows = []
for i, tr in enumerate(mine):
    ts = time.time()
    use = [g for g in sets[tr] if g in genes]
    res = scdrs.method.score_cell(adata, use, n_ctrl=N_CTRL, random_seed=SEED)
    col = "norm_score" if "norm_score" in res.columns else res.columns[0]
    s = res[col].to_numpy()
    v1, vA, vB = axis_vals(s, G["MAT"], G["PRE"]), axis_vals(s, G["AST"], G["IMM"]), axis_vals(s, G["NEU"], G["REST"])

    def pack(v):
        if len(v) < 3:
            return (np.nan, np.nan, len(v))
        m = float(np.median(v))
        try:
            p = float(wilcoxon(v).pvalue)
        except ValueError:
            p = np.nan
        return (m, p, len(v))

    d1, p1, n1 = pack(v1); dA, pA, nA = pack(vA); dB, pB, nB = pack(vB)
    rows.append({"trait": tr, "n_genes_hit": len(use), "d1": d1, "P1": p1, "n1": n1,
                 "dA": dA, "PA": pA, "nA": nA, "dB": dB, "PB": pB, "nB": nB,
                 "median_diff_d1": d1})
    log("[%d/%d] %s d1=%+.3f(n%d) dA=%+.3f(n%d) dB=%+.3f(n%d) (%.0fs)"
        % (i + 1, len(mine), tr, d1, n1, dA, nA, dB, nB, time.time() - ts))
    pd.DataFrame(rows).to_csv(os.path.join(TAB, "a7_delta2_shard%d.csv" % A.shard), index=False)
log("shard %d done" % A.shard)
