# -*- coding: utf-8 -*-
"""A4_GSE301908_delta2_full74.py -- three axes x all 74 traits for GSE301908 (human snRNA-seq,
293,970 nuclei released).

The previous version, `A4_GSE301908_delta2_mirror.py`, scored only 12 key traits (its own header said
"~2 h, 12 traits only").  This script runs all 74 and computes **all three axes in one scoring pass**,
so d1 / d2A / d2B share one internal scale.

Object construction reuses the memory-cheap reader validated in the previous version plus the stratified
lineage cap (the full 293,970 nuclei make scDRS' gene-level loess hit a near-singularity - two different
errors were observed before this was settled):
  float32 CSR cache -> genes detected in >= 20 nuclei -> nucleus QC (nFeature >= 200, pct.mt < 25)
  -> <= 12,000 nuclei per lineage
Axis definitions (same form as A2; majorCluster labels):
  d1  = median(Oligo) - median(OPC)
  d2A = median(Astro) - median(Micro + Lym)
  d2B = median(Neuron) - median(all non-neuronal)
Statistical unit = donor (the msNN prefix of orig.ident).  n_ctrl=200, random_seed=0, i.e. verbatim the
parameters of the delivered A4 d1 sweep, so Panel C's d1 ranks and Panel B's cell come from one basis.
(The d2 axes were scored with n_ctrl=100 in the earlier mirror script; in THIS script all three axes
come out of the single n_ctrl=200 pass above - the parameter that made the delivered numbers is the one
in the code.)

Basis decision 2026-11-04 - input transform: this script reads `GSE301908_full` (log-normalised,
float64), NOT `GSE301908_counts` (raw counts, float32).  Reason: every cohort that feeds Panels A/B/C/D
is given to scDRS as a log-normalised matrix (the GSE277430 source file is literally named
`GSE277430_normCount_plain.rds`; the share of integer-valued non-zero entries measures 0 %), and the
A4 #5 in the Panel B cell also comes from the log basis - switching to counts here would make A4 the
only cohort of the five in a different transform family.  The counts version survives only as a
delivered sensitivity analysis (S7g: rho(log, counts) = 0.831, MS is #5 on both, SCZ #68 vs #71).

Bug fixed against the earlier mirror: that script read this float64 file as **float32** (the values come
out as garbage of order -3.4e38).  Its output is no longer in tables/, so **no delivered number was
contaminated**.  This script reads float64 and keeps A4d's QC + <= 12,000-per-lineage stratified sampling.

Usage:  python A4_GSE301908_delta2_full74.py --shard 0 --nshards 3
        python A4_GSE301908_delta2_full74.py --merge --nshards 3
Output: tables/a4_delta2_shard{N}.csv  ->  tables/trait_sweep_delta2_GSE301908_74.csv
Runtime: scDRS scoring dominates; run under the scoring kernel
$DSCI_PY (scdrs 1.0.2).
"""
import argparse
import os
import time

import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.stats import wilcoxon

# Root of the author's analysis tree (FY/, figures/, routeA_probe/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
OUT = WORK + "/figures/routeA_out/A4_GSE301908_v1"
TAB = os.path.join(OUT, "tables")
IN = WORK + "/routeA_probe/GSE301908_full"
GS = WORK + "/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs"
N_CTRL, SEED = 200, 0   # verbatim the delivered A4 d1 sweep (n_ctrl=200, seed=0): C's d1 ranks stay same-source as B's cell
CAP = 12000
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
    d = pd.concat([pd.read_csv(os.path.join(TAB, "a4_delta2_shard%d.csv" % i))
                   for i in range(A.nshards)], ignore_index=True)
    d = d[d.trait != "TRAIT"].sort_values("d1", ascending=False).reset_index(drop=True)
    d["r1"] = d.d1.rank(ascending=False, method="min").astype(int)
    d["rA"] = d.dA.rank(ascending=False, method="min").astype(int)
    d["rB"] = d.dB.rank(ascending=False, method="min").astype(int)
    d.to_csv(os.path.join(TAB, "trait_sweep_delta2_GSE301908_74.csv"), index=False)
    log("merged %d rows -> tables/trait_sweep_delta2_GSE301908_74.csv" % len(d))
    old = pd.read_csv(os.path.join(TAB, "trait_sweep_74_GSE301908.csv"))
    for pat, tag in [("Multiple_sclerosis", "MS"), ("Schizophrenia", "SCZ")]:
        r = d[d.trait.str.contains(pat)].iloc[0]
        o = old[old.trait.str.contains(pat)].iloc[0]
        log("%s new: d1 %+.3f #%d | d2A %+.3f #%d | d2B %+.3f #%d   || old OL-only d1 rank #%s"
            % (tag, r.d1, r.r1, r.dA, r.rA, r.dB, r.rB, o["rank"]))
    raise SystemExit(0)

import anndata as ad
import scdrs

mine = [t for i, t in enumerate(traits) if i % A.nshards == A.shard]
log("shard %d/%d owns %d traits" % (A.shard, A.nshards, len(mine)))

genes = pd.read_csv(os.path.join(IN, "full_genes.tsv"), header=None)[0].astype(str).tolist()
barcodes = pd.read_csv(os.path.join(IN, "full_barcodes.tsv"), header=None)[0].astype(str).tolist()
meta = pd.read_csv(os.path.join(IN, "full_meta.tsv"), sep="\t")
n_cells, n_genes = len(barcodes), len(genes)
with open(os.path.join(IN, "full_csr.bin"), "rb") as fh:
    indptr = np.fromfile(fh, dtype=np.int32, count=n_cells + 1)
    nnz = int(indptr[-1])
    indices = np.fromfile(fh, dtype=np.int32, count=nnz)
    fh.seek((n_cells + 1) * 4 + nnz * 4)
    data = np.fromfile(fh, dtype=np.float64, count=nnz)      # the data section of full_csr.bin is float64 (log-normalised)
M = csc_matrix((data, indices, indptr), shape=(n_genes, n_cells)).T.tocsr()
log("loaded %d nuclei x %d genes (nnz=%d) in %.0fs" % (n_cells, n_genes, nnz, time.time() - t0))
del data, indices

ct = meta["majorCluster"].astype(str)
raw = meta["orig.ident"].astype(str)
meta["patient"] = raw.str.extract(r"^(ms\d+)")[0].fillna(raw)
log("%d donors" % meta["patient"].nunique())
adata = ad.AnnData(X=M, obs=meta.copy(), var=pd.DataFrame(index=genes))
adata.obs_names = barcodes
det = np.asarray((M > 0).sum(0)).ravel()
adata = adata[:, det >= 20].copy()
nf = pd.to_numeric(adata.obs["nFeature_RNA"], errors="coerce")
mt = pd.to_numeric(adata.obs["percent.mt"], errors="coerce")
adata = adata[(nf >= 200).values & (mt < 25).values].copy()
log("after gene/nucleus QC %s" % (adata.shape,))
ctc = adata.obs["majorCluster"].astype(str)
rng = np.random.default_rng(1)
sel = []
for c in ctc.unique():
    v = np.where((ctc == c).values)[0]
    sel.append(v if len(v) <= CAP else np.sort(rng.choice(v, CAP, replace=False)))
adata = adata[np.sort(np.concatenate(sel))].copy()
log("after stratified sampling (<=%d per lineage) %s" % (CAP, adata.shape))
scdrs.pp.preprocess(adata)
ct = adata.obs["majorCluster"].astype(str)
log("preprocess ok %.0fs | lineages %s" % (time.time() - t0, ct.value_counts().to_dict()))

G = {
    "MAT": (ct == "Oligo").values,
    "PRE": (ct == "OPC").values,
    "AST": (ct == "Astro").values,
    "IMM": ct.isin(["Micro", "Lym"]).values,
    "NEU": (ct == "Neuron").values,
    "REST": (ct != "Neuron").values,
}
pat = adata.obs["patient"].astype(str).values
donors = sorted(set(pat))


def axis_vals(s, a, b):
    out = []
    for dn in donors:
        m = pat == dn
        x, y = s[m & a], s[m & b]
        if len(x) >= MIN_CELL and len(y) >= MIN_CELL:
            out.append(float(np.median(x) - np.median(y)))
    return np.array(out)


rows = []
for i, tr in enumerate(mine):
    ts = time.time()
    use = [g for g in sets[tr] if g in adata.var_names]
    res = scdrs.method.score_cell(adata, use, n_ctrl=N_CTRL, random_seed=SEED)
    col = "norm_score" if "norm_score" in res.columns else res.columns[0]
    s = res[col].to_numpy()

    def pack(v):
        if len(v) < 3:
            return (np.nan, np.nan, len(v))
        m = float(np.median(v))
        try:
            p = float(wilcoxon(v).pvalue)
        except ValueError:
            p = np.nan
        return (m, p, len(v))

    v1, vA, vB = axis_vals(s, G["MAT"], G["PRE"]), axis_vals(s, G["AST"], G["IMM"]), axis_vals(s, G["NEU"], G["REST"])
    d1, p1, n1 = pack(v1); dA, pA, nA = pack(vA); dB, pB, nB = pack(vB)
    rows.append({"trait": tr, "n_genes_hit": len(use), "d1": d1, "P1": p1, "n1": n1,
                 "dA": dA, "PA": pA, "nA": nA, "dB": dB, "PB": pB, "nB": nB})
    log("[%d/%d] %s d1=%+.3f(n%d) dA=%+.3f(n%d) dB=%+.3f(n%d) (%.0fs)"
        % (i + 1, len(mine), tr, d1, n1, dA, nA, dB, nB, time.time() - ts))
    pd.DataFrame(rows).to_csv(os.path.join(TAB, "a4_delta2_shard%d.csv" % A.shard), index=False)
log("shard %d done" % A.shard)
