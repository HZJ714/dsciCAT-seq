#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s15b_allgene_pop_counts.py -- per-gene x per-population DETECTED / TOTAL nuclei counts for
   every gene in every released matrix (not just the 65 curated rows).

Why: 7Da of FIG6C_FIG7D_NARRATIVE_DESIGN_20261116 (their P0-b) asks whether the *defining gene set of
each subtype* is enriched for gold-standard drug targets. That needs a population-specific gene set
built from the whole tested universe, and the universe itself for the Fisher denominator.
The recipe is the one already committed (S8_bubble/stats_v9.py, reproduced by s10a): detection =
stored value > 0 for that nucleus; populations are the harmonised labels used everywhere else.
Read-only on all inputs.
"""
import os
import numpy as np
import pandas as pd
import gzip

# Root of the author's analysis tree (figures/, routeA_out/, routeA_probe/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
OUT = WORK + "/figures/SUPPL_ANALYSES_20261112/tables"
MAP = {
 "GSE118257": {**{f"Oligo{i}": "OL" for i in range(1, 7)}, "ImOlGs": "iOL", "COPs": "iOL", "OPCs": "OPC",
               "Astro": "AST", "Micro/Mac": "MG", "Immune": "LYM", "Endo": "VASC", "Pericyte/SMC": "VASC", "Neuron": "NEU"},
 "GSE301908": {"Oligo": "OL", "OPC": "OPC", "Astro": "AST", "Micro": "MG", "Lym": "LYM", "Endo": "VASC", "Neuron": "NEU"},
 "GSE231585": {"OLIGO": "OL", "OPC": "OPC", "AST": "AST", "MG": "MG", "LYM": "LYM", "END": "VASC", "NEU": "NEU"},
 "GSE277430": {"OLIGO": "OL", "OPC": "OPC"}}
POPS = ["OL", "iOL", "OPC", "AST", "MG", "LYM", "VASC", "NEU"]


def save(ds, genes, det, tot):
    df = pd.DataFrame(det, index=pd.Index(genes, name="gene"), columns=POPS[:det.shape[1]])
    df.to_csv(f"{OUT}/S15_det_{ds}.csv")
    pd.DataFrame(tot, index=POPS[:tot.shape[0]], columns=["n"]).to_csv(f"{OUT}/S15_n_{ds}.csv")
    print("   saved", ds, df.shape, "n=", tot)


# ---------------------------------------------------------------- A3 : text gz, all genes
def A3():
    P = WORK + "/routeA_probe/GSE118257_matrix.txt.gz"
    ANN = WORK + "/routeA_probe/GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz"
    POOL = {"Astro": ["Astrocytes", "Astrocytes2"], "Micro/Mac": ["Microglia_Macrophages", "Macrophages"],
            "Immune": ["Immune_cells"], "Endo": ["Endothelial_cells1", "Endothelial_cells2"],
            "Pericyte/SMC": ["Pericytes", "Vasc_smooth_muscle"], "Neuron": [f"Neuron{i}" for i in range(1, 6)]}
    OL = ["OPCs", "COPs", "ImOlGs", "Oligo1", "Oligo2", "Oligo3", "Oligo4", "Oligo5", "Oligo6"]
    with gzip.open(P, "rt") as f:
        cells = f.readline().rstrip("\n").split("\t")
    lab = pd.read_csv(ANN, sep="\t", index_col=0).loc[cells, "Celltypes"].values
    msk = {r: (lab == r) for r in OL}
    msk.update({r: np.isin(lab, pl) for r, pl in POOL.items()})
    gidx = np.full(len(cells), -1, np.int64)
    for k, r in enumerate(MAP["GSE118257"]):
        if r in msk:
            gidx[msk[r]] = POPS.index(MAP["GSE118257"][r])
    np_ = np.bincount(gidx[gidx >= 0], minlength=len(POPS)).astype(np.int64)
    with gzip.open(P, "rt") as f:
        f.readline()
        genes, rows = [], []
        for i, line in enumerate(f):
            p = line.rstrip("\n").split("\t", 1)
            v = np.frombuffer(bytearray(p[1].replace("\t", " ").encode("ascii", "ignore")), dtype="i1")
            # fast positive mask: any byte other than '0'/' ' means a non-zero token
            toks = p[1].split("\t")
            det = np.fromiter((t != "0" and t != "" for t in toks), dtype=bool, count=len(toks))
            genes.append(p[0])
            rows.append(np.bincount(gidx[det], minlength=len(POPS)).astype(np.int64))
            if i % 5000 == 0:
                print("   A3", i, flush=True)
    save("GSE118257", genes, np.array(rows), np_)


# ---------------------------------------------------------------- A4 / A2 : CSR dump
def csr(ds, D, gene_file, bin_file, meta_file, labcol):
    gn = [l.rstrip("\n") for l in open(D + gene_file)]
    mt = pd.read_csv(D + meta_file, sep="\t", dtype=str)
    cl = mt[labcol].map(MAP[ds]).values
    present = [p for p in POPS if p in set(cl[~pd.isna(cl)])]
    pid = {p: i for i, p in enumerate(present)}
    gidx = np.array([pid.get(x, -1) if not pd.isna(x) else -1 for x in cl], np.int64)
    nc = len(gidx)
    ip = np.fromfile(D + bin_file, dtype=np.int32, count=nc + 1)
    nz = int(ip[-1])
    dt = np.float32 if ds == "GSE301908" else np.float64
    ix = np.memmap(D + bin_file, dtype=np.int32, mode="r", offset=(nc + 1) * 4, shape=(nz,))
    da = np.memmap(D + bin_file, dtype=dt, mode="r", offset=(nc + 1) * 4 + nz * 4, shape=(nz,))
    np_ = np.bincount(gidx[gidx >= 0], minlength=len(present)).astype(np.int64)
    det = np.zeros((len(gn), len(present)), np.int64)
    for s in range(0, nz, 50_000_000):
        e = min(s + 50_000_000, nz)
        seg = np.asarray(ix[s:e]); val = np.asarray(da[s:e])
        pos = np.nonzero(val > 0)[0]
        cell = np.searchsorted(ip, pos + s, side="right") - 1
        gi = gidx[cell]
        keep = gi >= 0
        gene = seg[pos[keep]]
        flat = gene * len(present) + gi[keep]
        det += np.bincount(flat, minlength=det.size).reshape(det.shape)
        print("   ", ds, "%d/%d nonzero" % (e, nz), flush=True)
    df = pd.DataFrame(det, index=pd.Index(gn, name="gene"), columns=present)
    df.to_csv(f"{OUT}/S15_det_{ds}.csv")
    pd.DataFrame(np_, index=present, columns=["n"]).to_csv(f"{OUT}/S15_n_{ds}.csv")
    print("   saved", ds, df.shape)


# ---------------------------------------------------------------- A7 : h5ad CSC
def A7():
    import anndata as ad
    import scipy.sparse as sp
    R = WORK + "/figures/routeA_out/A7_GSE231585_rna_v1/derived/"
    a = ad.read_h5ad(R + "raw_qc.h5ad")
    s = ad.read_h5ad(R + "scored.h5ad")
    lab = s.obs["celltype"].reindex(a.obs_names).map(MAP["GSE231585"]).values
    present = [p for p in POPS if p in set(lab[~pd.isna(lab)])]
    pid = {p: i for i, p in enumerate(present)}
    gidx = np.array([pid.get(x, -1) if not pd.isna(x) else -1 for x in lab], np.int64)
    X = (a.X if sp.issparse(a.X) else sp.csr_matrix(a.X)).tocsc()
    np_ = np.bincount(gidx[gidx >= 0], minlength=len(present)).astype(np.int64)
    det = np.zeros((X.shape[1], len(present)), np.int64)
    for j in range(X.shape[1]):
        col = X[:, j]
        r = col.indices[col.data > 0]
        gi = gidx[r]
        k = gi >= 0
        if k.any():
            det[j] += np.bincount(gi[k], minlength=len(present))
        if j % 5000 == 0:
            print("   A7", j, flush=True)
    df = pd.DataFrame(det, index=pd.Index(np.asarray(a.var_names), name="gene"), columns=present)
    df.to_csv(f"{OUT}/S15_det_GSE231585.csv")
    pd.DataFrame(np_, index=present, columns=["n"]).to_csv(f"{OUT}/S15_n_GSE231585.csv")
    print("   saved A7", df.shape)


if __name__ == "__main__":
    import sys
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "A3"):
        A3()
    if which in ("all", "A4"):
        csr("GSE301908", WORK + "/routeA_probe/GSE301908_counts/", "full_genes.tsv", "full_csr.bin",
            "full_meta.tsv", "majorCluster")
    if which in ("all", "A2"):
        csr("GSE277430", WORK + "/routeA_probe/GSE277430/OL_subset/", "OL_genes.tsv", "OL_csr.bin",
            "OL_meta.tsv", "expertAnno.l1")
    if which in ("all", "A7"):
        A7()
