# -*- coding: utf-8 -*-
"""A4d-step1 (v2, corrected against what was measured): assemble the RAW COUNTS of GSE301908 - the 17
*_raw_feature_bc_matrix.h5 files the authors deposited separately - into one count CSR covering
exactly the same nuclei with exactly the same labels as the log-normalised version already in use.

Measured facts (10-30 15:10, verified on this machine):
  - h5 = Cell Ranger raw_feature_bc_matrix: ~2M barcodes per library, 36,601 Gene Expression
    features, nnz 42-63M and 100% integers (min 1), so these really are raw counts.
  - barcodes alone cannot identify a library (any two groups share 40-50% of them), but every
    orig.ident group in the rds is fully contained in its true library (ms1r1 12,373/12,373 all in
    MS1; ms1r2 17,504 all in MS3; ms4r1 in MS5; ms4r2 in MS6; ms6r0 in MS2), so requiring
    "overlap == group size" identifies the library uniquely and without ambiguity.
  - identifying the libraries also matches the GEO names (MS1..MS14 / Ctrl1..Ctrl3) to the rds donor
    labels (msN + region rX), so the count version carries both a library and a donor unit.

Output $DSCI_WORK_ROOT/routeA_probe/GSE301908_counts/: full_csr.bin (indptr/indices int32, data float32)
 + full_barcodes.tsv / full_genes.tsv / full_meta.tsv / library_group_map.csv
"""
import os
import sys
import time

import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csc_matrix, vstack

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

SRC = WORK + "/routeA_probe/GSE301908_raw"
IND = WORK + "/routeA_probe/GSE301908_full"
OUTD = os.environ.get("A4D_OUT", WORK + "/routeA_probe/GSE301908_counts")
MANIFEST = "/tmp/a4_manifest.tsv"
os.makedirs(OUTD, exist_ok=True)
t0 = time.time()


def log(*a):
    print(*a, flush=True)


man = pd.read_csv(MANIFEST, sep="\t", header=None, names=["f", "size"])
meta = pd.read_csv(os.path.join(IND, "full_meta.tsv"), sep="\t")
meta = meta.rename(columns={meta.columns[0]: "cell"})
meta["grp"] = meta["orig.ident"].astype(str)
meta["bc"] = [c.split("_", 1)[1] if "_" in c else c for c in meta["cell"].astype(str)]
# local mapping inside a group (barcodes repeat across libraries, so never build a global table)
grows = {g: {b: i for b, i in zip(d["bc"], d.index)} for g, d in meta.groupby("grp")}
gsize = {g: len(v) for g, v in grows.items()}
# whether the rds groups are contiguous (decides if a reorder is needed)
starts = {g: min(v.values()) for g, v in grows.items()}
contig = all(max(v.values()) - min(v.values()) + 1 == len(v) for v in grows.values())
log("rds %d nuclei / %d groups; the groups are%s contiguous in meta" % (len(meta), len(grows), "" if contig else " NOT"))

blocks, order_rows, gmap = [], [], []
h5_perfect = {}
for f, sz in man.values:
    p = os.path.join(SRC, f)
    if not (os.path.exists(p) and os.path.getsize(p) == int(sz)):
        log("  skipping the incomplete download %s" % f[:34])
        continue
    lib = f.split("_")[1]                                       # MS1 / Ctrl3 ...
    with h5py.File(p, "r") as h:
        g0 = h["matrix"] if "matrix" in h else h[list(h.keys())[0]]
        bc = [b.decode() if isinstance(b, bytes) else b for b in g0["barcodes"][:]]
        fg = g0["features"]
        if isinstance(fg, h5py.Group):
            nm = [x.decode() for x in fg["name"][:]]
            tp = [x.decode() for x in fg["feature_type"][:]]
        else:
            nm = [x[1].decode() for x in fg[:]]
            tp = [x[2].decode() for x in fg[:]]
        kg = np.array([t == "Gene Expression" for t in tp])
        genes = np.array(nm)[kg].tolist()
        indptr = g0["indptr"][:].astype(np.int64)
        indices = g0["indices"][:].astype(np.int64)
        data = g0["data"][:].astype(np.float32)
    if not (np.allclose(data, np.floor(data)) and data.min() >= 0):
        log("  !! %s has non-integer counts -> abort" % lib)
        raise SystemExit(3)
    if "genes_ref" in dir():
        assert genes == genes_ref, "%s gene list differs from the previous one" % lib
    else:
        genes_ref = genes
    bset = set(bc)
    hit = {g: len(set(grows[g]) & bset) for g in gsize}
    perfect = [g for g, n in hit.items() if n == gsize[g]]
    if len(perfect) == 0:
        log("  !! %s does not fully contain any group (best hits %s) -> abort" % (lib, sorted(hit.items(), key=lambda x: -x[1])[:3]))
        raise SystemExit(4)
    h5_perfect[lib] = perfect
    log("  %-6s 100%% contains %s" % (lib, perfect))
    del bset, hit, perfect

# every rds group must be fully contained in at least one library; when several libraries contain it
# (like the pooled MS9, which also holds ms15r1), prefer the dedicated library with the fewest groups;
group2lib = {}
for lib, pl in sorted(h5_perfect.items(), key=lambda kv: (len(kv[1]), kv[0])):
    for g in pl:
        if g not in group2lib:
            group2lib[g] = lib
        elif len(h5_perfect[group2lib[g]]) > len(pl):
            group2lib[g] = lib
for g, lib in sorted(group2lib.items()):
    extras = [l for l, pl in h5_perfect.items() if l != lib and g in pl]
    if extras:
        log("  WARNING group %s shared by %s (taking the dedicated library %s, counts from %s)" % (g, extras, lib, lib))
if len(group2lib) < len(gsize):
    missing = set(gsize) - set(group2lib)
    log("  !! %d/%d groups unassigned: %s -> abort" % (len(missing), len(gsize), sorted(missing)[:8]))
    raise SystemExit(5)

# take the columns of each group out of every library
blocks, order_rows, gmap = [], [], []
for f, sz in man.values:
    lib = f.split("_")[1]                                       # MS1 / Ctrl3 ...
    if lib not in h5_perfect:
        continue
    with h5py.File(os.path.join(SRC, f), "r") as h:
        g0 = h["matrix"] if "matrix" in h else h[list(h.keys())[0]]
        fg = g0["features"]
        if isinstance(fg, h5py.Group):
            nm = [x.decode() for x in fg["name"][:]]
            tp = [x.decode() for x in fg["feature_type"][:]]
        else:
            nm = [x[1].decode() for x in fg[:]]
            tp = [x[2].decode() for x in fg[:]]
        kg = np.array([t == "Gene Expression" for t in tp])
        genes = np.array(nm)[kg].tolist()
        indptr = g0["indptr"][:].astype(np.int64)
        indices = g0["indices"][:].astype(np.int64)
        data = g0["data"][:].astype(np.float32)
        bcs = [b.decode() if isinstance(b, bytes) else b for b in g0["barcodes"][:]]
    if not (np.allclose(data, np.floor(data)) and data.min() >= 0):
        log("  !! %s has non-integer counts -> abort" % lib)
        raise SystemExit(3)
    if "genes_ref" in dir():
        assert genes == genes_ref, "%s gene list differs from the previous one" % lib
    else:
        genes_ref = genes
    for g in h5_perfect[lib]:
        loc = grows[g]
        pos = [i for i, b in enumerate(bcs) if b in loc]
        sub = csc_matrix((data, indices, indptr), shape=(len(nm), len(bcs)))[np.where(kg)[0], :][:, pos]
        M = sub.T.tocsr()                       # cells x genes (fix: slicing the rows of a genes-x-cells matrix in cell order had taken the first n genes)
        rows = [loc[bcs[i]] for i in pos]
        order = np.argsort(rows)
        assert M.shape[0] == len(rows) and M.shape[1] == len(genes), "%s shape %s != (%d,%d)" % (lib, M.shape, len(rows), len(genes))
        blocks.append(M[order])
        order_rows.extend([rows[i] for i in order])
        donor = pd.Series([g]).str.extract(r"^(ms\d+|ctrl\d*)")[0].fillna(g).iat[0]
        gmap.append({"library_geo": lib, "rds_group": g, "donor": donor,
                     "n_cells_rds": len(rows), "h5_barcodes": len(bcs), "h5_nnz": int(M.nnz)})
        log("  %-6s -> %-7s %d nuclei (%.0f%% covered) nnz=%.1fM integers OK (%.0fs)"
            % (lib, g, len(rows), 100 * len(rows) / gsize[g], M.nnz / 1e6, time.time() - t0))
    del data, indices, sub, M

pd.DataFrame(gmap).to_csv(os.path.join(OUTD, "library_group_map.csv"), index=False)
if len(blocks) < len(gsize):
    log("only %d/%d groups available -> not assembling yet (re-run once the download finishes)" % (len(blocks), len(gsize)))
    raise SystemExit(2)

A = vstack(blocks, format="csr")
order_rows = np.array(order_rows)
if not contig or not np.array_equal(order_rows, np.arange(len(meta))):
    perm = np.empty(len(meta), dtype=np.int64)
    perm[order_rows] = np.arange(len(order_rows))
    A = A[perm].tocsr()
    log("reordered to the rds nucleus order")
log("final %d nuclei x %d genes, nnz=%.1fM (%.1f counts per nucleus), min=%d max=%d, integers %.0f%%"
    % (A.shape[0], A.shape[1], A.nnz, A.nnz / A.shape[0], A.data.min(), A.data.max(),
       100 * float(np.allclose(A.data, np.floor(A.data)))))
assert A.shape[0] == len(meta)

with open(os.path.join(OUTD, "full_csr.bin"), "wb") as fh:
    A.indptr.astype(np.int32).tofile(fh)
    A.indices.astype(np.int32).tofile(fh)
    A.data.astype(np.float32).tofile(fh)
pd.Series(meta["cell"].astype(str)).to_csv(os.path.join(OUTD, "full_barcodes.tsv"), header=False, index=False)
pd.Series(genes_ref).to_csv(os.path.join(OUTD, "full_genes.tsv"), header=False, index=False)
libmap = {r["library_geo"]: r["rds_group"] for r in gmap}
meta["library_geo"] = meta["grp"].map({v: k for k, v in libmap.items()})
meta.drop(columns=["bc"]).to_csv(os.path.join(OUTD, "full_meta.tsv"), sep="\t", index=False)
log("done in %.0fs -> %s" % (time.time() - t0, OUTD))
