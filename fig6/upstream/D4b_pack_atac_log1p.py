#!/usr/bin/env python3
"""D4b_pack_atac_log1p.py — build the log1p chromatin matrices that panel D should have scored.

Why this exists: D4_build_h5ad.py packed the three chromatin arms straight out of the ArchR
GeneScoreMatrix, which is a 1e4-per-nucleus-normalised but **un-logged** quantity (max 1388.6 own /
297.1 10x / 372.7 HT). scDRS was invoked with --flag-raw-count False, so it applied expm1() to those
values: inf in the own arm (31/31 run logs warn), 1e162-scale numbers in the HT arm (25/25 logs warn
in the multiply that follows), silently ~1e129 in the 10x arm. The expression arms were always correct
(log1p CPM, max <= 8.69), and so was panel C's chromatin channel (Signac gene_score data slot =
log1p, max 6.47). See A6_fig6_input_transform_audit.py for the whole-table audit.

This script does not touch the delivered files. It writes a parallel tree
`D_crosschemistry/h5ad_log1p/` holding X.data -> log1p(X.data), everything else copied verbatim, and
asserts the round trip (expm1 of the new X reproduces the old X to float tolerance) plus the obs/var
identity. Run preflight_h5ad.py on the result before scoring.

    python3 D4b_pack_atac_log1p.py                       # 3 chromatin arms
    OUT_ROOT=... python3 D4b_pack_atac_log1p.py
"""
import os
import sys

import anndata as ad
import numpy as np
import scipy.sparse as sp

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

SRC_ROOT = os.environ.get("SRC_ROOT", WORK + "/figures/routeA_out/D_crosschemistry/h5ad")
OUT_ROOT = os.environ.get("OUT_ROOT", WORK + "/figures/routeA_out/D_crosschemistry/h5ad_log1p")
JOBS = (os.environ.get("JOBS") or
        "dsciCAT_seq_atac 10x_atac HT_scCAT_seq_HT_atac").split()
TOL = 1e-9


def pack(job):
    src, dst = f"{SRC_ROOT}/{job}.h5ad", f"{OUT_ROOT}/{job}.h5ad"
    a = ad.read_h5ad(src)
    X = a.X
    if not sp.issparse(X):
        raise SystemExit(f"{job}: expected a sparse X, got {type(X)}")
    before = float(X.data.max())
    X.data = np.log1p(X.data)
    a.write(dst, compression="gzip")

    b = ad.read_h5ad(dst)
    o = ad.read_h5ad(src)                                  # re-read: what scDRS was fed before
    rt = b.X.copy(); rt.data = np.expm1(rt.data)            # what scDRS rebuilds internally
    d = float(abs(rt - o.X).max())                          # must be ~1e-13 (float64 rounding)
    same_obs = list(map(str, a.obs_names)) == list(map(str, b.obs_names))
    same_var = list(map(str, a.var_names)) == list(map(str, b.var_names))
    same_cov = os.path.exists(f"{SRC_ROOT}/{job}.cov.tsv")
    if same_cov:
        import shutil
        shutil.copyfile(f"{SRC_ROOT}/{job}.cov.tsv", f"{OUT_ROOT}/{job}.cov.tsv")
    rowsum = np.asarray(b.X.sum(1)).ravel()
    print(f"{job:22s} {b.shape} max {before:9.3f} -> {b.X.data.max():6.3f} | "
          f"expm1 round-trip abs err {d:.2e} | obs {same_obs} var {same_var} | "
          f"log1p per-nucleus sum median {np.median(rowsum):.1f}")
    assert d < 1e-8 * max(1.0, before), f"{job}: round trip off by {d}"
    assert same_obs and same_var and same_cov, f"{job}: metadata mismatch"
    assert float(b.X.data.max()) < 30.0, f"{job}: still not a log-scale matrix"
    assert np.isfinite(b.X.data).all(), f"{job}: non-finite entries"
    return dst


if __name__ == "__main__":
    os.makedirs(OUT_ROOT, exist_ok=True)
    files = [pack(j) for j in JOBS]
    print("\nwrote:", ", ".join(os.path.basename(f) for f in files), "->", OUT_ROOT)
    rc = os.system(f'{sys.executable} {os.path.join(os.path.dirname(os.path.abspath(__file__)),"preflight_h5ad.py")} '
                   + " ".join(files))
    sys.exit(0 if rc == 0 else 2)
