#!/usr/bin/env python3
"""preflight_h5ad.py — refuse to score a matrix whose scale convention doesn't match the flag.

scDRS `compute-score --flag-raw-count False` applies `expm1()` to the input matrix internally
(scdrs/pp.py:613), i.e. it promises the matrix is already on a log scale.  `--flag-raw-count True`
instead expects raw counts and normalises to 1e4 itself.  Either promise, broken silently, changes
the control-cell fit and therefore every score — and it does NOT raise: it prints a numpy
RuntimeWarning into the per-batch compute.log and carries on.

This gate is run by D4b/D5 before launching a batch, on the exact file that will be passed to
`--h5ad-file`.  Exit code 2 = do not run.

    python3 preflight_h5ad.py A.h5ad B.h5ad ...            # default: expect log scale
    python3 preflight_h5ad.py --raw-count A.h5ad           # expect raw counts
"""
import argparse
import sys
import numpy as np

LOG_MAX = 30.0          # log1p of anything realistic is < 30; a log-scale matrix cannot exceed it
RAW_MIN = 30.0          # raw counts routinely exceed 30 in at least one cell


def check(path, expect_raw):
    import anndata as ad
    import scipy.sparse as sp
    a = ad.read_h5ad(path)                 # full read: these files are <= 400 MB on disk
    X = a.X
    if sp.issparse(X):
        data = np.asarray(X.data, dtype=np.float64)
    else:
        data = np.asarray(X, dtype=np.float64).ravel()
        data = data[data != 0]
    del a, X
    mx = float(data.max())
    shape = "-"
    # what scDRS will actually see after its own transform
    seen = data if expect_raw else np.expm1(data)
    n_inf = int(np.isinf(seen).sum())
    n_neg = int((data < -1e-9).sum())
    ok = (mx >= RAW_MIN) if expect_raw else (mx <= LOG_MAX and n_neg == 0)
    verdict = "OK" if ok and n_inf == 0 else "REJECT"
    print(f"{verdict:7s} {path}\n"
          f"        nnz={len(data)}  max={mx:.3f}  expm1_max={float(seen.max()):.3g}  "
          f"n_inf={n_inf}  n_negative={n_neg}  "
          f"expected={'raw counts' if expect_raw else 'log scale'}", file=sys.stderr)
    return ok and n_inf == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5ad", nargs="+")
    ap.add_argument("--raw-count", action="store_true",
                    help="input is raw counts (--flag-raw-count True); default expects log scale")
    args = ap.parse_args()
    bad = [p for p in args.h5ad if not check(p, args.raw_count)]
    if bad:
        print(f"\npreflight FAILED for {len(bad)} file(s).", file=sys.stderr)
        print("An un-logged 1e4-normalised matrix fed with --flag-raw-count False makes scDRS compute\n"
              "expm1(value): 1388 -> inf, 373 -> 1e162. The mean-variance control fit then degenerates\n"
              "and the scores are meaningless, but scDRS only warns and still writes score files.\n"
              "Fix the input (D4b_pack_atac_log1p.py) or switch the flag, not the figure.", file=sys.stderr)
        return 2
    print(f"preflight: {len(args.h5ad)} file(s) pass", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
