#!/usr/bin/env python3
"""A6_fig6_input_transform_audit.py - audits the "input-transform spec" of every scDRS input channel of Fig 6 (2026-11-04)

Background: in its --flag-raw-count False branch, scDRS compute-score applies expm1() to the input matrix,
i.e. it requires the input to already be on the log scale. This script takes every scDRS output tree that Fig 6 uses,
traces it back to the .h5ad hard-coded in the launching script, measures that file's actual scale, and counts the
overflow/expm1 warning lines in the tree's compute.log. Together the two make the "is the input compliant" verdict, not just whether it errored.

Usage: /opt/conda/bin/python3 A6_fig6_input_transform_audit.py [--out CSV]
Read-only; writes no deliverables (except the CSV given by --out).
"""
import os, sys, glob, gzip, argparse
import numpy as np
import anndata as ad

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

D = WORK + "/figures/routeA_out/D_crosschemistry"
P = WORK + "/FY/oligo"
# which chromatin input counts as "current": on delivery day it was h5ad/ (not log-scaled, defective); after the fix it is h5ad_log1p/.
# with STRICT=1, any arm judged "not log-scaled" exits with code 2, used as a hard gate by make_all_current.sh.
DH5 = os.environ.get("D_H5_SUBDIR", "h5ad")
STRICT = os.environ.get("STRICT", "0") == "1"

# (panel, tree name, output-tree dir, input h5ad, input-h5ad origin = the launching script)
TREES = [
    ("A/B MS single-trait (RNA)",  "scdrs_ms_out",        f"{P}/scdrs_ms_out",
     f"{P}/3k_dsci_oligo_subset.h5ad",                 "scripts/scdrs_ms_pipeline.py:101 H5AD_FILE"),
    ("C    RNA 74-trait",          "full74",              f"{P}/scdrs_traits_draft/full74",
     f"{P}/3k_dsci_oligo_subset.h5ad",                 "scripts/MS_GWAS/run_full74.sh:5 H5AD"),
    ("C    chromatin 74-trait",    "out_gs74",            f"{P}/fig6_rerun_20261019/dualomics_probe/out_gs74",
     f"{P}/fig6_rerun_20261019/dualomics_probe/genescore_oligo_5977.h5ad",
     "FY/.../dualomics_probe/run_full74_genescore.sh:6 H5AD"),
    ("C corrob. same-nucleus dual-channel MS", "scdrs_genescore_MS", f"{P}/fig6_rerun_20261019/dualomics_probe/scdrs_genescore_MS",
     f"{P}/fig6_rerun_20261019/dualomics_probe/genescore_oligo_5977.h5ad", "same as above"),
]
for job in ("dsciCAT_seq_rna", "dsciCAT_seq_atac", "10x_rna", "10x_atac",
            "HT_scCAT_seq_HT_rna", "HT_scCAT_seq_HT_atac"):
    TREES.append(("D    " + job, job, f"{D}/scdrs_log1p_20261104/{job}" if DH5 == "h5ad_log1p" and "atac" in job else f"{D}/scdrs/{job}",
                  f"{D}/{DH5}/{job}.h5ad", "scripts/routeA/D5_run_scdrs.sh:31 H5AD=$H5/${job}.h5ad"))
TREES.append(("D (abandoned GIM)", "ABANDONED_GIM", f"{D}/scdrs/ABANDONED_dsciCAT_seq_rna_GIM",
              f"{D}/derived/GENESCORE_FROM_GIM.h5ad", "judged abandoned by V1_rna_source_scale.py"))


def count_warn(tree_dir):
    """lines of expm1 / overflow warnings in compute.log and the number of jobs involved"""
    nline = njob = njobtot = 0
    for f in glob.glob(os.path.join(tree_dir, "**", "compute.log"), recursive=True):
        njobtot += 1
        try:
            txt = open(f, errors="ignore").read()
        except Exception:
            continue
        hits = sum(1 for ln in txt.splitlines()
                   if "overflow" in ln.lower() or "expm1" in ln.lower())
        if hits:
            njob += 1
            nline += hits
    return nline, njob, njobtot


def scale_of(path):
    if not os.path.exists(path):
        return None
    a = ad.read_h5ad(path)
    X = a.X
    from scipy import sparse
    if sparse.issparse(X):
        mx = float(X.data.max()); nz = int(X.nnz)
        ps = np.asarray(X.sum(1)).ravel()                    # per-nucleus sum (raw)
        Xe = X.copy(); Xe.data = np.expm1(np.asarray(Xe.data, dtype=np.float64))
        eps = np.asarray(Xe.sum(1)).ravel()                  # per-nucleus sum (after expm1 = what scDRS actually sees)
        del Xe
    else:
        mx = float(np.nanmax(X)); nz = int((X != 0).sum())
        ps = np.nansum(X, axis=1); eps = np.nansum(np.expm1(X), axis=1)
    return dict(shape=str(a.shape), nnz=nz, max=mx,
                psum_med=float(np.median(ps)), epsum_med=float(np.median(eps)),
                epsum_max=float(np.max(eps)),
                n_over_709=int((np.ravel(X.data if sparse.issparse(X) else X) > 709).sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{D}/tables/fig6_input_transform_audit_20261104"
                    f"{'' if DH5 == 'h5ad' else '_log1p'}.csv")
    a = ap.parse_args()
    rows = []
    for panel, tree, tdir, h5, src in TREES:
        nl, nj, nt = (0, 0, 0)
        if os.path.isdir(tdir):
            nl, nj, nt = count_warn(tdir)
        s = scale_of(h5)
        if s is None:
            verdict = "input file not on disk"
        else:
            logged = s["max"] < 30
            normed = abs(s["epsum_med"] - 1e4) < 500 or abs(s["psum_med"] - 1e4) < 500
            if logged:
                verdict = "log-scale OK (honours --flag-raw-count False)"
            else:
                verdict = "not log-scaled BAD (scDRS will expm1 it again)"
        r = dict(panel=panel, tree=tree, active_input_subdir=DH5, jobs=nt, warn_lines=nl, warn_jobs=nj,
                 h5ad=os.path.basename(h5), src=src, verdict=verdict, **(s or {}))
        rows.append(r)
        print("%-26s jobs=%-4d warn=%-3d/%-3d  max=%-9s expm1 per-nucleus-sum median=%-9s  %s" % (
            panel, nt, nl, nj, ("-" if s is None else f'{s["max"]:.2f}'),
            ("-" if s is None else f'{s["epsum_med"]:.1f}'), verdict))
    import pandas as pd
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    pd.DataFrame(rows).to_csv(a.out, index=False)
    print("\n->", a.out)
    bad = [r for r in rows if r["verdict"].startswith("not log-scaled")]
    if bad and STRICT:
        print("STRICT: %d current inputs non-compliant -> %s" % (len(bad), [r["panel"] for r in bad]), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
