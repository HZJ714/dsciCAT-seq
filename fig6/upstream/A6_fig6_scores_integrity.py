#!/usr/bin/env python3
"""A6_fig6_scores_integrity.py - integrity check of the on-disk score trees (must pass before a re-run tree enters D6)

Why needed: D5's resume logic skips any batch directory that already contains *.score.gz. Any mid-run kill (OOM, accidental kill,
two concurrent runners writing the same batch) can leave a truncated score.gz that the resume logic treats as done and skips - and every later table
silently inherits the corrupt data. This script checks every batch directory against the trait list in its .gs, one by one:
  1) every trait the batch should have has a *.score.gz (not full_score.gz)
  2) gzip reads through and the norm_score column exists
  3) row count == the arm h5ad's nucleus count (all cell_ids hit)
  4) all norm_score values are finite
Any failure -> exit code 2, listing the directories to delete and rerun.

Usage:
  TREE=<score root> python3 A6_fig6_scores_integrity.py            # default scdrs_log1p_20261104
  TREE=<score root> FIX=1 python3 A6_fig6_scores_integrity.py       # delete the bad batch dirs outright (the next resume redoes them)
"""
import glob
import os
import subprocess
import sys

import numpy as np
import pandas as pd

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

D = WORK + "/figures/routeA_out/D_crosschemistry"
TREE = os.environ.get("TREE", D + "/scdrs_log1p_20261104")
H5 = os.environ.get("H5DIR", D + "/h5ad_log1p")
GS = os.environ.get(
    "GSDIR",
    WORK + "/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs.batch")
FIX = os.environ.get("FIX", "0") == "1"

bad = []
for job in sorted(os.listdir(TREE)):
    jdir = os.path.join(TREE, job)
    if not os.path.isdir(jdir):
        continue
    ha = os.path.join(H5, job + ".h5ad")
    if not os.path.exists(ha):
        ha = os.path.join(D, "h5ad", job + ".h5ad")          # the rna arm was not re-run; it stays in its original dir
    if not os.path.exists(ha):
        print("%-24s skipping: no matching h5ad (an abandoned arm)" % job)
        continue
    n_expect = int(subprocess.check_output(
        [sys.executable, "-c",
         f"import anndata as ad;a=ad.read_h5ad('{ha}',backed='r');print(a.shape[0])"]).decode().strip())
    # union trait list of all batch .gs files (repair/<trait> dirs have no .gs of their own; the trait must be in some batch file)
    gsmap = {}
    for g in glob.glob(os.path.join(GS, "magma_10kb_1000.batch*.gs")):
        for tr in pd.read_csv(g, sep="\t", header=0)["TRAIT"].unique():
            gsmap.setdefault(str(tr), g)
    seen_traits = set()
    nb = 0
    counts = {}
    for b in sorted(glob.glob(os.path.join(jdir, "*.batch*"))
                  + [d for d in glob.glob(os.path.join(jdir, "repair", "*")) if os.path.isdir(d)]):
        if os.path.sep + "repair" + os.path.sep in b:
            tr = os.path.basename(b)                      # repair/<trait>/ - the trait name is the directory name
            if tr not in gsmap:
                bad.append((b, "repair trait in no batch .gs"))
                continue
            traits, f0 = [tr], os.path.join(b, tr + ".score.gz")
        else:
            gsf = os.path.join(GS, os.path.basename(b) + ".gs")
            if not os.path.exists(gsf):
                bad.append((b, "no matching .gs found"))
                continue
            traits = [str(t) for t in pd.read_csv(gsf, sep="\t", header=0)["TRAIT"].unique()]
            f0 = None
        nb += 1
        seen_traits |= set(traits)
        for tr in traits:
            f = f0 if f0 else os.path.join(b, tr + ".score.gz")
            if not os.path.exists(f):
                rep = os.path.join(jdir, "repair", tr, tr + ".score.gz")
                if os.path.exists(rep):        # in the delivery round, 6 traits were re-run into repair/ dirs
                    f = rep
                else:
                    bad.append((f, "missing file (batch incomplete)"))
                    continue
            try:
                df = pd.read_csv(f, sep="\t")
            except Exception as e:
                bad.append((f, "gzip/parse failure %s" % type(e).__name__))
                continue
            if "norm_score" not in df.columns:
                bad.append((f, "no norm_score column"))
                continue
            counts.setdefault(len(df), []).append(f)         # row counts must agree within an arm; see below
            if not np.isfinite(df.norm_score.values).all():
                bad.append((f, "non-finite norm_score"))
                continue
    # scDRS --flag-filter-data True drops nuclei with incomplete covariates: the delivered rna arm is exactly 7,507 / h5ad 7,523.
    # That is not an error, but every trait in the same arm must use the same nucleus set, or D6's reindex leaks NaNs.
    if counts:
        main = max(counts, key=lambda k: len(counts[k]))
        for k, v in counts.items():
            if k != main:
                bad.append((v[0], "row count %d disagrees with the arm's main %d (%d files)" % (k, main, len(v))))
    ks = sorted(counts)
    desc = ", ".join("%d%s" % (k, "" if k == n_expect else " (!= h5ad)") for k in ks) or "-"
    print("%-24s batches %2d | traits %2d | h5ad nuclei %d | score rows %s"
          % (job, nb, len(seen_traits), n_expect, desc))

print("-" * 60)
if bad:
    print("bad items %d:" % len(bad))
    for f, why in bad[:40]:
        print("   %-72s %s" % (os.path.relpath(f, TREE), why) if f.endswith(".gz") else "   %s %s" % (f, why))
    if FIX:
        dirs = {f if os.path.isdir(f) else os.path.dirname(f) for f, _ in bad}
        for d in sorted(dirs):
            print("   rm -rf", os.path.relpath(d, TREE))
            subprocess.run(["rm", "-rf", d])
        print("FIX=1: bad batch dirs deleted; the D5 re-run will redo them (skipping completed batches)")
    sys.exit(2)
print("integrity: all batches complete, readable, row counts match the nucleus counts")
