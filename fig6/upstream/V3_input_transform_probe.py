#!/usr/bin/env python3
"""V3 — is the panel-D chromatin channel fed to scDRS in the right transform?

scDRS `compute-score --flag-raw-count False` applies expm1() to the input matrix internally, so the
input must be log-transformed.  The three panel-D *expression* matrices are (Seurat `data` slot,
expm1 sums to exactly 1e4 per nucleus).  The three panel-D *chromatin* matrices are the ArchR
GeneScoreMatrix, 1e4-normalised but NOT logged (max 1388.6 / 297.1 / 372.7), so expm1 overflows to inf
in the own arm and distorts the mean-variance trend in the other two.

This script has two halves.

  # 1. write the corrected inputs (log1p of the same matrices) and print the scdrs commands
  python3 V3_input_transform_probe.py prepare  <out_dir>
  #    -> then run the printed scdrs calls (MAGMA batch 5 = MS / PBC / ReactionTime, ~6 min total
  #       for all three chemistries; ~3.5 min per trait for the own arm)

  # 2. compare the re-scored batch against the delivered scores
  python3 V3_input_transform_probe.py compare  <delivered_scdrs_root> <probe_root1> <probe_root2> ...
  #    -> prints and returns the per-nucleus Spearman and the shift of delta1 = med(Oligo)-med(OPC)

  Whole-panel fix (once decided): rerun 74 traits x 3 chemistries = 222 chromatin runs, then rerun
  D6 -> the notebook -> compose.  Expression arms, panels A/B/C and all page gates are unaffected.

Measured result of the batch-5 probe (see V3_input_transform_probe.csv):

  dsciCAT-seq(own) MS  rho_nuc 0.464  delta1 +0.262 -> -0.151   (sign flip)
  10x Multiome     MS  rho_nuc 0.914  delta1 -0.430 -> -0.474
  HT-scCAT-seq     MS  rho_nuc 0.373  delta1 -0.297 -> +0.064   (sign flip)
"""
import glob
import os
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

ROOT = os.environ.get("DSCI_FIG_ROOT", WORK + "/figures/routeA_out")
DC = os.path.join(ROOT, "D_crosschemistry")
ARMS = ["dsciCAT_seq_atac", "10x_atac", "HT_scCAT_seq_HT_atac"]
LABELS = {"dsciCAT_seq_atac": "cells_dsciCAT_seq.csv",
          "10x_atac": "cells_10x.csv",
          "HT_scCAT_seq_HT_atac": "cells_HT_scCAT_seq_HT.csv"}
GS_BATCH = os.path.join(os.environ.get("DSCI_FY_ROOT", WORK + "/FY/oligo"),
                        "scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs.batch/"
                        "magma_10kb_1000.batch5.gs")
SCDRS = os.environ.get("SCDRS_BIN", "scdrs")
MAT, OPC = "Oligo", "OPC"


def prepare(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    cmds = []
    for j in ARMS:
        src = os.path.join(DC, "h5ad", j + ".h5ad")
        a = ad.read_h5ad(src)
        X = sp.csr_matrix(a.X).astype(np.float64)
        raw_max, raw_sum = X.data.max(), np.median(np.asarray(X.sum(1)).ravel())
        X.data = np.log1p(X.data)
        h5 = os.path.join(out_dir, "%s_log1p.h5ad" % j)
        ad.AnnData(X=X, obs=a.obs.copy(), var=a.var.copy(), uns=a.uns).write_h5ad(h5)
        cov = os.path.join(out_dir, "%s_log1p.cov.tsv" % j)
        subprocess.run(["cp", os.path.join(DC, "h5ad", j + ".cov.tsv"), cov], check=True)
        print("[V3] %-24s un-logged max %8.1f  per-nucleus sum %8.1f  -> log1p max %.3f"
              % (j, raw_max, raw_sum, X.data.max()))
        cmds.append('%s compute-score --h5ad-file %s --h5ad-species mouse --gs-file %s \\\n'
                    '    --gs-species human --out-folder %s/out_%s --flag-filter-data True \\\n'
                    '    --flag-raw-count False --n-ctrl 1000 --flag-return-ctrl-raw-score False \\\n'
                    '    --flag-return-ctrl-norm-score True --cov-file %s'
                    % (SCDRS, h5, GS_BATCH, out_dir, j, cov))
    print("\n# run these:\n" + "\n".join(cmds))
    return cmds


def _load(d):
    res = {}
    for g in sorted(glob.glob(d + "/**/*.score.gz", recursive=True)):
        if "full_score" in g:
            continue
        t = pd.read_csv(g, sep="\t", index_col=0)
        if "norm_score" not in t:
            continue
        res[os.path.basename(g).replace(".score.gz", "")] = (
            np.array([str(i) for i in t.index]), t["norm_score"].values)
    return res


def compare(delivered_root, probe_roots):
    rows = []
    for j in ARMS:
        old = _load(os.path.join(delivered_root, j))
        new = None
        for p in probe_roots:
            cand = {}
            for cand_dir in (p if os.path.basename(p).startswith("out_") else "",
                               os.path.join(p, j), os.path.join(p, "out_" + j)):
                cand = _load(cand_dir) if cand_dir else {}
                if cand:
                    break
            if cand:
                new = cand
                break
        if not new:
            print("[V3] no probe output for", j)
            continue
        lab = pd.read_csv(os.path.join(DC, "tables", LABELS[j])).set_index("cell")["label"]
        for tr in sorted(new):
            if tr not in old:
                continue
            ca, sa = old[tr]
            cb, sb = new[tr]
            m = pd.concat([pd.Series(sb, index=cb).rename("new"),
                           pd.Series(sa, index=ca).rename("old")], axis=1, join="inner")
            cls = lab.reindex(m.index).values
            x, y = m["new"].values, m["old"].values
            o, p = cls == MAT, cls == OPC
            d_new = np.median(x[o]) - np.median(x[p])
            d_old = np.median(y[o]) - np.median(y[p])
            rows.append(dict(tech=j, trait=tr, n=int((o | p).sum()),
                             rho_nuc=round(spearmanr(x, y).statistic, 4),
                             delta1_delivered=round(d_old, 4), delta1_log1p_input=round(d_new, 4),
                             shift=round(d_new - d_old, 4),
                             sign_flip=bool(np.sign(d_old) != np.sign(d_new))))
    d = pd.DataFrame(rows)
    print(d.to_string(index=False))
    return d


if __name__ == "__main__":
    if sys.argv[1] == "prepare":
        prepare(sys.argv[2])
    elif sys.argv[1] == "compare":
        compare(sys.argv[2], sys.argv[3:])
    else:
        sys.exit(__doc__)
