#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S6v16_export_caches.py — build the small input caches that Figure_S6.ipynb reads.

Why these files exist
---------------------
Four S6 panels start from objects that are far too large to ship with a figure
repository (a GENCODE GTF, a 5,977-nucleus h5ad, and 6 x 74-trait scDRS score
sweeps totalling ~6 GB).  Each panel, however, only ever touches a handful of
columns of those objects.  This script extracts exactly those columns once, from
the authoritative analysis tree, and writes them as flat CSVs:

  s6a_gene_coords.csv        chrom, tss, name   for every gene on the four chromosomes
                             that host the S6A loci (PHLDB1/chr11, PHGDH/chr1,
                             BCAS1/chr20, PLLP/chr16).  Panel A only ever asks
                             "which TSS falls inside lead-SNP +/- 500 kb", so the
                             chromosome list, not the window, is the filter — a
                             changed window therefore cannot silently miss genes.
  s6c_umap_score.csv         UMAP_1, UMAP_2, cell_type, PASS_Multiple_sclerosis_rna
                             for the nuclei that have both (panel C; identical
                             join to the delivered Fig 6B chain).
  s6f_<job>_<chan>.csv       cell, norm_score (MS only), label — the two columns
                             panel F uses, for the three chemistries x two channels.
  fig6_cellscore.csv         cell_type, UMAP_1, UMAP_2, norm_score (Fig 6 A/B).

Discipline
----------
* Nothing is recomputed statistically here: no medians, no correlations, no
  p-values.  This is a column projection of delivered objects, so a cached value
  and a freshly computed value are the same number by construction (asserted
  below for the F arms against the delivered MWU ledger is *not* done here — the
  panel script does that assertion every run).
* Paths come from environment variables only (DSCI_WORK_ROOT, ...); the script
  never writes outside the cache directory it is pointed at.
* Re-running is idempotent and cheap; the notebooks read the caches if present
  and fall back to the full objects otherwise.

Run:
  DSCI_WORK_ROOT=<work> \
  S6_CACHE=<repo>/figS6/figS6_out/caches \
  python3 S6v16_export_caches.py
"""
import os
import sys

import numpy as np
import pandas as pd

WORK = os.environ.get('DSCI_WORK_ROOT', '')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                      # ms_locus_overlap.py lives next to this file
sys.path.insert(0, os.path.join(WORK, "scripts"))

CACHE = os.environ.get("S6_CACHE", os.path.join(HERE, "..", "figS6_out", "caches"))
os.makedirs(CACHE, exist_ok=True)

SCMS = os.environ.get("S6_SCMS", os.path.join(WORK, "FY", "oligo", "scdrs_ms_out"))
R = os.environ.get("DSCI_FIG_ROOT", os.path.join(WORK, "figures", "routeA_out"))
LOCUS = os.environ.get("S6_LOCUS_TABLE", os.path.join(SCMS, "ms_locus_overlap_table.csv"))
GTF = os.environ.get("S6_GTF", os.path.join(SCMS, "gencode.v45.basic.annotation.gtf.gz"))
H5 = os.environ.get("S6_H5AD", os.path.join(WORK, "FY", "oligo", "3k_dsci_oligo_subset.h5ad"))
CELLTAB = os.environ.get("S6_CELL_LEVEL",
                         os.path.join(R, "O1_own_multiome_paired", "tables", "cell_level_table.csv"))
DCHEM = os.path.join(R, "D_crosschemistry")
CELLS = {"dsciCAT_seq": "cells_dsciCAT_seq.csv", "10x": "cells_10x.csv",
         "HT_scCAT_seq_HT": "cells_HT_scCAT_seq_HT.csv"}
JOBS = {"dsciCAT_seq": ("dsciCAT_seq_rna", "dsciCAT_seq_atac"),
        "10x": ("10x_rna", "10x_atac"),
        "HT_scCAT_seq_HT": ("HT_scCAT_seq_HT_rna", "HT_scCAT_seq_HT_atac")}
MS = "PASS_Multiple_sclerosis"


def log(*a):
    print("[caches]", *a, flush=True)


# ---------------------------------------------------------------- 1. S6A gene coordinates
HIT = pd.read_csv(LOCUS)
SEL4 = ["PHLDB1", "PHGDH", "BCAS1", "PLLP"]
chroms = sorted({str(HIT[HIT.gene == g].iloc[0]["chr"]) for g in SEL4})
log("S6A loci on chromosomes %s" % chroms)

import ms_locus_overlap as mlo                      # noqa: E402  (delivered window logic)
GENES = mlo.load_genes(GTF)
keep = [(c, t, n) for (c, t, n) in GENES if str(c) in chroms]
pd.DataFrame(keep, columns=["chrom", "tss", "name"]).to_csv(
    os.path.join(CACHE, "s6a_gene_coords.csv"), index=False)
log("s6a_gene_coords.csv  %d genes (of %d parsed)" % (len(keep), len(GENES)))
# every window must actually contain genes, otherwise panel A would draw an empty track
for g in SEL4:
    h = HIT[HIT.gene == g].iloc[0]
    lo, hi = int(h.snp_pos) - 500_000, int(h.snp_pos) + 500_000
    n = sum(1 for c, t, nm in keep if str(c) == str(h["chr"]) and lo <= t <= hi)
    assert n > 0, "locus %s has no cached genes in +/-500 kb" % g
    log("  %-7s chr%s %d-%d : %d neighbouring genes" % (g, h["chr"], lo, hi, n))

# ---------------------------------------------------------------- 2. S6C / Fig6 UMAP + score
import anndata as ad                                # noqa: E402
obs = ad.read_h5ad(H5, backed="r").obs
CT = pd.read_csv(CELLTAB, low_memory=False)
dU = obs[["UMAP_1", "UMAP_2", "cell_type"]].join(
    CT.set_index("cell_id")[[MS + "_rna"]], how="inner")
dU = dU.astype({"UMAP_1": float, "UMAP_2": float, MS + "_rna": float})
dU["cell_type"] = dU["cell_type"].astype(str)
dU.to_csv(os.path.join(CACHE, "s6c_umap_score.csv"))
log("s6c_umap_score.csv   %d nuclei" % len(dU))

ms_score = pd.read_csv(os.path.join(SCMS, MS + ".score.gz"), sep="\t", index_col=0)["norm_score"]
df6 = obs[["cell_type", "UMAP_1", "UMAP_2"]].join(ms_score.to_frame("norm_score")).dropna()
df6["cell_type"] = df6.cell_type.astype(str)
df6.to_csv(os.path.join(CACHE, "fig6_cellscore.csv"))
log("fig6_cellscore.csv   %d nuclei (Fig 6 A/B)" % len(df6))
assert len(df6) == 5962, "Fig 6 A/B nucleus count moved: %d" % len(df6)

# ---------------------------------------------------------------- 3. S6F per-nucleus MS scores
for tag, (rna_job, atac_job) in JOBS.items():
    lab = pd.read_csv(os.path.join(DCHEM, "tables_log1p_20261104", CELLS[tag])).set_index("cell")["label"]
    for chan, job in (("rna", rna_job), ("atac", atac_job)):
        # rna arm = delivered tree; atac arm = the log1p-corrected rerun (2026-11-04)
        base = os.path.join(DCHEM, "scdrs" if chan == "rna" else "scdrs_log1p_20261104", job)
        fs = [os.path.join(d, f) for d, _, fl in os.walk(base) for f in fl
              if f == MS + ".full_score.gz"]
        assert fs, "no score files under %s" % base
        df = pd.concat([pd.read_csv(f, sep="\t", index_col=0, low_memory=False) for f in fs])
        df = df[df.index.isin(lab.index)].copy()
        df["label"] = lab.loc[df.index].values
        df[["norm_score", "label"]].to_csv(os.path.join(CACHE, "s6f_%s.csv" % job))
        log("s6f_%-22s %6d nuclei  mean %+0.3f" % (job + ".csv", len(df), df.norm_score.mean()))
        # 10x (854) and HT (571) really are small libraries; the count is a
        # property of the deposited data, not of this extraction.
        assert len(df) > 100, job

log("done -> %s" % CACHE)
