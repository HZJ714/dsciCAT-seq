#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D2f: pack the "measured RNA channel" into chemistry-comparable npz files (genes x cells, values = log1p(CP10k)).

Sources (all the Seurat v4 RNA$data layer = measured expression after LogNormalize, not a CCA-imputed matrix):
  dsciCAT-seq   <- D2e: the RNA$data of the FY13/FY14 3k_*_multiome.RDS, columns cut to the OL nuclei
                   (FY13 30,830 genes / FY14 30,992 genes -> the two libraries intersect first)
  HT-scCAT-seq  <- HT_harmony_expr.mtx of fig5's export_tech_expr.R (the same RNA$data extraction)
  10x Multiome  <- tenx_expr.mtx of fig5's export_tech_expr.R

The cross-chemistry main analysis is restricted to the INTERSECTION of all sources' genes - so "one chemistry measures more genes" is no longer
a confounder of the cross-chemistry rho (handover sheet §2's n_matched_genes requirement); the ATAC channel is already one and the same
24,333-gene GeneScoreMatrix, consistent by construction. The intersection size and the number of dropped genes go into the log and the table.

Column order is unified to the rows of tables/cells_{tag}.csv, pairing nucleus-by-nucleus with gs_{tag}.npz (the ATAC channel).
The ArchR GeneIntegrationMatrix that V1 rejected was renamed rnagim_* in D2d (evidence kept, not the main analysis).
"""
import os, numpy as np, pandas as pd
from scipy.io import mmread
from scipy.sparse import save_npz, csc_matrix, hstack

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

OUT = os.path.join(os.environ.get('DSCI_FIG_ROOT', WORK + "/figures"), 'routeA_out/D_crosschemistry')
DER, TAB = OUT + "/derived", OUT + "/tables"
U = os.environ.get("DSCI_ONLINE_ROOT", ".")
FIG5 = U + "/GitHub/fig5/fig5_out"
SAY = lambda fmt, *a: print("[%s] %s" % (pd.Timestamp.now().strftime("%H:%M:%S"), fmt % a if a else fmt), flush=True)


def mmread_any(path):
    """fig5's *_expr.mtx.gz is actually an uncompressed mtx (the name lies); decide by the magic bytes"""
    import gzip as _gz
    with open(path, "rb") as fh:
        magic = fh.read(2)
    with (_gz.open(path, "rb") if magic == b"\x1f\x8b" else open(path, "rb")) as fh:
        return mmread(fh).tocsc()

# 0) rename the GIM version down (skip if already done)
for tag in ["dsciCAT_seq", "10x", "HT_scCAT_seq_HT"]:
    for src, dst in [("rna_%s.npz" % tag, "rnagim_%s.npz" % tag), ("rnalin_%s.npz" % tag, "rnagimlin_%s.npz" % tag),
                     ("rna_%s_genes.txt" % tag, "rnagim_%s_genes.txt" % tag)]:
        if os.path.exists(os.path.join(DER, src)) and not os.path.exists(os.path.join(DER, dst)):
            os.replace(os.path.join(DER, src), os.path.join(DER, dst))
SAY("GIM version = derived/rnagim_*/rnagimlin_* (rejected by V1; evidence only)")


# 1) read the three chemistries (the two dsciCAT libraries read separately)
sources = []          # (name, matrix, genes, matched_cells_in_order)
for s_ in ["FY13", "FY14"]:
    M = mmread_any(os.path.join(DER, "rnaSeurat_dsciCAT_%s.mtx.gz" % s_))
    genes = [l.strip() for l in open(os.path.join(DER, "rnaSeurat_dsciCAT_%s.mtx_genes.txt" % s_)) if l.strip()]
    cells = pd.read_csv(os.path.join(DER, "rnaSeurat_dsciCAT_%s.mtx_meta.csv" % s_)).cell.tolist()
    SAY("  dsciCAT %s: %s matrix | %d nuclei", s_, M.shape, len(cells))
    sources.append(("dsciCAT_" + s_, M, genes, cells))
want = {}
for tag in ["dsciCAT_seq", "10x", "HT_scCAT_seq_HT"]:
    want[tag] = pd.read_csv(os.path.join(TAB, "cells_%s.csv" % tag)).cell.tolist()
for tag, mtx, meta in [("HT_scCAT_seq_HT", "HT_harmony_expr.mtx.gz", "HT_harmony_meta.csv"),
                       ("10x", "tenx_expr.mtx.gz", "tenx_meta.csv")]:
    M = mmread_any(os.path.join(FIG5, mtx))
    genes = [l.strip() for l in open(os.path.join(FIG5, mtx.replace("_expr.mtx.gz", "_genes.txt"))) if l.strip()]
    cells = pd.read_csv(os.path.join(FIG5, meta), index_col=0).index.tolist()
    assert len(genes) == M.shape[0] and len(cells) == M.shape[1]
    SAY("  %s: %s matrix | %d nuclei", tag, M.shape, len(cells))
    sources.append((tag, M, genes, cells))

# 2) the all-source gene intersection (both dsciCAT libraries take part -> one gene set guaranteed)
inter = sorted(set.intersection(*[set(sg[2]) for sg in sources]))
SAY("RNA intersection genes = %d | per source: %s", len(inter), {s[0]: len(s[2]) for s in sources})
pd.DataFrame({"gene": inter}).to_csv(os.path.join(TAB, "rna_common_genes.csv"), index=False)

# 3) cut per chemistry to (intersection genes x OL nuclei x cells-table order)
rows = []
for tag in ["dsciCAT_seq", "10x", "HT_scCAT_seq_HT"]:
    ws = want[tag]
    if tag == "dsciCAT_seq":
        parts = []
        for nm, M, genes, cells in sources:
            if not nm.startswith("dsciCAT_"):
                continue
            gi = np.array([genes.index(g) for g in inter]); ci = [cells.index(c) for c in ws if c in set(cells)]
            parts.append((M[np.ix_(gi, ci)], [c for c in ws if c in set(cells)]))
        Mall = hstack([p[0] for p in parts], format="csc")
        allc = [c for _, hs in parts for c in hs]
        order = [allc.index(c) for c in ws]
        assert np.array(allc)[order].tolist() == ws, "dsciCAT column alignment failed"
        M = Mall[:, order]
    else:
        nm, M, genes, cells = [s for s in sources if s[0] == tag][0]
        gi = np.array([genes.index(g) for g in inter]); ci = [cells.index(c) for c in ws]
        M = M[np.ix_(gi, ci)]
    assert M.shape == (len(inter), len(ws)), (M.shape, len(inter), len(ws))
    nz = np.asarray((M > 0).sum(axis=0)).ravel(); cs = np.asarray(M.sum(axis=0)).ravel()
    save_npz(os.path.join(DER, "rna_%s.npz" % tag), M)
    with open(os.path.join(DER, "rna_%s_genes.txt" % tag), "w") as fh:
        fh.write("\n".join(inter) + "\n")
    pd.DataFrame({"cell": ws, "nGene_RNA_measured": nz, "logsum_RNA": cs}).to_csv(
        os.path.join(TAB, "rnaseurat_qc_%s.csv" % tag), index=False)
    rows.append(dict(tech=tag, n_nuclei=len(ws), n_genes=len(inter), nnz=int(M.nnz),
                     density_pct=round(100 * M.nnz / (M.shape[0] * M.shape[1]), 3),
                     vmin=round(float(M.data.min()), 4), vmax=round(float(M.data.max()), 4),
                     logsum_median=round(float(np.median(cs)), 1), nGene_median=float(np.median(nz))))
    SAY("  %-18s RNA(measured) %d x %d | nnz=%d | density %.1f%% | range %.3f..%.3f | per-nucleus log-sum median %.0f",
        tag, M.shape[0], M.shape[1], M.nnz, 100 * M.nnz / np.prod(M.shape), M.data.min(), M.data.max(), np.median(cs))
    del M
pd.DataFrame(rows).to_csv(os.path.join(TAB, "rna_seurat_sources.csv"), index=False)
SAY("ALL_FINISHED | unified RNA gene table = %d genes" % len(inter))
