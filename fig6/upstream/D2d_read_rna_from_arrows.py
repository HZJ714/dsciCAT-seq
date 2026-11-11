#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D2d: take the GeneIntegrationMatrix (RNA channel) directly from the ArchR ArrowFiles.

Why not load the whole thing in R: handover sheet §2 requires the RNA to come from the same ArchR project too
(GeneIntegrationMatrix: Info/Units="NormCounts", Params=cca, per-nucleus sum exactly 1e4,
 the three chemistries share one §4/§5 pipeline = the tangible evidence of an "identical pipeline").
But getMatrixFromProject sucks the full 22,112 x 70,643 matrix into R; anon measured at 28 GB
(cgroup 32 GB) and gets OOM-killed - so it is done per arrow file, per chromosome, keeping only the OL-lineage columns.

Scales: two outputs are written
  rna_{tag}.mtx.gz     = log1p(CP10k)   the same scale as C's RNA channel (Seurat LogNormalize data) -> to scDRS
  rnalin_{tag}.mtx.gz  = linear CP10k   for the V1 scale check only

Built-in check: the same reader rebuilds the 10x OL subset of the GeneScoreMatrix
and compares it element-wise against the gs_10x.mtx.gz that D1 exported via the ArchR R interface; it must be exactly identical.
"""
import os, gzip
import numpy as np, pandas as pd, h5py
from scipy.sparse import csc_matrix, coo_matrix, load_npz, save_npz
from scipy.io import mmread, mmwrite

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

PRJ = os.path.join(os.environ.get("DSCI_ONLINE_ROOT", "."), "FY/archr-subset-copnfol")
ARR = os.path.join(PRJ, "ArrowFiles")
OUT = os.path.join(os.environ.get('DSCI_FIG_ROOT', WORK + "/figures"), 'routeA_out/D_crosschemistry')
DER, TAB, CHK = OUT + "/derived", OUT + "/tables", OUT + "/derived/_check"
os.makedirs(CHK, exist_ok=True)
# ArchR geneDF / arrow block row order = chr name first (alphabetical, matching the group order inside the arrow), then FeatureDF$idx
CHR_ORDER = sorted(["chr" + str(c) for c in list(range(1, 20)) + ["X"]])
TECH2TAG = {"dsciCAT-seq": "dsciCAT_seq", "HT-scCAT-seq-HT": "HT_scCAT_seq_HT", "10x": "10x"}
ARROWS = {"dsciCAT-seq": ["FY13.arrow", "FY14.arrow"],
          "HT-scCAT-seq-HT": ["EL13.arrow", "EL14.arrow", "EL15.arrow", "EL16.arrow"],
          "10x": ["10x.arrow"]}
def SAY(fmt, *a):
    msg = fmt % a if a else fmt
    print("[%s] %s" % (pd.Timestamp.now().strftime("%H:%M:%S"), msg), flush=True)


def arrow_cells(path, mat):
    with h5py.File(path, "r") as f:
        return [c.decode() if isinstance(c, bytes) else str(c) for c in np.asarray(f[mat]["Info"]["CellNames"])]


def read_matrix(path, mat, want_cells):
    """-> (csc matrix genes x len(want_cells), gene-name list). Reads only the needed columns, chunked per chr."""
    with h5py.File(path, "r") as f:
        info = f[mat]["Info"]
        fd = np.asarray(info["FeatureDF"])
        names = [x.decode() if isinstance(x, bytes) else str(x) for x in fd["name"]]
        seqn = [x.decode() if isinstance(x, bytes) else str(x) for x in fd["seqnames"]]
        start = np.asarray(fd["start"]); idxg = np.asarray(fd["idx"])
        cn = [x.decode() if isinstance(x, bytes) else str(x) for x in np.asarray(info["CellNames"])]
        pos = {c: k for k, c in enumerate(cn)}
        miss = [c for c in want_cells if c not in pos]
        if miss:
            raise KeyError("%s: %d nuclei not in this sample (e.g. %s)" % (os.path.basename(path), len(miss), miss[:3]))
        cols0 = np.array([pos[c] for c in want_cells], dtype=np.int64)
        col_out = {int(v): k for k, v in enumerate(cols0)}      # this sample's column numbers -> output-matrix column numbers
        bychr = {}
        for k, g in enumerate(seqn):
            bychr.setdefault(g, []).append(k)
        # output row order = the ArchR geneDF order: chr names alphabetically first, idx within a chr - exactly as in the D1/ArchR R interface
        gseq = np.concatenate([np.array(sorted(bychr[c], key=lambda k: idxg[k]), dtype=np.int64)
                               for c in CHR_ORDER if c in bychr])
        rowpos = np.empty(len(idxg), dtype=np.int64); rowpos[gseq] = np.arange(len(gseq))
        names_out = [names[k] for k in gseq]
        RO, CO, VA = [], [], []
        for c in CHR_ORDER:
            if c not in bychr:
                continue
            order = np.array(sorted(bychr[c], key=lambda k: idxg[k]), dtype=np.int64)  # rows within the block = that chr's idx ascending
            g = f[mat][c]
            jL = np.asarray(g["jLengths"][0], dtype=np.int64)
            jV = np.asarray(g["jValues"][0], dtype=np.int64) - 1   # the block stores only non-empty columns: jV = the column's 0-based global number
            nrow = len(np.asarray(g["rowSums"][0]))
            if nrow != len(order):
                raise RuntimeError("%s %s: block rows %d != that chr's gene count %d" % (os.path.basename(path), c, nrow, len(order)))
            p = np.concatenate([[0], np.cumsum(jL)])
            colpos = {int(v): k for k, v in enumerate(jV)}
            x = np.asarray(g["x"][0], dtype=np.float64)
            i = np.asarray(g["i"][0], dtype=np.int64)
            for ci in cols0:                                  # column by column: the OL nuclei are few enough
                k = colpos.get(int(ci))
                if k is None:
                    continue
                s, e = p[k], p[k + 1]
                if e <= s:
                    continue
                RO.append(rowpos[order[i[s:e] - 1]]); CO.append(np.full(e - s, col_out[int(ci)])); VA.append(x[s:e])
        if not RO:
            raise RuntimeError("empty matrix " + path)
        R = np.concatenate(RO); C = np.concatenate(CO); V = np.concatenate(VA)
        return csc_matrix((V, (R, C)), shape=(len(names_out), len(cols0))), names_out


def write_mtx(path, M):
    """writes binary npz (fast, compact); the mtx text version is exported separately only when manual checking needs it."""
    from scipy.sparse import save_npz
    save_npz(path, M.tocsc())


def main():
    meta = pd.read_csv(os.path.join(TAB, "fig2_light_meta.csv"))
    SAY("light meta %s | chemistries %s", meta.shape, meta.technique.value_counts().to_dict())

    # ---------- check: rebuild the GeneScoreMatrix (10x OL 854 nuclei) and compare element-wise with the D1 export ----------
    ref = pd.read_csv(os.path.join(TAB, "cells_10x.csv"))
    G, gn = read_matrix(os.path.join(ARR, "10x.arrow"), "GeneScoreMatrix", ref.cell.tolist())
    d1 = mmread(os.path.join(DER, "gs_10x.mtx.gz"))
    d1g = [l.strip() for l in open(os.path.join(DER, "gs_10x_genes.txt")) if l.strip()]
    SAY("[check] rebuilt %s (nnz %d) vs D1 %s | gene names agree %s",
        G.shape, G.nnz, d1.shape, gn == d1g)
    assert gn == d1g, "gene-name order disagrees"
    d1 = d1.tocsc()
    assert d1.shape == G.shape, "shape differs: D1 %s vs rebuilt %s" % (d1.shape, G.shape)
    dif = abs(d1 - G).max()
    SAY("[check] max absolute difference = %s -> the reader and the ArchR R interface %s", dif, "exactly the same source" if dif == 0 else "NOT identical")
    assert dif == 0.0

    # ---------- export the RNA channel per chemistry ----------
    for tech, tag in TECH2TAG.items():
        tab = pd.read_csv(os.path.join(TAB, "cells_%s.csv" % tag))
        cells = tab.cell.tolist()
        mats = []
        for a in ARROWS[tech]:
            path = os.path.join(ARR, a)
            cn = set(arrow_cells(path, "GeneIntegrationMatrix"))
            loc = [c for c in cells if c in cn]
            if not loc:
                SAY("  %s: no target nuclei", a); continue
            M, gnames = read_matrix(path, "GeneIntegrationMatrix", loc)
            cs = np.asarray(M.sum(axis=0)).ravel()
            SAY("  %-11s %d genes x %d nuclei | nnz=%d | per-nucleus sum 1e4: %s (min %.1f max %.1f)",
                a, M.shape[0], M.shape[1], M.nnz, bool(np.all(np.abs(cs - 1e4) < 1e-6)), cs.min(), cs.max())
            if 'gnames_ref' in dir():
                assert gnames == gnames_ref, "gene tables disagree across samples"
            gnames_ref = gnames
            mats.append((M, loc))
        # after concatenating samples, reorder to the rows of cells_{tag}.csv (exactly the same column order as gs_{tag}.mtx.gz -> same-nucleus pairing)
        from scipy.sparse import hstack
        Mlin = hstack([m for m, _ in mats], format="csc")
        allcols = [c for _, loc in mats for c in loc]
        order = np.array([allcols.index(c) for c in cells])
        assert np.array(allcols)[order].tolist() == cells
        Mlin = Mlin[:, order]
        Mlog = Mlin.copy(); Mlog.data = np.log1p(Mlog.data)
        fp_l = os.path.join(DER, "rna_%s.npz" % tag)
        fp_r = os.path.join(DER, "rnalin_%s.npz" % tag)
        write_mtx(fp_l, Mlog); write_mtx(fp_r, Mlin)
        with open(os.path.join(DER, "rna_%s_genes.txt" % tag), "w") as fh:
            fh.write("\n".join(gnames_ref) + "\n")
        SAY("  wrote %s + %s | %d genes x %d nuclei | log nnz=%d | log range %.4f..%.3f | per-nucleus log-sum median %.0f",
            os.path.basename(fp_l), os.path.basename(fp_r), len(gnames_ref), Mlin.shape[1], Mlog.nnz,
            Mlog.data.min(), Mlog.data.max(), np.median(np.asarray(Mlog.sum(axis=0)).ravel()))
        # per-nucleus detected gene count (backup/cross-check covariate)
        pd.DataFrame({"cell": cells, "tech": tech, "label": tab.label,
                      "nGene_GIM": np.asarray((Mlin > 0).sum(axis=0)).ravel(),
                      "umi_sum_GIM": np.asarray(Mlin.sum(axis=0)).ravel()}).to_csv(
            os.path.join(TAB, "gim_qc_%s.csv" % tag), index=False)
        del Mlin, Mlog
    SAY("ALL_FINISHED")


if __name__ == "__main__":
    main()
