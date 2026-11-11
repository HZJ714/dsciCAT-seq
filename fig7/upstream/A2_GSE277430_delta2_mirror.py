# -*- coding: utf-8 -*-
"""A2_GSE277430_delta2_mirror.py -- the mirror control: the 74 traits projected onto a SECOND
identity axis, delta2 (i.e. not OL maturation).  Review receipt §3, priority P0.

  delta2A = donor-level median(AST) - median(IMM)
  delta2B = donor-level median(NEU) - median(GLIA)

Reading intended: MS should peak on delta1 (the OL-maturation side) while SCZ / IQ traits should peak
here (the neuronal side).  That contrast is what makes the MS peak specific to the OL-maturation axis
instead of being an artefact of "high scDRS score anywhere".

Outputs (written under OUT/tables):
  trait_sweep_delta2.csv    the two mirror axes for all 74 traits
  trait_sweep_dual_axis.csv trait_sweep_delta2.csv merged onto the delivered delta1 sweep
                            (trait_sweep_74.csv); this merged file is the Panel C input for GSE277430

Core parameters, frozen: n_ctrl=100, random_seed=0, wilcoxon-free (donor medians only).  n_ctrl here
is 100, not the 200 the delta1 sweeps use - it is the value the delivered table was produced with and
is recorded as such (house spec: state the parameter that made the delivered number).

Runtime: the scDRS scoring dominates (~40 s/trait on one core -> ~50 min for 74 traits).
Run under the scoring kernel $DSCI_PY (scdrs 1.0.2).
"""
import os
import time
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.stats import wilcoxon
import anndata as ad
import scdrs

# Root of the author's analysis tree (FY/, figures/, routeA_probe/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
OUT = WORK + "/figures/routeA_out/A2_GSE277430_v1"
IN = WORK + "/routeA_probe/GSE277430/full"
GS = WORK + "/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs"
os.makedirs(os.path.join(OUT, "tables"), exist_ok=True)
t0 = time.time()

def log(*a):
    print(*a, flush=True)

genes = pd.read_csv(os.path.join(IN, "full_genes.tsv"), header=None)[0].astype(str).tolist()
barcodes = pd.read_csv(os.path.join(IN, "full_barcodes.tsv"), header=None)[0].astype(str).tolist()
meta = pd.read_csv(os.path.join(IN, "full_meta.tsv"), sep="\t")
n_cells = len(barcodes); n_genes = len(genes)
with open(os.path.join(IN, "full_csr.bin"), "rb") as fh:
    indptr = np.fromfile(fh, dtype=np.int32, count=n_cells + 1)
    indices = np.fromfile(fh, dtype=np.int32, count=indptr[-1])
    data = np.fromfile(fh, dtype=np.float64, count=indptr[-1])
M = csc_matrix((data, indices, indptr), shape=(n_genes, n_cells)).T.tocsr()
log("whole lineage %d x %d (nnz=%d)" % (n_genes, n_cells, len(data)))

adata = ad.AnnData(X=M, obs=meta.set_index("barcodes").loc[barcodes],
                   var=pd.DataFrame(index=genes))
adata.obs_names = barcodes
scdrs.pp.preprocess(adata)
log("preprocess ok %.0fs" % (time.time() - t0))

meta["patient"] = meta["patient"].astype(str)
l1 = meta["expertAnno.l1"].astype(str)
groups = {
    "AST": l1 == "AST",
    "IMM": l1 == "IMM",
    "NEU": l1.isin(["EXC NEU", "INH NEU"]),
    "GLIA": l1.isin(["AST", "IMM", "OLIGO", "OPC", "VAS"]),
}
log("nuclei per axis group: AST %d, IMM %d, NEU %d, GLIA %d" % (
    groups["AST"].sum(), groups["IMM"].sum(), groups["NEU"].sum(), groups["GLIA"].sum()))

with open(GS) as fh:
    traits = {}
    for line in fh:
        sp = line.rstrip("\n").split("\t")
        traits[sp[0]] = [g for g in sp[1].split(",") if g in genes]

def donor_med(score, mask, gname):
    s = pd.Series(score[mask], index=meta["patient"][mask].values)
    return s.groupby(level=0).median().rename(gname)

rows = []
for i, (trait, glist) in enumerate(traits.items()):
    ts = time.time()
    res = scdrs.method.score_cell(adata, gene_list=glist, n_ctrl=100,
                                  random_seed=0, verbose=False)
    sc = res["norm_score"].values
    ast = donor_med(sc, groups["AST"], "AST")
    imm = donor_med(sc, groups["IMM"], "IMM")
    neu = donor_med(sc, groups["NEU"], "NEU")
    gli = donor_med(sc, groups["GLIA"], "GLIA")
    c1 = ast.index.intersection(imm.index)
    c2 = neu.index.intersection(gli.index)
    dA = ast[c1].values - imm[c1].values
    dB = neu[c2].values - gli[c2].values
    rows.append({"trait": trait,
                 "dA_AST_minus_IMM": np.median(dA), "nA_donors": len(dA),
                 "dB_NEU_minus_GLIA": np.median(dB), "nB_donors": len(dB)})
    log("[%d/74] %s dA=%.3f dB=%.3f (%.0fs)" % (
        i + 1, trait, rows[-1]["dA_AST_minus_IMM"], rows[-1]["dB_NEU_minus_GLIA"],
        time.time() - ts))

tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(OUT, "tables", "trait_sweep_delta2.csv"), index=False)
# merge onto delta1 when the delivered sweep table is present
p1 = os.path.join(OUT, "tables", "trait_sweep_74.csv")
if os.path.exists(p1):
    d1 = pd.read_csv(p1)
    merged = d1[["trait", "rank", "median_diff", "P_wilcoxon"]].merge(tab, on="trait")
    merged.to_csv(os.path.join(OUT, "tables", "trait_sweep_dual_axis.csv"), index=False)
    log("dual-axis table written: tables/trait_sweep_dual_axis.csv")
log("done %.0fs" % (time.time() - t0))
