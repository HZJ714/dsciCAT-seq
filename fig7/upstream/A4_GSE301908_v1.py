# -*- coding: utf-8 -*-
"""A4 GSE301908 (template, metadata columns auto-detected): scDRS MS scoring, then donor-level statistics.
Primary criterion: mature (OLIGO) vs precursor (OPC), donor-level paired Wilcoxon
Lesion gradient: between-group MWU whenever metadata carries a pathology / lesion / domain column
Specificity: 5 unrelated traits + random-set percentile (same panel as A2)
"""
import os
import time
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.stats import mannwhitneyu, wilcoxon
import anndata as ad
import scdrs

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

IN = WORK + "/routeA_probe/GSE301908"
OUT = WORK + "/figures/routeA_out/A4_GSE301908_v1"
GS = WORK + "/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs"
os.makedirs(os.path.join(OUT, "tables"), exist_ok=True)
os.makedirs(os.path.join(OUT, "pdf"), exist_ok=True)
t0 = time.time()

def log(*a):
    print(*a, flush=True)

# ---- read the OL subset ----
genes = pd.read_csv(os.path.join(IN, "OL_genes.tsv"), header=None)[0].astype(str).tolist()
barcodes = pd.read_csv(os.path.join(IN, "OL_barcodes.tsv"), header=None)[0].astype(str).tolist()
meta = pd.read_csv(os.path.join(IN, "OL_meta.tsv"), sep="\t")
n_cells = len(barcodes); n_genes = len(genes)
with open(os.path.join(IN, "OL_csr.bin"), "rb") as fh:
    indptr = np.fromfile(fh, dtype=np.int32, count=n_cells + 1)
    indices = np.fromfile(fh, dtype=np.int32, count=indptr[-1])
    data = np.fromfile(fh, dtype=np.float64, count=indptr[-1])
M = csc_matrix((data, indices, indptr), shape=(n_genes, n_cells)).T.tocsr()
log("OL: %d genes x %d nuclei (nnz=%d)" % (n_genes, n_cells, len(data)))

# detect the metadata column names
log("meta columns:", list(meta.columns))

# ---- detect the mature / precursor columns ----
def find_col(keys, exclusions=()):
    for k in keys:
        for c in meta.columns:
            if c.lower() == k.lower() and c not in exclusions:
                return c
    return None

ctcol = find_col(["expertAnno.l1", "majorCluster", "celltype", "cell_type", "cell.type", "annotation",
                  "celltype_annotation", "cellclass", "lineage", "cluster", "predicted.celltype"])
if ctcol is None:
    log("!! no celltype column found; check the meta column names and set it by hand"); sys.exit(1)
ct = meta[ctcol].astype(str)
log("celltype column = %s; value counts:" % ctcol)
print(ct.value_counts().to_string())
mature_kw = ["oligo", "oligodendrocyte"]
prec_kw = ["opc", "precursor", "polydendrocyte"]
is_mat = ct.str.lower().str.contains("|".join(mature_kw), regex=True).values
is_pre = ct.str.lower().str.contains("|".join(prec_kw), regex=True).values
log("mature nuclei %d, precursor nuclei %d" % (is_mat.sum(), is_pre.sum()))
if is_mat.sum() < 50 or is_pre.sum() < 50:
    log("!! too few mature/precursor nuclei, check the annotation"); sys.exit(1)

# ---- donor column (the ms{num} prefix of orig.ident is the donor; r0/1/2 are regions) ----
patcol = find_col(["patient", "donor", "subject", "orig.ident", "Sample", "ID"])
if patcol is None:
    log("!! no donor column found"); sys.exit(1)
raw_pat = meta[patcol].astype(str)
pref = raw_pat.str.extract(r"^(ms\d+)")
meta["patient"] = pref[0].fillna(raw_pat)
log("donor column = %s; %d raw values -> %d donors" % (patcol, raw_pat.nunique(), meta["patient"].nunique()))
log("donor list: %s" % sorted(meta["patient"].unique()))

adata = ad.AnnData(X=M, obs=meta.copy(), var=pd.DataFrame(index=genes))
adata.obs_names = barcodes
scdrs.pp.preprocess(adata)
log("preprocess ok %.0fs" % (time.time() - t0))

with open(GS) as fh:
    gs = {}
    for line in fh:
        sp = line.rstrip("\n").split("\t")
        gs[sp[0]] = [g for g in sp[1].split(",") if g in genes]
MS = gs["PASS_Multiple_sclerosis"]
log("MS set hits %d/%d" % (len(MS), 1000))

res = scdrs.method.score_cell(adata, gene_list=MS, n_ctrl=200, random_seed=0, verbose=False)
score = res["norm_score"].values
meta["MS_norm_score"] = score

mat_med = pd.Series(score[is_mat]).groupby(meta["patient"][is_mat].values).median()
pre_med = pd.Series(score[is_pre]).groupby(meta["patient"][is_pre].values).median()
common = mat_med.index.intersection(pre_med.index)
d = mat_med[common].values - pre_med[common].values
w = wilcoxon(d)
log("\n[primary criterion] mature vs precursor, donor level: n=%d P=%.3g median difference %.3f, same sign %d/%d"
    % (len(d), w.pvalue, np.median(d), int((np.sign(d) == np.sign(np.median(d))).sum()), len(d)))
pd.DataFrame({"patient": common, "mature_minus_precursor": d}).to_csv(
    os.path.join(OUT, "tables", "patient_mature_minus_precursor.csv"), index=False)

meta_out = meta[["MS_norm_score", "patient"] + list(meta.columns)].drop(columns=["MS_norm_score"], errors="ignore")
meta_out.insert(0, "MS_norm_score", meta["MS_norm_score"])
meta_out.to_csv(os.path.join(OUT, "tables", "OL_meta_with_MS_score.csv"), index=False)

# ---- lesion gradient (if a pathology / lesion / domain column exists) ----
for k in ["pathology_class", "pathology", "lesion", "lesion_stage", "domain", "region", "tissue", "group"]:
    c = find_col([k], exclusions={patcol})
    if c:
        grp = meta[c].astype(str)
        log("\nlesion / region column = %s, distribution:" % c)
        print(grp.value_counts().to_string())
        rows = []
        for g in grp.unique():
            med = pd.Series(score[grp == g]).groupby(meta["patient"][grp == g].values).median()
            rows.append({"group": g, "n_donors": len(med), "median_score": med.median()})
        pd.DataFrame(rows).to_csv(os.path.join(OUT, "tables", "group_summary.csv"), index=False)
        break

# ---- specificity: 5 traits + random sets ----
spec_traits = ["PASS_Schizophrenia_Pardinas2018", "PASS_Alzheimers_Jansen2019",
               "PASS_Rheumatoid_Arthritis", "PASS_Type_2_Diabetes",
               "PASS_Intelligence_SavageJansen2018"]
srows = []
for t in spec_traits:
    r = scdrs.method.score_cell(adata, gene_list=gs[t], n_ctrl=200, random_seed=0, verbose=False)
    sc = r["norm_score"].values
    mm = pd.Series(sc[is_mat]).groupby(meta["patient"][is_mat].values).median()
    pp = pd.Series(sc[is_pre]).groupby(meta["patient"][is_pre].values).median()
    cc = mm.index.intersection(pp.index)
    dd = mm[cc].values - pp[cc].values
    ww = wilcoxon(dd)
    srows.append({"trait": t, "n_genes": len(gs[t]), "median_diff": np.median(dd),
                  "P_wilcoxon": ww.pvalue, "n_patients": len(dd),
                  "n_same_sign": int((np.sign(dd) == np.sign(np.median(dd))).sum())})
    log("specificity: %s delta=%.3f P=%.2g same sign %d/%d" % (t, np.median(dd), ww.pvalue, srows[-1]["n_same_sign"], len(dd)))
pd.DataFrame(srows).to_csv(os.path.join(OUT, "tables", "specificity_controls.csv"), index=False)

rng = np.random.default_rng(0)
rnd = []
for i in range(50):
    geneset = list(rng.choice(genes, size=len(MS), replace=False))
    r = scdrs.method.score_cell(adata, gene_list=geneset, n_ctrl=100, random_seed=0, verbose=False)
    sc = r["norm_score"].values
    mm = pd.Series(sc[is_mat]).groupby(meta["patient"][is_mat].values).median()
    pp = pd.Series(sc[is_pre]).groupby(meta["patient"][is_pre].values).median()
    cc = mm.index.intersection(pp.index)
    rnd.append(np.median(mm[cc].values - pp[cc].values))
rnd = np.array(rnd)
pct = (rnd < np.median(d)).mean() * 100
log("\nrandom sets (n=50): median %.3f p90 %.3f max %.3f -> MS=%.3f, percentile %.1f%%"
    % (np.median(rnd), np.percentile(rnd, 90), rnd.max(), np.median(d), pct))
pd.DataFrame({"random_median_diff": rnd}).to_csv(
    os.path.join(OUT, "tables", "random_sets_diffs.csv"), index=False)
log("done in %.0fs" % (time.time() - t0))
