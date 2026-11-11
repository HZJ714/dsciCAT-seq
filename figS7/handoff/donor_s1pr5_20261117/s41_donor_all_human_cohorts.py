"""S1PR5 detection re-tested at the DONOR unit in every human cohort the substrate allows.

Fig 7D draws the detection block from GSE301908 only, so the first donor-level re-test (s1pr5_per_donor_
detection.py) inherited that single cohort. This script asks the same question in all four human cohorts
and is explicit about which ones can answer it:

  GSE301908  full_csr.bin, integer counts, grp -> msNN subject        -> 11 donors   YES
  GSE118257  released count matrix (integers), Sample -> donor        ->  9 donors   YES
  GSE231585  raw_qc.h5ad (integer UMID counts), sample -> donor        ->  8 donors   YES
  GSE277430  the released matrix is log-normalised (values 0-7.2, non-integral), so ">= 1 read" is not
             defined on it; the OL_subset that the audit used holds only OLIGO + OPC, i.e. no comparator
             population exists at all.                                                        -> NOT TESTABLE

Rule, identical in every cohort and identical to the one that draws the red box in Fig 7D:
  margin = OL detection% - strongest non-lineage detection%   (comparators AST MG LYM VASC NEU;
  lineage groups OPC/iOL are excluded by construction, they are the other arm of the contrast)
  rule met  <=>  margin >= 1.0 pp  AND  OL% >= 1.5 x strongest%
A population only enters the comparison if it has >= MIN_POP nuclei in that donor; a donor needs
>= MIN_POP OL nuclei to be evaluable at all.

Every cohort's pooled totals are asserted against S15_det_*.csv / S15_n_*.csv, the tables Fig 7D itself
was drawn from, so a per-donor number cannot drift away from the pooled one it is supposed to explain.

Outputs: S41_det_donor_gene_pop.csv, S41_summary_by_cohort.csv, S41_S1PR5_per_donor.csv  (tables only).
Read-only on all inputs.
"""
import os
import numpy as np
import pandas as pd

ROOT = os.environ.get("DSCI_WORK_ROOT", ".")   # export it; do not hardcode a machine
PROBE = os.path.join(ROOT, "routeA_probe")
A7 = os.path.join(ROOT, "figures", "routeA_out", "A7_GSE231585_rna_v1", "derived")
TAB = os.path.join(ROOT, "figures", "SUPPL_ANALYSES_20261112", "tables")
OUT = os.environ.get("S41_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "tables"))
os.makedirs(OUT, exist_ok=True)

POPS = ["OL", "iOL", "OPC", "AST", "MG", "LYM", "VASC", "NEU"]
NONLIN = ["AST", "MG", "LYM", "VASC", "NEU"]
MIN_POP = 30
GENE = "S1PR5"
MARGIN, FOLD = 1.0, 1.5


def from_csr(D, binf, genef, metaf, labcol, MAP, donorf):
    """the house csr.bin layout: int32 indptr(n+1) | int32 indices(nnz) | float32|float64 data(nnz)"""
    genes = [l.rstrip("\n") for l in open(os.path.join(D, genef))]
    meta = pd.read_csv(os.path.join(D, metaf), sep="\t", dtype=str)
    n = len(meta)
    pop = meta[labcol].map(MAP).values
    donor = donorf(meta)
    with open(os.path.join(D, binf), "rb") as fh:
        ip = np.fromfile(fh, dtype=np.int32, count=n + 1)
    nnz = int(ip[-1])
    size = os.path.getsize(os.path.join(D, binf))
    dt = np.float64 if abs(size - ((n + 1) * 4 + nnz * 4 + nnz * 8)) < 16 else np.float32
    gi = genes.index(GENE)
    det = np.zeros(n, np.int64)
    CH = 40_000_000
    for s in range(0, nnz, CH):
        e = min(s + CH, nnz)
        ix = np.memmap(os.path.join(D, binf), dtype=np.int32, mode="r",
                       offset=(n + 1) * 4, shape=(nnz,))[s:e]
        va = np.memmap(os.path.join(D, binf), dtype=dt, mode="r",
                       offset=(n + 1) * 4 + nnz * 4, shape=(nnz,))[s:e]
        hit = np.nonzero((np.asarray(ix) == gi) & (np.asarray(va) > 0))[0]
        det[np.searchsorted(ip, hit + s, side="right") - 1] += 1
    return donor, pop, (det > 0).astype(int), dt


def a3_matrix():
    """GSE118257: the released gene-by-cell text matrix; only the S1PR5 line is read"""
    D = os.path.join(PROBE, "GSE118257_matrix.txt.gz")
    import gzip
    with gzip.open(D, "rt") as fh:
        bc = fh.readline().rstrip("\n").split("\t")
        vals = None
        for line in fh:
            f = line.split("\t", 1)
            if f[0] == GENE:
                vals = np.fromstring(f[1].rstrip("\n"), sep="\t", dtype=np.int64)
                break
    assert vals is not None and len(vals) == len(bc), (None if vals is None else len(vals), len(bc))
    ann = pd.read_csv(os.path.join(PROBE, "GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz"),
                      sep="\t").drop_duplicates("Detected").set_index("Detected").reindex(bc)
    assert ann.notna().all().all()
    MAP = {**{("Oligo%d" % i): "OL" for i in range(1, 7)}, "ImOlGs": "iOL", "COPs": "iOL", "OPCs": "OPC",
           "Astrocytes": "AST", "Astrocytes2": "AST", "Microglia_Macrophages": "MG", "Macrophages": "MG",
           "Immune_cells": "LYM", "Endothelial_cells1": "VASC", "Endothelial_cells2": "VASC",
           "Pericytes": "VASC", "Vasc_smooth_muscle": "VASC",
           **{("Neuron%d" % i): "NEU" for i in range(1, 6)}}
    donor = ann.Sample.str.extract(r"^(CO\d+|SD48/16|MS\d+)")[0].values
    return donor, ann.Celltypes.map(MAP).values, (vals > 0).astype(int), "counts"


def a7_h5ad():
    import anndata as ad
    import scipy.sparse as sp
    a = ad.read_h5ad(os.path.join(A7, "raw_qc.h5ad"))
    s = ad.read_h5ad(os.path.join(A7, "scored.h5ad"))
    lab = s.obs["celltype"].reindex(a.obs_names).map(
        {"OLIGO": "OL", "OPC": "OPC", "AST": "AST", "MG": "MG", "LYM": "LYM", "END": "VASC", "NEU": "NEU"}).values
    donor = pd.Series(s.obs["sample"].reindex(a.obs_names).values).str.extract(r"^(C\d+|MS\d+)")[0].values
    v = np.asarray(sp.csr_matrix(a.X)[:, list(a.var_names).index(GENE)].todense()).ravel()
    assert np.all(v == np.round(v)), "GSE231585 raw_qc.h5ad is not integer counts"
    return donor, lab, (v > 0).astype(int), "counts"


def pooled_check(name, pop, det):
    """the pooled arm must reproduce the table Fig 7D was drawn from"""
    ref = pd.read_csv(os.path.join(TAB, "S15_det_%s.csv" % name), index_col=0).loc[GENE]
    nref = pd.read_csv(os.path.join(TAB, "S15_n_%s.csv" % name), index_col=0)
    g = pd.DataFrame({"p": pop, "d": det}).groupby("p").d.agg(["sum", "count"])
    for p in ref.index:
        if p in g.index:
            assert g.loc[p, "sum"] == ref[p] and g.loc[p, "count"] == nref.loc[p, "n"], (name, p, g.loc[p].tolist(),
                                                                                         (ref[p], nref.loc[p, "n"]))
    print("  pooled totals reproduce S15_det_%s.csv for every population" % name)


def per_donor(name, donor, pop, det):
    T = pd.DataFrame({"donor": donor, "pop": pop, "det": det})
    rows = []
    for d, sub in T.groupby("donor"):
        g = sub.groupby("pop").det.agg(["sum", "count"])
        nol = int(g.loc["OL", "count"]) if "OL" in g.index else 0
        olp = 100.0 * g.loc["OL", "sum"] / nol if nol else np.nan
        cand = {}
        for p in NONLIN:
            cand[p] = (100.0 * g.loc[p, "sum"] / g.loc[p, "count"]
                       if p in g.index and g.loc[p, "count"] >= MIN_POP else np.nan)
        fin = {p: v for p, v in cand.items() if np.isfinite(v)}
        ev = nol >= MIN_POP and len(fin) > 0
        # the stricter convention used by s1pr5_per_donor_detection.py: every comparator arm has to be
        # evaluable in that donor, otherwise the "strongest other population" is chosen from a subset
        ev_strict = nol >= MIN_POP and len(fin) == len(NONLIN)
        best = max(fin, key=fin.get) if fin else None
        bm = fin.get(best, np.nan)
        marg = olp - bm if ev else np.nan
        fold = olp / bm if ev and np.isfinite(bm) else np.nan
        rows.append(dict(cohort=name, donor=d, n_ol=nol, OL_pct=olp, binding_pop=best, n_binding=(
            int(g.loc[best, "count"]) if best is not None else 0), binding_pct=bm, margin_pp=marg,
            fold=fold, evaluable=ev, evaluable_strict=ev_strict, ol_highest=bool(ev and marg > 0),
            rule=bool(ev and marg >= MARGIN and (fold >= FOLD or bm == 0)), n_comparators=len(fin)))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    frames = []
    print("GSE301908")
    donor, pop, det, dt = from_csr(
        os.path.join(PROBE, "GSE301908_counts"), "full_csr.bin", "full_genes.tsv", "full_meta.tsv",
        "majorCluster", {"Oligo": "OL", "OPC": "OPC", "Astro": "AST", "Micro": "MG", "Lym": "LYM",
                         "Endo": "VASC", "Neuron": "NEU"},
        lambda m: pd.Series(m.grp.str.extract(r"^(ms\d+)")[0].values).fillna(m.grp).values)
    pooled_check("GSE301908", pop, det); frames.append(per_donor("GSE301908", donor, pop, det))
    print("GSE118257")
    donor, pop, det, dt = a3_matrix()
    pooled_check("GSE118257", pop, det); frames.append(per_donor("GSE118257", donor, pop, det))
    print("GSE231585")
    donor, pop, det, dt = a7_h5ad()
    pooled_check("GSE231585", pop, det); frames.append(per_donor("GSE231585", donor, pop, det))

    D = pd.concat(frames, ignore_index=True)
    D.to_csv(os.path.join(OUT, "S41_S1PR5_per_donor.csv"), index=False)
    ev = D[D.evaluable]
    S = ev.groupby("cohort").agg(donors_evaluable=("donor", "size"), ol_highest=("ol_highest", "sum"),
                                 rule_met=("rule", "sum"), median_margin_pp=("margin_pp", "median"),
                                 median_OL_pct=("OL_pct", "median"), median_binding_pct=("binding_pct", "median"))
    # keyed by name: the groupby sorts, the frames do not
    by = {f.cohort.iloc[0]: f for f in frames}
    S["donors_all_convention"] = [len(by[c][by[c].evaluable_strict]) for c in S.index]
    S["rule_met_all_convention"] = [int(by[c][by[c].evaluable_strict].rule.sum()) for c in S.index]
    S["total_donors_in_cohort"] = [len(by[c]) for c in S.index]
    S.loc["ALL"] = [ev.donor.size, ev.ol_highest.sum(), ev.rule.sum(), ev.margin_pp.median(),
                    ev.OL_pct.median(), ev.binding_pct.median(),
                    int(D.evaluable_strict.sum()), int(D[D.evaluable_strict].rule.sum()), D.donor.nunique()]
    S.to_csv(os.path.join(OUT, "S41_summary_by_cohort.csv"))
    from scipy.stats import binomtest
    k, n = int(ev.ol_highest.sum()), len(ev)
    print("\n%s" % S.to_string())
    print("\ndirection across cohorts: OL is the highest-detecting population in %d/%d evaluable donors, "
          "sign test p = %.2g" % (k, n, binomtest(k, n, 1 / 6.0, alternative="greater").pvalue))
    print("full red-box rule met in %d/%d evaluable donors" % (int(ev.rule.sum()), n))
    print("\nnot testable: GSE277430 (released matrix is log-normalised, and the OL_subset the audit used "
          "carries no non-lineage population)")
