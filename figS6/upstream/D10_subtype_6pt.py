#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D10: fig6D_subtype_6pt.csv - six-subtype six-point cross-channel concordance recomputed with D's new scores (the handover-sheet §3 deliverable).

Why recompute: the six-point table of C/O1 used C's scores (the 5,962-nucleus, 31,263-gene h5ad).
D's dsciCAT arm now uses per-library measured RNA + 30,559 shared genes + the 7,523-nucleus control pool,
so we must confirm the negative conclusion "MS is discordant at the six-point level" still holds under the new scores, not a matrix-swap artefact.

Label source: C's h5ad obs["cell_type"] (the author's six subtypes MOL6/MOL2/MOL3/6/MOL3/MFOL/COP-NFOL),
taking only the nuclei present in both batches. The six-point test beyond Δ1 = for each trait, treat the six subtype mean norm_scores as six points,
cross-channel Spearman (n=6). Per-nucleus cross-channel rho (the same convention as O1's L4) is given alongside as a check.
"""
import os, glob, numpy as np, pandas as pd, anndata as ad
from scipy.stats import spearmanr

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

OUT = WORK + "/figures/routeA_out/D_crosschemistry"
SCDRS_SUB = os.environ.get("SCDRS_SUB", "scdrs")          # rna arm: the delivered tree (not re-run)
SCDRS_SUB_ATAC = os.environ.get("SCDRS_SUB_ATAC", "scdrs_log1p_20261104")  # atac arm: the 2026-11-04 corrected re-run
TAB = os.environ.get("TAB", OUT + "/tables")              # writes to tables_log1p_20261104 after the fix
H5C = WORK + "/FY/oligo/3k_dsci_oligo_subset.h5ad"
O1 = WORK + "/figures/routeA_out/O1_own_multiome_paired/tables"
SAY = lambda fmt, *a: print("[%s] %s" % (pd.Timestamp.now().strftime("%H:%M:%S"), fmt % a if a else fmt), flush=True)
ORDER = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]


def load_scores(job):
    sub = SCDRS_SUB_ATAC if job.endswith("_atac") else SCDRS_SUB
    out = {}
    for f in glob.glob(os.path.join(OUT, sub, job, "**", "*.score.gz"), recursive=True):
        if f.endswith(".full_score.gz"):
            continue
        tr = os.path.basename(f).replace(".score.gz", "")
        if tr in out:
            continue
        out[tr] = pd.read_csv(f, sep="\t", index_col=0)["norm_score"]
    return pd.DataFrame(out)


Sr, Sa = load_scores("dsciCAT_seq_rna"), load_scores("dsciCAT_seq_atac")
assert len(Sr) and len(Sa), "scores not on disk yet"
Oc = ad.read_h5ad(H5C, backed="r").obs
both = sorted(set(Sr.index) & set(Sa.index) & set(Oc.index))
traits = sorted(set(Sr.columns) & set(Sa.columns))
lab = Oc.loc[both, "cell_type"].astype(str).values
S, A = Sr.loc[both, traits], Sa.loc[both, traits]
SAY("nuclei %d | traits %d | subtypes %s", len(both), len(traits), dict(pd.Series(lab).value_counts()))

rows = []
for t in traits:
    m6 = np.array([S[t].values[lab == ct].mean() for ct in ORDER])
    a6 = np.array([A[t].values[lab == ct].mean() for ct in ORDER])
    r, p = spearmanr(m6, a6)
    rows.append(dict(trait=t, rho_6pt=r, P_6pt=p, n_points=6,
                     mean_rna=m6.mean(), mean_atac=a6.mean(),
                     **{"rna_%s" % c: v for c, v in zip(ORDER, m6)}, **{"atac_%s" % c: v for c, v in zip(ORDER, a6)}))
T = pd.DataFrame(rows)
T["sig_0.05"] = T.P_6pt < 0.05
SAY("\nkey three rows (review receipt §10.2 expectations: MS -0.143 / SCZ +0.771 / HEIGHT +0.943)")
key = {"PASS_Multiple_sclerosis": "MS", "PASS_Schizophrenia_Pardinas2018": "SCZ", "UKB_460K.body_HEIGHTz": "HEIGHT"}
print(T[T.trait.isin(key)].assign(name=lambda d: d.trait.map(key))[["name", "rho_6pt", "P_6pt"]].round(4).to_string(index=False))
SAY("all %d traits: median rho=%.3f | fraction rho>0 %.3f | fraction P<0.05 %.3f" %
    (len(T), T.rho_6pt.median(), (T.rho_6pt > 0).mean(), (T.P_6pt < 0.05).mean()))

# per-nucleus cross-channel rho (the same convention as O1's L4: each nucleus's vector correlation across the 74 traits)
def per_cell_rho(X, Y):
    rx = np.argsort(np.argsort(X, axis=1), axis=1).astype(float)
    ry = np.argsort(np.argsort(Y, axis=1), axis=1).astype(float)
    a = rx - rx.mean(1, keepdims=True); b = ry - ry.mean(1, keepdims=True)
    return (a * b).sum(1) / np.maximum(np.sqrt((a ** 2).sum(1) * (b ** 2).sum(1)), 1e-12)
pc = per_cell_rho(S.values, A.values)
null = np.array([per_cell_rho(S.values, A.values[np.random.default_rng(20261103).permutation(len(A))]).mean() for _ in range(200)])
SAY("per-nucleus cross-channel rho: median %.4f mean %.4f (O1 delivered median 0.0756 / mean 0.0790) | shuffled-null mean %.4f",
    np.median(pc), pc.mean(), null.mean())
T["per_nucleus_rho_median_new"] = np.median(pc)
T["per_nucleus_rho_mean_new"] = pc.mean()
T["perm_null_mean"] = null.mean()
T["n_nuclei"] = len(both)
T.round(5).to_csv(TAB + "/fig6D_subtype_6pt.csv", index=False)
SAY("written %s/fig6D_subtype_6pt.csv", TAB)
