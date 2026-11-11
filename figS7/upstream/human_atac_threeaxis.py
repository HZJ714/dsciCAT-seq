# -*- coding: utf-8 -*-
"""human_atac_threeaxis.py — chromatin (ATAC) counterpart of the RNA three-axis ledger: per-library d1 (GSE231584, 9 blind libraries / 10,870 nuclei), 74-trait chromatin ranks, cross-axis pairs.

M4_atac_supplement.py - the ATAC (chromatin) counterpart of the Fig S7 numbers (2026-11-08).

chromatin rank bands, the Fig7C chromatin cross-axis. This script computes what had never been computed;

human (GSE231584 snATAC, 9 blinded libraries / 10,870 nuclei):
  A. Fig7A version: per-library chromatin d1 for the MS trait (OL_mature - OPC), with library-level case/ctrl labels.
  B. Fig7B version: chromatin d1 rank band for the 74 traits (the same-donor table already carries
     rank_CHR, this only assembles it).  C. Fig7C version: cross-axis (d1 OL / d2A Astro-[Micro+Lymph] /
     d2B Neuron-non-neuronal) over 74 traits, robust z per axis, MS rank, between-axis rho - new computation.
mouse (GSE250589, m3 3 mice):  B. 74-trait ATAC d1 rank band (M3_paired3 has it, only assembled), plus
     rho(ln n_peak, ATAC d1).  C. cross-axis: not computable at the data ceiling (all non-OL nuclei in the
     3 mice ~ 86: Astro 62 / Micro 17 / Neuron 1) - recorded.  Also: the same-nucleus paired 74-trait
     version waits for the RNA arm (M3_full74_paired3.py refresh).

Output: SUPPORT_v1/tables/S7_atac_supplement.csv (+ *_crossaxis_human.csv) + a printed summary."""
import os, numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

R = WORK + "/figures/routeA_out"
OUT = R + "/SUPPORT_v1/tables"
os.makedirs(OUT, exist_ok=True)
MS = "PASS_Multiple_sclerosis"
nice = lambda f, *a: print("[%s] %s" % (pd.Timestamp.now().strftime("%H:%M:%S"), f % a if a else f), flush=True)

# ---------------------------------------------------------------- human
A6 = R + "/A6_GSE231584_atac_v1"
ann = pd.read_pickle(A6 + "/derived/annot.pkl")          # barcode/lib/celltype
scores = np.load(A6 + "/tables/trait_cell_scores.npz", allow_pickle=True)
traits = sorted(scores.files)
lib_labels = pd.read_csv(A6 + "/tables/atac_library_labels_SEXVALIDATED.csv")
lab = dict(zip(lib_labels.suffix.astype(int), zip(lib_labels.donor, lib_labels.domain)))

def d_axis(ct_hi, ct_lo):
    """delta = med(score|ct_hi) - med(score|ct_lo), pooled; returns one Series per trait."""
    out = {}
    for t in traits:
        s = scores[t]
        m = ann.celltype.isin(ct_hi).values
        n = ann.celltype.isin(ct_lo).values
        if m.sum() < 10 or n.sum() < 10:
            out[t] = np.nan
        else:
            out[t] = float(np.median(s[m]) - np.median(s[n]))
    return pd.Series(out)

OL = ["OL_mature"]; OPC = ["OPC"]; AST = ["Astro"]; IMM = ["Micro", "Lymph"]; NEU = ["Neuron"]
d1 = d_axis(OL, OPC); d2a = d_axis(AST, IMM); d2b = d_axis(NEU, ["Astro", "Micro", "Lymph", "OL_mature", "OPC", "Endo", "Epen"])

# per-library MS d1 (Fig7A version)
ms_s = scores[MS]
rows = []
for lib, g in ann.groupby("lib"):
    a = g.celltype.eq("OL_mature").values; b = g.celltype.eq("OPC").values
    if a.sum() < 10 or b.sum() < 10:
        continue
    rows.append(dict(lib=int(lib), n_ol=int(a.sum()), n_opc=int(b.sum()),
                     d1_ms=float(np.median(ms_s[g.index[a]]) - np.median(ms_s[g.index[b]]))))
HU = pd.DataFrame(rows)
labmap = HU.lib.map(lambda x: lab.get(int(x), ("NA", "NA")))
HU["donor"], HU["domain"] = zip(*labmap)
HU = HU.sort_values("lib")

def robust_z(s):
    m, mad = s.median(), (s - s.median()).abs().median()
    return (s - m) / (1.4826 * mad)

RZ = pd.DataFrame({"d1": d1, "d2a": d2a, "d2b": d2b}).dropna()
rz = RZ.apply(robust_z)
ms_rank = {k: int((rz[k] > rz[k].loc[MS]).sum()) + 1 for k in rz.columns}
r12 = spearmanr(RZ.d1, RZ.d2a, nan_policy="omit"); r13 = spearmanr(RZ.d1, RZ.d2b, nan_policy="omit")
r23 = spearmanr(RZ.d2a, RZ.d2b, nan_policy="omit")
CA = pd.DataFrame({
    "axis": ["OL maturation (Δ1)", "astro − immune (Δ2A)", "neuron − glia (Δ2B)"],
    "n_trait": [RZ.d1.notna().sum(), RZ.d2a.notna().sum(), RZ.d2b.notna().sum()],
    "MS_robust_z": [round(rz.d1.loc[MS], 3), round(rz.d2a.loc[MS], 3), round(rz.d2b.loc[MS], 3)],
    "MS_rank": [ms_rank["d1"], ms_rank["d2a"], ms_rank["d2b"]]})
CA.to_csv(OUT + "/S7_atac_crossaxis_human.csv", index=False)
CAX = pd.concat([RZ, rz], axis=1, keys=["Δ", "robust z"])
CAX.to_csv(OUT + "/S7_atac_crossaxis_human_full.csv")

nice("human cross-axis: rho(d1,d2A)=%.3f rho(d1,d2B)=%.3f rho(d2A,d2B)=%.3f (n=%d)" % (r12.statistic, r13.statistic, r23.statistic, len(RZ)))
nice("human cross-axis MS: OL rank %d/%d, AST-IMM rank %d, NEU rank %d" %
     (ms_rank["d1"], len(RZ), ms_rank["d2a"], ms_rank["d2b"]))
nice("human per-library MS d1 (ATAC): %s" % HU[["lib", "donor", "d1_ms"]].to_string(index=False).replace("\n", " | "))
if HU.donor.nunique() >= 2 and HU.groupby("donor").size().min() >= 1:
    # keep only donors with a definite label and a pairing
    pass
# case/ctrl test (only when there are >= 2 controls)
cc = HU[HU.donor.str.startswith("MS")].d1_ms.values
ctr = HU[HU.donor.str.startswith("C")].d1_ms.values
nice("human ATAC MS d1 case/ctrl: n_case=%d n_ctrl=%d  median %.3f vs %.3f" %
     (len(cc), len(ctr), np.median(cc) if len(cc) else np.nan, np.median(ctr) if len(ctr) else np.nan))

# ---------------------------------------------------------------- assemble rank bands (Fig7B version)
sd = pd.read_csv(A6 + "/tables/same_donor_RNA_vs_CHROMATIN_delta1_74traits.csv")
m3 = pd.read_csv(WORK + "/routeA_probe/m3_probe/tables_full74/M3_paired3_atac_vs_rna_74.csv")
band = sd[["trait", "rank_RNA", "rank_CHR"]].merge(
    m3[["trait", "rna_d1_3mouse", "atac_d1"]], on="trait", how="left")
band = band.rename(columns={"rank_RNA": "human_RNA_rank", "rank_CHR": "human_CHR_rank"})
band["mouse_RNA_rank_3mice"] = band.rna_d1_3mouse.rank(ascending=False).astype(int)
band["mouse_ATAC_rank_3mice"] = band.atac_d1.rank(ascending=False).astype(int)
band = band.drop(columns=["rna_d1_3mouse", "atac_d1"])
band.to_csv(OUT + "/S7_atac_rank_bands.csv", index=False)
nice("rank bands assembled: human RNA/CHR + mouse RNA/ATAC (3 mice); MS = %s" %
     band[band.trait == MS].iloc[0].to_dict())

# mouse chromatin n_peak confound (Fig7B version footnote)
lnp = m3.set_index("trait").ln_n_peak
r_m = spearmanr(m3.atac_d1, m3.ln_n_peak, nan_policy="omit")
nice("mouse ATAC d1 x ln n_peak rho = %.3f (n=%d) | human chromatin d1 x ln n_peak rho = -0.712 (recorded in D3)"
     % (r_m.statistic, m3.atac_d1.notna().sum()))

# ---------------------------------------------------------------- mouse cross-axis not computable (recorded)
A = pd.read_pickle(WORK + "/routeA_probe/mouse_EAE/anno_min.pkl")
A["bc"] = A.bc16.astype(str).str.split("-").str[0]
import glob
cnt = {}
for lib in ("LateFemale4", "EarlyCtrl4", "Late5"):
    f = sorted(glob.glob(WORK + "/routeA_probe/m3_probe/scdrs/%s/atacgs/*/*.score.gz" % lib))[0]
    d = pd.read_csv(f, sep="\t", index_col=0)
    bc = [str(x).split("-")[0] for x in d.index]
    cnt[lib] = A[A.bc.isin(set(bc))].cellType.value_counts().to_dict()
tot = {k: sum(v.values()) for k, v in cnt.items()}
nice("non-OL nuclei covered by the 3 mouse ATAC libraries: %s -> cross-axis not computable (Astro/Micro/Neuron together < 100)" % cnt)

# ---------------------------------------------------------------- summary rows
summary = pd.DataFrame([
    dict(item="human GSE231584 cross-axis rho(d1,d2A)", value="%.3f" % r12.statistic),
    dict(item="human GSE231584 cross-axis rho(d1,d2B)", value="%.3f" % r13.statistic),
    dict(item="human GSE231584 MS chromatin rank (same_donor)", value="%d/74" % sd[sd.trait == MS].rank_CHR.iloc[0]),
    dict(item="mouse GSE250589 MS ATAC rank (3 mice)", value="%d/74" % band[band.trait == MS].mouse_ATAC_rank_3mice.iloc[0]),
    dict(item="mouse ATAC d1 x ln n_peak rho", value="%.3f" % r_m.statistic),
    dict(item="human ATAC per-library MS d1 usable libraries", value=str(len(HU)))])
summary.to_csv(OUT + "/S7_atac_supplement_summary.csv", index=False)
nice("done. output -> SUPPORT_v1/tables/S7_atac_*.csv")
