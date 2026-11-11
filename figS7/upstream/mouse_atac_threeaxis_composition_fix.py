#!/usr/bin/env python3
"""mouse_atac_threeaxis_composition_fix.py — aggregation fix for the mouse ATAC run: the wide score table lacked the cellType/model columns.

M6_p2_fix.py - repair of the M6 P2 aggregation (the wide score tables lack a cellType column)

M6 P1 finished completely (25 libraries x 74 traits), but score_lib wrote the wide table without
cellType/model, so P2's d_axis hit AttributeError on S.cellType. This script:
  (1) reads each m6_scores/*.scores.csv.gz, refills cellType/model from anno_min.pkl by obs="gsm|bc16"
      (bc16 is unique inside a suffix), overwrites in place - after that a re-run of M6 P1 skips the
      existing libraries and P2 works;  (2) replicates M6's P2 exactly: d1/d2a/d2b aggregation -> robust z -> ranks -> FULL/SUM/split-half/summary json."""
import os, re, glob, json, time
import numpy as np, pandas as pd

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

MOUSE_EAE = WORK + "/routeA_probe/mouse_EAE"
SUP = WORK + "/figures/routeA_out/SUPPORT_v1/tables"
OUT = WORK + "/routeA_probe/m3_probe/m6_scores"
SEED = 20261030
MOL = {"MOL1", "MOL2", "MOL56"}; OPC_ = {"OPC", "COP"}
KEEP = MOL | OPC_ | {"Astro", "Micro", "Neuron", "Epen"}
MS = "PASS_Multiple_sclerosis"
t0 = time.time()
def log(*a): print("[%4.0fs]" % (time.time() - t0), *a, flush=True)

# ---------- label mapping (bc16 is unique inside a suffix but repeats across suffixes -> map gsm->suffix) ----------
A = pd.read_pickle(MOUSE_EAE + "/anno_min.pkl")
A["bc16"] = A["bc16"].astype(str); A["suffix"] = A.suffix.astype(str)
assert A.groupby("suffix").bc16.apply(lambda s: s.is_unique).all(), "bc16 not unique inside a suffix"
LAB = {sf: g.set_index("bc16")[["cellType", "model"]] for sf, g in A.groupby("suffix")}
MAP = pd.read_csv(SUP + "/S7_atac_mouse_gsm_suffix_map.csv").set_index("gsm")
def _sfx(v):
    try:
        return str(int(float(v)))
    except Exception:
        return "NA"
MAP["suffix"] = MAP.suffix.map(_sfx)
log("anno loaded %d rows / %d suffixes; gsm->suffix table %d rows" % (len(A), len(LAB), len(MAP)))

# ---------- (1) refill the columns in place ----------
files = sorted(glob.glob(OUT + "/*.scores.csv.gz"))
n_fix, miss = 0, []
for f in files:
    W = pd.read_csv(f, index_col=0)
    if "cellType" in W.columns:
        continue
    gsm = f.split("/")[-1].split(".")[0]
    sfx = MAP.loc[gsm, "suffix"]
    lab = LAB[sfx]
    bc = pd.Index([i.split("|")[1] for i in W.index])
    ct = lab["cellType"].reindex(bc).values
    md = lab["model"].reindex(bc).values
    if pd.isna(ct).any():
        miss.append((f, int(pd.isna(ct).sum())))
    W["cellType"] = ct; W["model"] = md
    W.to_csv(f, compression="gzip")
    n_fix += 1
log("refilled %d files; missing mapping %s" % (n_fix, miss or "none"))

# ---------- (2) P2 aggregation (replicates M6) ----------
S = pd.concat([pd.read_csv(f, index_col=0) for f in sorted(glob.glob(OUT + "/*.scores.csv.gz"))])
S = S[~S.index.duplicated()]
S["gsm"] = S.index.map(lambda i: i.split("|")[0])
log("per-nucleus scores pooled %d nuclei / %d libraries" % (len(S), S.gsm.nunique()))
OTHER_OF_NEU = set(KEEP) - {"Neuron"}

def d_axis(hi, lo):
    o = {}
    for tr in TRAITS:
        a = S.loc[S.cellType.isin(hi).values, tr]; b = S.loc[S.cellType.isin(lo).values, tr]
        o[tr] = float(a.median() - b.median()) if min(a.notna().sum(), b.notna().sum()) >= 10 else np.nan
    return pd.Series(o)

import scipy.stats as st
GS_FILE = WORK + "/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs"
import scdrs
GA = pd.read_csv(SUP + "/S7_atac_fullstats_GENE_STATS.csv", index_col=0)
dict_gs = scdrs.util.load_gs(GS_FILE, src_species="human", dst_species="mouse", to_intersect=set(GA.index))
TRAITS = sorted(dict_gs.keys())

AX = pd.DataFrame({"d1": d_axis(MOL, OPC_), "d2a": d_axis({"Astro"}, {"Micro"}),
                   "d2b": d_axis({"Neuron"}, OTHER_OF_NEU)})
def rz(s):
    m, mad = s.median(), (s - s.median()).abs().median()
    return (s - m) / (1.4826 * mad)
OK = AX.dropna()
RZ = OK.apply(rz); RANK = OK.rank(ascending=False).astype(int)
r12 = st.spearmanr(OK.d1, OK.d2a); r13 = st.spearmanr(OK.d1, OK.d2b)
FULL = pd.concat([AX.add_prefix("Δ"), RZ.add_prefix("z_"), RANK.add_prefix("rk_")], axis=1)
FULL.index.name = "trait"; FULL.round(4).to_csv(SUP + "/S7_atac_crossaxis_mouse_full.csv")

n_by = S.cellType.value_counts().to_dict()
NAME = {"d1": "OL maturation (Δ1)", "d2a": "astro − immune (Δ2A)", "d2b": "neuron − non-neuron (Δ2B)"}
GN = {"d1": "MOL %d vs OPC/COP %d" % (sum(n_by.get(k, 0) for k in MOL), n_by.get("OPC", 0) + n_by.get("COP", 0)),
      "d2a": "Astro %d vs Micro %d" % (n_by.get("Astro", 0), n_by.get("Micro", 0)),
      "d2b": "Neuron %d vs the rest %d" % (n_by.get("Neuron", 0), len(S) - n_by.get("Neuron", 0))}
SUM = pd.DataFrame([dict(axis=NAME[k], n_trait=int(OK[k].notna().sum()), MS_robust_z=round(float(RZ.loc[MS, k]), 3),
                         MS_rank=int(RANK.loc[MS, k]), group_n=GN[k]) for k in ("d1", "d2a", "d2b")])
SUM.to_csv(SUP + "/S7_atac_crossaxis_mouse.csv", index=False)

# split-half
mol = S[S.cellType.isin(MOL).values]
rng2 = np.random.default_rng(SEED + 1)
half = pd.Series(rng2.permutation(len(mol)) < len(mol) // 2, index=mol.index)
dsh = {}
for s in (True, False):
    sub = mol[half == s]
    oth = S[S.cellType.isin(OPC_).values]
    dsh[s] = pd.Series({tr: sub[tr].median() - oth[tr].median() for tr in TRAITS})
rk_h = pd.DataFrame(dsh).rank(ascending=False).astype(int)
stab = dict(ms_rank_halves=[int(rk_h.loc[MS, True]), int(rk_h.loc[MS, False])],
            rho_halves=round(float(st.spearmanr(dsh[True], dsh[False]).statistic), 3))

MAP = pd.read_csv(SUP + "/S7_atac_mouse_gsm_suffix_map.csv")
amb = MAP[MAP.overlap_frac < 0.5]
json.dump(dict(n_lib=int(S.gsm.nunique()), n_nuclei_scored=int(len(S)), n_trait=len(TRAITS),
               celltype_n=n_by, group_n=GN,
               rho_d1_d2a=round(float(r12.statistic), 3), rho_d1_d2b=round(float(r13.statistic), 3),
               p_d1_d2a=float(r12.pvalue), p_d1_d2b=float(r13.pvalue),
               ms_rank=dict(d1=int(RANK.loc[MS, "d1"]), d2a=int(RANK.loc[MS, "d2a"]), d2b=int(RANK.loc[MS, "d2b"])),
               ms_robust_z={k: round(float(RZ.loc[MS, k]), 3) for k in ("d1", "d2a", "d2b")},
               sd_axis={k: round(float(AX[k].std()), 3) for k in ("d1", "d2a", "d2b")},
               split_half=stab, MOL_CAP=150, PAR=12, seed=SEED, dropped_libs=amb.gsm.tolist()),
          open(SUP + "/S7_atac_crossaxis_mouse_summary.json", "w"), indent=1)
log(SUM.to_string(index=False))
log("rho(d1,d2A)=%.3f P=%.2g | rho(d1,d2B)=%.3f P=%.2g | split-half d1 rank %s rho=%.3f" % (
    r12.statistic, r12.pvalue, r13.statistic, r13.pvalue, stab["ms_rank_halves"], stab["rho_halves"]))
log("done")
