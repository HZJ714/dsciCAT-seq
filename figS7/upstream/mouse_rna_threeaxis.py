"""mouse_rna_threeaxis.py — mouse RNA three-axis run (GSE250589, 25 libraries), supplying the d2A/d2B columns that the RNA side had never computed; identical nuclei, labels, gene sets and seed as the ATAC run.

M7_rna_mouse_threeaxis.py - mouse RNA (GSE250589 atlas, 25 libraries, same nuclei/labels/seed as M6) three axes.

  Only d1 had ever been delivered for mouse RNA (M1/M3); d2A/d2B were never computed. This script fills
  them in with the M6 basis strictly matched:
  (1) the same GSM<->sample map (M6 writes S7_atac_mouse_gsm_suffix_map.csv, overlap>=0.5);
  (2) the same nucleus selection (anno_min.pkl, same SEED=20261030, same MOL_CAP/CAP ceilings);
  (3) the same scorer scdrs.method.score_cell(n_ctrl=200, random_seed=SEED), RNA channel = official
      Gene Expression -> CP10K -> log1p (as in the M1 mouse RNA pipeline; no gene_set_stats passed, scDRS
      computes them);  (4) the same aggregation: d1=MOL-OPC/COP, d2A=Astro-Micro, d2B=Neuron-non-Neuron,
      per-library median -> robust z -> ranks. Output: S7_rna_crossaxis_mouse_full.csv / S7_rna_crossaxis_mouse.csv / S7_rna_mouse_threeaxis_summary.json (all new file names)."""
import os, re, sys, glob, json, time, gzip
import numpy as np, pandas as pd

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
os.environ.setdefault("TQDM_DISABLE", "1")
import h5py

MOUSE_EAE = WORK + "/routeA_probe/mouse_EAE"
PROBE = WORK + "/routeA_probe/m3_probe"
SUP = WORK + "/figures/routeA_out/SUPPORT_v1/tables"
OUT = PROBE + "/m7_scores"
os.makedirs(OUT, exist_ok=True)
GS_FILE = WORK + "/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs"
MS = "PASS_Multiple_sclerosis"
SEED = 20261030
MOL_CAP = int(os.environ.get("M7_MOL_CAP", os.environ.get("M6_MOL_CAP", "150")))
CAP = {"Astro": 150, "Micro": 120, "OPC": 120, "COP": 120, "Epen": 80, "Neuron": 10_000}
PAR = int(os.environ.get("M7_PAR", "4"))
DRY = os.environ.get("M7_DRY") == "1"
SMOKE = os.environ.get("M7_SMOKE") == "1"
RNA_GS_FILE = SUP + "/S7_rna_fullstats_GENE_STATS.csv"
if not os.path.exists(RNA_GS_FILE):                        # wait for M7_gene_stats (RNA full-data moments)
    import time as _t
    for _ in range(180):
        if os.path.exists(RNA_GS_FILE):
            break
        _t.sleep(30)
    if not os.path.exists(RNA_GS_FILE):
        raise SystemExit("RNA full-data GENE_STATS %s missing (M7_gene_stats.py not finished)" % RNA_GS_FILE)
t0 = time.time()
def log(*a): print("[%4.0fs]" % (time.time() - t0), *a, flush=True)

GS = pd.read_csv(RNA_GS_FILE, index_col=0)
UNION = set(GS.index)
log("RNA full-data GENE_STATS loaded: %d genes (M7_gene_stats, all nuclei of 28 libraries)" % len(UNION))

MOL = {"MOL1", "MOL2", "MOL56"}; OPC_ = {"OPC", "COP"}
KEEP = MOL | OPC_ | {"Astro", "Micro", "Neuron", "Epen"}
OTHER_OF_NEU = set(KEEP) - {"Neuron"}

# ---------- 1. reuse the M6 GSM<->sample map (never re-guess it) ----------
MF = SUP + "/S7_atac_mouse_gsm_suffix_map.csv"
if not os.path.exists(MF):
    raise SystemExit("M6 map %s missing (M6 did not reach that step)" % MF)
MAP = pd.read_csv(MF).set_index("gsm")
def _sfx(v):                                        # the CSV wrote suffix as a float ("1.0"); normalise back to the anno string key
    try:
        return str(int(float(v)))
    except Exception:
        return "NA"
MAP["suffix"] = MAP.suffix.map(_sfx)
GOOD = MAP[MAP.overlap_frac >= 0.5]
log("map reuse: %d/%d libraries usable" % (len(GOOD), len(MAP)))

# ---------- 2. nucleus selection with exactly the M6 rule ----------
A = pd.read_pickle(MOUSE_EAE + "/anno_min.pkl")
A["bc16"] = A["bc16"].astype(str); A["suffix"] = A.suffix.astype(str)
anno = A[A.cellType.isin(KEEP)][["bc16", "suffix", "cellType", "orig.ident", "model"]]
rngm = np.random.default_rng(SEED)
sel = []
for gsm, r in GOOD.iterrows():
    g = anno[anno.suffix == r["suffix"]]
    if not len(g):
        continue
    parts = []
    for ct, sub in g.groupby("cellType"):
        cap = MOL_CAP if ct in MOL else CAP.get(ct, 200)
        parts.append(sub if len(sub) <= cap else sub.iloc[rngm.permutation(len(sub))[:cap]])
    sel.append(pd.concat(parts).assign(gsm=gsm))
SEL = pd.concat(sel, ignore_index=True)
vc = SEL.cellType.value_counts().to_dict()
log("selected %d nuclei / %d libraries | %s" % (len(SEL), SEL.gsm.nunique(), vc))
EXPECT = {"MOL56": 3750, "MOL2": 3660, "Astro": 2124, "Micro": 1583, "OPC": 1552, "Epen": 995,
          "COP": 673, "MOL1": 513, "Neuron": 78}          # measured from the M6 log (14,928 nuclei / 25 libraries)
assert len(SEL) == sum(EXPECT.values()) and all(vc.get(k) == v for k, v in EXPECT.items()), \
    "nucleus selection differs from M6 (expected %d nuclei, per class %s), stopping" % (sum(EXPECT.values()), EXPECT)

import anndata as ad, scdrs
dict_gs = scdrs.util.load_gs(GS_FILE, src_species="human", dst_species="mouse", to_intersect=UNION)
TRAITS = sorted(dict_gs.keys())
if SMOKE:
    TRAITS = [MS] + [t for t in TRAITS if t != MS][:3]
log("gs loaded %d traits%s (to_intersect=RNA UNION, same as M6)" % (len(TRAITS), " (SMOKE)" if SMOKE else ""))

# ---------- 3. per-library RNA channel + 74-trait scoring (checkpoint per trait) ----------
def rna_matrix(gsm, barcodes):
    with h5py.File("%s/%s_gex.h5" % (MOUSE_EAE, gsm), "r") as h:
        g = h["matrix"]
        gt = np.array([x.decode() if isinstance(x, bytes) else x for x in g["features"]["feature_type"][:]])
        nm = np.array([x.decode() if isinstance(x, bytes) else x for x in g["features"]["name"][:]])
        bc_all = np.array([x.decode().split("-")[0] for x in g["barcodes"][:]])
        full = sp_csc(g)
    gi = np.where(gt == "Gene Expression")[0]
    want = np.isin(bc_all, barcodes)
    X = full[gi][:, want].T.toarray().astype(np.float32)
    names = np.unique(nm[gi])
    # official symbol de-duplication: keep the maximum row per duplicated name (equivalent to the M1 var_names_make_unique preprocessing, but reproducible)
    idx = {}
    for j, s in enumerate(nm[gi]):
        s = str(s)
        if s not in idx or X[:, j].max() > X[:, idx[s]].max():
            idx[s] = j
    keep = sorted(idx.values()); sym = [str(nm[gi][j]) for j in keep]
    return X[:, keep], sym, bc_all[want]

def sp_csc(g):
    from scipy.sparse import csc_matrix
    return csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=(int(g["shape"][0]), int(g["shape"][1])))

def _pick(a):
    g = SEL[SEL.gsm == a.uns["gsm"]]
    lab = pd.Series(g.cellType.values, index=pd.Index(g.bc16)).groupby(level=0).first()
    present = [c for c in a.obs_names if c.split("|")[1] in set(lab.index)]
    grp = {}
    for ct in sorted(set(lab.values)):
        bcs = set(lab[lab == ct].index)
        ids = [i for i in a.obs_names if i.split("|")[1] in bcs]
        if ids:
            grp[ct] = ids
    sel = sorted({i for v in grp.values() for i in v})
    return sel, grp

def score_lib(gsm):
    dst = "%s/%s.scores.csv.gz" % (OUT, gsm)
    if os.path.exists(dst):
        return gsm, ("already present, skipped", 0)
    mark_ = "%s/%s.done" % (OUT, gsm)
    long_ = "%s/%s.long.csv" % (OUT, gsm)
    done_tr = []
    if os.path.exists(mark_):
        done_tr = [x.strip() for x in open(mark_) if x.strip() in TRAITS]
    import pandas as _pd
    g = SEL[SEL.gsm == gsm]
    X, sym, bc = rna_matrix(gsm, g.bc16.values)
    colmask = np.array([x in UNION for x in sym])          # keep only genes covered by the RNA GENE_STATS (same as M6)
    X = X[:, colmask]; sym = [x for x, k in zip(sym, colmask) if k]
    Xn = np.log1p(X * (1e4 / np.maximum(X.sum(1), 1)[:, None])).astype(np.float32)
    a = ad.AnnData(X=Xn)
    a.var_names = pd.Index([sym[i] for i in range(len(sym))])
    a.var_names_make_unique()
    a.obs_names = pd.Index(["%s|%s" % (gsm, b) for b in bc])
    a.uns["gsm"] = gsm
    a.uns["SCDRS_PARAM"] = {"FLAG_SPARSE": False, "FLAG_COV": False, "FLAG_ADJ_PROP": False,
                            "GENE_STATS": GS.loc[a.var_names].copy()}
    sel, grp = _pick(a)
    out = pd.DataFrame(index=a.obs_names)
    lab = pd.Series(g.cellType.values, index=pd.Index(g.bc16)).groupby(level=0).first()
    out["cellType"] = lab.reindex(pd.Index([i.split("|")[1] for i in out.index])).values
    out["model"] = pd.Series(g["model"].values, index=pd.Index(g.bc16)).groupby(level=0).first().reindex(
        pd.Index([i.split("|")[1] for i in out.index])).values
    assert out.cellType.notna().all(), gsm
    tlib = time.time()
    if done_tr:
        log("    %s resumed: %d/%d traits already done" % (gsm, len(done_tr), len(TRAITS)))
    for k, tr in enumerate(TRAITS):
        if tr in done_tr:
            continue
        if k % 10 == 0:
            log("    %s trait %2d/%d (library elapsed %4.0fs)" % (gsm, k, len(TRAITS), time.time() - tlib))
        use = [x for x in dict_gs[tr][0] if x in set(a.var_names)]
        if len(use) < 5:
            v = np.full(len(a.obs_names), np.nan)
        else:
            v = scdrs.method.score_cell(a, use, n_ctrl=200, random_seed=SEED, verbose=False)["norm_score"].values
        _pd.DataFrame({"obs": a.obs_names, "trait": tr, "score": v}).to_csv(
            long_, mode="a", header=not os.path.exists(long_), index=False)
        with open(mark_, "a") as fh:
            fh.write(tr + "\n")
    L = _pd.read_csv(long_)
    W = L[L.trait.isin(set(open(mark_).read().split()))].pivot(index="obs", columns="trait", values="score")
    for tr in TRAITS:
        if tr not in W.columns:
            W[tr] = np.nan
    W = W[TRAITS]
    W["cellType"] = out.cellType.values; W["model"] = out.model.values
    W.to_csv(dst, compression="gzip")
    return gsm, W.shape

libs = sorted(SEL.gsm.unique())
if SMOKE:
    libs = libs[:1]
if DRY:
    libs = libs[:2]; log("DRY mode: only %s" % libs)
log("P1 start %d libraries x %d traits (PAR=%d, MOL_CAP=%d)" % (len(libs), len(TRAITS), PAR if not DRY else 2, MOL_CAP))
done = 0
from concurrent.futures import ProcessPoolExecutor
with ProcessPoolExecutor(max_workers=PAR if not DRY else 2) as ex:
    for gsm, shape in ex.map(score_lib, libs):
        done += 1; log("  %d/%d %s → %s" % (done, len(libs), gsm, shape))
if DRY:
    log("DRY done, no aggregation"); raise SystemExit(0)

# ---------- 4. P2 aggregation (same functions, same basis as M6) ----------
S = pd.concat([pd.read_csv(f, index_col=0) for f in sorted(glob.glob(OUT + "/*.scores.csv.gz"))])
S = S[~S.index.duplicated()]
log("per-nucleus scores pooled %d nuclei / %d libraries" % (len(S), len(set(i.split("|")[0] for i in S.index))))

def d_axis(hi, lo):
    o = {}
    for tr in TRAITS:
        a = S.loc[S.cellType.isin(hi).values, tr]; b = S.loc[S.cellType.isin(lo).values, tr]
        o[tr] = float(a.median() - b.median()) if min(a.notna().sum(), b.notna().sum()) >= 10 else np.nan
    return pd.Series(o)

AX = pd.DataFrame({"d1": d_axis(MOL, OPC_), "d2a": d_axis({"Astro"}, {"Micro"}),
                   "d2b": d_axis({"Neuron"}, OTHER_OF_NEU)})
def rz(s):
    m, mad = s.median(), (s - s.median()).abs().median()
    return (s - m) / (1.4826 * mad)
OK = AX.dropna()
RZ = OK.apply(rz); RANK = OK.rank(ascending=False).astype(int)
import scipy.stats as st
r12 = st.spearmanr(OK.d1, OK.d2a); r13 = st.spearmanr(OK.d1, OK.d2b)
FULL = pd.concat([AX.add_prefix("Δ"), RZ.add_prefix("z_"), RANK.add_prefix("rk_")], axis=1)
FULL.index.name = "trait"; FULL.round(4).to_csv(SUP + "/S7_rna_crossaxis_mouse_full.csv")

NAME = {"d1": "OL maturation (Δ1)", "d2a": "astrocyte − microglia (Δ2A)", "d2b": "neuron − non-neuron (Δ2B)"}
rows = []
for ax in ("d1", "d2a", "d2b"):
    v = AX[ax].dropna()
    rows.append(dict(axis=NAME[ax], n_trait=int(v.size), ms_value=float(v[MS]),
                     ms_z=float(rz(v)[MS]), ms_rank=int((v > v[MS]).sum()) + 1,
                     top3=", ".join(v.sort_values(ascending=False).index[:3])))
SUM = pd.DataFrame(rows)
SUM.to_csv(SUP + "/S7_rna_crossaxis_mouse.csv", index=False)
json.dump(dict(seed=SEED, mol_cap=MOL_CAP, n_nuclei=int(len(S)), n_libs=int(S.gsm.nunique()) if "gsm" in S else len(libs),
               celltype_counts=S.cellType.value_counts().to_dict(),
               rho_d1_d2a=[round(float(r12.statistic), 4), float(r12.pvalue)],
               rho_d1_d2b=[round(float(r13.statistic), 4), float(r13.pvalue)],
               ms_rank={NAME[a]: int(SUM[SUM.axis == NAME[a]].ms_rank.iloc[0]) for a in ("d1", "d2a", "d2b")},
               note="same nuclei/labels/seed as M6; RNA channel = official Gene Expression, CP10K+log1p, scdrs score_cell n_ctrl=200"),
          open(SUP + "/S7_rna_mouse_threeaxis_summary.json", "w"), indent=1, ensure_ascii=False)
log("done: d1 MS z %.3f rank %d | d2A %.3f/%d | d2B %.3f/%d | rho(d1,d2a) %.3f rho(d1,d2b) %.3f" % (
    SUM.ms_z.iloc[0], SUM.ms_rank.iloc[0], SUM.ms_z.iloc[1], SUM.ms_rank.iloc[1],
    SUM.ms_z.iloc[2], SUM.ms_rank.iloc[2], r12.statistic, r13.statistic))
