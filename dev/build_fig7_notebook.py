#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build Figure_7.ipynb — panel-by-panel reproduction of the journal Figure 7
(external-cohort replication of the OL-maturation contrast), organised exactly
like Figure_6.ipynb. Run:  python3 build_fig7_notebook.py
"""
import os
import nbformat as nbf

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "fig7")


def md(s):
    return nbf.v4.new_markdown_cell(s)


def code(s):
    return nbf.v4.new_code_cell(s)


cells = []

# =============================================================================
# Title
# =============================================================================
cells.append(md("""# Figure 7 — External replication of the mature-OL vs OPC scDRS contrast

This notebook reproduces the three panels of **Figure 7** from the analysis objects upstream of it:
per-nucleus scDRS scores of four independent human MS cohorts and of the mouse EAE cohort
GSE250589, plus the 74-trait sweep tables that the panels read. Journal assembly (single-column
sheet, panel letters, final cropping) is done outside the notebook.

| Panel | File | Size | Content |
|---|---|---|---|
| A | `fig7v13_A_donors.pdf` | 7.92 × 4.44 cm | unit-level Δ1 = median(mature OL) − median(OPC) of the MS scDRS score: four human cohorts (donor- or library-level dots) + mouse GSE250589 EAE libraries split by disease stage (Early/Peak/Late); Mann-Whitney *P* above each column |
| B | `fig7v13_B_ranks.pdf` | 7.71 × 5.99 cm (delivered sheet) | 74 GWAS traits × 4 human cohorts: colour = raw rank of cohort-level Δ1 (1 = strongest OL contrast), row order = mean rank after gene-set-size adjustment; rows labelled **MS / AD / SCZ** — the same short forms used in panel C and in Figure 6, full names and citations given in the legend at first use |
| C | `fig7v13_C_threeaxes.pdf` | 7.90 × 5.71 cm | same four cohorts, 2 × 2: Δ on three cell-identity axes (OL maturation / astro−immune / neuron−glia), robust *z* within cohort × axis; MS and SCZ traits marked with their OL-axis rank; `*` on GSE231585 = the neuron−glia axis rests on 717 paired nuclei only |

The deconfounding diagnostics (confound ρ strip, raw→adjusted rank transitions, matched-window
controls) are supplementary (S7b) and are *not* drawn here."""))

# =============================================================================
# Section 0 — environment, paths, versions
# =============================================================================
cells.append(md("""---

## 0. Environment, paths and versions

Three roots, all overridable, default to the author's layout:
`FIG7_ROOT` (per-cohort analysis objects), `FIG7_PROBE` (downloaded cohort matrices),
`FIG7_OUT` (panel PDFs). Nothing in the table below is produced by this notebook; the
`upstream/` directory contains the scripts that produced it (see `upstream/README.md`)."""))

cells.append(code("""import os, glob, time
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr, gaussian_kde
import statsmodels.api as sm

CM = 2.54
ROOT = os.environ.get("FIG7_ROOT", os.path.join(os.environ.get("DSCI_WORK_ROOT", "."), "figures", "routeA_out"))
PROBE = os.environ.get("FIG7_PROBE", os.path.join(os.environ.get("DSCI_WORK_ROOT", "."), "routeA_probe"))
OUT = os.environ.get("FIG7_OUT", "fig7_out")
os.makedirs(OUT, exist_ok=True)
MS = "PASS_Multiple_sclerosis"
SCZ = "PASS_Schizophrenia_Pardinas2018"
AD = "PASS_Alzheimers_Jansen2019"
"""))

cells.append(md("""### 0.1 Where the inputs come from

| cohort | matrix (scale as delivered) | annotation | scDRS run | unit of Δ1 |
|---|---|---|---|---|
| A2 GSE277430 | `GSE277430/OL_subset/OL_csr.bin` — **log1p(CP10K)** dump | `expertAnno.l1` (OLIGO / OPC) + `disease` | `score_cell(n_ctrl=100, random_seed=0)` (dual-axis table) | donor |
| A3 GSE118257 | `GSE118257_matrix.txt.gz` — **raw counts** | `GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz` | `score_cell(n_ctrl=200, random_seed=20261029)` | donor |
| A4 GSE301908 | `GSE301908/OL_csr.bin` — **log1p(CP10K)** dump | `majorCluster` (Oligo / OPC) + `patient_info` | `score_cell(n_ctrl=200, random_seed=0)` | library |
| A7 GSE231585 | `A7_GSE231585_rna_v1/derived/scored.h5ad` — **raw counts** | author-panel cell types | `score_cell(n_ctrl=200, random_seed=20261030)` | sample (5 MS / 3 Ctrl) |
| M GSE250589 | mouse EAE spinal cord, 19 libraries (13 EAE / 6 Ctrl) | author annotation, lesion phase | per-library Δ1 table `M1b_authoranno_v1/tables/per_lib_delta1_author_vs_self.csv` | library, split Early/Peak/Late |

All four human sweeps score the same gene sets: `magma_10kb_1000.74_traits.gs`
(74 MAGMA gene sets, 1000 genes each). `scdrs.pp.preprocess` is applied in every cohort script
before `score_cell`, as run (frozen; see 0.3). The gene-set peak-size covariate
(`n_peak_raw`, used only for the B-panel row order) comes from the authors' own CNP0010263
ARC peak set: deduplicated ATAC peaks within TSS ± 100 kb of the gene-set genes
(`A6_GSE231584_atac_v1/tables/trait_chromatin_sizeMatched.csv`)."""))

cells.append(md("""### 0.2 Software versions"""))

cells.append(code("""import sys, scipy, statsmodels, anndata, matplotlib
print("python", sys.version.split()[0])
for m in (np, pd, scipy, sm, anndata, matplotlib):
    print(m.__name__, m.__version__)
try:                              # scoring itself runs upstream; this notebook only reads the
    import scdrs                  # delivered score tables, so scdrs is optional here
    print("scdrs", scdrs.__version__, "(this kernel)")
except ImportError:
    print("scdrs 1.0.2 - not in this kernel; the upstream sweeps ran in")
    print("  the scoring environment: $DSCI_PY (python 3.9.23, scdrs 1.0.2), see upstream/README.md")"""))

cells.append(md("""### 0.3 Frozen criteria (fixed before the analysis was run, not chosen afterwards)

| | |
|---|---|
| contrast | Δ1 = median(score of mature OL) − median(score of OPC), per unit (donor / library / sample) |
| cohort-level Δ1 | median over units of the per-unit Δ1 (Panel B colour input) |
| rank | descending rank over the 74 traits, `method="min"`; rank 1 = strongest OL contrast |
| gene-set-size adjustment | OLS residual of Δ1 ~ ln(n_peak_raw), re-ranked; reported as the row order (mean adjusted rank over the four cohorts). The raw ranks remain the cell colour. Confound strength per cohort (Spearman ρ of Δ1 with ln n_peak): −0.84 / −0.69 / −0.78 / −0.87 — shown in S7b, not in the main figure |
| matched-window control | traits with n_peak within ±15 % of the MS set (S7b) |
| robust z (Panel C) | (Δ − median) / (1.4826 × MAD), computed within cohort × axis, so the four panels share one y-axis |
| GSE231585 caveat | the neuron−glia axis (Δ2B) is computed on 717 paired nuclei after thinning → marked `*` |
| scDRS | `score_cell(n_ctrl=100 for A2, 200 for A3/A4/A7, random_seed fixed per cohort)`; input scale per cohort as in 0.1 (two cohorts log1p, two raw; `scdrs.pp.preprocess` as run) |"""))

cells.append(code("""def rank_desc(x):
    # Descending rank over traits; 1 = strongest (largest) delta.
    return (-pd.Series(np.asarray(x, float), index=x.index if hasattr(x, "index") else None)
            ).rank(method="min").astype(int).values

def ols_resid(y, cov):
    # Residual of y on [1, cov] - removes the smooth gene-set-size effect.
    X = sm.add_constant(np.asarray(cov, float))
    return np.asarray(sm.OLS(np.asarray(y, float), X).fit().resid)

def robust_z(v):
    v = np.asarray(v, float); med = np.nanmedian(v); mad = 1.4826 * np.nanmedian(np.abs(v - med))
    return (v - med) / (mad if mad > 1e-12 else 1.0)

def mwu_p(x, y):
    return float(mannwhitneyu(x, y, alternative="two-sided").pvalue)

# Arial: register the shipped TTF when the OS font is absent (cosmetic only; the
# delivered figure is assembled with the same font). Overridable via ARIAL_TTF.
_A = next((p for p in [os.environ.get("ARIAL_TTF"),
                   os.path.join(os.environ.get("DSCI_WORK_ROOT", "."), "fonts", "ARIAL.TTF"),
                   "fonts/ARIAL.TTF",
                   "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"]
                 if p and os.path.exists(p)), "")
if os.path.exists(_A):
    import matplotlib.font_manager as fm
    fm.fontManager.addfont(_A)
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "axes.linewidth": 0.6})"""))

# =============================================================================
# Panel A
# =============================================================================
cells.append(md("""---

## 1. Panel A — unit-level Δ1 across four human cohorts and mouse EAE stages

**Data chain**: per-nucleus MS scDRS scores + cell-type/donor annotation (A2/A4), or the
delivered donor-level tables (A3/A7/M). For every unit we take
`median(MS score | mature OL) − median(MS score | OPC)`; the panel plots one column per cohort
(human) or lesion stage (mouse), with a Mann-Whitney *P* per column (human: MS vs Ctrl donors;
mouse: stage vs Ctrl libraries)."""))

cells.append(code("""t0 = time.time()
rows = []
# A2: per-nucleus table -> donor medians
a2 = pd.read_csv(f"{ROOT}/A2_GSE277430_v1/tables/OL_meta_with_MS_score.csv",
                 usecols=["expertAnno.l1", "patient", "disease", "MS_norm_score"])
dm = a2[a2["expertAnno.l1"] == "OLIGO"].groupby("patient")["MS_norm_score"].median()
dp = a2[a2["expertAnno.l1"] == "OPC"].groupby("patient")["MS_norm_score"].median()
dis = a2.groupby("patient")["disease"].first()
ix = dm.index.intersection(dp.index)
for p, v in (dm[ix] - dp[ix]).items():
    rows.append(dict(study="A2 GSE277430", unit="donor", id=str(p), cond=str(dis[p]), d1=float(v)))
# A3: delivered donor-level pairs
a3 = pd.read_csv(f"{ROOT}/A3_GSE118257_v1/tables/patient_pair_scores.csv")
for _, r in a3.iterrows():
    rows.append(dict(study="A3 GSE118257", unit="donor", id=str(r.patient), cond=str(r.Condition), d1=float(r["diff"])))
# A4: per-nucleus table -> library medians (one patient per library)
a4 = pd.read_csv(f"{ROOT}/A4_GSE301908_v1/tables/OL_meta_with_MS_score.csv")
a4["patient"] = a4["orig.ident"].astype(str)
dm4 = a4[a4.majorCluster == "Oligo"].groupby("patient")["MS_norm_score"].median()
dp4 = a4[a4.majorCluster == "OPC"].groupby("patient")["MS_norm_score"].median()
st4 = a4.groupby("patient")["patient_info"].first()
ix4 = dm4.index.intersection(dp4.index)
for p, v in (dm4[ix4] - dp4[ix4]).items():
    rows.append(dict(study="A4 GSE301908", unit="library", id=str(p),
                     cond={"patient": "MS", "ctrl": "Ctrl"}[str(st4[p])], d1=float(v)))
# A7: author-panel donor table (donor id encodes condition: MS* vs others)
a7 = pd.read_csv(f"{ROOT}/A7_GSE231585_rna_v1/tables/delta1_MS_authorpanel_donors.csv")
for _, r in a7.iterrows():
    rows.append(dict(study="A7 GSE231585", unit="donor", id=str(r.donor),
                     cond="MS" if str(r.donor).startswith("MS") else "Ctrl", d1=float(r.d1_new)))
# M: mouse GSE250589, per-library delta1 of the MS gene set (author annotation)
mm = pd.read_csv(f"{ROOT}/M1b_authoranno_v1/tables/per_lib_delta1_author_vs_self.csv")
mm = mm[(mm.trait == MS) & (mm.anno == "author")]
for _, r in mm.iterrows():
    rows.append(dict(study="M GSE250589", unit="library", id=str(r.lib),
                     cond="Ctrl" if r.phase == "Ctrl" else "EAE", d1=float(r.delta1), phase=r.phase))
A = pd.DataFrame(rows)
A["cond"] = A.cond.replace({"MS_Mye": "MS", "MS_Demye": "MS", "CTRL": "Ctrl", "CONTROL": "Ctrl"})
print("A: %d units | %s" % (len(A), A.groupby(["study", "cond"]).size().to_dict()), "%.0fs" % (time.time() - t0))"""))

cells.append(code("""# between-unit P per human cohort (recomputed here; identical to the delivered table)
PP = {}
for st, g in A.groupby("study"):
    x, y = g[g.cond.isin(["MS", "EAE"])].d1.values, g[g.cond == "Ctrl"].d1.values
    PP[st] = mwu_p(x, y)
print({k: round(v, 4) for k, v in PP.items()})

COLC = {"MS": "#D51F26", "Ctrl": "#272E6A", "EAE": "#F28E2B"}
STUDIES = ["A2 GSE277430", "A3 GSE118257", "A4 GSE301908", "A7 GSE231585"]
TPS = ["Early", "Peak", "Late"]
_mg = A[A.study == "M GSE250589"]          # library ids carry the stage, e.g. Early_Ctrl_r3 / Late_EAE_r2
COLG = [A[A.study == s] for s in STUDIES] + [_mg[_mg.id.str.startswith(t)] for t in TPS]
HEADS = [s.split(" ")[1] for s in STUDIES] + TPS

rng = np.random.default_rng(0)
fig, ax = plt.subplots(figsize=(7.82 / CM, 5.40 / CM))
fig.subplots_adjust(left=0.13, right=0.985, top=0.86, bottom=0.20)
for i, g in enumerate(COLG):
    for cond, sub in g.groupby("cond"):
        ax.scatter(i + rng.uniform(-0.16, 0.16, len(sub)), sub.d1, s=15,
                   color=COLC.get(cond, "0.4"), alpha=0.9, edgecolor="white", lw=0.35, zorder=4)
    for cond in ["Ctrl", "MS", "EAE"]:
        s_ = g[g.cond == cond].d1.values
        if len(s_):
            xm = i + (0.22 if cond == "Ctrl" else -0.22)
            ax.plot([xm - 0.13, xm + 0.13], [np.median(s_)] * 2, color=COLC[cond], lw=1.6, zorder=5)
    pv = PP[STUDIES[i]] if i < 4 else mwu_p(g[g.cond == "EAE"].d1.values, g[g.cond == "Ctrl"].d1.values)
    ax.text(i, 2.62, "P=%.3f" % pv, fontsize=6.0, va="top", ha="center", color="0.15", zorder=6)
ax.axhline(0, color="0.42", ls=":", lw=0.7)
ax.set_xlim(-0.5, 6.5); ax.set_ylim(-0.45, 3.42); ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
ax.set_xticks(range(7)); ax.set_xticklabels([""] * 7)
for i, hd in enumerate(HEADS):
    ax.text(i, -0.62, hd.replace("GSE", "GSE\\n"), fontsize=6.0, ha="center", va="top", linespacing=1.32)
ax.set_ylabel("mature OL − OPC\\n(median score)", fontsize=8.0, linespacing=1.35)
ax.tick_params(labelsize=8.0, width=0.5, length=2.5)
handles = [plt.Line2D([], [], marker="o", ls="", ms=5, color=c, label=l)
           for l, c in [("MS", COLC["MS"]), ("EAE", COLC["EAE"]), ("Ctrl", COLC["Ctrl"])]]
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.155), ncol=3,
          fontsize=6.0, frameon=False)
ax.axvspan(-0.5, 3.5, color="#274B6D", alpha=0.06)
ax.axvspan(3.5, 6.5, color="#B2691E", alpha=0.07)
ax.text(1.5, 3.22, "Human MS (brain)", fontsize=6.0, fontweight="bold", ha="center", color="#1F3E6E")
ax.text(5.0, 3.22, "Mouse EAE (spinal cord)", fontsize=6.0, fontweight="bold", ha="center", color="#8A4B12")
# synced with the delivered assembly (author 2026-11-05): below the light-orange brace, write the shared key GSE250589 of the three mouse columns,
# otherwise a reader sees "Mouse EAE (spinal cord)" without knowing which accession.
ax.text(5.0, 3.00, "GSE250589", fontsize=6.0, ha="center", va="top", color="#8A4B12", zorder=6)
fig.savefig(f"{OUT}/fig7v13_A_donors.pdf", bbox_inches="tight", pad_inches=0.03)
plt.close(fig)
import fitz
d = fitz.open(f"{OUT}/fig7v13_A_donors.pdf")[0]
print("fig7v13_A_donors.pdf  %.2f x %.2f cm" % (d.rect.width / 28.3465, d.rect.height / 28.3465))"""))

# =============================================================================
# Panel B
# =============================================================================
cells.append(md("""---

## 2. Panel B — 74 GWAS traits × 4 cohorts: rank heatmap of Δ1

**Data chain**: the four cohort sweep tables (one Δ1 per trait per cohort, computed by the
`upstream/` sweep scripts with the frozen scDRS parameters). Here we rebuild the Δ1 matrix,
its raw ranks (cell colour) and the gene-set-size-adjusted ranks (row order) — the same
quantities as the delivered `seven4_confound_ALL_cohorts.csv`.

Why adjust: across all four cohorts Δ1 correlates strongly with the size of the gene set's
regulatory domain (Spearman ρ of Δ1 with ln n_peak: −0.84 / −0.69 / −0.78 / −0.87, shown in S7b).
A larger peak domain inflates Δ1 mechanically, so the *raw* rank mixes biology with gene-set
size. Residualising Δ1 on ln n_peak (OLS) and re-ranking gives the adjusted rank; the row order
of the heatmap is the mean adjusted rank over the four cohorts. MS moves to adjusted rank
1 / 1 / 3 / 1 (mean 1.5)."""))

cells.append(code("""TILES = {  # cohort -> (sweep table, delta1 column)
    "A2": ("A2_GSE277430_v1/tables/trait_sweep_dual_axis.csv", "d1_OLmat"),
    "A3": ("A3_GSE118257_v1/tables/trait_sweep_delta2_GSE118257_74.csv", "d1_OLIGO_minus_OPC"),
    "A4": ("A4_GSE301908_v1/tables/trait_sweep_delta2_GSE301908_74.csv", "d1"),
    "A7": ("A7_GSE231585_rna_v1/tables/trait_sweep_delta2_GSE231585_74.csv", "d1"),
}
D = {}
for k, (f, c) in TILES.items():
    t = pd.read_csv(f"{ROOT}/{f}")
    D[k] = t.set_index("trait")[c].dropna()
M = pd.DataFrame(D)                      # 74 traits x 4 cohorts, cohort-level delta1
assert M.shape == (74, 4), M.shape
# gene-set peak size (own CNP0010263 ARC reference, TSS +/- 100 kb, deduplicated)
npeak = pd.read_csv(f"{ROOT}/A6_GSE231584_atac_v1/tables/trait_chromatin_sizeMatched.csv") \\
        .set_index("trait").n_peak_raw.reindex(M.index).astype(float)
lnp = np.log(npeak.values)
R_raw = pd.DataFrame({k: rank_desc(M[k].values) for k in M.columns}, index=M.index)
R_adj = pd.DataFrame({k: rank_desc(ols_resid(M[k].values, lnp)) for k in M.columns}, index=M.index)
order = R_adj.mean(axis=1).sort_values(ascending=False).index   # best (rank ~1) on top
print("MS raw ranks:", list(R_raw.loc[MS]), "| adjusted:", list(R_adj.loc[MS]),
      "| mean adj %.2f" % R_adj.loc[MS].mean())
print("confound rho (ln n_peak vs delta1):",
      {k: round(spearmanr(M[k], lnp).statistic, 3) for k in M.columns})"""))

cells.append(code("""Mr = R_raw.loc[order].values                    # rows in display order
fig, ax = plt.subplots(figsize=(8.60 / CM, 5.50 / CM))
# author 2026-11-05: with row names shortened (MS/AD/SCZ) the left margin narrows in sync, 0.215 -> 0.165, matching the delivered assembly script
fig.subplots_adjust(left=0.165, right=0.815, top=0.97, bottom=0.235)
im = ax.imshow(Mr, aspect="auto", cmap="viridis_r", vmin=0.5, vmax=74.5)
ax.set_xticks(range(4))
ax.set_xticklabels(["GSE277430", "GSE118257", "GSE301908", "GSE231585"],
                   fontsize=8.0, rotation=90, ha="right")
ax.tick_params(axis="x", pad=1.5)
for _t in ax.get_xticklabels():
    _t.set_fontweight("normal")
ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)
MARK = [(MS, "MS", "#D51F26"), (AD, "AD", "#2E7D32"), (SCZ, "SCZ", "#272E6A")]   # short forms, as in panel C / Fig 6
for t_, nm_, col_ in MARK:
    ry = list(order).index(t_)
    ax.text(-0.018, ry, nm_, transform=ax.get_yaxis_transform(), ha="right", va="center",
            fontsize=6.0, color=col_, zorder=6)
    ax.plot([-0.004, 0.028], [ry, ry], transform=ax.get_yaxis_transform(),
            color=col_, lw=0.9, clip_on=False, zorder=6)
cax = fig.add_axes([0.875, 0.335, 0.014, 0.44])
cb = fig.colorbar(im, cax=cax, orientation="vertical", ticks=[1, 20, 40, 60, 74])
cb.outline.set_linewidth(0.4)
cb.ax.tick_params(labelsize=8.0, length=1.6, pad=1.5)
fig.text(0.882, 0.815, "rank", fontsize=6.0, ha="center", va="bottom")
fig.savefig(f"{OUT}/fig7v13_B_ranks.pdf", bbox_inches="tight", pad_inches=0.03)
plt.close(fig)
d = fitz.open(f"{OUT}/fig7v13_B_ranks.pdf")[0]
print("fig7v13_B_ranks.pdf  %.2f x %.2f cm" % (d.rect.width / 28.3465, d.rect.height / 28.3465))"""))

# =============================================================================
# Panel C
# =============================================================================
cells.append(md("""---

## 3. Panel C — the same cells, three cell-identity axes

**Data chain**: the same four sweep tables, which also carry Δ on two further axes computed
from the same per-nucleus scDRS MS scores: astro−immune (median AST − median IMM) and
neuron−glia (median NEU − median GLIA). Within each cohort × axis the 74 trait-level Δ values
are converted to robust *z* (median/MAD), so the four panels share one y-axis; the OL-maturation
axis is Δ1. Dots mark MS (red) and SCZ (blue) with their rank on the plotted axis
(`#rank` on the OL axis = the same ranks coloured in Panel B).

Caveat printed as `*`: in GSE231585 the neuron−glia axis rests on 717 paired nuclei
(after thinning to equal neuron/glia depth); the other two axes use the full sample."""))

cells.append(code("""AX = {  # cohort -> (d1, dA, dB) column names in the sweep table
    "A2": ("d1_OLmat", "d2_AST_IM", "d2_NEU_GL"),
    "A3": ("d1_OLIGO_minus_OPC", "dA_AST_minus_IMM", "dB_NEU_minus_GLIA"),
    "A4": ("d1", "dA", "dB"),
    "A7": ("d1", "dA", "dB"),
}
GSE = {"A2": "GSE277430", "A3": "GSE118257", "A4": "GSE301908", "A7": "GSE231585"}
AXCOL = ["#4E7A3A", "#A6651E", "#5B5688"]
AXL = ["OL maturation", "astro − immune", "neuron − glia"]
YLO, YHI = -4.7, 5.5

def violin_clipped(ax, v, x, half, col):
    # KDE violins sampled only inside the visible window (keeps PDF paths in-axis).
    g0, g1 = max(v.min(), YLO + 0.06), min(v.max(), YHI - 0.06)
    if g1 <= g0:
        return
    grid = np.linspace(g0, g1, 220)
    den = gaussian_kde(v, bw_method=0.35)(grid)
    den = den / den.max() * half
    ax.fill_betweenx(grid, x - den, x + den, facecolor=col, edgecolor=col, lw=0.5,
                     alpha=0.30, zorder=3)

fig, axes = plt.subplots(2, 2, figsize=(7.90 / CM, 5.40 / CM))
fig.subplots_adjust(left=0.105, right=0.985, top=0.905, bottom=0.215, wspace=0.13, hspace=0.34)
for k, (co, (c1, cA, cB)) in enumerate(AX.items()):
    ax = axes[k // 2][k % 2]
    t = pd.read_csv(f"{ROOT}/{TILES[co][0]}")
    Z = pd.DataFrame({c: robust_z(t[c].values) for c in (c1, cA, cB)}, index=t.trait)
    R = pd.DataFrame({c: rank_desc(t[c].values) for c in (c1, cA, cB)}, index=t.trait)
    for j, c in enumerate((c1, cA, cB)):
        v = Z.loc[Z[c].between(YLO, YHI), c].values
        violin_clipped(ax, v, j, 0.16, AXCOL[j])
        q = np.nanpercentile(v, [25, 50, 75])
        ax.plot([j - 0.05, j + 0.05], [q[1], q[1]], color="black", lw=0.8, zorder=6)
        ax.add_patch(plt.Rectangle((j - 0.05, q[0]), 0.10, q[2] - q[0], fill=False,
                                   edgecolor="black", lw=0.6, zorder=6))
    for tr, col in [(MS, "#D51F26"), (SCZ, "#272E6A")]:
        if tr not in Z.index:
            continue
        for j, c in enumerate((c1, cA, cB)):
            if not np.isnan(Z.loc[tr, c]):
                ax.scatter(j + 0.10, Z.loc[tr, c], s=16, color=col, zorder=7,
                           edgecolor="white", lw=0.3)
                ax.text(j + 0.17, Z.loc[tr, c], "#%d" % R.loc[tr, c], fontsize=6.0,
                        va="center", color=col)
    ax.axhline(0, color="0.42", lw=0.5, ls=":")
    ax.set_xlim(-0.5, 2.5); ax.set_ylim(YLO, YHI)
    ax.set_xticks(range(3))
    for j in range(3):
        ax.text(j, YLO - 0.22, AXL[j], ha="right", va="top", fontsize=6.0,
                rotation=60, color="black")
    ax.set_xticklabels([])
    ax.set_title(GSE[co] + (" *" if co == "A7" else ""), fontsize=8.0, loc="left")
    ax.tick_params(labelsize=8.0, width=0.5, length=2.5)
    if k % 2 == 0:
        ax.set_ylabel("robust z (Δ)", fontsize=8.0)
    else:
        ax.set_yticklabels([])
h = [plt.Line2D([], [], marker="o", ls="", ms=5, color="#D51F26", label="MS"),
     plt.Line2D([], [], marker="o", ls="", ms=5, color="#272E6A", label="SCZ")]
axes[0][1].legend(handles=h, loc="upper right", bbox_to_anchor=(0.985, 0.945), ncol=2,
                  fontsize=6.0, frameon=False)
fig.savefig(f"{OUT}/fig7v13_C_threeaxes.pdf", bbox_inches="tight", pad_inches=0.03)
plt.close(fig)
d = fitz.open(f"{OUT}/fig7v13_C_threeaxes.pdf")[0]
print("fig7v13_C_threeaxes.pdf  %.2f x %.2f cm" % (d.rect.width / 28.3465, d.rect.height / 28.3465))"""))

nb = nbf.v4.new_notebook(metadata={
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.8"}})
nb["cells"] = cells
out = os.path.join(OUT, "Figure_7.ipynb")
nbf.write(nb, out)
print("wrote", out, len(cells), "cells")
