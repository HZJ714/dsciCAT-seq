# -*- coding: utf-8 -*-
"""mouse_atac_threeaxis.py — mouse ATAC three-axis run on GSE250589 (25 libraries, 14,928 nuclei); same nuclei, cell-type labels, gene sets and seed as the mouse RNA run.

=========================================================================================
Structural changes against v16 (the numeric basis is unchanged, everything live-computed + asserted):
 1) **no titles beside the panel letters** (all four v16 row-head sentences deleted, meaning to the legend);
 2) **no explanatory sentences on the face** ("raw = rank before stratification...", "staircase = size drives
    raw ranks", "15 nearest neighbours...", "rows: mean adjusted d1 rank...", the rho(d1,d2A)/rho(d1,d2B)
    heads all deleted and moved to the legend);
 3) **no (i)(ii)(iii)**: consecutive panel letters A-J, explained letter by letter in the legend;
 4) **de-duplicated**: v16 a2 deleted (stratified rank heatmap = the d1 column of the new ledger J), the d2
    violins deleted (= duplicate of main Fig 7C), the d4 text ledger deleted (goes to Table S7 / legend);
 5) **low-information chart types demoted to small tables**: the c1 four-point "obs vs size-predicted"
    (forest-like) merged into two columns of the C table; the c4 three-point "11 trait-level rho pairs"
    became the G small table (median [range]);
 6) **figS7b merged in**: new panel J = three-axis ledger (74 traits x 6 groups x d1/d2A/d2B), groups =
    four human RNA cohorts + **human ATAC (GSE231584)** + **mouse ATAC (GSE250589, atlas level, new M6
    run)**; axis titles written once (under the first group), vertical rules between groups;
 7) **one shared legend only**: the rank heatmaps (A / H / J) share one viridis_r 1-74 scale, drawn once;
    the tight pool (B) keeps its own 1-16 bar;

While M6 has no results, M6_WAIT=1 checks the layout only (the three mouse ATAC columns stay grey); a
formal figure must wait for the table. Output: SUPPORT_v1/pdf/figS7_v17_M6.pdf; face-number sensor table tables/figS7_v17_M6_face_numbers.csv"""
import sys, os

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
import os
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
os.environ.setdefault("TQDM_DISABLE", "1")
sys.path.insert(0, WORK + "/scripts/routeA")
from support_common import *
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.stats import spearmanr, rankdata

GRN, RED, BLU = "#0B8457", "#D51F26", "#1F77B4"
COH = ["A2", "A3", "A4", "A7"]
GSE = {"A2": "GSE277430", "A3": "GSE118257", "A4": "GSE301908", "A7": "GSE231585"}
import s7b_panel as SP
MS, AD, SCZ = SP.MS, SP.AD, SP.SCZ
PL = "UKB_460K.blood_PLATELET_COUNT"
LET, FS = 8.0, 6.0                     # two tiers: 8.0 only for panel letters, everything else 6.0
W, H = 17.0, 16.6
CMv = 2.54
NPOOL = 15
M6_WAIT = os.environ.get("M6_WAIT") == "1"

def ax_cm(x0, y0, w, h):
    return fig.add_axes([x0 / W, y0 / H, w / W, h / H])

def letter(x, y, ch):
    fig.text(x / W, y / H, ch, fontsize=LET, ha="left", va="bottom", color="0")

# ---------------------------------------------------------------- data sources
S4 = read("FIG7_v1/tables/seven4_confound_ALL_cohorts.csv")
MW = read("SUPPORT_v1/tables/fig7Sb_matched_window.csv")
CP = read("SUPPORT_v1/tables/figS7c_crosschannel_pairs.csv")
PK = read("A6_GSE231584_atac_v1/tables/delta1_vs_peakdomain_allchannels.csv")
T4 = read("O1_own_multiome_paired/tables/T4_delta1_74traits_both_channels.csv")
MO = read("M2_mouse_multiome_paired/tables/cross_channel_delta1_same_nuclei.csv")
HUM = S4.merge(PK[["trait", "median_groupdiff"]], on="trait")
SD = read("A6_GSE231584_atac_v1/tables/same_donor_RNA_vs_CHROMATIN_delta1_74traits.csv")
P3 = pd.read_csv(WORK + "/routeA_probe/m3_probe/tables_full74/M3_paired3_atac_vs_rna_74.csv")
FS_ = pd.read_csv("%s/SUPPORT_v1/tables/S7_atac_rescore_fullstats.csv" % R)[["trait", "atac_d1_fullstats"]]
P3 = P3.merge(FS_, on="trait", how="left")
P3["atac_d1_cli"] = P3["atac_d1"]; P3["atac_d1"] = P3["atac_d1_fullstats"]
P3["rna_d1_m1"] = P3["rna_d1_3mouse"]
assert P3.rna_d1_m1.notna().all() and P3.atac_d1.notna().all()

S = S4.sort_values(["n_peak_raw", "trait"], ascending=[False, True]).reset_index(drop=True)
assert len(S) == 74 and S.n_peak_raw.nunique() == 73
iMS, iPL = int(S.index[S.trait == MS][0]), int(S.index[S.trait == PL][0])
assert (iMS + 1, iPL + 1) == (37, 56)
LO, HI = float(MW.window_lo.iloc[0]), float(MW.window_hi.iloc[0])
inw = ((S.n_peak_raw >= LO) & (S.n_peak_raw <= HI)).values
w0, w1 = int(np.where(inw)[0][0]), int(np.where(inw)[0][-1])
assert inw[w0:w1 + 1].all() and w1 - w0 + 1 == int(MW.n_in_window.iloc[0]) == 40

def rz(v):
    v = np.asarray(v, float); m = np.median(v)
    return (v - m) / (1.4826 * np.median(np.abs(v - m)))

RHO_RAW = {c: float(MW[MW.cohort == c].rho_lnp_vs_d1.iloc[0]) for c in COH}
RHO_ADJ = {c: float(spearmanr(S.lnp.values, rz(S[c].values) - np.polyval(np.polyfit(S.lnp, rz(S[c]), 1), S.lnp)).statistic) for c in COH}
assert [round(abs(RHO_RAW[c]), 2) for c in COH] == [0.84, 0.69, 0.78, 0.87]
OBS, EXP = {}, {}
for c in COH:
    z = rz(S[c].values); b, a0 = np.polyfit(S.lnp.values, z, 1)
    OBS[c] = z[iPL] - z[iMS]; EXP[c] = b * (S.lnp.values[iPL] - S.lnp.values[iMS])

# tight pool
lnp = S.lnp.values
d = np.abs(lnp - lnp[iMS]); d[iMS] = np.inf
POOL = np.sort(np.append(np.sort(np.argsort(d)[:NPOOL]), iMS))
span_pool = float(np.exp(lnp[POOL].max() - lnp[POOL].min()))
assert len(POOL) == 16 and iPL not in set(POOL.tolist()) and span_pool < 1.20
POOLRANK = np.full((74, 4), np.nan); PR = {}
for k, c in enumerate(COH):
    r = rankdata(-S[c].values[POOL].astype(float), method="min")
    PR[c] = int(r[list(POOL).index(iMS)]); POOLRANK[POOL, k] = r
assert [PR[c] for c in COH] == [1, 1, 1, 1]
idxw = np.where(inw)[0]; LOC = {}
for c in COH:
    yv = S[c].values.astype(float); b1, b0 = np.polyfit(lnp[inw], yv[inw], 1)
    res = -(yv - (b1 * lnp + b0))
    rr = rankdata(res[idxw], method="min")
    LOC[c] = (int(rr[list(idxw).index(iMS)]), int(rr[list(idxw).index(iPL)]))
assert [LOC[c][0] for c in COH] == [1, 2, 4, 1] and [LOC[c][1] for c in COH] == [2, 1, 3, 2]

def colof(trait, key, c):
    return int(S4[S4.trait == trait].iloc[0]["%s_%s" % (key, c)])
TAB = {}
for tr, nm in [(MS, "MS"), (PL, "PLATELET"), (AD, "AD"), (SCZ, "SCZ")]:
    TAB[nm] = {"raw": [colof(tr, "rank_raw", c) for c in COH], "adj": [colof(tr, "rank_adj", c) for c in COH],
               "mat": ([LOC[c][0] for c in COH] if nm == "MS" else
                       [LOC[c][1] for c in COH] if nm == "PLATELET" else [None] * 4)}
assert TAB["MS"]["raw"] == [2, 4, 7, 2] and TAB["MS"]["adj"] == [1, 1, 3, 1] and TAB["MS"]["mat"] == [1, 2, 4, 1]
assert TAB["PLATELET"]["raw"] == [1, 2, 1, 1] and TAB["PLATELET"]["adj"] == [11, 5, 8, 7] and TAB["PLATELET"]["mat"] == [2, 1, 3, 2]
assert TAB["AD"]["raw"] == [25, 41, 26, 9] and TAB["AD"]["adj"] == [18, 49, 27, 12]
assert TAB["SCZ"]["raw"] == [67, 73, 63, 67] and TAB["SCZ"]["adj"] == [68, 74, 42, 71]

# derived quantities for the D row (now the H/I panels)
rank_ms = lambda df, col: int((df[col] > df[df.trait == MS].iloc[0][col]).sum()) + 1
D = {}
D["human_RNA"], D["human_CHR"] = rank_ms(SD, "RNA_delta1"), rank_ms(SD, "CHR_delta1")
D["mouse_RNA"], D["mouse_ATAC"] = rank_ms(P3, "rna_d1_m1"), rank_ms(P3, "atac_d1")
assert (D["human_RNA"], D["human_CHR"]) == (2, 21) and (D["mouse_RNA"], D["mouse_ATAC"]) == (2, 28)
D["rho_mouse74"] = float(spearmanr(P3.rna_d1_m1, P3.atac_d1, nan_policy="omit").statistic)
D["P_mouse74"] = float(spearmanr(P3.rna_d1_m1, P3.atac_d1, nan_policy="omit").pvalue)
r_lnp_m = spearmanr(P3.atac_d1, P3.ln_n_peak, nan_policy="omit").statistic
assert abs(r_lnp_m - 0.050) < 0.02

# ---------------------------------------------------------------- panel J: three-axis ledger (both ATAC channels)
RM = read("SUPPORT_v1/tables/fig7Sc_rank_matrix_4cohorts.csv").set_index("trait")
AXRHO = read("SUPPORT_v1/tables/fig7Sc_axis_pair_rho.csv")
CAXF = pd.read_csv(R + "/SUPPORT_v1/tables/S7_atac_crossaxis_human_full.csv", index_col=0, header=[0, 1])
CAXF.index = [i if isinstance(i, str) else i for i in CAXF.index]
GRP_J = []                                            # (column-group name, [three axis-rank Series], (rho12, rho13))
for k, c in enumerate(COH):
    q = AXRHO[AXRHO.cohort == c]
    GRP_J.append((GSE[c], [RM["c%d_1" % k], RM["c%d_A" % k], RM["c%d_B" % k]],
                  (float(q[q.pair == "d1 vs d2A"].rho.iloc[0]), float(q[q.pair == "d1 vs d2B"].rho.iloc[0]))))
hATAC = {}
for ax in ("d1", "d2a", "d2b"):
    v = CAXF[("\u0394", ax)]
    hATAC[ax] = v.rank(ascending=False, method="min")
GRP_J.append(("GSE231584", [hATAC["d1"], hATAC["d2a"], hATAC["d2b"]], None))
MOUSE_J = None
MF = "%s/SUPPORT_v1/tables/S7_atac_crossaxis_mouse_full.csv" % R
if os.path.exists(MF):
    MX = pd.read_csv(MF, index_col=0)
    _r12 = float(MX["\u0394d1"].corr(MX["\u0394d2a"], method="spearman"))
    _r13 = float(MX["\u0394d1"].corr(MX["\u0394d2b"], method="spearman"))
    MOUSE_J = [("GSE250589", [MX["\u0394d1"].rank(ascending=False, method="min"),
                              MX["\u0394d2a"].rank(ascending=False, method="min"),
                              MX["\u0394d2b"].rank(ascending=False, method="min")], (_r12, _r13))]
elif M6_WAIT:
    say("M6_WAIT=1: the three mouse ATAC axis columns stay empty (layout check only, no formal figure)")
else:
    raise SystemExit("S7_atac_crossaxis_mouse_full.csv missing (M6 not finished), no formal figure")

if MOUSE_J:
    GRP_J += MOUSE_J
rho_hu = [float(CAXF[("\u0394", "d1")].corr(CAXF[("\u0394", a)], method="spearman")) for a in ("d2a", "d2b")]
GRP_J[4] = (GRP_J[4][0], GRP_J[4][1], (rho_hu[0], rho_hu[1]))      # human ATAC group gains rho(d1,d2A)/rho(d1,d2B)
# row order = mean adjusted d1 rank over the four cohorts (same key as Fig 7B)
keyadj = RM[[c for c in RM.columns if c.endswith("_1")]].mean(axis=1)
JROWS = keyadj.sort_values().index.tolist()
assert len(JROWS) == 74
ZJ = np.column_stack([rms[si].reindex(JROWS).values.astype(float) for _, rms, _ in GRP_J for si in range(3)])
assert ZJ.shape == (74, 3 * len(GRP_J)), ZJ.shape
MS_J = JROWS.index(MS) + 1
RANK_HU = {ax: int(hATAC[ax].loc[MS]) for ax in ("d1", "d2a", "d2b")}
assert (RANK_HU["d1"], RANK_HU["d2a"], RANK_HU["d2b"]) == (31, 65, 43), RANK_HU
rho21, rho22 = rho_hu
assert abs(rho_hu[0] + 0.659) < 0.01 and abs(rho_hu[1] + 0.681) < 0.01

# ---------------------------------------------------------------- canvas
fig = plt.figure(figsize=(W / CMv, H / CMv))
plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.weight"] = "normal"
plt.rcParams["text.color"] = "0"; plt.rcParams["axes.labelcolor"] = "0"
plt.rcParams["xtick.color"] = "0"; plt.rcParams["ytick.color"] = "0"

cmap_rank = plt.get_cmap("viridis_r")
im_shared = None

def ladder(ax, key, rho, with_rho=True):
    global im_shared
    Z = np.array([S["%s_%s" % (key, c)].values.astype(float) for c in COH]).T
    im = ax.imshow(Z, aspect="auto", cmap=cmap_rank, vmin=1, vmax=74)
    ax.set_xlim(-0.5, 3.5); ax.set_ylim(73.5, -0.5)
    ax.set_xticks([]); ax.set_yticks([]); ax.tick_params(length=0)
    ax.add_patch(plt.Rectangle((-0.5, w0 - 0.5), 4.0, w1 - w0 + 1, facecolor="none",
                               edgecolor="0", lw=0.6, ls=(0, (2, 1.5)), zorder=6))
    if with_rho:
        for k, c in enumerate(COH):
            ax.text(k, -2.9, "%+.2f" % rho[c], ha="center", va="center", fontsize=FS,
                    color="#8B1A1A" if abs(rho[c]) > 0.5 else "0", clip_on=False)
    im_shared = im_shared or im
    return im

# ================================================================ row 1: A rank ladder / B tight pool / C exact rank table
# four-row grid (top to bottom): row1 12.20-15.60 | row2 8.50-11.40 | row3 5.00-7.70 | row4 1.00-4.20, letters at panel top-left +0.12
yA, hA = 12.20, 3.40
letter(1.80, yA + hA + 0.12, "A")
a1 = ax_cm(1.95, yA, 2.20, hA)
im1 = ladder(a1, "rank_raw", RHO_RAW)
letter(4.70, yA + hA + 0.12, "B")
cmap_pool = ListedColormap(cmap_rank(np.linspace(0, 1, 17))[:16]); cmap_pool.set_bad("0.88")
a3 = ax_cm(4.85, yA, 2.20, hA)
im3 = a3.imshow(np.ma.masked_invalid(POOLRANK), aspect="auto", cmap=cmap_pool, vmin=1, vmax=16)
a3.set_xlim(-0.5, 3.5); a3.set_ylim(73.5, -0.5); a3.set_xticks([]); a3.set_yticks([]); a3.tick_params(length=0)
a3.set_facecolor("0.88")
for k, c in enumerate(COH):
    a3.text(k, -2.9, "1/16" if PR[c] == 1 else "%d/16" % PR[c], ha="center", va="center",
            fontsize=FS, color=RED, clip_on=False)
for axx in (a1, a3):
    axx.axhline(iMS, color=RED, lw=0.9, zorder=7); axx.axhline(iPL, color="0", lw=0.9, zorder=7)
fig.text(1.88 / W, (yA + hA * (73.5 - iMS) / 74.0) / H, "MS", fontsize=FS, color=RED, ha="right", va="center")
fig.text(1.88 / W, (yA + hA * (73.5 - iPL) / 74.0) / H, "PLATELET", fontsize=FS, color="0", ha="right", va="center")
cbA = fig.colorbar(im1, cax=ax_cm(7.35, yA + hA - 2.95, 0.22, 2.95), orientation="vertical")
cbA.set_ticks([1, 20, 40, 60, 74]); cbA.ax.tick_params(labelsize=FS, length=1.6)
cbA.outline.set_linewidth(0.4); cbA.ax.set_title("rank", fontsize=FS, pad=2)
# tight-pool colour bar deleted: every rank-type panel (A/H/J) shares one rank scale; the within-pool rank is already printed on the face as "1/16"

# ---- C: exact rank table (two-line cohort heads to avoid mutual pressure) + the two mechanism rows as a separate small table
letter(8.75, yA + hA + 0.12, "C")
x_lab, x_num0 = 9.45, 9.70
cw, gapg = 0.40, 0.34
y_grp, y_key, y_row0, rh = 15.50, 14.98, 14.60, 0.44
cols = []
for gi, c in enumerate(COH):
    for si, key in enumerate(("raw", "adj", "mat")):
        cols.append((c, key, x_num0 + gi * (3 * cw + gapg) + si * cw))
for gi, c in enumerate(COH):
    xc = x_num0 + gi * (3 * cw + gapg) + 1.5 * cw
    fig.text(xc / W, y_grp / H, "GSE\n%s" % GSE[c][3:], fontsize=FS, ha="center", va="center",
             color="0", linespacing=1.30)
    for si, key in enumerate(("raw", "adj", "mat")):
        fig.text((x_num0 + gi * (3 * cw + gapg) + si * cw + cw / 2) / W, y_key / H,
                 "str" if key == "adj" else key, fontsize=FS, ha="center", va="center", color="0")
ROWCOL = {"MS": RED, "PLATELET": "0", "AD": GRN, "SCZ": BLU}
for ri, nm in enumerate(("MS", "PLATELET", "AD", "SCZ")):
    yc = y_row0 - ri * rh
    fig.text(x_lab / W, yc / H, nm, fontsize=FS, ha="right", va="center", color=ROWCOL[nm])
    for (c, key, x) in cols:
        v = TAB[nm][key][COH.index(c)]
        fig.text((x + cw - 0.05) / W, yc / H, "\u2014" if v is None else str(v),
                 fontsize=FS, ha="right", va="center", color="0")
    if ri:
        yl = yc + rh / 2 - 0.02
        fig.add_artist(plt.Line2D([x_num0 / W, (cols[-1][2] + cw) / W], [yl / H] * 2, lw=0.4, color="0", zorder=1))
# two mechanism rows (the forest-like dot plot demoted to a small table): row names right-aligned to x_lab, values centred on each cohort group
y_mech = y_row0 - 4 * rh - 0.24
fig.text((x_num0 + 1.5 * cw) / W, (y_mech + 0.30) / H, "PLATELET \u2212 MS", fontsize=FS,
         ha="center", va="center", color="0")
for si, (kk, src) in enumerate((("obs", OBS), ("size-pred", EXP))):
    yy = y_mech - si * 0.42
    fig.text((x_num0 - 0.05) / W, yy / H, kk, fontsize=FS, ha="right", va="center", color="0")
    for gi, c in enumerate(COH):
        xc = x_num0 + gi * (3 * cw + gapg) + 1.5 * cw
        fig.text(xc / W, yy / H, "%+.2f" % src[c], fontsize=FS, ha="center", va="center", color="0")

# ================================================================ row 2: D/E/F cross-channel scatters + G rho small table
yC, hC = 8.75, 2.90
def pairrow(sub):
    return CP[CP.pair.str.contains(sub, regex=False)].iloc[0]

def chan(ax, x, y, col, mspt=None, xl="", yl=""):
    ax.scatter(x, y, s=6, color=col, lw=0, alpha=0.85, zorder=3)
    ax.axhline(0, color="0", lw=0.6, zorder=2); ax.axvline(0, color="0", lw=0.6, zorder=2)
    if mspt:
        ax.scatter([mspt[0]], [mspt[1]], s=26, facecolor=RED, edgecolor="black", lw=0.5, zorder=5)
    ax.tick_params(labelsize=FS, length=1.6, colors="0")
    ax.set_xlabel(xl, fontsize=FS); ax.set_ylabel(yl, fontsize=FS)
    ax.spines[["top", "right"]].set_visible(False)

def note(ax, lines, fx=0.22, fy=0.18, corner="tl"):
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    ax.set_xlim(x0, x1 + fx * (x1 - x0))
    if corner == "br":
        ax.set_ylim(y0 - fy * (y1 - y0), y1)
    else:
        ax.set_ylim(y0, y1 + fy * (y1 - y0))
    xy = (0.985, 0.035, "right", "bottom") if corner == "br" else (0.035, 0.965, "left", "top")
    ax.text(xy[0], xy[1], "\n".join(lines), transform=ax.transAxes, fontsize=FS, ha=xy[2],
            va=xy[3], color="0", linespacing=1.35, zorder=6)

letter(1.80, yC + hC + 0.12, "D")
d0 = ax_cm(1.95, yC, 2.60, hC)
o = pairrow("O1 \u81ea\u5bb6"); msr = T4[T4.trait == MS].iloc[0]
chan(d0, T4.d1_rna, T4.d1_atac, BLU, mspt=(msr.d1_rna, msr.d1_atac), xl="RNA \u03941 (scDRS)", yl="chromatin \u03941 (gene score)")
note(d0, ["\u03c1 %.3f \u2192 %.3f" % (o.rho_raw, o.rho_adj), "MS rank 1 \u2192 18"])
letter(4.70, yC + hC + 0.12, "E")
d2 = ax_cm(5.60, yC, 2.60, hC)
u = pairrow("GSE231585 \u2194 GSE231584"); msu = HUM[HUM.trait == MS].iloc[0]
chan(d2, HUM["A7"], HUM["median_groupdiff"], RED, mspt=(msu["A7"], msu["median_groupdiff"]),
     xl="RNA \u03941 (GSE231585)", yl="chromatin \u03941 (snATAC)")
mu_r = int(HUM[HUM.trait == MS].rank_raw_A7.iloc[0]); mu_c = int(HUM[HUM.trait == MS].rank_raw_chr.iloc[0])
note(d2, ["\u03c1 %.3f \u2192 %.3f" % (u.rho_raw, u.rho_adj), "P %.2f" % u.P_adj, "MS rank %d \u2192 %d" % (mu_r, mu_c)])
letter(8.30, yC + hC + 0.12, "F")
d3 = ax_cm(9.20, yC, 2.30, hC)
mm = pairrow("Mouse GSE250589")
chan(d3, MO.rna_delta1, MO.chr_delta1_MSpk * 1e3, GRN, xl="RNA \u03941 (MS scDRS)", yl="chromatin \u03941 (MS peaks, \u00d71000)")
d3.set_yticks([2.5, 0.0, -2.5])
note(d3, ["\u03c1 %.3f, P %.3f" % (mm.rho_raw, mm.P_raw)], fx=0.30, fy=0.30)

# ---- G: 11 trait-level cross-channel rho pairs (three-point plot -> small table)
letter(11.70, yC + hC + 0.12, "G")
KIND = [("RNA\u2194RNA", "RNA\u2194RNA"), ("unpaired", "RNA\u2194chromatin(\u5206\u5f00\u5efa\u5e93)"),
        ("same nucleus", "RNA\u2194chromatin(\u540c\u6838\u914d\u5bf9)"), ("probe", "RNA\u2194chromatin(\u63a2\u9488)")]
gx, gy0, gdy = 12.05, 11.05, 0.58
fig.text(gx / W, (gy0 + 0.52) / H, "pair type", fontsize=FS, ha="left", va="center", color="0")
fig.text((gx + 3.10) / W, (gy0 + 0.52) / H, "\u03c1 raw", fontsize=FS, ha="right", va="center", color="0")
fig.text((gx + 4.20) / W, (gy0 + 0.52) / H, "\u03c1 strat.", fontsize=FS, ha="right", va="center", color="0")
fig.text((gx + 4.85) / W, (gy0 + 0.52) / H, "n", fontsize=FS, ha="right", va="center", color="0")
GROWS = []
for i, (lab, kind) in enumerate(KIND):
    sub = CP[CP.kind == kind]
    raw = sub.rho_raw.dropna().values; adj = sub.rho_adj.dropna().values
    y = gy0 - i * gdy
    fig.text(gx / W, y / H, lab, fontsize=FS, ha="left", va="center", color="0")
    fig.text((gx + 3.10) / W, y / H, "%.2f\u2013%.2f" % (raw.min(), raw.max()), fontsize=FS, ha="right", va="center", color="0")
    fig.text((gx + 4.20) / W, y / H, "\u2014" if len(adj) == 0 else "%.2f\u2013%.2f" % (adj.min(), adj.max()),
             fontsize=FS, ha="right", va="center", color="0")
    fig.text((gx + 4.85) / W, y / H, "%d" % len(sub), fontsize=FS, ha="right", va="center", color="0")
    GROWS.append((lab, float(raw.min()), float(raw.max()), adj))

# ================================================================ row 3: H two-channel rank bands / I mouse 74-trait scatter
yD, hD = 5.05, 2.60
letter(1.80, yD + hD + 0.12, "H")
rk = {"human RNA": dict(zip(SD.trait, SD.rank_RNA)), "human ATAC": dict(zip(SD.trait, SD.rank_CHR)),
      "mouse RNA": dict(zip(P3.trait, P3.rna_d1_m1.rank(ascending=False).astype(int))),
      "mouse ATAC": dict(zip(P3.trait, P3.atac_d1.rank(ascending=False).astype(int)))}
cols4 = ["human RNA", "human ATAC", "mouse RNA", "mouse ATAC"]
Z4 = np.array([[rk[c][t] for t in S.trait] for c in cols4]).T
xy = [1.95, 3.05, 4.15, 5.25]
for j, c in enumerate(cols4):
    ax = ax_cm(xy[j], yD + 0.40, 0.95, hD - 0.80)
    ax.imshow(Z4[:, j:j + 1], aspect="auto", cmap=cmap_rank, vmin=1, vmax=74)
    ax.set_xlim(-0.5, 0.5); ax.set_ylim(73.5, -0.5); ax.set_xticks([]); ax.set_yticks([]); ax.tick_params(length=0)
    ax.axhline(iMS, color=RED, lw=0.9, zorder=7)
    fig.text((xy[j] + 0.475) / W, (yD + hD - 0.22) / H, c.replace(" ", "\n"), fontsize=FS, ha="center", va="bottom",
             color="0", linespacing=1.55)
    fig.text((xy[j] + 0.475) / W, (yD + 0.28) / H, "%d/74" % int(rk[c][MS]), fontsize=FS, ha="center", va="center", color=RED)
letter(6.55, yD + hD + 0.12, "I")
i1 = ax_cm(7.45, yD, 2.70, hD)
ms3 = P3[P3.trait == MS].iloc[0]
chan(i1, P3.rna_d1_m1, P3.atac_d1, GRN, mspt=(ms3.rna_d1_m1, ms3.atac_d1),
     xl="RNA \u03941 (3 mice, M1)", yl="ATAC \u03941 (gene score)")
note(i1, ["\u03c1 %.3f, P %.2g" % (D["rho_mouse74"], D["P_mouse74"]),
          "MS rank %d \u2192 %d" % (D["mouse_RNA"], D["mouse_ATAC"])], fx=0.28, fy=0.30)

# ---- K: small table of chromatin-channel three-axis MS ranks (fills the row-3 right margin, puts the ATAC axis numbers on the face)
letter(10.60, yD + hD + 0.12, "K")
kk_x, kk_y0, kk_dy = 10.85, yD + hD - 0.45, 0.52
fig.text((kk_x + 2.50) / W, (kk_y0 + 0.95) / H, "MS rank of 74", fontsize=FS, ha="center", va="center", color="0")
fig.text((kk_x + 1.90) / W, (kk_y0 + 0.46) / H, "human", fontsize=FS, ha="center", va="center", color="0")
fig.text((kk_x + 3.10) / W, (kk_y0 + 0.46) / H, "mouse", fontsize=FS, ha="center", va="center", color="0")
KVAL = []
for si, axl in enumerate(("\u03941", "\u03942A", "\u03942B")):
    yy = kk_y0 - si * kk_dy
    fig.text(kk_x / W, yy / H, axl, fontsize=FS, ha="left", va="center", color="0")
    vh = int(GRP_J[4][1][si].loc[MS])
    fig.text((kk_x + 1.90) / W, yy / H, str(vh), fontsize=FS, ha="center", va="center", color="0")
    vm = None
    if len(GRP_J) > 5:
        vm = int(GRP_J[5][1][si].loc[MS])
        fig.text((kk_x + 3.10) / W, yy / H, str(vm), fontsize=FS, ha="center", va="center", color="0")
    else:
        fig.text((kk_x + 3.10) / W, yy / H, "\u2014", fontsize=FS, ha="center", va="center", color="0")
    KVAL.append((axl, vh, vm))

# ================================================================ row 4: J three-axis ledger (cross-channel)
yJ, hJ = 1.00, 3.20
letter(1.80, yJ + hJ + 0.12, "J")
j1 = ax_cm(1.95, yJ, 12.60, hJ - 0.50)
imJ = j1.imshow(ZJ, aspect="auto", cmap=cmap_rank, vmin=1, vmax=74)
j1.set_xlim(-0.5, ZJ.shape[1] - 0.5); j1.set_ylim(73.5, -0.5)
j1.set_xticks([]); j1.set_yticks([]); j1.tick_params(length=0)
for gi in range(1, len(GRP_J)):
    j1.axvline(3 * gi - 0.5, color="0", lw=0.5, zorder=6)
j1.axhline(JROWS.index(MS), color=RED, lw=0.9, zorder=7)
for nm, col in ((AD, GRN), (SCZ, BLU)):
    if nm in JROWS:
        j1.axhline(JROWS.index(nm), color=col, lw=0.7, zorder=7)
for gi, (gname, _, rho_pair) in enumerate(GRP_J):
    j1.text(3 * gi + 1, -3.5, gname, ha="center", va="center", fontsize=FS, color="0", clip_on=False)
    if rho_pair is not None:
        j1.text(3 * gi + 1, -11.0, "\u03c1 %.2f / %.2f" % rho_pair, ha="center", va="center",
                fontsize=FS, color="0", clip_on=False)
for si, lab in enumerate(["\u03941", "\u03942A", "\u03942B"]):
    j1.text(si + 0.5, 75.6, lab, ha="center", va="center", fontsize=FS, color="#8B1A1A" if si == 1 else "0", clip_on=False)
fig.text(1.88 / W, (yJ + (hJ - 0.50) * (73.5 - JROWS.index(MS)) / 74.0) / H, "MS", fontsize=FS,
         color=RED, ha="right", va="center")

# ================================================================ face-number sensor table
face = []
for c in COH:
    face += [("A \u03c1(ln n_peak,\u03941 raw)", GSE[c], RHO_RAW[c], "%+.2f" % RHO_RAW[c], "fig7Sb_matched_window.csv"),
             ("B MS tight-pool rank", GSE[c], PR[c], "%d of 16" % PR[c], "seven4, 15 nearest |dln n_peak|"),
             ("C obs lead", GSE[c], OBS[c], "%+.2f" % OBS[c], "seven4 robust z PLATELET \u2212 MS"),
             ("C size-pred lead", GSE[c], EXP[c], "%+.2f" % EXP[c], "fit slope x dln n_peak")]
for nm in ("MS", "PLATELET", "AD", "SCZ"):
    for key in ("raw", "adj", "mat"):
        for c, v in zip(COH, TAB[nm][key]):
            face.append(("C table %s %s" % (nm, key), GSE[c], -1 if v is None else v,
                         "\u2014" if v is None else str(v), "seven4_confound_ALL_cohorts.csv"))
for tag, key in [("D own", "O1 \u81ea\u5bb6"), ("E unpaired", "GSE231585 \u2194 GSE231584"), ("F mouse", "Mouse GSE250589")]:
    q = pairrow(key)
    face.append((tag, "rho_raw", float(q.rho_raw), "%.3f" % q.rho_raw, "figS7c_crosschannel_pairs.csv"))
    if not pd.isna(q.rho_adj):
        face.append((tag, "rho_adj", float(q.rho_adj), "%.3f" % q.rho_adj, "figS7c_crosschannel_pairs.csv"))
for lab, lo_, hi_, adj in GROWS:
    face.append(("G pair-type rho range", lab, float(lo_), "%.2f\u2013%.2f" % (lo_, hi_), "figS7c_crosschannel_pairs.csv"))
    if len(adj):
        face.append(("G pair-type rho strat range", lab, float(adj.min()), "%.2f\u2013%.2f" % (adj.min(), adj.max()), "figS7c_crosschannel_pairs.csv"))
for nm, v in [("H MS rank human RNA", D["human_RNA"]), ("H MS rank human ATAC", D["human_CHR"]),
              ("H MS rank mouse RNA", D["mouse_RNA"]), ("H MS rank mouse ATAC", D["mouse_ATAC"])]:
    face.append((nm, "MS", v, "%d/74" % v, "same_donor / M3_paired3"))
face += [("I rho mouse 74", "mouse", D["rho_mouse74"], "%.3f" % D["rho_mouse74"], "M3_paired3"),
         ("I P mouse 74", "mouse", D["P_mouse74"], "%.2g" % D["P_mouse74"], "M3_paired3")]
for axl, vh, vm in KVAL:
    face.append(("K ATAC ledger human MS rank %s" % axl, "MS", vh, str(vh), "S7_atac_crossaxis_human_full.csv"))
    if vm is not None:
        face.append(("K ATAC ledger mouse MS rank %s" % axl, "MS", vm, str(vm), "S7_atac_crossaxis_mouse_full.csv"))
for gi, (gname, rms, rho_pair) in enumerate(GRP_J):
    if rho_pair is not None:
        face += [("J ledger rho %s d1-d2A" % gname, "all", rho_pair[0], "\u03c1 %.2f / %.2f" % rho_pair, "fig7Sc_axis_pair_rho.csv"),
                 ("J ledger rho %s d1-d2B" % gname, "all", rho_pair[1], "%.2f" % rho_pair[1],
                  "fig7Sc_axis_pair_rho.csv" if gi < 4 else "S7_atac_crossaxis_mouse_full.csv")]
    for si, axl in enumerate(("\u03941", "\u03942A", "\u03942B")):
        face.append(("J ledger MS rank %s %s (colour cell, not face text)" % (gname, axl), "MS",
                     int(rms[si].loc[MS]), "cell", "see table source"))
pd.DataFrame(face, columns=["tile", "item", "value", "print", "source"]).to_csv(
    os.path.join(TABDIR, "figS7_v17_M6_face_numbers.csv"), index=False)
say("face sensor table %d rows -> figS7_v17_M6_face_numbers.csv" % len(face))
say("J ledger: MS human ATAC %d/%d/%d | mouse ATAC %s | rho(d1,d2A) %.2f rho(d1,d2B) %.2f" % (
    RANK_HU["d1"], RANK_HU["d2a"], RANK_HU["d2b"],
    "pending" not in str(MOUSE_J), rho_hu[0], rho_hu[1]))
save(fig, "figS7_v17_M6", "v17: no titles beside letters / sentences to legend / de-duplicated / three-axis ledger with human + mouse ATAC")
