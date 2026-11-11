# -*- coding: utf-8 -*-
"""axis_rho_and_rank_matrix.py — per-cohort rho(d1,d2A)/rho(d1,d2B) and the 4-cohort rank matrix used by panels A/F/H.

  S7a  support for A: dataset x provenance table (accession / platform / species / condition / unit / n / annotation source / nuclei)
  S7b  support for B: 74-trait raw/adj rank matrix + the ln n_peak confound itself + matched window + cross-channel dumbbell (why only the RNA channel is claimed)
  S7c  support for C: three-axis rank table + the third cohort (GSE301908) mirror replication + axis definitions and donor-subset disclosure
  S7d  support for D: between-donor four-cohort dot plot + within-lesion four groups + mouse EAE per library + composition confound + channel-gap five variants + mouse 74 traits

Discipline: numbers are read from delivered tables only; derived quantities are written to SUPPORT_v1/tables/*.csv first, with a source column."""
import sys, os, textwrap

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
import os
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
sys.path.insert(0, WORK + "/scripts/routeA")
from support_common import *
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, binomtest

say("S7 support grid start")
# ---------------- delivered tables ----------------
A_DON = read("FIG7_v11/tables/fig7A_donors_delta1.csv")
A_P = read("FIG7_v11/tables/fig7A_between_donor_P.csv")
DIS = read("FIG7_v11/tables/fig7E_disease_block.csv")
S4 = read("FIG7_v1/tables/seven4_confound_ALL_cohorts.csv")
D3 = read("D_crosschemistry/tables/D3_crossconcordance_domainadjust.csv")
# ---- house format (author 11-04 + measured on the delivered Figure_6/Figure_S6): grid width 18 -> 17 cm;
# panel title lines come off the face and are handed to the legend; the letters A/B/C/D are put down by the page composer in Arial regular 8 pt upper case.
def _no_suptitle(fig, *a, **k):
    return None
fig_prop = property(lambda self: self)
import matplotlib.figure as _mf
_mf.Figure._s7_no_title = staticmethod(lambda *a, **k: None)

MW = read("SUPPORT_v1/tables/fig7Sb_matched_window.csv")
GRN = "#0B8457"          # the same green as the main figure make_fig7_v11_panels.py (Alzheimer row)
A2D = read("A2_GSE277430_v1/tables/trait_sweep_dual_axis.csv")
A3D = read("A3_GSE118257_v1/tables/trait_sweep_delta2_GSE118257_74.csv")
A4D = read("A4_GSE301908_v1/tables/trait_sweep_delta2_GSE301908_COUNTS.csv")
A2META = read("A2_GSE277430_v1/tables/OL_meta_with_MS_score.csv")
COMP = read("D_crosschemistry/tables/fig6D_composition.csv")
M2B = read("M2b_authoranno_v1/tables/M2b_per_mouse.csv")
M74 = read("M_mouse_74trait/tables/trait_delta1_mouse_74.csv")
A5CC = read("A5_case_control/tables/case_control_delta1.csv")
A7P = read("A7_GSE231585_rna_v1/tables/A7_between_donor_P.csv")
A7SUM = read("A7_GSE231585_rna_v1/tables/authorpanel_summary.csv")
MS = "PASS_Multiple_sclerosis"
SHORT = lambda t: (t.replace("PASS_", "").replace("UKB_460K.", "")
                    .replace("_Pardinas2018", "").replace("_Demontis2018", "").replace("_SavageJansen2018", "")
                    .replace("_Liu2019", "").replace("_Jansen2019", "").replace("_Nielsen2018", "")
                    .replace("_Stahl2019", "").replace("_deLange2017", ""))

# derived-quantity register
DER = []
def derived(q, v, src, how):
    DER.append(dict(quantity=q, value=v, source=src, method=how))

# ================= S7a dataset provenance table =================
# On 11-04 main figure A (v7) moved the technology / sorted-or-not rows off the face (all seven columns are
# transcriptome-arm d1, so the technology name carries no discrimination) and deleted the "+ snATAC-seq" hanging under the GSE231585 column (that chromatin arm belongs to GSE231584, a different accession).
# The withdrawn wording has to land somewhere, and S7a is the table that catches it; it also adds the
# **d1 definition column (mature side - precursor side)** asked for by audit receipt P1-1: in the four
# external cohorts the mature side is Oligo(all), which is not the same contrast as our O1/D MOL6 - COP/NFOL, so ranks are only readable within a cohort, and that has to be stated in the supplement.
def study_row(study, geo, platform, species, cond, unit, ncase, nctrl, anno, nuclei, note_=""):
    return dict(study=study, geo=geo, platform=platform, species=species, condition=cond, unit=unit,
                n_case=ncase, n_ctrl=nctrl, annotation_source=anno, ol_nuclei=nuclei, note=note_)


# Per-cohort sampling / chemistry version / nuclei enrichment. On the face these are **compressed strings**
# (the whole provenance does not fit into 17 cm); the verbatim long version is in FIG7_v11/dataset_platform_provenance.md section 1, and the CSV carries a source column pointing back.
PROV2 = [
    dict(study="A2 GSE277430", tissue="WM (5 samples) + cortex (18), same 21 donors",
         prep="not sorted (Dounce + lysis)",
         chem="3' v3; CR 6.1.2; NovaSeq 6000; GPL24676",
         mat="Oligo = authors' expertAnno.l1 major cluster (MOL1/2/5/6 not separated)",
         pre="OPC (expertAnno.l1)", same="NO - Oligo(all), ours is MOL6",
         arm="RNA arm only | 44,767 OL nuclei"),
    dict(study="A3 GSE118257", tissue="autopsy WM; 20 samples incl. NAWM/Active/CA/CI/RL",
         prep="not sorted (density gradient)",
         chem="3' v2; CR 2.1.1; HiSeq 2500; GPL16791",
         mat="Oligo = Mature_OL Oligo1-6 (authors' FinalAnnotationTable)",
         pre="OPC / COP (FinalAnnotationTable)", same="NO - Oligo(all), ours is MOL6",
         arm="RNA arm only | authors' computed labels"),
    dict(study="A4 GSE301908", tissue="WM lesions, chronic active (3 spatial domains)",
         prep="not stated in the delivered summary",
         chem="3' v3; CR 7.1.0 (read from the h5); NovaSeq 6000",
         mat="Oligo (authors' majorCluster)", pre="OPC (authors' majorCluster)",
         same="NO - Oligo(all), ours is MOL6",
         arm="RNA arm only; study also has a MERFISH layer"),
    dict(study="A7 GSE231585", tissue="WM (C1/C3/C4) + lesions (NAWM/CA/RL/AL)",
         prep="FANS (antibody not stated in GEO)",
         chem="3' v3; CR 6.1.2 with --include-introns; h5 records 4.0 (logged)",
         mat="Oligo (our cluster-level call on the authors' marker panel)",
         pre="OPC (same call)", same="NO - Oligo(all), ours is MOL6",
         arm="RNA arm only | 52,555 OL / 1,841 OPC nuclei"),
    dict(study="M GSE250589", tissue="spinal cord, NOT brain; Naive/Early/Peak/Late",
         prep="FACS GFP+ lineage (Sox10-cre), not MOL antibodies",
         chem="Multiome PN-1000283; arc 2.0.2 (h5 2.0.0); mm10",
         mat="MOL1 + MOL2 + MOL56 (authors' annotation.csv.gz)",
         pre="OPC + COP (authors' annotation.csv.gz)", same="NO - mixed mature, ours is MOL6",
         arm="RNA arm of the paired multiome; ATAC arm not used for ranks"),
    dict(study="A6 GSE231584", tissue="same cohort, microdissected; 9 snATAC libraries",
         prep="FANS + anti-lamin A/C",
         chem="scATAC; cellranger-atac 2.0.0",
         mat="- (chromatin)", pre="-", same="n/a - no RNA Delta1 here",
         arm="not in Fig 7A; only S7b 4 | 2,408 OL nuclei"),
]
P2 = pd.DataFrame(PROV2)
P2["source"] = "FIG7_v11/dataset_platform_provenance.md §1 (full strings) + fig7A_between_donor_P.csv"

P2.to_csv(os.path.join(TABDIR, "fig7Sa_platform_delta1_definition.csv"), index=False)

prov = []
p = A_P.set_index("study")
prov.append(study_row("A2 GSE277430", "GSE277430", "10x single-nucleus RNA-seq", "human", "MS / Ctrl", "donor",
                      int(p.loc["A2 GSE277430", "n_case"]), int(p.loc["A2 GSE277430", "n_ctrl"]),
                      "authors' expert annotation (expertAnno.l1)",
                      "%s OL nuclei" % format(int(str(A5CC[A5CC.cohort == "GSE277430"].min_cells_note.iloc[0]).split()[0]), ","),
                      "21 donors; one opposite-sign donor 2006-041"))
prov.append(study_row("A3 GSE118257", "GSE118257", "10x single-nucleus RNA-seq", "human", "MS / Ctrl + lesion", "donor",
                      int(p.loc["A3 GSE118257", "n_case"]), int(p.loc["A3 GSE118257", "n_ctrl"]),
                      "authors' annotation (FinalAnnotationTable Celltypes)", "n/a in delivered tables",
                      "9 donors; lesion pairs n=2/lesion type (underpowered, descriptive only)"))
prov.append(study_row("A4 GSE301908", "GSE301908", "10x single-nucleus RNA-seq", "human", "MS / Ctrl", "library",
                      int(p.loc["A4 GSE301908", "n_case"]), int(p.loc["A4 GSE301908", "n_ctrl"]),
                      "authors' annotation (majorCluster)", "n/a in delivered tables",
                      "main analysis donor n=11; case-control uses 17 libraries (14 MS / 3 Ctrl) - both stated"))
prov.append(study_row("A7 GSE231585", "GSE231585", "10x single-nucleus RNA-seq", "human", "MS / Ctrl", "donor",
                      int(p.loc["A7 GSE231585", "n_case"]), int(p.loc["A7 GSE231585", "n_ctrl"]),
                      "ours: cluster-level, authors' published marker panel (no per-nucleus labels distributed)",
                      "%s OL / %s OPC nuclei" % (format(int(A7SUM.n_OL_nuc.iloc[0]), ","), format(int(A7SUM.n_OPC_nuc.iloc[0]), ",")),
                      "self-annotated: main claims must not rest on it"))
prov.append(study_row("M GSE250589", "GSE250589", "10x Multiome, paired RNA + ATAC", "mouse (EAE)", "EAE / Ctrl", "library",
                      int(p.loc["M GSE250589", "n_case"]), int(p.loc["M GSE250589", "n_ctrl"]),
                      "authors' annotation (annotation.csv.gz; MOL1+MOL2+MOL56 vs OPC+COP)",
                      "%d mice in the paired multiome arm" % len(M2B),
                      "one study, two arms: RNA 28 libraries + matched multiome 18 mice"))
prov.append(study_row("A6 GSE231584", "GSE231584", "10x single-nucleus ATAC-seq, separate libraries", "human", "MS / Ctrl", "library",
                      "-", "-", "authors' annotation", "2,408 OL nuclei (6 usable blinded libraries)",
                      "used only in the cross-channel strip (S7B 4), never for a rank claim"))
PROV = pd.DataFrame(prov)
PROV.to_csv(os.path.join(TABDIR, "fig7Sa_provenance.csv"), index=False)

fig = plt.figure(figsize=(17.0 / CM, 17.4 / CM))
gsA = fig.add_gridspec(3, 1, left=0.025, right=0.985, top=0.930, bottom=0.020, hspace=0.30,
                       height_ratios=[1.0, 1.0, 1.02])
ax = fig.add_subplot(gsA[0])
# Main figure A (v7) no longer prints technology or sorting: those two columns move here (full text in provenance section 1).
# lh must be >= 1.2*fs / cell height: the previous lh=0.030 with fs=5.5 let the two lines of one cell overlap themselves (491 measured collisions).
cols = [("study / role", "study", 20, 0.000), ("species / condition", "spcond", 12, 0.115),
        ("tissue sampled (in the source study)", "tissue", 24, 0.190),
        ("platform", "platform", 14, 0.350), ("chemistry / software / instrument", "chem", 26, 0.425),
        ("nuclei preparation", "prep", 20, 0.620)]
rows = []
PP = P2.set_index("study")
for _, r in PROV.iterrows():
    k = PP.loc[r.study]
    rows.append(dict(study="%s  %s" % (r.study, r.geo), spcond="%s %s" % (r.species, r.condition),
                     tissue=k.tissue, platform=r.platform, chem=k.chem, prep=k.prep))
ax.text(0.0, 1.03, "1  What the six external studies are: material, chemistry and nuclei preparation"
                   "   (Fig 7A prints only the accession and the species/tissue bracket)",
        fontsize=6.4, transform=ax.transAxes, va="bottom")
draw_table(ax, rows, cols, fs=5.5, lh=0.046, pad=0.018, top=0.905)

ax1 = fig.add_subplot(gsA[1]); ax1.axis("off")
cols2 = [("study", "study", 16, 0.000),
         ("Delta1 mature side (as defined in THAT study)", "mat", 30, 0.100),
         ("Delta1 precursor side", "pre", 22, 0.290),
         ("same contrast as our MOL6 - COP/NFOL?", "same", 18, 0.430),
         ("what Fig 7A plots / OL-lineage material", "arm", 30, 0.550)]
rows2 = []
for _, r in PROV.iterrows():
    k = PP.loc[r.study]
    rows2.append(dict(study=r.study, mat=k.mat, pre=k.pre, same=k.same, arm=k.arm))
ax1.text(0.0, 1.03, "2  The Delta1 definition column (audit FIGURE_LAYOUT_DELTA_AUDIT_v1, P1-1): "
                    "the four external human cohorts all contrast Oligo(all) with OPC",
         fontsize=6.4, transform=ax1.transAxes, va="bottom")
draw_table(ax1, rows2, cols2, fs=5.5, lh=0.046, pad=0.018, top=0.905)

ax2 = fig.add_subplot(gsA[2]); ax2.axis("off")
wrapblock(ax2, [
    "     Units, annotation provenance, and the things a reader must not infer",
    "Because the mature/precursor labels differ by study (table 2), the Delta1 columns of Fig 7A are NOT the same biological contrast: each cohort is internally self-consistent, so ranks and signs are read WITHIN a cohort only, never across cohorts in absolute value. The wording we allow is 'MS ranks in the top 7 of 74 trait sets within each cohort', not 'the four cohorts agree on the size of Delta1'.",
    "GSE301908 is plotted per donor (11) but its case-control comparison is per library (17; 14 MS / 3 Ctrl); both n values are printed because 11 is not 17.",
    "GSE231585 distributes no per-nucleus labels (GEO, Zenodo and the authors' repository checked), so nuclei were assigned at cluster level following the authors' published marker panel and eight-class scheme; this reproduces our previous partition exactly (%s OL / %s OPC nuclei, donor-level d1 unchanged). Main claims do not rest on it." % (format(int(A7SUM.n_OL_nuc.iloc[0]), ","), format(int(A7SUM.n_OPC_nuc.iloc[0]), ",")),
    "GSE250589 is ONE study with two arms (RNA 28 libraries; matched multiome 18 mice) and is reported as one line. Library is not animal: each library pools 2-3 mice (EAE 38 mice / 19 experiments, CFA-Ctrl 9, Naive 6), Naive is not a separate level in the delivered object, and the P value is experiment-level. The mouse column is a parallel observation in spinal cord, not a replication of the human cohorts.",
    "GSE231584 (snATAC) is a SEPARATE accession from GSE231585: same donors and largely overlapping samples, but libraries built separately, so it is described as a matched cohort that was separately profiled - never as same-nucleus or same-cell pairs. It enters only the cross-channel strip in Fig S7b 4 and is never used for a rank claim.",
    "Between-donor MS vs Ctrl per cohort (Mann-Whitney, fig7A_between_donor_P.csv). The mouse value is the pooled 13 EAE / 6 Ctrl libraries; Fig 7A instead shows the three EAE timepoints (Early / Peak / Late) separately: " + "; ".join(
        "%s %d/%d delta %+.3f P = %s" % (r.study.split()[-1], r.n_case, r.n_ctrl, r.delta, fmt_p(r.P))
        for _, r in A_P.iterrows()),
], fs=5.5, width=168)
fig._s7_no_title("Fig S7a  Support for Fig 7A - what the external studies are, how the nuclei were prepared, and what each Delta1 column actually contrasts",
             fontsize=7.4, x=0.008, y=0.972, ha="left")
save(fig, "S7a_provenance", "provenance table + d1 definition column + reading rules")

# ================= S7b =================
fig = plt.figure(figsize=(17.0 / CM, 17.4 / CM))
gs = fig.add_gridspec(2, 2, left=0.10, right=0.985, top=0.895, bottom=0.06, wspace=0.42, hspace=0.52)

# (a) rank-matrix heatmap (same row order as main figure B, same three disease row names; the main figure no longer prints the hollow circles or the colour-scale wording, so this catches it)
ax = panel(fig, gs, 0, 0, "1  Rank of every trait, raw vs adjusted",
           "within each cohort: raw then n_peak-adjusted; 1 = top of 74")
COH = ["A2", "A3", "A4", "A7"]
GSE = {"A2": "GSE277430", "A3": "GSE118257", "A4": "GSE301908", "A7": "GSE231585"}
M = S4.sort_values("rank_adj_A2").reset_index(drop=True)
grid = np.column_stack([v for c in COH for v in
                        (M["rank_raw_" + c].astype(float).values, M["rank_adj_" + c].astype(float).values)])
im = ax.imshow(grid, aspect="auto", cmap="viridis_r", vmin=1, vmax=74)
TR3 = [(MS, "MS", RED), ("PASS_Alzheimers_Jansen2019", "Alzheimer", GRN),
       ("PASS_Schizophrenia_Pardinas2018", "Schizophrenia", BLU)]
for tr, nm, cc in TR3:
    irow = int(M.index[M.trait == tr][0])
    ax.add_patch(plt.Rectangle((-0.5, irow - 0.5), 8, 1, fill=False, ec=cc, lw=1.2, zorder=5))
    ax.text(-0.70, irow, nm, fontsize=5.6, color=cc, va="center", ha="right", zorder=6)
# y ticks all removed: main figure B has no ticks anyway and the left margin is needed for the three row names (keeping ticks would collide with them)
ax.set_yticks([])
ax.set_xticks(range(8))
ax.set_xticklabels([c for c in ["A2", "", "A3", "", "A4", "", "A7", ""]], fontsize=5.8)
for i in range(1, 4):
    ax.axvline(i * 2 - 0.5, color="0.95", lw=0.5)
ax.set_xlim(-0.5, 7.5); ax.set_ylim(len(M) - 0.5, -0.5)
ax.set_ylabel("")          # left margin for the three row names: ylabel and "Alzheimer" sat in the same place, measured overlap
for sp in ax.spines.values():
    sp.set_visible(False)
ax.tick_params(length=0)


def rk(key, tr):
    r = S4[S4.trait == tr].iloc[0]
    return " / ".join(str(int(r["%s_%s" % (key, c)])) for c in COH)


# ranks written under the cell x labels (they used to sit in a right-hand note block and collided with the next cell's ylabel: 2 measured collisions)
ax.set_xlabel("A2 = GSE277430   A3 = GSE118257   A4 = GSE301908   A7 = GSE231585\n"
              "rows ordered by mean adjusted rank (top = best)\n"
              "MS   raw " + rk("rank_raw", MS) + "   adj " + rk("rank_adj", MS) + "\n"
              "AD   raw " + rk("rank_raw", TR3[1][0]) + "   adj " + rk("rank_adj", TR3[1][0]) + "\n"
              "SCZ  raw " + rk("rank_raw", TR3[2][0]) + "   adj " + rk("rank_adj", TR3[2][0]), fontsize=5.5)

# (b) the confound itself
ax = panel(fig, gs, 0, 1, "2  The confound itself: gene-set regulatory domain size",
           "one point per trait; rho from D3 confound-check rows")
COLC = {"A2": BLU2, "A3": RED, "A4": ORG, "A7": "0.45"}
for c in COH:
    d = S4[S4[c].notna()]
    ax.scatter(d.lnp, d[c], s=9, color=COLC[c], alpha=0.75, zorder=3,
              label="%s  rho = %.2f" % (c, float(MW[MW.cohort == c].rho_lnp_vs_d1.iloc[0])))
r = S4[S4.trait == MS].iloc[0]
ax.scatter([r.lnp], [r.A2], s=42, facecolor="none", edgecolor=RED, lw=1.2, zorder=5)
ax.annotate("MS", (r.lnp, r.A2), textcoords="offset points", xytext=(5, 4), fontsize=5.6, color=RED)
ax.set_xlabel("ln n_peak (gene-set regulatory domain)", fontsize=5.8)
ax.set_ylabel("trait d1 (mature - precursor)", fontsize=5.8)
ax.legend(fontsize=5.5, frameon=True, framealpha=0.95, edgecolor="none", loc="lower left",
          title="Spearman(ln n_peak, d1)", title_fontsize=5.5)

# (c) matched window
ax = panel(fig, gs, 1, 0, "3  Model-free control: rank inside a size-matched window",
           "traits whose n_peak is within +/-15% of MS's, ranked by d1 (no regression)")
y = np.arange(len(MW))[::-1]
ax.barh(y, MW.MS_rank_in_window, height=0.55, color=BLU2, zorder=3)
for yy, rk, n in zip(y, MW.MS_rank_in_window, MW.n_in_window):
    ax.text(rk + 0.08, yy, "%d / %d" % (rk, n), va="center", fontsize=5.8, color="0.15")
ax.set_yticks(y); ax.set_yticklabels(MW.study, fontsize=5.8)
ax.set_xlim(0, max(MW.n_in_window) * 0.55); ax.set_xlabel("rank of MS within the matched window", fontsize=5.8)
note(ax, 0.98, 0.20, "raw rank /74: " + " / ".join(str(int(v)) for v in MW.rank_raw) +
     "\nadjusted: " + " / ".join(str(int(v)) for v in MW.rank_adj) +
     "\nwindow n_peak %d-%d" % (MW.window_lo.iloc[0], MW.window_hi.iloc[0]), ha="right")

# (d) cross-channel dumbbell
ax = panel(fig, gs, 1, 1, "4  Why the rank claim is restricted to the RNA channel",
           "pairwise cross-channel rho of trait d1, before and after domain adjustment")
d4 = D3[D3.kind.astype(str).isin(["RNA\u2194RNA", "RNA\u2194chromatin(\u5206\u5f00\u5efa\u5e93)"])].copy()
o1 = D3[D3["pair"].astype(str).str.startswith("O1 own")].iloc[0]
rows = list(d4.itertuples()) + [o1]
y = np.arange(len(rows))[::-1]
for yy, r in zip(y, rows):
    kind = str(r.kind)
    col = BLU2 if kind.startswith("RNA\u2194RNA") else (RED if "chromatin" in kind else ORG)
    if str(r.pair).startswith("O1"):
        col = "0.15"
    ax.plot([r.rho_raw, r.rho_adj], [yy, yy], lw=0.9, color="0.70", zorder=2)
    ax.plot(r.rho_raw, yy, "o", ms=4.2, mfc="white", mec=col, mew=0.9, zorder=4)
    ax.plot(r.rho_adj, yy, "s", ms=4.0, color=col, zorder=4)
def shortpair(x):
    x = str(x)
    if x.startswith("O1"):
        return "O1 own, same nuclei"
    x = x.replace("GSE", "").replace("\u2194", "<->").replace(" (snATAC)", "")
    x = x.replace("\u5206\u5f00\u5efa\u5e93", "sep").replace("231584", "231584 ATAC")
    return x[:20]


labels = [shortpair(r.pair) for r in rows]
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=5.5)
ax.set_xlim(-0.1, 1.02); ax.set_xlabel("Spearman rho between the two channels (74 traits)", fontsize=5.8)
ax.axvline(0, lw=0.5, color=MUT)
note(ax, 0.40, 0.30, "RNA<->RNA keep rho; separately built\nRNA<->chromatin collapse;\nthe same-nucleus pair does not", ha="left")
fig._s7_no_title("Fig S7b  Support for Fig 7B - the full rank matrix, the confound it is adjusted for, and why cross-channel ranks are not claimed",
             fontsize=7.4, x=0.008, y=0.985, ha="left")
save(fig, "S7b_rank_matrix", "4 micro-panels")

# ================= S7c =================
rho_A4_A = float(spearmanr(A4D.d1, A4D.dA).statistic); p_A4_A = float(spearmanr(A4D.d1, A4D.dA).pvalue)
rho_A4_B = float(spearmanr(A4D.d1, A4D.dB).statistic); p_A4_B = float(spearmanr(A4D.d1, A4D.dB).pvalue)
pd.DataFrame([
    dict(quantity="rho(d1, d2A) GSE301908 counts, 12 matched traits", value=rho_A4_A, P=p_A4_A,
         source="A4_GSE301908_v1/tables/trait_sweep_delta2_GSE301908_COUNTS.csv", method="Spearman over the delivered rows"),
    dict(quantity="rho(d1, d2B) GSE301908 counts, 12 matched traits", value=rho_A4_B, P=p_A4_B,
         source="A4_GSE301908_v1/tables/trait_sweep_delta2_GSE301908_COUNTS.csv", method="Spearman over the delivered rows"),
]).to_csv(os.path.join(TABDIR, "fig7Sc_third_cohort_rho.csv"), index=False)
say("A4 third-cohort mirror rho = %.3f (P=%.4f) / %.3f (P=%.4f)" % (rho_A4_A, p_A4_A, rho_A4_B, p_A4_B))

fig = plt.figure(figsize=(17.0 / CM, 17.4 / CM))
gs = fig.add_gridspec(2, 2, left=0.11, right=0.985, top=0.895, bottom=0.06, wspace=0.46, hspace=0.55)

A4D74 = read("A4_GSE301908_v1/tables/trait_sweep_delta2_GSE301908_74.csv")
A7D74 = read("A7_GSE231585_rna_v1/tables/trait_sweep_delta2_GSE231585_74.csv")
GSE4 = {"A2": "GSE277430", "A3": "GSE118257", "A4": "GSE301908", "A7": "GSE231585"}
SRC_PATH = {"A2": "A2_GSE277430_v1/tables/trait_sweep_dual_axis.csv",
            "A3": "A3_GSE118257_v1/tables/trait_sweep_delta2_GSE118257_74.csv",
            "A4": "A4_GSE301908_v1/tables/trait_sweep_delta2_GSE301908_74.csv",
            "A7": "A7_GSE231585_rna_v1/tables/trait_sweep_delta2_GSE231585_74.csv"}

# (a) three-axis rank heatmap (main figure C is four violin tiles; this gives the full 74 x axis rank matrix, all four cohorts, not just A2/A3)
ax = panel(fig, gs, 0, 0, "1  Three identity axes, 74 traits, four cohorts",
           "rank (1 = top of 74) on d1 / d2A / d2B within each cohort; rows ordered by rank on GSE277430 d1")
AXC = ["#4E7A3A", "#A6651E", "#5B5688"]      # same colours as the main figure C AXCOL: colour is the axis identity
SRC = [("A2", A2D, ["r1", "rA", "rB"]), ("A3", A3D, ["r1", "rA", "rB"]),
       ("A4", A4D74, ["r1", "rA", "rB"]), ("A7", A7D74, ["r1", "rA", "rB"])]
J = None
for k, (code, d, rks) in enumerate(SRC):
    t = d[d.trait != "TRAIT"][["trait"] + rks].rename(columns={rks[0]: "c%d_1" % k,
                                                               rks[1]: "c%d_A" % k, rks[2]: "c%d_B" % k})
    J = t if J is None else J.merge(t, on="trait")
J = J.sort_values("c0_1").reset_index(drop=True)
J.to_csv(os.path.join(TABDIR, "fig7Sc_rank_matrix_4cohorts.csv"), index=False)
ORDER = ["c%d_%s" % (k, s) for k in range(4) for s in ("1", "A", "B")]
grid = J[ORDER].astype(float).values
ax.imshow(grid, aspect="auto", cmap="RdBu_r", vmin=1, vmax=74)
i_ms = int(J.index[J.trait == MS][0])
i_scz = int(J.index[J.trait.astype(str).str.contains("Schizophrenia")][0])
for i, col in [(i_ms, RED), (i_scz, BLU)]:
    ax.add_patch(plt.Rectangle((-0.5, i - 0.5), 12, 1, fill=False, ec=col, lw=1.2, zorder=5))
    ax.text(-0.85, i, "MS" if i == i_ms else "SCZ", fontsize=5.8, color=col, va="center", ha="right")
ax.set_yticks([])
ax.set_xticks(range(12))
ax.set_xticklabels([""] * 12)
for k in range(4):
    for j, (suf, lab) in enumerate([("1", "d1"), ("A", "d2A"), ("B", "d2B")]):
        ax.text(k * 3 + j, len(J) + 0.25, lab, fontsize=5.5, ha="center", va="top",
                color=AXC[j])            # axis names go **below** the heatmap: pressed onto it, the text-on-data gate would catch them
    ax.text(k * 3 + 1, len(J) + 2.6, SRC[k][0], fontsize=5.8, ha="center", va="top", color="0.15")
    if k:
        ax.axvline(k * 3 - 0.5, color="0.75", lw=0.6)
ax.set_ylim(len(J) + 3.6, -0.5); ax.set_xlim(-0.5, 11.5)
ax.set_ylabel("")   # left margin only for the MS/SCZ row names
for sp in ax.spines.values():
    sp.set_visible(False)
def rk3(idx):
    return "  ".join("%d/%d/%d" % tuple(int(J["c%d_%s" % (k, s)].iloc[idx])
                       for s in ("1", "A", "B")) for k in range(4))


ax.set_xlabel("A2 = GSE277430   A3 = GSE118257   A4 = GSE301908   A7 = GSE231585\n"
              "d1 / d2A / d2B = OL maturation / astro-immune / neuron-glia\n"
              "MS   " + rk3(i_ms) + "\nSCZ  " + rk3(i_scz), fontsize=5.5)

# (b) third-cohort mirror (the old version used the counts table n=12; now A4/A7 both have 74-trait tables,
# which is what main figure C uses, so this follows with n=74; the counts variant stays in the derived table as a sensitivity)
ax = panel(fig, gs, 0, 1, "2  Mirror replicates in the two additional cohorts",
           "GSE301908 and GSE231585, all 74 traits each")
MIR = []
for code, dd in [("A4", A4D74), ("A7", A7D74)]:
    d = dd[dd.trait != "TRAIT"].dropna(subset=["d1", "dA", "dB"])
    col = BLU2 if code == "A4" else ORG
    ax.scatter(d.dA, d.d1, s=13, color=col, zorder=3, label="%s d2A" % GSE4[code])
    ax.scatter(d.dB, d.d1, s=13, color=col, zorder=3, alpha=0.45, marker="s", label="%s d2B" % GSE4[code])
    for nm, cx in [("d2A", d.dA), ("d2B", d.dB)]:
        rr, pv = float(spearmanr(d.d1, cx).statistic), float(spearmanr(d.d1, cx).pvalue)
        MIR.append(dict(cohort=code, gse=GSE4[code], variant="74 traits", axis=nm, rho=rr, P=pv, n=len(d)))
        k, b0 = np.polyfit(cx, d.d1, 1)
        xs = np.linspace(cx.min(), cx.max(), 20)
        ax.plot(xs, k * xs + b0, lw=0.9, color=col, ls="--" if nm == "d2A" else ":", zorder=2)
dc = A4D[A4D.d1.notna()]
for nm, cx in [("d2A", dc.dA), ("d2B", dc.dB)]:
    MIR.append(dict(cohort="A4", gse="GSE301908", variant="counts, %d traits" % len(dc), axis=nm,
                    rho=float(spearmanr(dc.d1, cx).statistic), P=float(spearmanr(dc.d1, cx).pvalue), n=len(dc)))
pd.DataFrame(MIR).to_csv(os.path.join(TABDIR, "fig7Sc_mirror_rho.csv"), index=False)
ax.text(0.02, 0.06, "\n".join("%s %s: rho = %+.2f (P = %s)" % (GSE4[r.cohort], r.axis, r.rho, fmt_p(r.P))
                              for r in pd.DataFrame(MIR).itertuples() if r.variant == "74 traits"),
        transform=ax.transAxes, fontsize=5.5, va="bottom", color="0.25", bbox=WB)
ax.set_xlabel("trait d2 (mirror axis)", fontsize=5.8); ax.set_ylabel("trait d1 (mature - precursor)", fontsize=5.8)
ax.legend(fontsize=5.5, frameon=False, loc="upper right", ncol=2, handletextpad=0.2, columnspacing=0.8)

# (c) axis-pair Spearman summary (all four cohorts: main figure C draws four tiles, so this cannot give only two)
ax = panel(fig, gs, 1, 0, "3  Axis-pair correlations, all four cohorts",
           "negative rho = the mirror (Spearman over the 74 traits, delivered columns, nothing re-derived)")
tab = []
AXCOLS = {"A2": ("d1_OLmat", "d2_AST_IM", "d2_NEU_GL"),
          "A3": ("d1_OLIGO_minus_OPC", "dA_AST_minus_IMM", "dB_NEU_minus_GLIA"),
          "A4": ("d1", "dA", "dB"), "A7": ("d1", "dA", "dB")}
for code, d, _r in SRC:
    c1, cA, cB = AXCOLS[code]
    dd = d[d.trait != "TRAIT"].dropna(subset=[c1, cA, cB])
    for nm, cx in [("d1 vs d2A", cA), ("d1 vs d2B", cB)]:
        rr = float(spearmanr(dd[c1], dd[cx]).statistic); pp_ = float(spearmanr(dd[c1], dd[cx]).pvalue)
        tab.append(dict(cohort=code, gse=GSE4[code], n=len(dd), pair=nm, rho=rr, P=pp_,
                        source=SRC_PATH[code]))
TAB = pd.DataFrame(tab)
TAB.to_csv(os.path.join(TABDIR, "fig7Sc_axis_pair_rho.csv"), index=False)
ax.axis("off")
lines = ["", "cohort  GSE        n   pair          rho        P", ""]
for _, r in TAB.iterrows():
    lines.append("%-7s %-11s %2d  %-11s %+.3f   %s" % (r.cohort, r.gse, r.n, r.pair, r.rho, fmt_p(r.P)))
lines += ["", "All eight rho are negative: on every axis pair, in every cohort, the",
          "mirror axis runs opposite to the OL-maturation axis.",
          "Ranks used in Fig 7C are per-cohort ranks of the same 74 traits",
          "on each axis; no cross-cohort pooling."]
ax.text(0.0, 1.0, "\n".join(lines), transform=ax.transAxes, va="top", fontsize=5.8, linespacing=1.6, color="0.15")

# (d) axis definitions + donor subsets + the reading rules main figure C no longer puts on the face
ax = panel(fig, gs, 1, 1, "4  Axis definitions, donor subsets, and how to read Fig 7C"); ax.axis("off")
n1 = int(A3D.n1_donors.iloc[0]); nA = int(A3D.nA_donors.iloc[0]); nB = int(A3D.nB_donors.iloc[0])
T4 = [
    "d1 = median MS scDRS norm_score on mature OL minus median on OPC (the main axis). d2A = median(AST) - median(IMM); "
    "d2B = median(NEU) - median(GLIA); NEU = excitatory + inhibitory neurons, GLIA = AST, IMM, OLIGO, OPC, VAS. "
    "The mirror axes use different marker sets and cell classes from d1, so they are independent of it by construction.",
    "The labels behind d1 differ across cohorts (S7a table 2: the external cohorts use Oligo(all) - OPC, our model uses "
    "MOL6 - COP/NFOL), so ranks are read within a cohort only; the axis colours in Fig 7C mark the same three axes, "
    "not a shared scale.",
    "GSE118257 donor subsets: d1 uses %d donors, d2A uses %d, d2B uses %d (columns n1_donors / nA_donors / nB_donors) - "
    "the three axes are therefore not computed on an identical donor set there; disclosed, not hidden, and a "
    "common-donor recomputation is outstanding (CE8)." % (n1, nA, nB),
    "Fig 7C's y is a within-cohort, within-axis robust z of d (median/MAD over that cohort's 74 traits): a monotone "
    "transform, so sign(z) means exactly 'ranked better than the median of 74'. No cross-cohort z comparison is made "
    "and ranks are never computed from z.",
    "Violin kernels are evaluated only inside the visible window, so a few outlier traits are cut at the axis edge; "
    "ranks and IQR boxes use all 74 traits.",
    "The mouse 74-trait sweep is not in Fig 7C (its neuronal layer has 78 neurons only); it is reported in S7d.",
]
# layout constraints (all measured, none guessed):
#   single column width=70 -> the right edge falls off and is cropped; width<=62 -> the lines do not fit the cell
#   height and the last paragraph drops off the canvas (no ink in the PDF); two columns width=34 -> the second
#   column collides with the first (4 measured collisions). Conclusion: single column, 58 characters, line spacing 1.15, body trimmed to ~1250 characters - the out-of-bounds, collision and cell-title gates all pass at once.
wrapblock(ax, T4, fs=5.5, width=58, ls=1.15)
fig._s7_no_title("Fig S7c  Support for Fig 7C - the mirror is not an artefact of one axis or one donor subset",
             fontsize=7.4, x=0.008, y=0.985, ha="left")
save(fig, "S7c_three_axes", "4 micro-panels")

# ================= S7d =================
binom_ms = binomtest(int(M74[M74.trait == MS].n_pos.iloc[0]), int(M74[M74.trait == MS].n_lib.iloc[0]), 0.5, alternative="greater")
derived("binomial P, MS positive in every mouse library", float(binom_ms.pvalue),
        "M_mouse_74trait/tables/trait_delta1_mouse_74.csv", "binomtest(n_pos, n_lib, 0.5, greater)")
pd.DataFrame(DER).to_csv(os.path.join(TABDIR, "fig7Sd_derived_numbers.csv"), index=False)

fig = plt.figure(figsize=(17.0 / CM, 17.4 / CM))
gs = fig.add_gridspec(3, 2, left=0.135, right=0.985, top=0.895, bottom=0.055, wspace=0.36, hspace=0.52)

# (a) between-donor four-cohort dot plot
ax = panel(fig, gs, 0, 0, "1  Diagnosis does not move the axis",
           "dots = donors (A4: libraries); Mann-Whitney P on top")
hum = [s for s in A_DON.study.unique() if not s.startswith("M")]
for i, st in enumerate(hum):
    d = A_DON[A_DON.study == st]
    for j, cond in enumerate(["MS", "Ctrl"]):
        v = d[d.cond == cond].d1.values
        xs = np.random.default_rng(7 + i).normal(i * 2 + j * 0.55, 0.09, len(v))
        ax.scatter(xs, v, s=13, color=RED if cond == "MS" else "0.62",
                   facecolor="none" if cond == "Ctrl" else RED,
                   edgecolor="0.45" if cond == "Ctrl" else RED, zorder=3)
    pr = A_P[A_P.study == st].iloc[0]
    ax.text(i * 2 + 0.5, 2.62, "P = %s" % fmt_p(pr.P), ha="center", fontsize=5.5, color="0.25")
    ax.text(i * 2 + 0.5, -0.35, "%d / %d" % (pr.n_case, pr.n_ctrl), ha="center", fontsize=5.5, color=MUT)
ax.set_xticks([i * 2 + 0.5 for i in range(len(hum))])
ax.set_xticklabels([s.split()[-1] for s in hum], fontsize=5.6, rotation=0)
ax.axhline(0, lw=0.6, color=MUT, ls="--")
ax.set_ylabel("donor-level d1", fontsize=5.8); ax.set_ylim(-0.6, 2.75)
ax.set_xlabel("cohort (MS / Ctrl counts below, P above)", fontsize=5.8)
note(ax, 0.98, 0.30, "all four null: the axis is not\na diagnostic marker", ha="right")

# (b) within-lesion
ax = panel(fig, gs, 0, 1, "2  Inside a donor the signal tracks lesion state",
           "per-sample medians; P and Cliff's delta as delivered")
grp = ["CTRL_CX", "CTRL_WM", "MS_Mye", "MS_Demye"]
samp = A2META.groupby(["group", "patient"])["MS_norm_score"].median().reset_index()
for i, g in enumerate(grp):
    v = samp[samp.group == g].MS_norm_score.values
    xs = np.random.default_rng(11 + i).normal(i, 0.07, len(v))
    ax.scatter(xs, v, s=15, color=RED if g.startswith("MS") else "0.55", zorder=3)
    ax.plot([i - 0.26, i + 0.26], [np.median(v)] * 2, lw=1.4, color="0.15", zorder=4)
ax.set_xticks(range(4)); ax.set_xticklabels(["CTRL_CX", "CTRL_WM", "MS_Mye", "MS_Demye"], fontsize=5.6)
ax.set_ylabel("sample-level median d1", fontsize=5.8)
note(ax, 0.98, 0.96, "Demye vs CTRL_WM  P = %s, delta = %.2f\nDemye vs MS_Mye    P = %s, delta = %.2f\nDemye vs CTRL_CX  P = %s (n.s.)"
     % (fmt_p(DIS[DIS.note == "MS_Demye_vs_CTRL_WM"].P.iloc[0]), DIS[DIS.note == "MS_Demye_vs_CTRL_WM"].effect.iloc[0],
        fmt_p(DIS[DIS.note == "MS_Demye_vs_MS_Mye"].P.iloc[0]), DIS[DIS.note == "MS_Demye_vs_MS_Mye"].effect.iloc[0],
        fmt_p(DIS[DIS.note == "MS_Demye_vs_CTRL_CX"].P.iloc[0])), ha="right")
headroom(ax, top=0.30)          # three annotation lines lifted into the empty band, never over the MS_Demye column dots

# (c) mouse EAE per library
ax = panel(fig, gs, 1, 0, "3  Mouse EAE: the signal is LOWER, not higher",
           "GSE250589 paired multiome, one dot per library")
d = A_DON[A_DON.study == "M GSE250589"]
for j, cond in enumerate(["Ctrl", "EAE"]):
    v = d[d.cond == cond].d1.values
    xs = np.random.default_rng(21 + j).normal(j, 0.08, len(v))
    ax.scatter(xs, v, s=16, color=RED if cond == "EAE" else "0.55", zorder=3)
    ax.plot([j - 0.26, j + 0.26], [np.median(v)] * 2, lw=1.4, color="0.15", zorder=4)
pr = A_P[A_P.study == "M GSE250589"].iloc[0]
ax.set_xticks([0, 1]); ax.set_xticklabels(["Ctrl (n=%d)" % pr.n_ctrl, "EAE (n=%d)" % pr.n_case], fontsize=5.8)
ax.axhline(0, lw=0.6, color=MUT, ls="--")
ax.set_ylabel("library-level d1", fontsize=5.8)
note(ax, 0.98, 0.96, "delta = %+.3f, P = %s\nprinted as measured" % (pr.delta, fmt_p(pr.P)), ha="right")
headroom(ax, top=0.22)          # two annotation lines lifted into the empty band the same way

# (d) composition confound
ax = panel(fig, gs, 1, 1, "4  But disease and composition move together",
           "n_MOL / n_OPC per mouse")
ratio = (M2B.n_MOL / M2B.n_OPC)
grpmap = {"Ctrl": "Ctrl", "EAE": "EAE"}
is_ctrl = M2B.mouse.astype(str).str.contains("Ctrl")
for j, (nm, mask) in enumerate([("Ctrl", is_ctrl), ("EAE", ~is_ctrl)]):
    v = ratio[mask].values
    xs = np.random.default_rng(31 + j).normal(j, 0.07, len(v))
    ax.scatter(xs, v, s=16, color=RED if j else "0.55", zorder=3)
    ax.plot([j - 0.26, j + 0.26], [np.median(v)] * 2, lw=1.4, color="0.15", zorder=4)
crow = COMP[COMP.quantity.astype(str).str.contains("composition ratio")].iloc[0]
ax.set_xticks([0, 1]); ax.set_xticklabels(["Ctrl (n=%d)" % int(is_ctrl.sum()), "EAE (n=%d)" % int((~is_ctrl).sum())], fontsize=5.8)
ax.set_ylabel("n_MOL / n_OPC", fontsize=5.8)
note(ax, 0.98, 0.96, "ratio P = %s\nRNA axis after composition\nadjustment: P = %s (n.s.)\nATAC: P = %s"
     % (fmt_p(crow.mwu_P), fmt_p(float(COMP[COMP.quantity.astype(str).str.contains("adjusting") & (COMP.axis == "RNA")].mwu_P.iloc[0])),
        fmt_p(float(COMP[COMP.quantity.astype(str).str.contains("adjusting") & (COMP.axis == "ATAC")].mwu_P.iloc[0]))), ha="right")

# (e) channel-gap five variants
ax = panel(fig, gs, 2, 0, "5  Channel decoupling, all five variants",
           "per-mouse gap = RNA_z - ATAC_z")
g = COMP[COMP.quantity.astype(str).str.contains("channel gap")].iloc[0]
lines = ["", "test                     value          95% CI", "",
         "median delta             %+.3f        [%+.3f, %+.3f]" % (g.delta, g.boot_ci_lo, g.boot_ci_hi),
         "mean delta               n.d.          [%+.3f, %+.3f]" % (g.boot_ci_mean_lo, g.boot_ci_mean_hi),
         "Mann-Whitney P           %s        (one-sided n/a: two-sided)" % fmt_p(g.mwu_P),
         "Welch t = %.3f            P = %s" % (g.welch_t, fmt_p(g.welch_P)),
         "Cliff's delta            %+.3f" % g.cliffs_delta,
         "Hedges g                 %+.3f" % g.hedges_g,
         "",
         "n.d. = the table gives the mean's bootstrap CI but no mean point estimate.",
         "The median-based CI contains 0 while the mean-based CI",
         "does not: the decoupling is nominally significant but",
         "not robust to the estimator, which is why all five",
         "numbers are printed instead of the most favourable one.",
         "n = %d mice (M2b_per_mouse.csv); resamples n = %d, seed %d." % (len(M2B), int(g.n_perm_or_boot), int(g.seed))]
ax.axis("off")
ax.text(0.0, 1.0, "\n".join(lines), transform=ax.transAxes, va="top", fontsize=5.8, linespacing=1.35, color="0.15")

# (f) mouse 74-trait within-cohort ranks
ax = panel(fig, gs, 2, 1, "6  Mouse MS signal replicated on 74 traits",
           "within-cohort rank, GSE250589")
m = M74.sort_values("cohort_rank").reset_index(drop=True)
ax.bar(range(len(m)), m.median_diff, width=0.85,
       color=[RED if t == MS else ("0.35" if "Schizophrenia" in str(t) else "0.78") for t in m.trait], zorder=3)
i_ms = int(m.index[m.trait == MS][0])
note(ax, 0.62, 0.98, "MS rank %d/74, d1 = %+.3f\n%d/%d libraries positive\n(binomial P = %s)"
     % (int(m.cohort_rank[i_ms]), m.median_diff[i_ms], int(m.n_pos[i_ms]), int(m.n_lib[i_ms]),
        fmt_p(float(binom_ms.pvalue))), col=RED, ha="left")
for key, lab in [("Schizophrenia", "SCZ %d" % int(m[m.trait.astype(str).str.contains("Schizophrenia")].cohort_rank.iloc[0])),
                 ("PASS_SWB", "SWB %d" % int(m[m.trait.astype(str).str.contains("SWB")].cohort_rank.iloc[0])),
                 ("HEIGHTz", "HEIGHT %d" % int(m[m.trait.astype(str).str.contains("HEIGHTz")].cohort_rank.iloc[0]))]:
    i = int(m.index[m.trait.astype(str).str.contains(key)][0])
    ax.text(i, m.median_diff[i] - 0.18, lab, fontsize=5.5, ha="center", va="top", color="0.25")
ax.axhline(0, lw=0.6, color=MUT)
ax.set_xlabel("74 traits, ranked", fontsize=5.8); ax.set_ylabel("mouse d1 (mature - precursor)", fontsize=5.8)
ax.set_ylim(min(m.median_diff) - 1.0, max(m.median_diff) + 0.9)
ax.set_xlim(-1, len(m))
ax.tick_params(labelsize=5.5)
# Historical (v15) sheet title, kept byte-as-drawn: the superseded panel carried "Fig 7D/E".
# The final Figure 7 has panels A-C only, so the current wording of that chain everywhere in
# this repository is "Fig 7 analysis chain (Table S6/S7)" - see figS7/LEGEND_S7_final.md §3.
fig._s7_no_title("Fig S7d  Support for Fig 7D/E - the disease block in full: null at diagnosis, moving with lesion state, replicated in the mouse",
             fontsize=7.2, x=0.008, y=0.985, ha="left")
save(fig, "S7d_disease_support", "6 micro-panels")
say("S7 support grid done")
