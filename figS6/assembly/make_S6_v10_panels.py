# -*- coding: utf-8 -*-
"""Fig S6 **v10**（2026-11-04 夜，作者 20:5x 布局指令）——四行七格瘦身版。

行1 = A(位点轨) + B(窗口置换)   —— 原样保留（A 的逐基因命中清单进 Table S6）
行2 = C(分辨率 bin UMAP 三联+对照) —— 对应 Fig 6B 的分箱稳健性，原生重建（v6 素材字号 5.0 不过门）
行3 = D(亚型两级一致性，3 分面 + 顶部横排亚型图例) + E(置换零分布推断；跨化学 rho 6 值已进 Table S6)
行4 = F(单元 bootstrap 小提琴，梯级名用交付表 orig.ident：dsciCAT-seq-50k/100k) +
      G(逐核箱线 2行×3列 = 小点(每核) + 箱体，右下角双侧 Mann-Whitney U 的 P；
        Δ1 与 rho_adj/PASS 全部退回 Table S6 与图注)
F(双轨)/G(小鼠) 两格出图进表（Table S6 已有全部数值）。

纪律：面板字母由 compose 层统一 Arial 常规 8pt 大写；**图外不留灰色小字段落**（全部进
figS6_v10_legend_draft.md）；跨模态面板不画恒等线（审计 A1）；数字只从交付表/tables_v8 读，
每个上图数字登记 figS6v10_number_index.csv；可见字 >=5.5pt；Arial only。
"""
import os, sys, glob

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, WORK + "/scripts")
from support_common import *
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from arial_style import use_arial, FONTS_PDF
import fitz

R_OUT = R
SUP2 = R_OUT + "/SUPPORT_v2"
PDF2 = SUP2 + "/pdf"
V8 = SUP2 + "/tables_v8"
PNGDIR = os.environ.get("PNGOUT", "/tmp/sup2_preview"); os.makedirs(PNGDIR, exist_ok=True)
CMp = 2.54
use_arial({"axes.titlesize": 6.0, "axes.labelsize": 5.8,
           "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "font.size": 5.8})
say("S6 v10 开始（四行七格）")

STATE_NAMES = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
SUBPAL = {"COP/NFOL": "#D51F26", "MFOL": "#272E6A", "MOL2": "#208A42",
          "MOL3": "#89288F", "MOL3/6": "#F47D2B", "MOL6": "#FEE500"}
TECHPAL = {"10x Multiome": "#D55E00", "HT-scCAT-seq": "#009E73", "dsciCAT-seq (own)": "0.30"}
MS, SCZ, HEIGHT = ("PASS_Multiple_sclerosis", "PASS_Schizophrenia_Pardinas2018",
                   "UKB_460K.body_HEIGHTz")
NAMES = {MS: ("MS", "o"), SCZ: ("SCZ", "s"), HEIGHT: ("HEIGHT", "^")}
BLUERED = LinearSegmentedColormap.from_list("bluered", ["#9ECAE1", "#FFFFFF", "#E60000"])

SIX = pd.read_csv(V8 + "/figS6v8_sixpoint_all74.csv")
MED = pd.read_csv(V8 + "/figS6v8_state_median_74x6x2.csv")
LVAR = pd.read_csv(V8 + "/figS6v8_ladder_variants.csv")
UNIT = pd.read_csv(V8 + "/figS6v8_unit_ladder.csv")
DREP = pd.read_csv(V8 + "/figS6v8_D_replication.csv")
WPERM = pd.read_csv(V8 + "/figS6v8_window_permutation.csv").set_index("stat")
WNULL = pd.read_csv(V8 + "/figS6v8_window_permutation_null.csv")
BOOT = pd.read_csv(V8 + "/figS6v10_bootstrap_units.csv")
SUM = pd.read_csv(V8 + "/figS6v10_unit_summary.csv").set_index("unit")
NULLV = pd.read_csv(V8 + "/figS6v10_permnull_o1.csv")
NULS = pd.read_csv(V8 + "/figS6v10_permnull_summary.csv").set_index("stat")
CELL = read("O1_own_multiome_paired/tables/cell_level_table.csv", low_memory=False)
T2 = read("O1_own_multiome_paired/tables/T2_MS_per_nucleus_rna_vs_atac.csv")
L1 = read("O1_own_multiome_paired/tables/L1_by_library.csv").set_index("lib")
PERC = read("D_crosschemistry/tables_log1p_20261104/fig6D_delta1_74traits_per_tech.csv")
LAD = read("D_crosschemistry/tables_log1p_20261104/fig6D_crosschem_ladder.csv")
HIT = pd.read_csv(WORK + "/FY/oligo/scdrs_ms_out/ms_locus_overlap_table.csv")

def drep(fam, ch, pair, stat):
    q = DREP[(DREP.family == fam) & (DREP.channel == ch) & (DREP.pair == pair) & (DREP.stat == stat)]
    return float(q.iloc[0].value)

ARMS = ["dsciCAT-seq (own)", "10x Multiome", "HT-scCAT-seq"]
piv = {a: PERC[PERC.tech == a].set_index("trait") for a in ARMS}
RHO_XA = {(a, b, ch): drep("跨制库化学 同通道（同一性状 Δ1 两两一致）", ch.upper(),
                           "%s vs %s" % (a, b), "Spearman rho of Δ1")
          for a, b, ch in [("dsciCAT-seq (own)", "10x Multiome", "rna"), ("dsciCAT-seq (own)", "HT-scCAT-seq", "rna"),
                           ("dsciCAT-seq (own)", "10x Multiome", "atac"), ("dsciCAT-seq (own)", "HT-scCAT-seq", "atac"),
                           ("10x Multiome", "HT-scCAT-seq", "rna"), ("10x Multiome", "HT-scCAT-seq", "atac")]}
LADJ = {r.tech: (float(r.rho_adj), bool(r.ladder_reproduced), float(r.noise_ceiling))
        for _, r in LAD.iterrows()}
RHO_O1 = float(LVAR[(LVAR.family == "O1 同核跨通道梯子") & (LVAR.variant == "raw")].rho.iloc[0])
RHO_O1_WIDE = float(LVAR[LVAR.variant.astype(str).str.contains("sens")].rho.iloc[0])
RHO_DOWN = float(LVAR[LVAR.variant.astype(str).str.contains("rho_raw") &
                      (LVAR.family.astype(str).str.contains("D"))].rho.iloc[0])
RHO_NUC = float(T2[T2.trait.astype(str).str.contains("Multiple_sclerosis")].rho_raw.iloc[0])
RHO_MOUSE = float(UNIT[UNIT.level == "mouse"].rho.iloc[0])
LIBRHO = {"FY13": float(L1.loc["FY13", "rho_xchan"]), "FY14": float(L1.loc["FY14", "rho_xchan"])}
rho_of = {t: float(SIX[SIX.trait == t].rho_6pt.iloc[0]) for t in [MS, SCZ, HEIGHT]}
P_of = {t: float(SIX[SIX.trait == t].P_6pt.iloc[0]) for t in [MS, SCZ, HEIGHT]}
med_rho = float(SIX.rho_6pt.median())
MSR = {a: piv[a].loc[MS] for a in ARMS}

txt_rows = []
def dtxt(key, val, source, note=""):
    txt_rows.append(dict(key=key, value=val, source=source, note=note)); return val
dtxt("hexbin gridsize (display parameter)", "15 / 25", "matplotlib hexbin default comparison (Fig 6B uses 25)")
dtxt("bootstrap B / seed", "1000 / 20261104", V8 + "/figS6v10_unit_summary.csv")
dtxt("G 右下角统计量", "two-sided Mann-Whitney U on per-cell scDRS norm_score (OPC vs Oligo)",
     V8 + "/figS6v10_G_opc_oligo_mwu.csv")
dtxt("F 文库梯级身份", "FY13 = dsciCAT-seq-50k; FY14 = dsciCAT-seq-100k (交付表 orig.ident)",
     R_OUT + "/D_crosschemistry/tables_log1p_20261104/cells_dsciCAT_seq.csv")
IDX = []
def put(pdf, key, val, printed, source):
    IDX.append(dict(pdf=pdf, key=key, value=val, printed=printed, source=source))

def newcell(w, h):
    fig = plt.figure(figsize=(w / CMp, h / CMp))
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    fig._s6v10 = [w, h]
    return fig

def savecell(fig, name):
    p = os.path.join(PDF2, name + ".pdf"); fig.savefig(p); plt.close(fig)
    d = fitz.open(p); w = d[0].rect.width * CMp / 72; h = d[0].rect.height * CMp / 72; d.close()
    say("[cell] %-26s %5.2f x %5.2f cm" % (name, w, h)); return p

def axcm(fig, x, y, w, h):
    W, H = fig._s6v10
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])

def clean(ax, xlab=None, ylab=None, tick=5.5, lab=5.8):
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
    for sp in ["left", "bottom"]: ax.spines[sp].set_linewidth(0.5)
    ax.tick_params(width=0.4, length=1.8, pad=0.9, labelsize=tick)
    if xlab: ax.set_xlabel(xlab, fontsize=lab, labelpad=0.8)
    if ylab: ax.set_ylabel(ylab, fontsize=lab, labelpad=0.8)

CLIP = 1.5
def diag(ax, lo, hi):
    ax.plot([lo, hi], [lo, hi], color="0.55", lw=0.6, ls=(0, (2, 2)), zorder=1)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
def xmod(ax, lo, hi, x=0.03, y=0.97):
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.text(x, y, "cross-modality:\nrho only; amplitudes\nnot comparable",
            transform=ax.transAxes, fontsize=5.5, va="top", ha="left", color="0.45")

# =============== 行1-A：位点轨（同 v8）===============
import ms_locus_overlap as mlo
GENES = mlo.load_genes(WORK + "/FY/oligo/scdrs_ms_out/gencode.v45.basic.annotation.gtf.gz")
SEL4 = ["PHLDB1", "PHGDH", "BCAS1", "PLLP"]
fig = newcell(11.48, 3.00); PA = "S6a_loci_tracks_v10.pdf"
gA = fig.add_gridspec(1, 4, left=0.048, right=0.995, top=1 - 0.42 / 3.00, bottom=0.275, wspace=0.13)
for k, gname in enumerate(SEL4):
    h = HIT[HIT.gene == gname].iloc[0]
    ax = fig.add_subplot(gA[0, k])
    ch = str(h["chr"]); lo, hi = int(h.snp_pos) - 500_000, int(h.snp_pos) + 500_000
    for c2, t, nm in GENES:
        if c2 == ch and lo <= t <= hi and nm.upper() != gname:
            ax.plot([t, t], [0.42, 0.58], color="grey", lw=0.8)
    ax.plot([h.gene_tss, h.gene_tss], [0.42, 0.58], color=RED, lw=2.0)
    ax.plot([h.snp_pos, h.snp_pos], [0.70, 0.86], color="black", lw=1.2)
    ha = "left" if h.gene_tss > hi - 0.18 * (hi - lo) else ("right" if h.gene_tss < lo + 0.18 * (hi - lo) else "center")
    ax.plot([h.gene_tss, h.gene_tss], [0.40, 0.28], color=RED, lw=0.5)
    ax.text(h.gene_tss, 0.17, gname, fontsize=6.2, color=RED, fontweight="bold", ha=ha, va="center")
    ax.text(0.02, 0.98, "GWAS lead-SNP P = %.0e" % h.pvalue, transform=ax.transAxes,
            ha="left", va="top", fontsize=5.7)
    ax.text(0.02, 0.78, "%+.0f kb from gene" % (h.distance_bp / 1000.0), transform=ax.transAxes,
            ha="left", va="top", fontsize=5.7, color="0.35")
    put(PA, "A.%s lead-SNP P" % gname, float(h.pvalue), "P = %.0e" % h.pvalue,
        "FY/oligo/scdrs_ms_out/ms_locus_overlap_table.csv")
    ax.set_xlim(lo, hi); ax.set_ylim(0, 1.0)
    ax.set_xticks([lo, hi]); ax.set_xticklabels(["%.1f" % (lo / 1e6), "%.1f" % (hi / 1e6)], fontsize=5.8)
    for lb, al in zip(ax.get_xticklabels(), ["left", "right"]): lb.set_ha(al)
    ax.set_yticks([]); ax.set_title("%s  chr%s" % (gname, ch), fontsize=6.8, pad=1.5)
    for sp in ["top", "right", "left"]: ax.spines[sp].set_visible(False)
    for sp in ax.spines.values(): sp.set_linewidth(0.5)
    ax.tick_params(width=0.5, length=2)
savecell(fig, "S6a_loci_tracks_v10")

# =============== 行1-B：窗口置换（同 v8）===============
fig = newcell(5.40, 3.00); PB = "S6b_permutation_strip_v10.pdf"
axb = axcm(fig, 1.00, 0.44, 4.28, 1.96)
obs_pct = float(WPERM.loc["observed (%)", "value"])
axb.hist(WNULL.null_pct.values, bins=np.arange(6, 17, 0.5), color=(0.69, 0.769, 0.871),
         edgecolor=(0.28, 0.40, 0.56), linewidth=0.4, zorder=3)
axb.axvline(obs_pct, color=RED, lw=1.1, zorder=4)
axb.set_xlim(5.4, 16.6); axb.set_xticks([6, 8, 10, 12, 14, 16])
clean(axb, xlab="% of the 331 platform genes in window", ylab="permutations")
ptxt = "obs 38/331 = %.1f%%\nn = 10,000;  P = %.3f" % (
    obs_pct, float(WPERM.loc["permutation P (one-sided, observed vs null)", "value"]))
axb.text(0.985, 0.96, ptxt, transform=axb.transAxes, fontsize=5.5, va="top", ha="right", color="0.2")
for key, lab in [("platform genes in window (observed count)", "B.obs count"), ("observed (%)", "B.obs pct"),
                 ("permutation P (one-sided, observed vs null)", "B.perm P"),
                 ("permutation draws", "B.n perm"), ("Fisher P", "B.Fisher P"),
                 ("background window rate (%)", "B.background")]:
    put(PB, lab, float(WPERM.loc[key, "value"]), key, V8 + "/figS6v8_window_permutation.csv")
savecell(fig, "S6b_permutation_strip_v10")

# =============== 行2-C：分辨率 bin UMAP（旧 figS6c_umap_robust 版式：RNA 四格 + 单色标；
# 亚型着色图在 Fig 4，不再重复；ATAC 块删除——本行只回答 Fig 6B 分箱稳健性）===============
import anndata as ad
fig = newcell(17.00, 4.40); PC = "S6c_umap_bins_v10.pdf"
obs = ad.read_h5ad(WORK + "/FY/oligo/3k_dsci_oligo_subset.h5ad", backed="r").obs
dU = obs[["UMAP_1", "UMAP_2", "cell_type"]].join(
    CELL.set_index("cell_id")[["PASS_Multiple_sclerosis_rna", "PASS_Multiple_sclerosis_atac"]],
    how="inner")
dU = dU.astype({"PASS_Multiple_sclerosis_rna": float, "PASS_Multiple_sclerosis_atac": float})
say("[C] UMAP 核对 %d 核" % len(dU))
# 色标范围 = 旧 robust 条同一规则：原始分对称 99 百分位（±~2.5），RNA 四格共用
xr = max(np.nanpercentile(dU.PASS_Multiple_sclerosis_rna, 99),
         abs(np.nanpercentile(dU.PASS_Multiple_sclerosis_rna, 1)))
TX = [0.62, 4.52, 8.42, 12.32]; TW, TY, TH = 3.55, 0.62, 3.20
def umap_tile(x, ttl):
    a = axcm(fig, x, TY, TW, TH); a.set_title(ttl, fontsize=5.8, pad=1.6, loc="left")
    a.set_xticks([]); a.set_yticks([])
    a.set_xlabel("UMAP_1", fontsize=5.5, labelpad=0.6)
    for sp in a.spines.values(): sp.set_linewidth(0.5)
    return a
sp1 = dU.UMAP_1.max() - dU.UMAP_1.min(); sp2 = dU.UMAP_2.max() - dU.UMAP_2.min()
mol = dU[dU.cell_type == "MOL6"][["UMAP_1", "UMAP_2"]].median()
cop = dU[dU.cell_type == "COP/NFOL"][["UMAP_1", "UMAP_2"]].median()
ax1 = umap_tile(TX[0], "single cells (no smoothing)")
ax1.scatter(dU.UMAP_1, dU.UMAP_2, s=0.35, c=dU.PASS_Multiple_sclerosis_rna, cmap=BLUERED,
            vmin=-xr, vmax=xr, linewidths=0, rasterized=True)
ax2 = umap_tile(TX[1], "gridsize 15")
hb2 = ax2.hexbin(dU.UMAP_1, dU.UMAP_2, C=dU.PASS_Multiple_sclerosis_rna, reduce_C_function=np.mean,
                 gridsize=15, mincnt=1, cmap=BLUERED, linewidths=0, edgecolors="none", zorder=2)
hb2.set_clim(-xr, xr); hb2.set_rasterized(True)
ax3 = umap_tile(TX[2], "gridsize 25 (Fig 6B)")
hb3 = ax3.hexbin(dU.UMAP_1, dU.UMAP_2, C=dU.PASS_Multiple_sclerosis_rna, reduce_C_function=np.mean,
                 gridsize=25, mincnt=1, cmap=BLUERED, linewidths=0, edgecolors="none", zorder=2)
hb3.set_clim(-xr, xr); hb3.set_rasterized(True)
ax4 = umap_tile(TX[3], "gridsize 40")
hb4 = ax4.hexbin(dU.UMAP_1, dU.UMAP_2, C=dU.PASS_Multiple_sclerosis_rna, reduce_C_function=np.mean,
                 gridsize=40, mincnt=1, cmap=BLUERED, linewidths=0, edgecolors="none", zorder=2)
hb4.set_clim(-xr, xr); hb4.set_rasterized(True)
for a in (ax1, ax2):   # 旧版式：前两格带 MOL6 / COP-NFOL 箭头注释
    a.annotate("MOL6", xy=(mol.UMAP_1, mol.UMAP_2), xytext=(mol.UMAP_1 + 0.10 * sp1, mol.UMAP_2 + 0.30 * sp2),
               fontsize=5.8, fontweight="bold", color=RED, ha="center",
               arrowprops=dict(arrowstyle="->", color=RED, lw=0.7, shrinkB=1))
    a.annotate("COP/NFOL", xy=(cop.UMAP_1, cop.UMAP_2), xytext=(cop.UMAP_1 + 0.55 * sp1, cop.UMAP_2 + 0.38 * sp2),
               fontsize=5.8, fontweight="bold", color=(0.15, 0.35, 0.75), ha="center",
               arrowprops=dict(arrowstyle="->", color=(0.15, 0.35, 0.75), lw=0.7, shrinkB=1))
cbax = fig.add_axes([16.05 / 17.00, 0.88 / 4.40, 0.30 / 17.00, 2.60 / 4.40])
cbR = matplotlib.colorbar.ColorbarBase(cbax, cmap=BLUERED, norm=matplotlib.colors.Normalize(-xr, xr))
cbR.outline.set_linewidth(0.4); cbR.ax.tick_params(labelsize=5.5, length=1.6, width=0.4, pad=0.8)
cbR.set_label("norm_score", fontsize=5.5, labelpad=1.2)
put(PC, "C.gridsize 15", 15, "gridsize 15", V8 + "/figS6v10_figure_text.csv")
put(PC, "C.gridsize 25", 25, "gridsize 25 (Fig 6B)", V8 + "/figS6v10_figure_text.csv")
put(PC, "C.gridsize 40", 40, "gridsize 40", V8 + "/figS6v10_figure_text.csv")
savecell(fig, "S6c_umap_bins_v10")

# =============== 行3-D：亚型两级一致性 → 3 分面（一性状一分面；图例/ρ-P 全部进图注与 Table S6）===============
fig = newcell(9.40, 4.60); PD = "S6d_subtype_profiles_v10.pdf"
_FAL = fitz.Font(fontfile=FONTS_PDF["normal"])          # 图例横排用：量字宽
axLG = axcm(fig, 0.62, 0.12, 8.66, 0.42); axLG.axis("off")
_x = 0.0                                                # 作者 2026-11-05：图例回到图面，三面之上横排一行
for st in STATE_NAMES:
    _w = _FAL.text_length(st, fontsize=6.0) * 2.54 / 72 / 8.66
    axLG.scatter([_x + 0.010], [0.45], s=22, marker="o", color=SUBPAL[st], edgecolor="0.35",
                 linewidths=0.3, transform=axLG.transAxes, clip_on=False, zorder=4)
    axLG.text(_x + 0.030, 0.45, st, transform=axLG.transAxes, fontsize=6.0, ha="left", va="center",
              color="0.15")
    _x += 0.030 + _w + 0.026
FX, FY, FW, FG = 1.00, 1.15, 2.50, 0.18   # 左 1.00（旋转 ylabel）+ 3×2.50 + 2×0.18 + 右 0.54 = 9.40
for i, t in enumerate([MS, SCZ, HEIGHT]):
    nm, mk = NAMES[t]
    axD = axcm(fig, FX + i * (FW + FG), FY, FW, 2.70)
    xr = MED[(MED.channel == "RNA") & (MED.trait == t)].iloc[0]
    xa = MED[(MED.channel == "ATAC(gene score)") & (MED.trait == t)].iloc[0]
    for st in STATE_NAMES:
        axD.scatter(float(xr["med_" + st]), float(xa["med_" + st]), s=26, marker="o",
                    color=SUBPAL[st], edgecolor="0.35", linewidths=0.3, zorder=4)
    lim = 0.92
    axD.set_xlim(-lim, lim); axD.set_ylim(-lim, lim)
    axD.set_xticks([-0.8, -0.4, 0, 0.4, 0.8]); axD.set_yticks([-0.8, -0.4, 0, 0.4, 0.8])
    if i:                                   # 作者 2026-11-05：纵轴刻度只留最左一面
        axD.set_yticklabels([])
    axD.set_title(nm, fontsize=6.8, pad=1.4, fontweight="bold")
    if i == 1:                              # 三面轴标题相同 → x 标题只写中间一面
        axD.set_xlabel("median norm_score,\nRNA channel", fontsize=6.0, labelpad=0.8)
    if i == 0:
        axD.set_ylabel("median norm_score, ATAC channel", fontsize=6.0, labelpad=0.8)
    for sp in ["top", "right"]: axD.spines[sp].set_visible(False)
    for sp in ["left", "bottom"]: axD.spines[sp].set_linewidth(0.5)
    axD.tick_params(width=0.4, length=1.8, pad=0.9, labelsize=6.0)
savecell(fig, "S6d_subtype_profiles_v10")

# =============== 行3-E：推断面板（置换零分布，不再重画 Fig 6C 散点；跨化学 6 值已进 Table S6）===============
fig = newcell(7.40, 4.60); PE = "S6e_ladder_onepanel_v10.pdf"
axE = axcm(fig, 1.00, 0.50, 6.00, 3.30)
obs_r = float(NULS.loc["observed cross-channel rho", "value"])
p_per = float(NULS.loc["permutation P (one-sided)", "value"])
nl975 = float(NULS.loc["null 97.5 pct", "value"])
nlmen = float(NULS.loc["null mean", "value"])
nper = int(NULS.loc["permutation draws", "value"])
axE.hist(NULLV.rho_perm.values, bins=np.arange(-0.5, 0.92, 0.04), color=(0.69, 0.769, 0.871),
         edgecolor=(0.28, 0.40, 0.56), linewidth=0.4, zorder=3)
axE.axvspan(nl975, 0.85, color="0.93", lw=0, zorder=0)   # 超出 null 97.5% 的尾部阴影（观测所在区）
axE.axvline(nl975, color="0.35", lw=0.8, ls="--", zorder=4)
axE.axvline(obs_r, color=RED, lw=1.1, zorder=5)
axE.text(0.975, 0.97, "observed\n%+.3f" % obs_r, transform=axE.transAxes, fontsize=6.2,
         fontweight="bold", ha="right", va="top", color=RED, linespacing=1.25)
axE.text(0.02, 0.97, "P < 0.0001\n(0 / 10,000 draws)", transform=axE.transAxes, fontsize=6.0,
         ha="left", va="top", color="0.15", linespacing=1.25)
axE.text(0.88, 0.30, "null mean %+.3f\nnull 97.5%% = %.3f" % (nlmen, nl975),
         transform=axE.transAxes, fontsize=5.8, ha="right", va="top", color="0.35", linespacing=1.25)
axE.set_xlim(-0.5, 0.9); axE.set_xticks([-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8])
axE.set_ylim(0, 1550); axE.set_yticks([0, 500, 1000])
clean(axE, xlab="cross-channel Spearman rho, trait pairing permuted", ylab="number of permutations",
      tick=6.0, lab=6.2)
put(PE, "E.obs rho", obs_r, "%+.3f" % obs_r, V8 + "/figS6v10_permnull_summary.csv")
put(PE, "E.perm P", p_per, "%.4f" % p_per, V8 + "/figS6v10_permnull_summary.csv")
put(PE, "E.null 975", nl975, "%.3f" % nl975, V8 + "/figS6v10_permnull_summary.csv")
put(PE, "E.null mean", nlmen, "%.3f" % nlmen, V8 + "/figS6v10_permnull_summary.csv")
put(PE, "E.n perm", nper, "{:,}".format(nper), V8 + "/figS6v10_permnull_summary.csv")
savecell(fig, "S6e_ladder_onepanel_v10")

# =============== 行4-F：单元 bootstrap 小提琴（替代点梯）===============
fig = newcell(8.60, 6.50); PF = "S6f_unit_violins_v10.pdf"
axF = axcm(fig, 1.05, 0.88, 7.30, 4.72)   # 底部留 0.88cm：mouse 梯级三行标签（GSE 号上图面）
# 作者 2026-11-05：文库梯级不许写内部库号 FY13/FY14，改用交付表 orig.ident 的技术名；
# “own” 一词一律不出现（D-own raw → trait / dsciCAT / raw）。
order = [("nucleus", "dsciCAT-seq\nnucleus\n(descriptive)", "0.75"),
         ("library_FY13", "dsciCAT-\nseq-50k", "0.20"), ("library_FY14", "dsciCAT-\nseq-100k", "0.20"),
         # 作者 2026-11-05：mouse 梯级是外部数据 GSE250589（Fig 7 分析链同一批，Table S6/S7），图面必须带号防误读成自家 dsciCAT
         ("mouse", "mouse EAE\nGSE250589\n(n = 18)", "0.20"), ("trait_O1", "trait\n(n = 74)", "0.20"),
         # 作者 2026-11-05：内部分析编号 O1 不上图面（与 FY13/FY14 同规），数据源由图顶注释行统一声明
         ("trait_O1_wide", "trait wide", RED), ("trait_D_own", "trait\ndsciCAT\nraw", RED)]
pts, labels, cols = [], [], []
for i, (u, lab, col) in enumerate(order):
    v = BOOT[BOOT.unit == u].rho.values
    vp = axF.violinplot(v, positions=[i], widths=0.72, showmedians=False, showextrema=False)
    fc = col if col in ("0.20", RED) else "0.75"
    for b in vp["bodies"]:
        b.set_facecolor(fc); b.set_alpha(0.35 if col != "0.75" else 0.9); b.set_edgecolor(fc)
    b = vp["bodies"][0]; b.set_linewidth(0.5)
    q = SUM.loc[u]
    axF.plot([i, i], [q.ci_lo, q.ci_hi], color=fc, lw=0.8, zorder=4)
    axF.plot([i - 0.10, i + 0.10], [q.ci_lo] * 2, color=fc, lw=0.7, zorder=4)
    axF.plot([i - 0.10, i + 0.10], [q.ci_hi] * 2, color=fc, lw=0.7, zorder=4)
    axF.scatter([i], [q.boot_median], s=12, color=fc, edgecolor="0.15", linewidths=0.5, zorder=5)
    axF.text(i, 1.045, "%.3f" % q.rho_point, ha="center", va="bottom", fontsize=5.5,
             color="0.15", fontweight="bold", rotation=0)
    pts.append(q.rho_point)
for u, lab, _ in order:
    put(PF, "F.%s rho" % u, float(SUM.loc[u].rho_point), "%.3f" % SUM.loc[u].rho_point,
        V8 + "/figS6v10_unit_summary.csv")
axF.axhline(0, color="0.80", lw=0.5, zorder=1)
axF.set_xlim(-0.6, len(order) - 0.4); axF.set_ylim(-0.08, 1.14)
axF.set_xticks(range(len(order))); axF.set_xticklabels([l for _, l, _ in order], fontsize=5.5)
axF.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
clean(axF, ylab="cross-channel Spearman rho\n(bootstrap, B = 1,000)", lab=5.6)
# 作者 2026-11-05：数据源两行灰字注释从图面去掉，只留在图注（F 段首句）
savecell(fig, "S6f_unit_violins_v10")

# =============== 行4-G：逐核箱线 2行×3列矩阵（v6-G 版式 + log1p 修正数字）===============
fig = newcell(8.20, 6.50); PG = "S6g_boxes_matrix_v10.pdf"
CELLF = {"dsciCAT-seq (own)": R_OUT + "/D_crosschemistry/tables_log1p_20261104/cells_dsciCAT_seq.csv",
         "10x Multiome": R_OUT + "/D_crosschemistry/tables_log1p_20261104/cells_10x.csv",
         "HT-scCAT-seq": R_OUT + "/D_crosschemistry/tables_log1p_20261104/cells_HT_scCAT_seq_HT.csv"}
DIRF = {"dsciCAT-seq (own)": ("dsciCAT_seq_rna", "dsciCAT_seq_atac"),
        "10x Multiome": ("10x_rna", "10x_atac"), "HT-scCAT-seq": ("HT_scCAT_seq_HT_rna", "HT_scCAT_seq_HT_atac")}
def ms_scores(tech, chan):
    rn, an = DIRF[tech]
    base = R_OUT + "/D_crosschemistry/scdrs/" + rn if chan == "rna" else R_OUT + "/D_crosschemistry/scdrs_log1p_20261104/" + an
    fs = sorted(glob.glob(base + "/*/PASS_Multiple_sclerosis.full_score.gz"))
    dfs = [pd.read_csv(f, sep="\t", index_col=0, low_memory=False) for f in fs]
    df = pd.concat(dfs); df.index.name = "cell"
    lb = pd.read_csv(CELLF[tech]).set_index("cell")["label"]
    df = df[df.index.isin(lb.index)].copy()
    df["label"] = lb.loc[df.index].values
    return df
gG = fig.add_gridspec(2, 3, left=0.13, right=0.985, top=1 - 0.42 / 6.50, bottom=0.42 / 6.50,
                      wspace=0.13, hspace=0.20)
# p_fig/p_label 必须按**字符串**读：pandas 默认会把 "3.1e-04" 解析成浮点，再 str() 就变成
# "0.00031"，台账与图面写法对不上（反向门靠数值容差仍会放过，但逐字对不上 = 台账失真）。
MWU = pd.read_csv(V8 + "/figS6v10_G_opc_oligo_mwu.csv", dtype={"p_fig": str, "p_label": str})
ARM_LAB = {"dsciCAT-seq (own)": "dsciCAT-seq", "10x Multiome": "10x Multiome",
           "HT-scCAT-seq": "HT-scCAT-seq"}     # 作者 2026-11-05：图面一律不出现 "(own)" 这类内部叫法
rngG = np.random.default_rng(20261105)
# 观感回退到旧版 v20_D_two（fig6v20_D_panel.py 手绘箱）：暗描边 + 粗中线 + 2–98 须 +
# 正态抖动的散点**垫在箱体后面**；作者 2026-11-05：dsci 点多，散点颜色/大小单独降密。
def _dark(c, thr=0.55):
    rgb = np.array(matplotlib.colors.to_rgb(c))
    lum = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
    return tuple(rgb * (thr / lum)) if lum > thr else tuple(rgb)
def _lighten(c, f):                        # 向白混色：雾化的密集点云
    rgb = np.array(matplotlib.colors.to_rgb(c))
    return tuple(rgb * f + (1 - f))
GPT = {"dsciCAT-seq (own)": dict(s=0.9, alpha=0.20, lighten=0.55),   # 7.3k 点/臂：更小更雾
       "10x Multiome": dict(s=1.6, alpha=0.30, lighten=0.0),          # ~0.8k：旧版原尺寸
       "HT-scCAT-seq": dict(s=1.6, alpha=0.30, lighten=0.0)}          # ~0.5k
for ri, ch in enumerate(["rna", "atac"]):
    x = 0
    for ci, tech in enumerate(ARMS):
        ax = fig.add_subplot(gG[ri, ci])
        df = ms_scores(tech, ch)
        mo = df[df.label == "Oligo"].norm_score.values; op = df[df.label == "OPC"].norm_score.values
        st = GPT[tech]
        for pos, arr, col in [(0, op, OPC), (0.6, mo, OLIGO)]:        # 每核小点（旧版观感）
            ax.scatter(pos + rngG.normal(0, 0.07, len(arr)), arr, s=st["s"],
                       color=_lighten(col, 1 - st["lighten"]) if st["lighten"] else col,
                       alpha=st["alpha"], linewidths=0, rasterized=True, zorder=1)
        bp = ax.boxplot([op, mo], positions=[0, 0.6], widths=0.5, patch_artist=True, showfliers=False,
                        manage_ticks=False, whis=[2, 98],
                        medianprops=dict(lw=1.3), boxprops=dict(lw=0.8),
                        whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7))
        for ln in bp["caps"]:                       # 窄帽（旧版 v20_D 的 ±0.16/±0.30 比例→此处 ±0.13）
            c = float(np.mean(ln.get_xdata())); y = float(ln.get_ydata()[0])
            ln.set_data([c - 0.13, c + 0.13], [y, y])
        for b, col in zip(bp["boxes"], [OPC, OLIGO]):
            b.set_facecolor(col); b.set_alpha(0.30); b.set_edgecolor(_dark(col)); b.set_zorder(2)
        for ln, col in zip(bp["medians"], [OPC, OLIGO]): ln.set_color(_dark(col)); ln.set_zorder(3)
        for key in ("whiskers", "caps"):                       # 须/帽按各自箱位着色（左 OPC 橙 / 右 OLIGO 蓝）
            for ln in bp[key]: ln.set_color(OLIGO if ln.get_xdata()[0] > 0.2 else OPC); ln.set_zorder(2)
        # 右下角：原 Δ1 → 双侧 Mann-Whitney U 的精确 P（Δ1 退回 Table S6；方法写图注）
        _r = MWU[(MWU.tech == ARM_LAB[tech]) & (MWU.channel == ch)].iloc[0]
        p_raw, p_txt = float(_r.p_two_sided), str(_r.p_label)
        ax.text(0.97, 0.05, p_txt, fontsize=5.8, ha="right", va="bottom", color="0.15",
                fontweight="bold", transform=ax.transAxes)
        put(PG, "G.MS MWU P %s %s" % (ch, ARM_LAB[tech]), p_raw, str(_r.p_fig),
            V8 + "/figS6v10_G_opc_oligo_mwu.csv")
        ax.set_ylim(-1.62, 1.4); ax.set_xlim(-0.5, 1.1)
        ax.axhline(0, color="0.80", lw=0.5, zorder=1)
        ax.set_xticks([0, 0.6]); ax.set_xticklabels(["OPC", "Oligo"], fontsize=5.5)
        ax.set_yticks([-1, 0, 1])
        if ci:                                   # 作者 2026-11-05：与 D 一致，刻度数字只留最左一列
            ax.set_yticklabels([])
        for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
        for sp in ["left", "bottom"]: ax.spines[sp].set_linewidth(0.5)
        ax.tick_params(width=0.4, length=1.8, pad=0.9, labelsize=5.5)
        if ri == 0: ax.set_title(ARM_LAB[tech], fontsize=6.0, pad=1.8)
        if ci == 0:
            # 旋转单行标签的长度就是它的横向占位（>面高就会出画），所以 (MS) 只写进图注
            ax.set_ylabel("per-cell norm_score, %s" % ("RNA" if ch == "rna" else "ATAC"),
                          fontsize=5.6, labelpad=1.0)
# 底部 rho_adj / PASS 文字行已按作者 2026-11-05 指令出图：三个 rho_adj 与判据进 Table S6 + 图注/正文。
savecell(fig, "S6g_boxes_matrix_v10")

# =============== 索引落表 + 拼页（字母 = compose 层 Arial 常规 8pt）===============
pd.DataFrame(txt_rows).to_csv(V8 + "/figS6v10_figure_text.csv", index=False)
pd.DataFrame(IDX).to_csv(V8 + "/figS6v10_number_index.csv", index=False)
say("[index] %d 个上图数字已登记 figS6v10_number_index.csv" % len(IDX))
import compose_support_S6S7 as CC
rows = [[(PDF2 + "/S6a_loci_tracks_v10.pdf", "A"), (PDF2 + "/S6b_permutation_strip_v10.pdf", "B")],
        [(PDF2 + "/S6c_umap_bins_v10.pdf", "C")],
        [(PDF2 + "/S6d_subtype_profiles_v10.pdf", "D"), (PDF2 + "/S6e_ladder_onepanel_v10.pdf", "E")],
        [(PDF2 + "/S6f_unit_violins_v10.pdf", "F"), (PDF2 + "/S6g_boxes_matrix_v10.pdf", "G")]]
CC.GAP = 0.12
CC.compose(rows, PDF2 + "/figS6_v10_draft.pdf", PDF2 + "/figS6_v10_draft.png")
# 作者 2026-11-05：右下角 "data: log1p-corrected 2026-11-04" 角标取消（像水印），数据版本只写图注。
say("S6 v10 完成：figS6_v10_draft.pdf")
