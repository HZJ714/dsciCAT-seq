# -*- coding: utf-8 -*-
"""Fig S6 **v11 候选**（2026-11-06）：只动 F 格——单元 bootstrap 梯新增 n_peak 调整级。

作者指令：fig6 主图不动；附图呈现去混杂的重要性。依据 DECONFOUND_TIERA_20261106 产物④§1。
F 格在 trait(n=74) 之后插入梯级 `trait n_peak adj`（0.708，CI 0.547–0.813，B=1,000
seed 20261104，残差化带截距）；其余六格一律复用 v10 交付 PDF（只读），几何逐字节不变。
新文件全部进 SUPPORT_v2/pdf_v11 与 tables_v11，不覆盖任何 v10 交付件。
"""
import os, sys

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, WORK + "/scripts")
from support_common import *
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from arial_style import use_arial
import fitz

R_OUT = R
SUP2 = R_OUT + "/SUPPORT_v2"
PDF2 = SUP2 + "/pdf"                    # v10 交付格：只读
PDF11 = SUP2 + "/pdf_v11"; os.makedirs(PDF11, exist_ok=True)
V8 = SUP2 + "/tables_v8"
V11 = SUP2 + "/tables_v11"
PNGDIR = os.environ.get("PNGOUT", "/tmp/sup2_preview"); os.makedirs(PNGDIR, exist_ok=True)
CMp = 2.54
RED = (0.79, 0.13, 0.18)
use_arial({"axes.titlesize": 6.0, "axes.labelsize": 5.8,
           "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "font.size": 5.8})
say("S6 v11 候选开始（仅 F 格重画）")

BOOT = pd.concat([pd.read_csv(V8 + "/figS6v10_bootstrap_units.csv"),
                  pd.read_csv(V11 + "/figS6v11_bootstrap_npeak.csv"),
                  pd.read_csv(V11 + "/figS6v11_bootstrap_nhit.csv")], ignore_index=True)
SUM = pd.read_csv(V11 + "/figS6v11_unit_summary_ext.csv").set_index("unit")
IDX = []
def put(pdf, key, val, printed, source):
    IDX.append(dict(pdf=pdf, key=key, value=val, printed=printed, source=source))

def newcell(w, h):
    fig = plt.figure(figsize=(w / CMp, h / CMp))
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    fig._s6v10 = [w, h]
    return fig

def savecell(fig, name):
    p = os.path.join(PDF11, name + ".pdf"); fig.savefig(p); plt.close(fig)
    d = fitz.open(p); w = d[0].rect.width * CMp / 72; h = d[0].rect.height * CMp / 72; d.close()
    say("[cell] %-26s %5.2f x %5.2f cm" % (name, w, h)); return p

def axcm(fig, x, y, w, h):
    W, H = fig._s6v10
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])

def clean(ax, xlab=None, ylab=None, tick=5.5, lab=5.8):
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
    for sp in ["left", "bottom"]: ax.spines[sp].set_linewidth(0.5)
    ax.tick_params(width=0.5, length=2.0, labelsize=tick)
    if xlab: ax.set_xlabel(xlab, fontsize=lab, labelpad=1.0)
    if ylab: ax.set_ylabel(ylab, fontsize=lab, labelpad=1.0)

# =============== 行4-F：单元 bootstrap 小提琴（v10 版式 + 新增 n_peak 梯级）===============
fig = newcell(8.60, 6.50); PF = "S6f_unit_violins_v11.pdf"
axF = axcm(fig, 1.00, 0.88, 7.50, 4.72)   # v11b：9 级，槽宽 0.833 cm
# 2026-11-06 作者复审：横轴标签改为"重抽单位"开头自解释（nucleus/library/mouse/trait），
# 新增 n_hit 级（深色 = 锁定规则家族），n_peak 级保持红色（敏感性口径）；
# GSE250589 独占一行且两侧邻居窄行（100k / raw），(n = 18)、(n = 74) 语义移交图注。
order = [("nucleus", "nucleus\n(descr.)", "0.75"),
         ("library_FY13", "library\n50k", "0.20"), ("library_FY14", "library\n100k", "0.20"),
         ("mouse", "mouse\nGSE250589", "0.20"),
         ("trait_O1", "trait\nraw", "0.20"),
         ("trait_O1_nhit", "trait\nn_hit\nadj", "0.20"),
         ("trait_O1_npeak", "trait\nn_peak\nadj", RED),
         ("trait_O1_wide", "trait\nwide\ncontrast", RED),
         ("trait_D_own", "trait\ndsciCAT\nraw", RED)]
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
for u, lab, _ in order:
    put(PF, "F.%s rho" % u, float(SUM.loc[u].rho_point), "%.3f" % SUM.loc[u].rho_point,
        V11 + "/figS6v11_unit_summary_ext.csv")
axF.axhline(0, color="0.80", lw=0.5, zorder=1)
axF.set_xlim(-0.6, len(order) - 0.4); axF.set_ylim(-0.08, 1.14)
axF.set_xticks(range(len(order))); axF.set_xticklabels([l for _, l, _ in order], fontsize=5.5)
axF.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
clean(axF, ylab="cross-channel Spearman rho\n(bootstrap, B = 1,000)", lab=5.6)
savecell(fig, "S6f_unit_violins_v11")

# =============== 拼页：A–E、G = v10 交付格（只读），F = v11 新格 ===============
import compose_support_S6S7 as CC
rows = [[(PDF2 + "/S6a_loci_tracks_v10.pdf", "A"), (PDF2 + "/S6b_permutation_strip_v10.pdf", "B")],
        [(PDF2 + "/S6c_umap_bins_v10.pdf", "C")],
        [(PDF2 + "/S6d_subtype_profiles_v10.pdf", "D"), (PDF2 + "/S6e_ladder_onepanel_v10.pdf", "E")],
        [(PDF11 + "/S6f_unit_violins_v11.pdf", "F"), (PDF2 + "/S6g_boxes_matrix_v10.pdf", "G")]]
CC.GAP = 0.12
CC.compose(rows, PDF11 + "/figS6_v11_candidate.pdf", PDF11 + "/figS6_v11_candidate.png")
pd.DataFrame(IDX).to_csv(V11 + "/figS6v11_F_number_index.csv", index=False)
d = fitz.open(PDF11 + "/figS6_v11_candidate.pdf")
say("v11 候选页幅 %.2f x %.2f cm（v10 交付 = 17.00 x 18.86 cm）"
    % (d[0].rect.width * CMp / 72, d[0].rect.height * CMp / 72)); d.close()
