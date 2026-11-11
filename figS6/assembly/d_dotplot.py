# -*- coding: utf-8 -*-
"""d_dotplot.py — Fig6D 身份验证气泡图 v3（回执 P0-8b 2，Seurat DotPlot 体例，作者定 17 行）

v2→v3（作者指定）：
  * 20 行 → 17 行，清单由作者定：Pdgfra Tcf7l2 | Bmp4 Gpr17 | Mal Plp1 Apod Aspa | Sox6 Nfix Mef2c Creb5 | 风险 5。
    直接动机：MFOL 列原本 7 行达峰、整列被占满；删 Mbp/Mog/Opalin 后只剩 4 行，六列各有归属。
    Pdgfra 被捐回（尽管核内只有 3.4% 阳性），其说明从 OMITTED 移到图注的绝对数字。
  * 行高 0.30 → 0.225 cm，图例并成一行，目标与 6C（高 5.35 cm）同高。点最大直径 1.90 mm（< 行高）。

v1→v2 的改动（作者要求：基因别太多，要塞得进主图）：
  * 27 行 → 20 行。髓鞘组 12→6，删 Mag/Mobp/Car2/Ndrg1/Qk/Myrf/S1pr1（共线、峰位与组名矛盾、或核内不可测）。
  * 新增 stage-TF 组 4 个，直接复现 Fig5a 的 early/mid/late 命名：
      Sox6 (early, COP 峰 15x) / Nfix (mid, MOL2) / Mef2c (mid, MOL3 峰, Fig5b anchor) / Creb5 (late, MOL6 方向)
    —— Mef2c 是作者点名要加的，与 Fig5 呼应；图注必须带 FigS5b 原话的限定语
       "Mef2c peaks in MOL3 at 0.285 (relative, not strong, specificity)"。
  * 髓鞘组补 Aspa：全表唯一真正在 **MOL6** 达峰的髓鞘相关基因（16.5x，pct 75%）——
    对角块右端靠它撑住。Plp1/Mbp 保留但它们是名气锚点（1.3x/1.5x，pct 96–100%，无区分力）。
  * 删 Pdgfra（核内 mean 0.047 / 3.4% 阳性，画出来是针尖），剔除理由写进 meta 与图注，不隐瞒。
  * 删 Tbx2：FigS5b 的 0.512 是 **eRegulon 活性**，不是基因表达（核内 0–1%），两者不能混在一张表达图上。
  * 美化：行高 0.30 cm；组间实线 + 左侧竖括号（§8.2 不放行头文字）；每行峰值单元格淡灰方框
    （对角块从"要看出"变"被指出"）；色阶改全表 1/99 分位定限；risk 组行序改回贡献降序，
    并在右侧加贡献条（Δ gap upon gene removal），使本图自带原 D 柱状图的信息。
  * 命名口径：lineage marker 用**小鼠斜体**（实测的就是小鼠核），risk 基因用**人类大写**
    （它们来自人 GWAS 基因集）；两套并存的理由写图注，见 meta["naming"]。

口径（验收硬条件，与 v1 相同）：
  * 细胞 = 与 Fig6A 同一批：scDRS RNA∩ATAC 分数文件同核交集 5,962
    （203/490/1016/688/991/2574，与 fig6_out/fig6A_stats.csv 逐一相同；h5ad 本身 5,977，用错即不同批）。
  * 表达 = 3k_dsci_oligo_subset.h5ad，实测 expm1(X) 每核和恒为 10,000 → log1p(CP10k)（Seurat RNA$data）。
  * 人→鼠一律 scDRS 内置同源表（SKILL §15.6）。RNA-only。点大小=阳性核比例，点色=行内跨亚型 z。
  * 组 5 = LOGO 贡献降序 ∩ 平均阳性率 ≥10%；禁用 "driver"。
"""
import os, sys, json
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad

WORK = os.environ.get('DSCI_WORK_ROOT', '')
sys.path.insert(0, os.path.join(WORK, "GitHub/figS6/assembly"))
os.environ.setdefault("DSCI_WORK_ROOT", WORK)
from arial_style import use_arial, use_two_tier           # noqa: E402
use_arial(); use_two_tier()
import matplotlib                                         # noqa: E402
import matplotlib.pyplot as plt                           # noqa: E402
import fitz                                               # noqa: E402

CM = 1.0 / 2.54
OUT = WORK + "/figures/routeA_out/DRAFT_v22"
TAB, PAN = OUT + "/tables", OUT + "/panels"
os.makedirs(TAB, exist_ok=True); os.makedirs(PAN, exist_ok=True)
SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
RNA_H5 = WORK + "/FY/oligo/3k_dsci_oligo_subset.h5ad"
RROOT = WORK + "/FY/oligo/scdrs_traits_draft/full74"
AROOT = WORK + "/FY/oligo/fig6_rerun_20261019/dualomics_probe/out_gs74"
TRAIT = "PASS_Multiple_sclerosis"
TARGET_W = 10.72          # Fig6 行2 右列（17.81 − 6C 6.85 − GAP 0.24）
PCT_MIN = 10.0
N_RISK = 5
RH = 0.225                # 行高 cm（作者要求：行距再小一点，尽量与 6C 同高）
DOT_MM = 1.90             # 100% 阳性时的点直径（必须 < 行高，否则上下行相碰）

# v3 基因清单（作者指定 17 个：删 Cspg4/Vcan/Opalin/Mbp/Mog，捞回 Pdgfra/Bmp4）
# 理由：MFOL 列原本 7 行达峰（Opalin/Mal/Mbp/Mog/Plp1/Creb5/HDAC7），视觉上被 MFOL 占满。
# 删后 MFOL 只剩 Mal/Plp1/Creb5/HDAC7 四行，六列各有归属。
MARKERS = [
    # 作者意见：COP/NFOL 严格说不是 OPC，只是共享部分 OPC 基因，所以 precursor 与
    # intermediate 合成一个模块，避免读者把第一列读成“这就是 OPC”
    ("OPC / intermediate", ["Pdgfra", "Tcf7l2", "Bmp4", "Gpr17"]),
    ("mature myelin",   ["Mal", "Plp1", "Apod", "Aspa"]),
    ("stage TF",        ["Sox6", "Nfix", "Mef2c", "Creb5"]),
]
OMITTED = {  # 主动剔除，理由进 meta + 图注（§14.6 不许把疑问留给读者）
    "S1pr1": "pharmacological target, not a myelin-structure gene; nuclear 0.02 / 2.5% (see Fig. 7D)",
    "Qk": "peak COP/NFOL, 1.6x, 100% positive -> no discriminating power inside the myelin group",
    "Mag": "collinear with Mog/Mobp (same MFOL peak)",
    "Mobp": "collinear with Mog (same MFOL peak)",
    "Car2": "same MOL2 peak as Apod, lower detection (52% vs 81%)",
    "Ndrg1": "low abundance (nuclear mean 0.33)",
    "Myrf": "peak COP/NFOL contradicts the 'mature myelin' row label",
    "Cspg4": "dropped by author for row budget; same COP/NFOL peak as Tcf7l2",
    "Vcan": "dropped by author for row budget; same COP/NFOL peak as Tcf7l2",
    "Opalin": "dropped by author: MFOL column was over-represented (7 peaks)",
    "Mbp": "dropped by author: MFOL over-representation + 1.5x only, 92-100% positive (no discrimination)",
    "Mog": "dropped by author: MFOL over-representation, collinear with Mal",
    "Tbx2": "Fig. S5b 0.512 is eRegulon ACTIVITY, not expression; nuclear RNA 0-1% -> cannot sit on an expression dot plot",
}

sys.path.insert(0, os.path.expanduser('"~/.local/lib/python3.9/site-packages'))
from scdrs.util import load_homolog_mapping               # noqa: E402


def same_nucleus_cells():
    def idx_of(root):
        import glob
        for f in sorted(glob.glob(root + "/*/" + TRAIT + ".full_score.gz")):
            try:
                return pd.read_csv(f, sep="\t", index_col=0).index
            except OSError:
                continue
        raise RuntimeError("找不到 " + TRAIT + ".full_score.gz in " + root)
    return idx_of(RROOT).intersection(idx_of(AROOT))


# ---------------------------------------------------------------- 取数
cells = same_nucleus_cells()
a = ad.read_h5ad(RNA_H5)
key = (pd.Series(a.obs["cell_type"].astype(str)).str.upper()
       .str.replace("/", "", regex=False).str.replace("-", "", regex=False)).values
mask = a.obs_names.isin(set(cells))
grp = {s: np.where(mask & (key == s.replace("/", "")))[0] for s in SUB}
pos = {str(v): i for i, v in enumerate(a.var_names)}
order = np.concatenate([grp[s] for s in SUB])
bnd = np.cumsum([0] + [len(grp[s]) for s in SUB])
cnt = "/".join(str(len(grp[s])) for s in SUB)
print("细胞口径: 5,962 同核 → %s（fig6A_stats.csv: 203/490/1016/688/991/2574）" % cnt)
assert cnt == "203/490/1016/688/991/2574", "细胞口径与 Fig6A 不同批！"

_cache = {}
def stat(mo):
    if mo in _cache:
        return _cache[mo]
    i = pos.get(mo)
    if i is None:
        return None
    col = a.X[:, i]
    col = (col.toarray().ravel() if sp.issparse(col) else np.asarray(col).ravel())[order]
    mn = [float(col[bnd[k]:bnd[k + 1]].mean()) for k in range(6)]
    pc = [float((col[bnd[k]:bnd[k + 1]] > 0).mean() * 100) for k in range(6)]
    _cache[mo] = (mn, pc)
    return _cache[mo]

rows = []
for gname, gl in MARKERS:
    for mo in gl:
        s = stat(mo)
        if s is None:
            print("  ! 矩阵缺 marker", mo); continue
        mn, pc = s
        rows.append(dict(group=gname, gene=mo, gene_human=mo, contribution=np.nan))
m = load_homolog_mapping("hsapiens", "mmusculus")
L = pd.read_csv(WORK + "/figures/routeA_out/GENELEVEL_v1/logo_delta1_contributions.csv")
risk, excl = [], []
for _, r in L.iterrows():
    mo = m.get(r.gene)
    if mo is None:
        excl.append(dict(gene=r.gene, why="no ortholog in scDRS built-in table")); continue
    s = stat(mo)
    if s is None:
        excl.append(dict(gene=r.gene, why="mouse symbol absent from RNA matrix")); continue
    mn, pc = s
    if np.mean(pc) < PCT_MIN:
        excl.append(dict(gene=r.gene, why="mean %.1f%% positive < %.0f%%" % (np.mean(pc), PCT_MIN))); continue
    risk.append(dict(group="MS risk contribution", gene=mo, gene_human=r.gene, contribution=float(r.contribution)))
    if len(risk) >= N_RISK:
        break
rows += risk
print("组 5 入选 %d: %s | 前 %d 名内被门槛/无同源剔除 %d"
      % (len(risk), ", ".join(x["gene_human"] for x in risk), len(risk) + len(excl), len(excl)))

M = pd.DataFrame(rows)
M["gorder"] = M.group.map({g[0]: i for i, g in enumerate(MARKERS)}).fillna(len(MARKERS)).astype(int)
# 组内排序：marker 组按峰值亚型（成对角块），risk 组按贡献降序（与右侧条一致）
M["peak_idx"] = [int(np.argmax(stat(g)[0])) for g in M.gene]
M = M.sort_values(["gorder", "peak_idx", "gene"]).reset_index(drop=True)
for gi in M.gorder.unique():                     # 组内重排
    sel = M.index[M.gorder == gi]
    if M.loc[sel, "group"].iloc[0] == "MS risk contribution":
        M.loc[sel, "ord2"] = -M.loc[sel, "contribution"]
    else:
        M.loc[sel, "ord2"] = M.loc[sel, "peak_idx"]
M = M.sort_values(["gorder", "ord2"]).reset_index(drop=True)

MEAN = np.array([stat(g)[0] for g in M.gene], float)
PCT = np.array([stat(g)[1] for g in M.gene], float)
Zc = (MEAN - MEAN.mean(1, keepdims=True)) / (MEAN.std(1, keepdims=True) + 1e-9)
VQ = (float(np.percentile(Zc, 1)), float(np.percentile(Zc, 99)))     # 全表 1/99 分位定限
print("色阶限 = %.2f / %.2f（全表 1/99 分位）" % VQ)

# ---------------------------------------------------------------- 画图
n = len(M)
TRANSPOSE = os.environ.get("DOT_TRANSPOSE", "0") == "1"
DOT_MM_EFF = DOT_MM if not TRANSPOSE else 2.60      # 横版行高 0.45 cm，点放到 2.6 mm
S_PER_PCT = np.pi * (DOT_MM_EFF / 25.4 * 72 / 2) ** 2 / 100.0
LEG = [20, 50, 85]
RISK_GI = len(MARKERS)
sel_risk = M.index[M.gorder == RISK_GI]

if not TRANSPOSE:
    # ---- 竖版：基因在纵轴 17 行，亚型在横轴 6 列
    BOT_CM = 2.05                     # 旋转亚型标签 1.30 + 圆点图例 0.60 + 间隙 0.15
    HH = RH * n + BOT_CM
    fig = plt.figure(figsize=(TARGET_W * CM, HH * CM))
    LM, AXW = 0.175, 0.700          # 贡献条移除后主图加宽，列间距 1.25 cm
    AX_BOT, AX_H = BOT_CM / HH, 1 - (BOT_CM + 0.40) / HH
    ax = fig.add_axes([LM, AX_BOT, AXW, AX_H])
    xs, ys = np.meshgrid(np.arange(6), np.arange(n))
    for i in range(n):                # 每行峰值单元格淡灰方框
        pk = int(np.argmax(Zc[i]))
        ax.add_patch(matplotlib.patches.Rectangle((pk - 0.34, i - 0.30), 0.68, 0.60,
                     fill=False, edgecolor="0.55", lw=0.40, zorder=2))
    ax.scatter(xs.ravel(), ys.ravel(), s=(np.clip(PCT, 0, None) * S_PER_PCT).ravel(),
               c=Zc.ravel(), cmap="coolwarm", vmin=VQ[0], vmax=VQ[1],
               edgecolors="0.25", linewidths=0.28, zorder=3)
    ax.set_xlim(-0.55, 5.55); ax.set_ylim(-0.62, n - 0.38); ax.invert_yaxis()
    ax.set_xticks(range(6)); ax.set_xticklabels(SUB, fontsize=6.0, rotation=45, ha="right")
    ax.set_yticks(range(n)); ax.set_yticklabels(M.gene_human, fontsize=6.0)
    for lb, g in zip(ax.get_yticklabels(), M.group):      # 小鼠斜体 / 人类正体
        lb.set_style("italic" if g != "MS risk contribution" else "normal")
    ax.tick_params(length=1.6, pad=1.2)
    ax.grid(False)
    for gi in sorted(M.gorder.unique())[1:]:              # 组间横线
        y = M.index[M.gorder == gi].min() - 0.5
        ax.plot([-0.55, 5.55], [y, y], color="0.35", lw=0.6, zorder=4, clip_on=False)
    for gi in sorted(M.gorder.unique()):                  # 左侧竖括号（图框内侧）
        s = M.index[M.gorder == gi]
        ax.plot([-0.505, -0.425, -0.425, -0.505],
                [s.min() - 0.34, s.min() - 0.34, s.max() + 0.34, s.max() + 0.34],
                color="0.3", lw=0.6, zorder=4)
    for s in ax.spines.values():
        s.set_linewidth(0.5)
    bcb = fig.add_axes([LM + AXW + 0.020, AX_BOT, 0.020, AX_H])   # 竖直 colorbar，不占底部
    mb = fig.colorbar(matplotlib.cm.ScalarMappable(matplotlib.colors.Normalize(*VQ), "coolwarm"), cax=bcb)
    mb.set_ticks([round(VQ[0]), 0, round(VQ[1])])
    mb.ax.tick_params(labelsize=6.0, length=1.6, pad=1.0); mb.outline.set_linewidth(0.5)
    mb.set_label("scaled mean expression\n(row-wise z)", fontsize=6.0, labelpad=1.5, linespacing=1.25)
    # 贡献条已按作者要求移除，量值改由附图面板承担（原 fig6v22_D LOGO 柱状图）
    lx = fig.add_axes([0.175, 0.012, 0.40, 0.078]); lx.axis("off")  # 底部一行圆点图例
    for k, pv in enumerate(LEG):
        cx = 0.10 + 0.33 * k
        lx.scatter([cx], [0.70], s=pv * S_PER_PCT, facecolors="none", edgecolors="0.25", linewidths=0.4)
        lx.text(cx, 0.24, "%d%%" % pv, ha="center", va="top", fontsize=6.0)
    PNAME = "fig6v23_D_identity_dotplot.pdf"
else:
    # ---- 横版（打横）：亚型在纵轴 6 行，基因在横轴 17 列
    # 成熟轴由「左→右」变「上→下」（作者批准转置，内容不变，§14.6 不颠倒顺序）
    # 底部模块标注：线宽不短于模块名，窄组下沉一行（阶梯），保证线与名一一对应
    RH2, TOP_CM, LAB_CM, GRP_CM = 0.48, 0.30, 1.41, 0.72
    BOT_CM = LAB_CM + GRP_CM
    HH = RH2 * 6 + TOP_CM + BOT_CM
    fig = plt.figure(figsize=(TARGET_W * CM, HH * CM))
    LM, AXW = 0.120, 0.655        # 让出右栏 2.4 cm 给图例
    AX_BOT, AX_H = BOT_CM / HH, 1 - (BOT_CM + TOP_CM) / HH
    ax = fig.add_axes([LM, AX_BOT, AXW, AX_H])
    RGT = (LM + AXW) * TARGET_W                       # 主图右边界（cm）

    NG = len(set(M.group))
    GAP = 0.55
    colw = AXW * TARGET_W / (n + GAP * (NG - 1))
    xcol, gspan, cur = [], {}, 0.0
    for gi, g in enumerate(sorted(M.gorder.unique())):
        if gi:
            cur += GAP
        start = cur
        for _ in range((M.gorder == gi).sum()):
            xcol.append(cur); cur += 1.0
        gspan[gi] = (start, cur - 1.0)
    xcol = np.array(xcol, float)
    XLO, XHI = xcol.min() - 0.60, xcol.max() + 0.60
    Zt, Pt = Zc.T, PCT.T
    xs, ys = np.meshgrid(xcol, np.arange(6))
    for j_ in range(n):
        pk = int(np.argmax(Zc[j_]))
        ax.add_patch(matplotlib.patches.Rectangle((xcol[j_] - 0.31, pk - 0.31), 0.62, 0.62,
                     fill=False, edgecolor="0.55", lw=0.40, zorder=2))
    ax.scatter(xs.ravel(), ys.ravel(), s=(np.clip(Pt, 0, None) * S_PER_PCT).ravel(),
               c=Zt.ravel(), cmap="coolwarm", vmin=VQ[0], vmax=VQ[1],
               edgecolors="0.25", linewidths=0.28, zorder=3)
    ax.set_xlim(XLO, XHI); ax.set_ylim(-0.62, 5.86); ax.invert_yaxis()
    # va="top" 是关键：va="center" 会把竖排 bbox 中心放在 tick 上，文字一半长进轴内（实测 2.16 cm），
    # 这才是“基因名和图形重叠”的根因，加大 pad / ylim 都治不了。
    ax.set_xticks(xcol)
    ax.set_xticklabels(M.gene_human, fontsize=6.0, rotation=90, ha="center", va="top")
    ax.set_yticks(range(6)); ax.set_yticklabels(SUB, fontsize=6.0)
    ax.tick_params(axis="x", length=0, pad=3.0)       # 刻度线关掉，否则戳在末字母上像下划线
    ax.tick_params(axis="y", length=1.6, pad=1.2)
    ax.grid(False)
    for s_ in ax.spines.values():
        s_.set_linewidth(0.5)
    for gi in sorted(M.gorder.unique())[1:]:
        xb = (xcol[M.index[M.gorder == gi].min()] + xcol[M.index[M.gorder == gi - 1].max()]) / 2
        ax.plot([xb, xb], [-0.62, 5.86], color="0.35", lw=0.6, zorder=4, clip_on=False)

    # 线严格等于组列跨度（2/2 列合并成 4 列后，OPC 段与 mature myelin 等长）；名两行居中于线
    GNAME = {0: "OPC /\nintermediate", 1: "mature\nmyelin", 2: "stage TF", 3: "MS risk\ncontribution"}
    LINE_Y, NAME_TOP = 0.60, 0.48
    for gi in sorted(M.gorder.unique()):
        x0, x1 = gspan[gi]
        a = LM * TARGET_W + x0 * colw - 0.5 * colw + 0.06
        b = LM * TARGET_W + x1 * colw + 0.5 * colw - 0.06
        fig.add_artist(matplotlib.lines.Line2D([a / TARGET_W, b / TARGET_W], [LINE_Y / HH, LINE_Y / HH],
                                                color="0.25", lw=0.5, clip_on=False))
        fig.text(((a + b) / 2) / TARGET_W, NAME_TOP / HH, GNAME[gi], fontsize=6.0,
                 ha="center", va="top", linespacing=1.25)

    bcb = fig.add_axes([LM + AXW + 0.025, AX_BOT, 0.020, AX_H])
    mb = fig.colorbar(matplotlib.cm.ScalarMappable(matplotlib.colors.Normalize(*VQ), "coolwarm"), cax=bcb)
    mb.set_ticks([round(VQ[0]), 0, round(VQ[1])])
    mb.ax.tick_params(labelsize=6.0, length=1.6, pad=1.0); mb.outline.set_linewidth(0.5)
    mb.set_label("scaled mean expression\n(row-wise z)", fontsize=6.0, labelpad=2.0,
                 rotation=90, linespacing=1.30)
    LXM = RGT + 0.06
    CY = [1.45, 1.02, 0.59]                           # 整体下移：上一版 Expressed (%) 压在 colorbar 的 -2 刻度上
    fig.text(LXM / TARGET_W, 1.88 / HH, "Expressed (%)", fontsize=6.0, ha="left", va="center")
    lx = fig.add_axes([(LXM + 0.02) / TARGET_W, (CY[2] - 0.20) / HH, 0.055, (CY[0] - CY[2] + 0.40) / HH])
    lx.axis("off"); lx.set_xlim(0, 1); lx.set_ylim(0, 1)
    for k, pv in enumerate(LEG):                     # 圈用 axes 画，文字用 figure 坐标（off-axes text 不入 PDF）
        lx.scatter([0.50], [(CY[k] - CY[2] + 0.20) / (CY[0] - CY[2] + 0.40)], s=pv * S_PER_PCT,
                   facecolors="none", edgecolors="0.25", linewidths=0.4)
        fig.text((LXM + 0.56) / TARGET_W, CY[k] / HH, "%d%%" % pv, fontsize=6.0, ha="left", va="center")
    PNAME = "fig6v23_C_identity_dotplot_T.pdf"

p = PAN + "/" + PNAME
fig.savefig(p, format="pdf"); plt.close(fig)
d = fitz.open(p); r0 = d[0].rect; d.close()
w_cm, h_cm = r0.width / 72 * 2.54, r0.height / 72 * 2.54
print("写出 %-42s %.2f x %.2f cm %s  （6C 高 5.35 cm）" % (PNAME, w_cm, h_cm,
      "OK" if abs(w_cm - TARGET_W) <= 0.06 else "! 宽 %.2f ≠ 目标 %.2f" % (w_cm, TARGET_W),))

# ---------------------------------------------------------------- 落表
T = pd.concat([M[["group", "gene", "gene_human", "contribution"]],
               pd.DataFrame(MEAN, columns=["mean_" + s for s in SUB]),
               pd.DataFrame(PCT, columns=["pct_" + s for s in SUB]),
               pd.DataFrame(Zc, columns=["z_" + s for s in SUB])], axis=1)
T.insert(0, "row", range(1, n + 1))
T["peak_subtype"] = [SUB[i] for i in M.peak_idx]
T.to_csv(TAB + "/panel6C_dotplot_markers.csv", index=False)
json.dump(dict(version="v6", panel_letter="C (was D; swapped with the cross-omics scatter, now D)",
               narrative_order="A/B/C are RNA-only evidence, D is the cross-omics comparison (author request)",
               label_anchor="xtick rotation=90 ha=center va=top; va=center was the real cause of the gene-name-overlaps-dots bug (half the rotated bbox sits inside the axis, measured 2.16 cm deep)", orientation="genes_on_x" if TRANSPOSE else "genes_on_y",
               panel_file=PNAME,
               contribution_bars="removed from this panel per author; magnitude now carried by the supplementary LOGO bar panel",
               dot_max_diameter_mm=DOT_MM_EFF, row_pitch_cm=RH if not TRANSPOSE else RH2,
               grouping_style="thin rule exactly = group column span, module name (2 lines) centred below; no colour band",
               module_labels=[v.replace("\n", " ") for v in GNAME.values()] if TRANSPOSE else None,
               xtick_pad_pt=6.0 if TRANSPOSE else None, cells=int(mask.sum()), nuclei_per_subtype={s: int(len(grp[s])) for s in SUB},
               reference_n_from_fig6A_stats="203/490/1016/688/991/2574 = 5962",
               expression="log1p(CP10k), Seurat RNA$data layer of 3k_dsci_oligo_subset.h5ad",
               ortholog_rule="scDRS built-in mouse_human_homologs.txt (16,466 pairs)",
               n_rows=int(n), rows_per_group={g: int((M.group == g).sum()) for g in M.group.unique()},
               colour_limits_rowwise_z=list(VQ),
               risk_selected=[r["gene_human"] for r in risk],
               risk_order="contribution descending (matches the right-hand bar)",
               risk_pct_threshold=PCT_MIN, risk_excluded_before_selection=excl,
               markers_omitted=OMITTED,
               naming=("mouse symbols italic for measured lineage markers; "
                       "human uppercase for risk-contribution genes (they enter the human GWAS gene set)"),
               channels="RNA only"), open(TAB + "/panel6C_dotplot_meta%s.json" % ("_T" if TRANSPOSE else ""), "w"), indent=1)
print(T[["row", "group", "gene_human", "peak_subtype", "contribution"]].round(3).to_string(index=False))
