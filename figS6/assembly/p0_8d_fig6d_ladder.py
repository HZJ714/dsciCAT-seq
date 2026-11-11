#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig6D（终局 4 面板版）— 74 性状双组学梯子回归（原 Fig6C 版式回归主图）

规格 FIG6_7_SYNTHESIS_RECEIPT_20261112.md / FIG7_DRUG_HEATMAP_RECEIPT_20261112.md：
  D = 74 trait 双组学（原 Fig6C 梯子回归；MS/SCZ/HEIGHT 锚点）
  面数 ρ = 0.819 + P = 6.0e-4（add-one (1+5)/10001，env §0 定案；旧 5.0e-4 作废）
  三化学全谱 + 逐核箱型 → S6D（不进本面板）
数据 GitHub/fig6/fig6_out/fig6C_delta1_own74.csv（74 行，同核双通道）+ fig6C_stats.csv
版式 = 原 fig6v21_C_contrast.pdf 逐点复现（§14.13），只改两处：
  1) 字号归两档（锚标签 5.6 → 6.0，§9.1）
  2) 角注加 ρ / P（§14.8：角注块只留整块统计量，身份标注仍就近）
D 槽尺寸 6.85 × 5.35 cm（与 c2_compose_drafts 行2 C=10.72 相加 = 行宽）。
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms
from matplotlib import text as mpltext

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from arial_style import use_arial, use_two_tier, CM, ten_pow_sup
use_arial(); use_two_tier()

WORK = os.environ.get('DSCI_WORK_ROOT', '')
SRC = WORK + "/GitHub/fig6/fig6_out"
PAN = WORK + "/figures/routeA_out/P0_8C/panels"
os.makedirs(PAN, exist_ok=True)

RED, BLU, MUT, GRID = "#D51F26", "#272E6A", "0.42", "0.85"
MARK = {"PASS_Multiple_sclerosis": "MS",
        "PASS_Schizophrenia_Pardinas2018": "SCZ",
        "UKB_460K.body_HEIGHTz": "HEIGHT"}

T = pd.read_csv(SRC + "/fig6C_delta1_own74.csv")
ST = pd.read_csv(SRC + "/fig6C_stats.csv").iloc[0]
RHO = float(ST.rho_delta1)                       # 0.819089
P_ADD = (1 + 5) / 10001                           # env §0 定案：add-one，6.0e-4
assert abs(RHO - 0.819) < 1e-3 and abs(P_ADD - 6.0e-4) < 1e-5, (RHO, P_ADD)

FS_LAB, FS_TITLE = 6.0, 8.0                       # 两档制
fig, axc = plt.subplots(figsize=(6.37 * CM, 5.5 * CM))   # arial_style: CM = inch/cm，乘不是除
fig.subplots_adjust(left=0.135, right=0.985, top=0.925, bottom=0.115)
axc.grid(True, lw=0.32, color=GRID, zorder=0); axc.set_axisbelow(True)
x, y = T.d1_atac.values, T.d1_rna.values
lo, hi = min(x.min(), y.min()) * 1.12 - 0.06, max(x.max(), y.max()) * 1.12 + 0.06
axc.scatter(x, y, s=11, facecolor="0.75", edgecolor="none", zorder=2)
axc.plot([lo, hi], [lo, hi], ls="--", lw=0.5, color=MUT, zorder=1)
axc.axhline(0, color=GRID, lw=0.5, zorder=1); axc.axvline(0, color=GRID, lw=0.5, zorder=1)
axc.set_xlim(lo, hi); axc.set_ylim(lo, hi)
axc.set_xlabel("Gene score (MOL6 \u2212 COP/NFOL)", fontsize=FS_TITLE)
axc.set_ylabel("Gene expression\n(MOL6 \u2212 COP/NFOL)", fontsize=FS_TITLE, linespacing=1.15)
axc.tick_params(labelsize=FS_LAB, length=1.8, pad=1.6, width=0.5)
for sp in axc.spines.values():
    sp.set_linewidth(0.5)

# 锚点身份标注（§14.8 就近）。v21 的锚坐标是用 5.6pt 搜的，归到两档制 6.0pt 后
# 文本盒变大、会压到别的点（原脚本注释即预言此事），故这里重跑网格搜索：
#   约束 = 不出 axes、不压任何其他点（留 pad）、不碰角注与已放标签、引线 ≤ MAX_LEAD
#   评分 = 引线长度 + 0.6 × 到 v21 原锚位的距离（尽量保持原观感）
ANCHOR0 = {"Multiple_sclerosis": (0.430, 0.960), "HEIGHTz": (0.660, 0.720),
           "Schizophrenia": (0.600, 0.570)}
PAD, MAX_LEAD = 4.0, 90.0            # display px @120dpi（4px≈2.4pt > 灰点半径 1.9pt；90px≈54pt 引线）

fig.canvas.draw(); _rnd = fig.canvas.get_renderer()
_ar = axc.get_window_extent()
_px = axc.transData.transform(np.c_[x, y])
_pts = [(T.trait.iloc[i], _px[i, 0], _px[i, 1]) for i in range(len(_px))]

def _box_of(txt, fx, fy):
    """用**真实临时 Text**测盒（与 annotate 同参数同坐标系）。
    v21 那种“独立 Text(0,0) + 手拼高宽”的测法量不准多行盒（实测漏判 HEIGHT 压点），弃用。"""
    _t = axc.text(fx, fy, txt, transform=axc.transAxes, fontsize=FS_LAB,
                  fontweight="bold", linespacing=1.25, va="top", ha="left")
    fig.canvas.draw()
    bb = matplotlib.transforms.Bbox(_t.get_window_extent(_rnd).get_points())
    _t.remove()
    return bb

# 角注先占位（右下空白三角，§14.8 角注只留 ρ/P）
_cor = [axc.text(0.975, 0.115, "P = %.0e" % P_ADD, transform=axc.transAxes, fontsize=FS_LAB,
                  color="0.15", ha="right", va="baseline"),
        axc.text(0.975, 0.195, "\u03c1 = %.3f" % RHO, transform=axc.transAxes, fontsize=FS_LAB,
                  color="0.15", ha="right", va="baseline")]
fig.canvas.draw()
_reserved = [t.get_window_extent(_rnd) for t in _cor]

_GRID = [(fx, fy) for fy in np.arange(0.98, 0.06, -0.02) for fx in np.arange(0.02, 0.93, 0.02)]
_lab = []
for key, col in [("Multiple_sclerosis", RED), ("Schizophrenia", BLU), ("HEIGHTz", MUT)]:
    r = T[T.trait.str.contains(key, regex=False)].iloc[0]
    nm = MARK[r.trait]
    axc.scatter([r.d1_atac], [r.d1_rna], s=26, facecolor=col, edgecolor="white", lw=0.4, zorder=4)
    ptx, pty = axc.transData.transform((r.d1_atac, r.d1_rna))
    _txt = "%s\nRNA #%d  ATAC #%d" % (nm, r.rank_rna_raw, r.rank_atac_raw)
    best = None
    for fx, fy in _GRID:
        bx = _box_of(_txt, fx, fy)
        if (bx.x0 < _ar.x0 + 0.5 or bx.x1 > _ar.x1 - 0.5
                or bx.y0 < _ar.y0 + 0.5 or bx.y1 > _ar.y1 - 0.5):
            continue
        if any(tr != r.trait and bx.x0 - PAD <= px <= bx.x1 + PAD
               and bx.y0 - PAD <= py <= bx.y1 + PAD for tr, px, py in _pts):
            continue
        if any(bx.overlaps(rb) for rb in _reserved):
            continue
        lead = min(np.hypot(bx.x0 - ptx, bx.y0 - pty), np.hypot(bx.x1 - ptx, bx.y0 - pty),
                   np.hypot(bx.x0 - ptx, bx.y1 - pty), np.hypot(bx.x1 - ptx, bx.y1 - pty))
        if lead > MAX_LEAD:
            continue
        pref = axc.transAxes.transform(ANCHOR0[key])
        score = lead + 0.6 * np.hypot(bx.x0 - pref[0], bx.y1 - pref[1])
        if best is None or score < best[0]:
            best = (score, fx, fy, bx)
    assert best is not None, "no collision-free anchor for %s" % nm
    _, fx, fy, bx = best
    _reserved.append(bx)
    _lab.append((nm, r.trait, fx, fy, _txt,
                 axc.annotate(_txt, xy=(r.d1_atac, r.d1_rna), xytext=(fx, fy),
                              textcoords="axes fraction", fontsize=FS_LAB, color=col,
                              fontweight="bold", linespacing=1.25, zorder=6, va="top", ha="left",
                              arrowprops=dict(arrowstyle="-", lw=0.4, color=col,
                                              shrinkA=1.0, shrinkB=1.5))))
    print("anchor %-8s -> (%.2f, %.2f)" % (nm, fx, fy))

# ---- 硬门：用同一把尺（真实 Text 盒）复核：不出 axes、不压任何其他点 ----
fig.canvas.draw(); _rnd = fig.canvas.get_renderer()
_ar = axc.get_window_extent()
_px = axc.transData.transform(np.c_[x, y])
_pts = [(T.trait.iloc[i], _px[i, 0], _px[i, 1]) for i in range(len(_px))]
_bad = []
for nm, own, fx, fy, _txt, _a in _lab:
    box = _box_of(_txt, fx, fy)
    if (box.x0 < _ar.x0 - 0.5 or box.x1 > _ar.x1 + 0.5
            or box.y0 < _ar.y0 - 0.5 or box.y1 > _ar.y1 + 0.5):
        _bad.append("%s label outside the axes" % nm)
    _hit = [tr for tr, _x, _y in _pts if tr != own and
            box.x0 - PAD <= _x <= box.x1 + PAD and box.y0 - PAD <= _y <= box.y1 + PAD]
    if _hit:
        _bad.append("%s covers %s" % (nm, _hit))
for _t in _cor:
    box = _t.get_window_extent(_rnd)
    _hit = [tr for tr, _x, _y in _pts
            if box.x0 - PAD <= _x <= box.x1 + PAD and box.y0 - PAD <= _y <= box.y1 + PAD]
    if _hit:
        _bad.append("corner stat covers %s" % _hit)
assert not _bad, "label gate failed: %s" % _bad
print("label gate: 0 outside, 0 covered")

out = PAN + "/fig6v23_D_dualomics_ladder.pdf"
fig.savefig(out, bbox_inches="tight", pad_inches=0.03)
plt.close(fig)
print("写出", out)
