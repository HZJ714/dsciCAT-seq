#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""FigS6 v17 —— D/E 两格**共用一条图例**（作者 2026-11-17：D/E 感觉可以共用一个图例）。

为什么现在能共用：两块都改成"每基因跨 6 个 OL 亚型标准化"的 z 之后，RNA 与 ATAC 不再
是两套单位（log2(CP10k+1) vs gene score），而是同一个无量纲刻度 ⇒ 一条色标就够。
原先每格各带两行图例 = 4 行图例、两条重复的渐变条，读者要在两格之间来回对齐。

几何：全宽 17.00 cm × 1.05 cm，单行，居中；面板由 p0_8c_fig6c_{top1000,landing}.py 以
PANEL_LEGEND=0 关掉自带图例。字号两档（8.0/6.0），面上文字一律纯黑（house-spec §2）。
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import fitz
import sys

WORK = os.environ.get('DSCI_WORK_ROOT', '')
# 字体必须走 arial_style：只设 rcParams["font.family"]="Arial" 在本机**不会**生效
# （Arial 由脚本从 $DSCI_WORK_ROOT/fonts 注册，不在系统字体库里），结果静默退回
# DejaVu Sans，font_gate --strict 直接判字体名违规。
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arial_style import use_arial, use_two_tier
PAN = WORK + "/figures/routeA_out/P0_8C/panels"
CM = 2.54
BLACK, GOLD, GREY, EDGE = "#000000", "#b8860b", "#d9d9d9", "#333333"
PW, PH = 17.00, float(os.environ.get("LEG_H", 1.05))
ZCLIP = float(os.environ.get("PANEL_ZCLIP", 2.0))
CMAP = LinearSegmentedColormap.from_list("by", ["#2166ac", "#f7f7f7", "#b2182b"])

use_arial(); use_two_tier()


def tw(t, fs=6.0):
    return fitz.get_text_length(t, fontname="helv", fontsize=fs) * CM / 72.0


fig = plt.figure(figsize=(PW / CM, PH / CM))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, PW); ax.set_ylim(0, PH); ax.axis("off")

# 与两格体表同一条水平中线：体表底边在面板 y=BODY_B(1.35)，共用条带放在其下 0.34 行距处
Y = 0.50
MKW = {"grad": 1.45, "sw": 0.26, "car": 0.30, "bar": 0.42, "bar2": 0.42}
items = [("grad", "z (per gene, across 6 OL subtypes)"), ("sw", "not detected"),
         ("car", "off scale"), ("bar", "tested subtype (MOL6)"),
         ("bar2", "reference (COP/NFOL)")]   # Δ 列删除后，对比的两列只能靠条子声明
gap = 0.55
wid = [MKW[k] + (0.0 if k == "grad" else 0.06 + tw(t)) for k, t in items]
tot = sum(wid) + gap * (len(items) - 1)
assert tot <= PW - 0.30, "共用图例超宽 %.2f > %.2f" % (tot, PW - 0.30)
x = (PW - tot) / 2.0
for (kind, t), w in zip(items, wid):
    if kind == "grad":
        # 方向：低→高、蓝→红从左到右。v31/v32 的面板自带图例写成 linspace(HI, LO)，
        # 色条与体表同色但方向相反（红在左＝+2），读者要倒着看；共用条带顺手改回常规方向。
        # 体表映射是 CMAP((v−LO)/(HI−LO)) ⇒ 红＝高，与此处一致。
        ax.imshow(np.linspace(-ZCLIP, ZCLIP, 64)[None, :],
                  extent=(x, x + MKW[kind], Y - 0.075, Y + 0.075), cmap=CMAP, aspect="auto", zorder=3)
        ax.add_patch(Rectangle((x, Y - 0.075), MKW[kind], 0.15, fc="none", ec=EDGE, lw=0.3))
        # 端点刻度以条端为中心；与下一条目的文字之间靠 gap 保证不连读（见 SKILL §17.1）
        for px, lab in ((x, "−%.1f" % ZCLIP), (x + MKW[kind] / 2, "0"), (x + MKW[kind], "+%.1f" % ZCLIP)):
            ax.text(px, Y - 0.115, lab, fontsize=6.0, ha="center", va="top", color=BLACK)
        ax.text(x + MKW[kind] / 2, Y + 0.115, t, fontsize=6.0, ha="center", va="bottom", color=BLACK)
    elif kind == "sw":
        ax.add_patch(Rectangle((x, Y - 0.075), MKW[kind], 0.15, fc=GREY, ec="#8c8c8c", lw=0.3))
        ax.text(x + MKW[kind] + 0.06, Y, t, fontsize=6.0, ha="left", va="center", color=BLACK)
    elif kind == "car":
        ax.plot(x + 0.07, Y, marker="<", ms=2.2, color="#2166ac", mew=0.3)
        ax.plot(x + 0.23, Y, marker=">", ms=2.2, color="#b2182b", mew=0.3)
        ax.text(x + MKW[kind] + 0.06, Y, t, fontsize=6.0, ha="left", va="center", color=BLACK)
    else:
        ax.add_patch(Rectangle((x, Y - 0.055), MKW[kind], 0.075,
                               fc=GOLD if kind == "bar" else BLACK, ec="none"))
        ax.text(x + MKW[kind] + 0.06, Y, t, fontsize=6.0, ha="left", va="center", color=BLACK)
    x += w + gap
print("共用图例占位 %.2f–%.2f / %.2f cm" % ((PW - tot) / 2.0, x - gap, PW))

out = PAN + "/figS6v17_DE_legend.pdf"
fig.savefig(out, dpi=600)

d = fitz.open(out); p = d[0]
ws = p.get_text("words")
sp = [s for b in p.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]]
print("  %-34s %.2f x %.2f cm | 字右 %.2f 底 %.2f | 字号 %s | 非黑白 %s"
      % (os.path.basename(out), p.rect.width * CM / 72, p.rect.height * CM / 72,
         max(w[2] * CM / 72 for w in ws), max(w[3] * CM / 72 for w in ws),
         sorted({round(s["size"], 1) for s in sp}),
         [{True: "黑", False: "白"}.get(s["color"] in (0, 16777215), "#%06x" % s["color"])
          for s in sp if s["color"] not in (0, 16777215)]))
d.close()
