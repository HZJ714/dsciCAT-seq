# -*- coding: utf-8 -*-
"""p0_8c_fig6c_spectrum.py —— Fig6C v34：MAGMA top1000（= Fig6A/B 打分集）基因层"着陆谱"

域与 A/B 一致（作者指出 v31 用 38/331 属偏航：那是 FigS6A/B/E 的 *cis*-window 位点解析）。
统计完全不碰 scDRS 分数 ⇒ 无 self-inclusion。

显示量 = 被检验量：
  Δ[g] = 0.2810·(mean[MOL3/6]−mean[COP/NFOL]) + 0.7190·(mean[MOL6]−mean[COP/NFOL])
  检验：风险集 Δ 均值 vs **二维检测率/表达量匹配**背景（5×5 层内抽对照，10,000 次置换）
  正对照：8 个髓鞘/OL 锚基因（证明检验有灵敏度，阴性不是因为没功率）
所有数字一律从 p08c_dual_setlevel_matched2d.csv 读（单一来源，图面不重算，避免两套实现）：
  RNA  Δ̄=−0.0851  匹配零 −0.1104±0.0199  z=+1.27  P=0.100   锚 +0.3572  P=0.012
  ATAC Δ̄=−0.0177  匹配零 −0.0125±0.0087  z=−0.59  P=0.729   锚 +0.1556  P=0.010
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import fitz

CM = 2.54 / 72.0
WORK = os.environ.get('DSCI_WORK_ROOT', '')
TAB = WORK + "/figures/routeA_out/P0_8C/tables"
PAN = WORK + "/figures/routeA_out/P0_8C/panels"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arial_style import use_arial, use_two_tier
use_arial(); use_two_tier()
BLACK, RED, GREY, GOLD = "#000000", "#b2182b", "#bdbdbd", "#333333"

M = pd.read_csv(TAB + "/p08c_dual_landing_genomewide.csv")
L = pd.read_csv(TAB + "/p08c_genelevel_all6_pool_allgenes.csv")
T = pd.read_csv(TAB + "/p08c_dual_setlevel_matched2d.csv")
ms875 = set(L.gene.astype(str)); ms223 = set(L[~L.near_zero].gene.astype(str))
anch = ["Mbp", "Plp1", "Mog", "Mag", "Cnp", "Mobp", "Sox10", "Olig2"]
SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]

PANELS = []
for tag, nm, xlo, xhi, ticks in (("rna", "RNA", -1.65, 1.25, [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0]),
                                 ("atac", "ATAC", -0.42, 0.42, [-0.4, -0.2, 0.0, 0.2, 0.4])):
    d, pct = M[tag + "_d"].values, M[tag + "_pct"].values
    mean = M[[tag + "_" + s for s in SUB]].mean(1).values
    pool = (pct >= 10) & np.isfinite(d) & np.isfinite(mean)
    bg = pool & ~M.gene.isin(ms875).values
    rk = pool & M.gene.isin(ms223).values
    ak = pool & M.gene.isin(anch).values
    r = T[(T.channel == tag) & (T["set"] == "top1000可测223")].iloc[0]
    a = T[(T.channel == tag) & (T["set"] == "髓鞘锚(正对照)")].iloc[0]
    obs = d[rk].mean()
    assert abs(obs - r.mean_d) < 1e-4 and abs(d[ak].mean() - a.mean_d) < 1e-4, (tag, obs, r.mean_d)
    assert abs((obs - r.null_mean) / r.null_sd - r.z) < 0.02
    PANELS.append(dict(tag=tag, nm=nm, d=d, bg=bg, rk=rk, ak=ak, xlo=xlo, xhi=xhi, ticks=ticks,
                       obs=obs, nul=r.null_mean, sd=r.null_sd, p=r.p_greater, n=int(rk.sum()),
                       nb=int(bg.sum()), aobs=d[ak].mean(), ap=a.p_greater, nan=int(ak.sum())))
    print("%s n=%d Δ̄%+.4f 零%+.4f±%.4f P=%.3f | 锚 n=%d Δ̄%+.4f P=%.3f"
          % (nm, PANELS[-1]["n"], obs, r.null_mean, r.null_sd, r.p_greater, ak.sum(), d[ak].mean(), a.p_greater))

fig = plt.figure(figsize=(10.72 / 2.54, 5.50 / 2.54))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 10.72); ax.set_ylim(0, 5.50); ax.axis("off")
ax.text(0.75, 5.23, "MS GWAS top-1000 gene set vs detection-matched background",
        fontsize=8.0, ha="left", va="baseline", color=BLACK)
# 纵向：主标题 5.23 / 通道名 4.95 / n 行 4.72 / P 行 4.49 / 轴顶 4.10 / 轴底 1.60 / 刻度 1.51 / 轴标签 1.14 / 图例 0.50
Y0, H, W, X0, GAP = 1.60, 2.50, 4.02, 1.50, 0.58
Y1 = Y0 + H
for i, P in enumerate(PANELS):
    x0 = X0 + i * (W + GAP)
    bins = np.linspace(P["xlo"], P["xhi"], 46)
    hb = np.histogram(P["d"][P["bg"]], bins=bins, density=True)[0]
    hr = np.histogram(P["d"][P["rk"]], bins=bins, density=True)[0]
    ha = np.histogram(P["d"][P["ak"]], bins=bins, density=True)[0]
    ymax = max(hb.max(), hr.max()) * 1.45
    sx = lambda v: x0 + (v - P["xlo"]) / (P["xhi"] - P["xlo"]) * W
    sy = lambda v: Y0 + v / ymax * H
    for b0, b1, v in zip(bins[:-1], bins[1:], hb):
        if v > 0:
            ax.add_patch(Rectangle((sx(b0), Y0), max(sx(b1) - sx(b0) - 0.010, 0.01), sy(v) - Y0, fc=GREY, ec="none"))
    for b0, b1, v in zip(bins[:-1], bins[1:], hr):
        if v > 0:
            ax.add_patch(Rectangle((sx(b0), Y0), max(sx(b1) - sx(b0) - 0.010, 0.01), sy(v) - Y0,
                                   fc=RED, ec="none", alpha=0.55))
    ax.add_patch(Rectangle((x0, Y0), W, H, fc="none", ec="#333333", lw=0.4))
    # 匹配零 ± SD 带：画在轴内顶部留白带里（ymax 已放大 1.45，最高的柱也够不着）
    ax.add_patch(Rectangle((sx(P["nul"] - P["sd"]), Y1 - 0.36), sx(P["nul"] + P["sd"]) - sx(P["nul"] - P["sd"]),
                           0.20, fc="#ffffff", ec=BLACK, lw=0.4))
    ax.plot([sx(P["nul"])] * 2, [Y1 - 0.40, Y1 - 0.12], ls=(0, (1.2, 1.0)), lw=0.6, color=BLACK)
    ax.plot([sx(P["aobs"]), sx(P["aobs"])], [Y0, Y0 + H], lw=0.8, color=GOLD, ls=(0, (3.0, 1.5)))
    ax.plot([sx(P["obs"]), sx(P["obs"])], [Y0, Y0 + H], lw=1.1, color=RED)
    ax.text(x0 + W / 2, 4.95, P["nm"], fontsize=8.0, ha="center", va="baseline", color=BLACK)
    ax.text(x0 + W / 2, 4.72, "n = %d risk / %d background" % (P["n"], P["nb"]),
            fontsize=6.0, ha="center", va="baseline", color=BLACK)
    ax.text(x0 + W / 2, 4.49, "risk P = %.3f    anchors P = %.3f" % (P["p"], P["ap"]),
            fontsize=6.0, ha="center", va="baseline", color=BLACK)
    for xv in P["ticks"]:
        ax.plot([sx(xv), sx(xv)], [Y0 - 0.06, Y0], lw=0.4, color=BLACK)
        ax.text(sx(xv), Y0 - 0.09, ("%+.1f" % xv).replace("+0.0", "0").replace("-0.0", "0")
                if abs(xv) > 1e-9 else "0", fontsize=6.0, ha="center", va="top", color=BLACK)
    ax.text(x0 - 0.16, Y0 + H / 2, "density", fontsize=6.0, rotation=90, ha="center", va="center", color=BLACK)
    ax.text(x0 + W / 2, Y0 - 0.46, "\u0394 enriched \u2212 COP/NFOL (%s)"
            % ("log2" if P["tag"] == "rna" else "score"), fontsize=6.0, ha="center", va="top", color=BLACK)


def tw(t, fs=6.0):
    return fitz.get_text_length(t, fontname="helv", fontsize=fs) * CM


def legend(x0, x1, y, items, gap=0.40):
    MKW = {"sw": 0.26, "ln": 0.45, "dash": 0.45, "band": 0.45}
    wid = [MKW[k] + 0.06 + tw(t) for k, t in items]
    tot = sum(wid) + gap * (len(items) - 1)
    if x0 + tot > x1:
        raise SystemExit("图例超宽 %.2f > %.2f" % (x0 + tot, x1))
    x = x0 + (x1 - x0 - tot) / 2.0
    for (kind, t), w in zip(items, wid):
        if kind == "sw":
            ax.add_patch(Rectangle((x, y - 0.075), MKW[kind], 0.15, fc=GREY, ec="#8c8c8c", lw=0.3))
        elif kind == "ln":
            ax.plot([x + 0.04, x + 0.41], [y, y], lw=1.1, color=RED)
        elif kind == "dash":
            ax.plot([x + 0.04, x + 0.41], [y, y], lw=0.8, color=GOLD, ls=(0, (3.0, 1.5)))
        else:
            ax.add_patch(Rectangle((x, y - 0.06), MKW[kind], 0.12, fc="#ffffff", ec=BLACK, lw=0.4))
            ax.plot([x + MKW[kind] / 2] * 2, [y - 0.09, y + 0.09], ls=(0, (1.2, 1.0)), lw=0.6, color=BLACK)
        ax.text(x + MKW[kind] + 0.06, y, t, fontsize=6.0, ha="left", va="center", color=BLACK)
        x += w + gap
    return x - gap


xe = legend(X0, X0 + 2 * W + GAP, 0.50,
            [("sw", "background"), ("ln", "risk mean"), ("dash", "anchors"),
             ("band", "matched null \u00b1 SD")])
print("图例 %.2f–%.2f（可用至 %.2f）" % (X0, xe, X0 + 2 * W + GAP))

out = PAN + "/fig6v34_C_top1000_spectrum.pdf"
fig.savefig(out, dpi=600)
d = fitz.open(out); p_ = d[0]
Wp, Hp = p_.rect.width * CM, p_.rect.height * CM
ws = p_.get_text("words")
sp_ = [s for b in p_.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]]
print("尺寸 %.2f × %.2f | 文字最右 %.2f（限 %.2f）| 最下 %.2f（限 %.2f）| 字号 %s | 非黑白字色 %s"
      % (Wp, Hp, max(w[2] * CM for w in ws), Wp - 0.20, max(w[3] * CM for w in ws), Hp - 0.10,
         sorted({round(s["size"], 2) for s in sp_}),
         sorted({"%06x" % s["color"] for s in sp_ if s["color"] not in (0, 0xffffff)})))
d.close()
