#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""target_stats_GSE118257_v7.py — target x cell-population detection statistics for panel I: percentage of nuclei with >=1 read, and mean log2(CP10k+1) among detecting nuclei.

S8 / S7-I bubble map v5 (two font tiers 6/8, single-line row names, drug names in the legend; three colour scales, one per channel)."""
"""S8 / S7-I bubble map v5 (two font tiers 6/8, single-line row names, drug names in the legend; three colour scales, one per channel)."""
import os, sys, numpy as np, pandas as pd

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
sys.path.insert(0, os.environ.get("PYLIBS", ""))   # optional local site-packages
import matplotlib as mpl
from matplotlib import font_manager as _fm
import matplotlib.text as mtext
mpl.use("Agg")
mpl.rcParams["pdf.fonttype"] = 42; mpl.rcParams["ps.fonttype"] = 42
_fm.fontManager.addfont(WORK + "/fonts/ARIAL.TTF")
mpl.rcParams["font.family"] = _fm.FontProperties(fname=WORK + "/fonts/ARIAL.TTF").get_name()
mpl.rcParams.update({"font.size": 6.0, "axes.linewidth": 0.5, "axes.spines.top": False,
                    "axes.spines.right": False})
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

OUT = WORK + "/figures/S8_bubble"; TAB = OUT + "/tables"
PT = 72.0 / 2.54                        # pt per cm
LET, FS = 8.0, 6.0
RED, ORG, MUT = "#D51F26", "#FF7F0E", "0.34"
CMAP = LinearSegmentedColormap.from_list("w2r", ["#FFFFFF", "#FFD9B3", ORG, RED])
CMAPB = LinearSegmentedColormap.from_list("w2b", ["#F2F2F2", "#A6C8E0", "#1F77B4", "#103A6B"])
S_MIN, S_MAX = 4.0, 56.0
V_T, V_M, V_BC = (2.0, 5.0), (4.0, 8.0), (1.5, 4.5)

A3 = pd.read_csv(TAB + "/A3_GSE118257_target_stats.csv")
A2 = pd.read_csv(TAB + "/A2_GSE277430_target_stats.csv").rename(columns={"mean_log": "mean_exp"})
MO = pd.read_csv(TAB + "/MO_dsciCAT_target_stats.csv").rename(columns={"mean_log": "mean_exp"})
MSR = pd.read_csv(WORK + "/figures/routeA_out/A3_GSE118257_v1/tables/celltype_MS_score.csv")
MS_BY = dict(zip(MSR["Celltypes"], MSR["mean"]))
NCELL = A3.groupby("group")["cells"].first(); MEDGENE = A3.groupby("group")["med_genes"].first()
COLS = ["OPCs", "COPs", "ImOlGs", "Oligo1", "Oligo2", "Oligo3", "Oligo4", "Oligo5", "Oligo6",
        "Astro", "Micro/Mac", "Immune", "Endo", "Pericyte/SMC", "Neuron"]
TG = ["S1PR5", "S1PR1", "S1PR3", "S1PR4", "CD52", "ITGA4", "NFE2L2", "DHODH", "BTK",
      "KCNA2", "KCNA1", "KCNA4"]
MG = ["PDGFRA", "PLP1", "GFAP", "P2RY12", "CLDN5", "RBFOX3"]
CORE = ["S1PR5", "S1PR1", "NFE2L2", "CD52"]
MOUSE_MAP = dict(zip(TG, ["S1pr5", "S1pr1", "S1pr3", "S1pr4", "Cd52", "Itga4", "Nfe2l2",
                          "Dhodh", "Btk", "Kcna2", "Kcna1", "Kcna4"]))
MCOLS = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
MN = MO.groupby("group")["cells"].first().to_dict()
BCOLS = ["OLIGO", "OPC"]; BN = {"OLIGO": 35836, "OPC": 8931}


def PP(df, rows, cols):
    return (df.pivot(index="gene", columns="group", values="pct").reindex(index=rows, columns=cols),
            df.pivot(index="gene", columns="group", values="mean_exp").reindex(index=rows, columns=cols))


def dots(ax, P, E, cols, px, py, vr, yoff=0.0, cmap=CMAP):
    n = len(P)
    for r, g in enumerate(P.index):
        for c, k in enumerate(cols):
            p, e = P.loc[g, k], E.loc[g, k]
            x, y = (c + .5) * px, yoff + (n - r - .5) * py
            if not np.isfinite(p):
                ax.text(x, y, "nd", fontsize=FS, ha="center", va="center", color="0.45")
                continue
            if p <= 0:
                ax.plot(x, y, marker="o", ms=1.4, mfc="none", mec="0.45", mew=0.4, zorder=3)
            else:
                ax.scatter(x, y, s=S_MIN + (S_MAX - S_MIN) * p / 100.0,
                           c=[min(max(float(e) if np.isfinite(e) else 0, vr[0]), vr[1])],
                           cmap=cmap, vmin=vr[0], vmax=vr[1], edgecolors="0.25", linewidths=0.25,
                           zorder=5)


DISP = {"Pericyte/SMC": "Peri/SMC", "Micro/Mac": "Myeloid"}


def wrap(k):
    return DISP.get(k, k)


def xlab(ax, cols, y, px, nlab=None, two=True):
    for c, k in enumerate(cols):
        n = int((nlab or NCELL.to_dict())[k])
        s = wrap(k) + ("\n%s" % format(n, ",") if two else "")
        ax.text((c + .5) * px, y, s, fontsize=FS, ha="center", va="top", linespacing=1.15, zorder=9)


def ylab(ax, rows, x, py):
    n = len(rows)
    for r, g in enumerate(rows):
        ax.text(x, (n - r - .5) * py, g, fontsize=FS, ha="right", va="center", zorder=9)


def cbar(fig, W, H, x, ytop, w, vr, lab, cmap=CMAP):
    a = fig.add_axes([x / W, (ytop - 0.24) / H, w / W, 0.24 / H])
    a.imshow(np.linspace(vr[0], vr[1], 256).reshape(1, -1), aspect="auto", cmap=cmap,
               vmin=vr[0], vmax=vr[1], extent=[vr[0], vr[1], 0, 1])
    a.set_xticks([vr[0], vr[1]]); a.set_xticklabels(["%g" % vr[0], "%g" % vr[1]])
    a.set_yticks([]); a.tick_params(labelsize=FS, length=1.3, width=0.4, pad=0.6)
    a.text(0, 1.30, lab, fontsize=FS, ha="left", va="bottom", color="0",
           transform=a.transAxes)
    for s in a.spines.values():
        s.set_linewidth(0.35)
    return a


def size_key(fig, W, H, x, ytop, step):
    a = fig.add_axes([x / W, (ytop - 0.80) / H, (step * 3) / W, 0.80 / H])
    a.set_xlim(0, step * 3); a.set_ylim(0.80, 0); a.axis("off")
    for i, p in enumerate([1, 10, 50]):
        a.scatter((i + .5) * step, 0.28, s=S_MIN + (S_MAX - S_MIN) * p / 100.0, c="0.55",
                  edgecolors="0.25", linewidths=0.25)
        a.text((i + .5) * step, 0.58, "%d%%" % p, fontsize=FS, ha="center", va="top", color=MUT)
    return a


def diag(fig, Wcm, Hcm, tag):
    r = fig.canvas.get_renderer(); sc = 2.54 / fig.dpi
    bb = []
    for o in fig.findobj(mtext.Text):
        t = o.get_text().strip()
        if not t or not o.get_visible():
            continue
        e = o.get_window_extent(renderer=r)
        bb.append((t.replace("\n", "|"), e.x0 * sc, e.x1 * sc, e.y0 * sc, e.y1 * sc))
    over = []
    for i in range(len(bb)):
        for j in range(i + 1, len(bb)):
            a, b = bb[i], bb[j]
            ox = min(a[2], b[2]) - max(a[1], b[1]); oy = min(a[4], b[4]) - max(a[3], b[3])
            if ox > 0.02 and oy > 0.02:
                over.append((a[0][:18], b[0][:18], round(ox, 3), round(oy, 3)))
    bleed = [b for b in bb if b[1] < 0.03 or b[2] > Wcm - 0.03 or b[3] < 0.02 or b[4] > Hcm - 0.02]
    print("### %-8s text=%d overlap=%d bleed=%d" % (tag, len(bb), len(over), len(bleed)))
    for o in over[:14]:
        print("    OVER ", o)
    for b in bleed[:10]:
        print("    BLEED", b[0][:24], "%.2f-%.2f x %.2f-%.2f" % (b[1], b[2], b[3], b[4]))


class Cur:
    def __init__(s, H): s.y = H
    def take(s, h, gap=0.0):
        s.y -= gap; top = s.y; s.y -= h; return top


# ===================== 1) full set =====================
W1 = 17.0; PADL = 0.40; LBW = 1.55; XD = PADL + LBW
PX, PY = 0.855, 0.345
hA = PY * (len(TG) + len(MG)); hB = PY * len(TG)
segs = [(0.30, 0.0), (0.30, 0.02), (0.34, 0.06), (hA, 0.10), (0.66, 0.0), (0.86, 0.22),
        (0.92, 0.10), (0.30, 0.30), (hB, 0.0), (0.66, 0.0), (0.30, 0.26), (0.86, 0.10)]
H1 = round(sum(h + g for h, g in segs) + 0.34, 2)
fig = plt.figure(figsize=(W1 / 2.54, H1 / 2.54), dpi=300)
c = Cur(H1)
t_marg = c.take(0.30); t_titleA = c.take(0.30, 0.02); t_strip = c.take(0.34, 0.06)
t_dataA = c.take(hA, 0.10); t_xlabA = c.take(0.66); t_depth = c.take(0.86, 0.22)
t_note1 = c.take(0.92, 0.10); t_titleB = c.take(0.30, 0.30); t_dataB = c.take(hB)
t_xlabB = c.take(0.66); t_note2 = c.take(0.30, 0.26); t_note3 = c.take(0.30, 0.14); t_keys = c.take(0.86, 0.10)


def mkax(top, h, left=XD, w=None):
    w = PX * len(COLS) if w is None else w
    a = fig.add_axes([left / W1, (top - h) / H1, w / W1, h / H1])
    a.set_xlim(0, w); a.set_ylim(0, h); a.axis("off"); return a


axA = mkax(t_dataA, hA)
PA, EA = PP(A3, TG + MG, COLS)
axA.axhline(PY * len(TG), color="0.5", lw=0.6, ls=(0, (1.4, 1.4)), zorder=2)
dots(axA, PA.iloc[:len(TG)], EA.iloc[:len(TG)], COLS, PX, PY, V_T, yoff=PY * len(MG))
dots(axA, PA.iloc[len(TG):], EA.iloc[len(TG):], COLS, PX, PY, V_M)
ylab(axA, TG + MG, -0.10, PY)
xlab(axA, COLS, -0.10, PX)
axS = mkax(t_strip, 0.34)
gv = [MS_BY[k] for k in COLS if k in MS_BY]; lo, hi = min(gv), max(gv)
for i, k in enumerate(COLS):
    if k not in MS_BY:
        continue
    t = (MS_BY[k] - lo) / (hi - lo)
    axS.add_patch(plt.Rectangle((i * PX + 0.03, 0), PX - 0.06, 0.34, fc=CMAPB(0.10 + .90 * t),
                                ec="none", zorder=3))
    axS.text((i + .5) * PX, 0.17, "%+.2f" % MS_BY[k], fontsize=FS, ha="center", va="center",
             color="0" if t < 0.66 else "white", zorder=5)
axS.text(-0.10, 0.17, "MS scDRS", fontsize=FS, ha="right", va="center", color=MUT)
axD = mkax(t_depth, 0.86)
axD.set_ylim(0, 2100)
for i, k in enumerate(COLS):
    axD.bar((i + .5) * PX, MEDGENE[k], width=PX - 0.30, color="0.72", linewidth=0, zorder=3)
axD.plot([0, PX * len(COLS)], [0, 0], color="0.35", lw=0.5, zorder=5)
axD.plot([0, PX * len(COLS)], [1000, 1000], color="0.45", lw=0.4, ls=(0, (1.2, 1.2)), zorder=4)
axD.text(-0.10, 1000, "1,000", fontsize=FS, ha="right", va="center", color=MUT)
axD.text(0, 2100, "median genes detected per nucleus", fontsize=FS, ha="left", va="bottom",
         color=MUT)
fig.text(PADL / W1, t_titleA / H1, "A", fontsize=LET, va="top")
fig.text(XD / W1, t_titleA / H1, "GSE118257 \u2014 human MS + control nuclei, %s nuclei in 15 populations" % format(int(NCELL[COLS].sum()), ","),
         fontsize=FS, va="top")
fig.text(XD / W1, t_note1 / H1,
         "absent from this matrix, not plotted: MS4A1, S1PR2, CD19          "
         "\u25cb open circle = 0 nuclei with a read", fontsize=FS, va="top", color=MUT)
fig.text(XD / W1, (t_note1 - 0.30) / H1,
         "pooled columns: Astro = Astrocytes + Astrocytes2;  Myeloid = Microglia_Macrophages + Macrophages;  "
         "Peri/SMC = Pericytes + VSMC\nEndo = Endothelial_cells1 + 2;  Neuron = Neuron1\u20135;  "
         "OPCs, COPs, ImOlGs and Oligo1\u20136 are the author labels as given",
         fontsize=FS, va="top", color=MUT, linespacing=1.20)
PB, EB = PP(A2, TG, BCOLS)
xbw = PX * 2 * 0.86
axB = mkax(t_dataB, hB, w=xbw)
dots(axB, PB, EB, BCOLS, xbw / 2, PY, V_BC)
ylab(axB, TG, -0.10, PY)
xlab(axB, BCOLS, -0.10, xbw / 2, nlab=BN)
fig.text(PADL / W1, t_titleB / H1, "B", fontsize=LET, va="top")
fig.text(XD / W1, t_titleB / H1, "GSE277430", fontsize=FS, va="top")
xC = XD + xbw + 1.12
PC, EC = PP(MO, [MOUSE_MAP[g] for g in TG], MCOLS); PC.index = TG; EC.index = TG
xcw = 0.95 * 6
axC = mkax(t_dataB, hB, left=xC, w=xcw)
dots(axC, PC, EC, MCOLS, xcw / 6, PY, V_BC)
xlab(axC, MCOLS, -0.10, xcw / 6, nlab=MN)
fig.text((xC - 0.52) / W1, t_titleB / H1, "C", fontsize=LET, va="top")
fig.text(xC / W1, t_titleB / H1, "mouse dsciCAT (6 subtypes)", fontsize=FS, va="top")
size_key(fig, W1, H1, XD, t_keys + 0.06, 0.90)
fig.text((XD - 0.10) / W1, (t_keys - 0.30) / H1, "% nuclei detected", fontsize=FS, ha="right",
         va="center", color=MUT)
xb = XD + 3.10
cbar(fig, W1, H1, xb, t_keys - 0.06, 2.30, V_T, "A targets: log2(CP10k+1)")
cbar(fig, W1, H1, xb + 2.62, t_keys - 0.06, 2.30, V_M, "A reference rows")
cbar(fig, W1, H1, xb + 5.44, t_keys - 0.06, 2.30, V_BC, "B, C: log-normalised count")
fig.text((xb + 8.15) / W1, (t_keys - 0.18) / H1,
         "colour scales differ by channel",
         fontsize=FS, ha="left", va="center", color=MUT)
diag(fig, W1, H1, "FULL")
fig.savefig(OUT + "/pdf/FigS8_v6_full.pdf")
fig.savefig(OUT + "/png/FigS8_v6_full.png", dpi=300)
plt.close(fig)

# ===================== 2) compact variant (the old S7 row-4 panel I) =====================
W2, H2 = 17.0, 2.95
PYc, PXc = 0.40, 0.800
LBL2 = 1.05; XPAD = 0.40; xd2 = XPAD + LBL2; wd2 = PXc * len(COLS)
fig2 = plt.figure(figsize=(W2 / 2.54, H2 / 2.54), dpi=300)
axf = fig2.add_axes([xd2 / W2, (H2 - 0.44 - PYc * 4) / H2, wd2 / W2, PYc * 4 / H2])
axf.set_xlim(0, wd2); axf.set_ylim(0, PYc * 4); axf.axis("off")
PF, EF = PP(A3, CORE, COLS)
dots(axf, PF, EF, COLS, PXc, PYc, (2.2, 4.4))
ylab(axf, CORE, -0.10, PYc)
xlab(axf, COLS, -0.10, PXc)
fig2.text(XPAD / W2, (H2 - 0.03) / H2, "I", fontsize=LET, va="top")
xk = xd2 + wd2 + 0.16
size_key(fig2, W2, H2, xk, 1.55, 0.55)
cbar(fig2, W2, H2, xk + 1.95, 1.30, 1.15, (2.2, 4.4), "")
diag(fig2, W2, H2, "COMPACT")
fig2.savefig(OUT + "/pdf/FigS7I_v7_compact.pdf")
fig2.savefig(OUT + "/png/FigS7I_v7_compact.png", dpi=300)
plt.close(fig2)
import fitz
for f in ["FigS8_v6_full", "FigS7I_v7_compact"]:
    d = fitz.open(OUT + "/pdf/%s.pdf" % f); r = d[0].rect; d.close()
    print("%-18s %.2f x %.2f cm" % (f, r.width / PT, r.height / PT))
