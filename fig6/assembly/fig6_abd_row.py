# -*- coding: utf-8 -*-
"""Fig6 single-row version (A / B / D only) - 2026-11-17

Background: the gene-level Fig6C was judged negative / untenable as a main-figure claim (receipt §16-17), so C left the main figure
and A/B/D moved from the old 2x2 into a SINGLE ROW.

Source: the panel code follows GitHub/fig6/build_fig6_notebook.py (the delivered v21 generator) line by line,
with only three changes -
  1) target frames: A/B/D re-render to equal height H, total width <= 17.9 cm (the delivered v21 width);
  2) A's subtype labels rotate 90° (in a single row each cell has only ~0.8 cm; horizontal labels would inevitably collide);
  3) D's column widths and chemistry names compress (12 columns -> ~0.29 cm each; short names dCAT/10x/HT,
     full names go into the caption).
STATS=full keeps A's four-line statistics table (MWU P / group P / group q / Cliff delta);
STATS=none moves the statistics into the caption, leaving A with boxes + jitter + delta marks (narrower and easier to read).

Usage: DSCI_WORK_ROOT=<work> python3 fig6_abd_row.py
"""
import os, sys, time
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import anndata as ad
from scipy.stats import mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("DSCI_WORK_ROOT", ".")   # export it; do not hardcode a machine
sys.path.insert(0, HERE)
from arial_style import use_arial, FONTS_PDF             # the font-whitelist convention (same source as the delivered)
CM = 2.54
OUT = os.environ.get("FIG6_OUT", os.path.join(WORK, "GitHub/fig6/fig6_out"))
DST = os.environ.get("FIG6_ROW_OUT", os.path.join(WORK, "figures/routeA_out/P0_8C"))
STATS = os.environ.get("STATS", "full")
FY = os.path.join(WORK, "FY")
ASPECT = os.environ.get("ASPECT", "auto")          # auto = fill the row height (vertical stretch); equal = keep the shape
H5_R = os.path.join(FY, "oligo", "3k_dsci_oligo_subset.h5ad")
MS_SCORE = os.path.join(FY, "oligo", "scdrs_ms_out", "PASS_Multiple_sclerosis.score.gz")
GS_GROUP = os.path.join(FY, "oligo", "scdrs_ms_out", "PASS_Multiple_sclerosis.scdrs_group.cell_type")
SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
STALLION = {"COP/NFOL": "#D51F26", "MFOL": "#272E6A", "MOL2": "#208A42",
            "MOL3": "#89288F", "MOL3/6": "#F47D2B", "MOL6": "#FEE500"}
TECH = {"dsciCAT_seq": "dsciCAT-seq (own)", "10x": "10x Multiome", "HT_scCAT_seq_HT": "HT-scCAT-seq"}
NAME = ["dCAT", "10x", "HT"]                    # short names; full names go into the caption
NAME_FULL = ["dsciCAT-seq", "10x Multiome", "HT-scCAT-seq"]
MARK = {"PASS_Multiple_sclerosis": "MS", "PASS_Schizophrenia_Pardinas2018": "SCZ",
        "UKB_460K.body_HEIGHTz": "HEIGHT"}
TRAIT_COLOR = {"MS": "#D51F26", "SCZ": "#272E6A", "HEIGHT": "#7E7E7E"}
use_arial({"font.size": 8.0, "axes.linewidth": 0.7, "xtick.labelsize": 8.0,
           "ytick.labelsize": 8.0, "savefig.dpi": 600})


def cliffs(x, y):
    n1, n2 = len(x), len(y)
    gt = sum((a > y).sum() for a in x)
    lt = sum((a < y).sum() for a in x)
    return (gt - lt) / float(n1 * n2)


def bh(p):
    """Word-identical to the delivered build_fig6_notebook.py (NaN excluded from the family)."""
    p = np.asarray(p, float); q = np.full(len(p), np.nan)
    idx = np.where(~np.isnan(p))[0]; m = len(idx)
    if m == 0:
        return q
    o = idx[np.argsort(p[idx])]
    qs = p[o] * m / np.arange(1, m + 1)
    q[o] = np.minimum(np.minimum.accumulate(qs[::-1])[::-1], 1.0)
    return q


def fmtp(p):
    return "%.3f" % p if p >= 0.001 else "%.1fe%d" % (p / 10 ** int(np.floor(np.log10(p))),
                                                      int(np.floor(np.log10(p))))


def dark(c, thr=0.55):
    rgb = np.array(mplc.to_rgb(c)); lum = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
    return tuple(rgb * (thr / lum)) if lum > thr else tuple(rgb)


# ---------------------------------------------------------------- data (the same chain as v21)
t0 = time.time()
obs = ad.read_h5ad(H5_R, backed="r").obs
ms = pd.read_csv(MS_SCORE, sep="\t", index_col=0)["norm_score"]
df = obs[["cell_type", "UMAP_1", "UMAP_2"]].join(ms.to_frame("norm_score")).dropna()
df["cell_type"] = df.cell_type.astype(str)
grp = pd.read_csv(GS_GROUP, sep=None, engine="python", index_col="group")
ref = df.loc[df.cell_type == "COP/NFOL", "norm_score"].values
rows = []
for ct in SUB:
    v = df.loc[df.cell_type == ct, "norm_score"].values
    rows.append(dict(cell_type=ct, n=len(v), median=float(np.median(v)),
                     mean=v.mean(), sem=v.std(ddof=1) / np.sqrt(len(v)),
                     mwu_p_vs_COPNFOL=np.nan if ct == "COP/NFOL"
                     else float(mannwhitneyu(v, ref, alternative="two-sided").pvalue),
                     group_P=float(grp.loc[ct, "assoc_mcp"]),
                     cliffs_delta=0.0 if ct == "COP/NFOL" else cliffs(v, ref)))
A = pd.DataFrame(rows)
A["mwu_q_BH5"] = bh(A.mwu_p_vs_COPNFOL.values)
A["group_q_BH5"] = bh(A.group_P.where(A.cell_type != "COP/NFOL").values)
OLD = pd.read_csv(os.path.join(OUT, "fig6A_stats.csv"))
assert np.allclose(A["n"].values, OLD["n"].values) and np.allclose(A["median"].values, OLD["median"].values), \
    "A's statistics must exactly match the delivered fig6A_stats.csv"
print("nuclei %d | statistics match the delivered | %.0fs" % (len(df), time.time() - t0))
data = [df.loc[df.cell_type == ct, "norm_score"].values for ct in SUB]
WB = dict(fc="white", ec="none", alpha=0.9, pad=0.25)


# ---------------------------------------------------------------- panel A
def render_A(w, h, path, ROT_A=90):
    fig, ax = plt.subplots(figsize=(w / CM, h / CM))
    fig.subplots_adjust(left=0.150, right=0.985, top=0.975, bottom=(0.34 if ROT_A == 90 else 0.20))
    bp = ax.boxplot(data, showfliers=False, whis=[5, 95], patch_artist=True, widths=0.55,
                    medianprops=dict(color="black", lw=1.0), whiskerprops=dict(lw=0.7),
                    capprops=dict(lw=0.7))
    for i, (patch, ct) in enumerate(zip(bp["boxes"], SUB)):
        patch.set_facecolor(STALLION[ct]); patch.set_alpha(0.55); patch.set_linewidth(0.6)
        patch.set_edgecolor(dark(STALLION[ct], 0.45))
        xj = np.random.RandomState(7).normal(i + 1, 0.07, len(data[i]))
        ax.scatter(xj, data[i], s=1.5, color=dark(STALLION[ct]), alpha=0.35, lw=0, zorder=3)
    ax.axhline(0, color="grey", ls="--", lw=0.7)
    wmax = max(l.get_ydata().max() for l in bp["whiskers"])
    R1, R2, R3, R4 = wmax + 0.16, wmax + 0.60, wmax + 1.04, wmax + 1.48
    for i, ct in enumerate(SUB):
        x = i + 1
        if ct == "COP/NFOL":
            if STATS == "full":
                for y, lab in ((R1, "MWU P"), (R2, "group P"), (R3, "group q"), (R4, "Cliff δ")):
                    ax.text(x, y, lab, ha="center", va="bottom", fontsize=5.5, color="0.15",
                            bbox=WB, zorder=8)
            continue
        r = A.loc[A.cell_type == ct].iloc[0]
        if STATS == "full":
            for y, txt, fs, col in ((R1, fmtp(r.mwu_p_vs_COPNFOL), 5.6, "0.2"),
                                    (R2, "%.3f" % r.group_P, 5.6, "0.35"),
                                    (R3, "%.3f" % r.group_q_BH5 + ("†" if r.group_q_BH5 < 0.05 else ""),
                                     5.6, "black" if r.group_q_BH5 < 0.05 else "0.45"),
                                    (R4, "%+.2f" % r.cliffs_delta, 5.6, "0.3")):
                ax.text(x, y, txt, ha="center", va="bottom", fontsize=fs, color=col, bbox=WB, zorder=8)
        else:
            ax.text(x, R1, "%+.2f" % r.cliffs_delta + ("†" if r.group_q_BH5 < 0.05 else ""),
                    ha="center", va="bottom", fontsize=6.0, color="black", zorder=8)
        v = data[i]
        ax.errorbar(x, v.mean(), yerr=v.std(ddof=1) / np.sqrt(len(v)), fmt="o", ms=2.4, mfc="white",
                    mec="black", mew=0.6, capsize=1.5, capthick=0.6, lw=0.6, ecolor="black", zorder=6)
    ax.set_xticks(range(1, len(SUB) + 1))
    sep = "\n" if ROT_A == 0 else "  "          # horizontal: two lines (name/sample count); vertical: one line
    ax.set_xticklabels[None] if False else None
    ax.set_xticklabels([ct + sep + "n=%s" % format(int(A.loc[A.cell_type == ct, "n"].iloc[0]), ",")
                        for ct in SUB], fontsize=6.0, rotation=ROT_A, ha="center")
    for lbl, ct in zip(ax.get_xticklabels(), SUB):
        lbl.set_color(dark(STALLION[ct])); lbl.set_fontweight("bold")
    top = R4 + 0.30 if STATS == "full" else wmax + 0.9
    ax.set_ylim(min(-2.9, df.norm_score.min() - 0.1), top)
    ax.set_yticks([t for t in [-2, -1, 0, 1, 2, 3] if t < top])
    ax.set_ylabel("scDRS norm_score\n(MS GWAS)", linespacing=1.35, fontsize=8.0)   # 8.0 tier + relaxed spacing so the two-line box never overlaps
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.02); plt.close(fig)


# ---------------------------------------------------------------- panel B
BLUERED = LinearSegmentedColormap.from_list("bluered", ["#9ECAE1", "#FFFFFF", "#E60000"])


def render_B(w, h, path):
    fig = plt.figure(figsize=(w / CM, h / CM))
    gs = fig.add_gridspec(1, 2, width_ratios=[22, 0.75], wspace=0.045, left=0.10, right=0.88,
                          top=0.92, bottom=0.10)
    ax, axc = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    hb = ax.hexbin(df["UMAP_1"], df["UMAP_2"], C=df["norm_score"], reduce_C_function=np.mean,
                   gridsize=25, mincnt=5, cmap=BLUERED, linewidths=0, edgecolors="none")
    hb.set_rasterized(True)
    vals = hb.get_array()
    v = max(np.nanpercentile(vals, 99), abs(np.nanpercentile(vals, 1)))
    hb.set_clim(-v, v)
    ax.set_xticks([]); ax.set_yticks([])
    if ASPECT == "equal":
        ax.set_aspect("equal")                     # no vertical stretch; the price is blank above and below
    xl, yl = ax.get_xlim(), ax.get_ylim()          # pin the data range first so annotations never stretch the axes
    fig.canvas.draw()
    bb = ax.get_window_extent()
    ar_body = bb.width / bb.height
    ar_data = (xl[1] - xl[0]) / (yl[1] - yl[0])
    print("  B aspect: body %.2f/%.2f = %.2f, data %.2f => vertical %s x %.2f"
          % (bb.width, bb.height, ar_body, ar_data, "compressed" if ar_body > ar_data else "stretched",
             (ar_data / ar_body) if ar_body < ar_data else (ar_body / ar_data)))
    for nm, du, dv, ha, col in [("COP/NFOL", 0.10, 0.16, "left", "#084594"),
                                 ("MOL6", -0.06, 0.10, "right", "#B30000"),
                                 ("MOL2", -0.10, -0.12, "right", "#B30000")]:
        s = df[df.cell_type == nm]
        px, py = s.UMAP_1.median(), s.UMAP_2.median()
        u, vv = (px - xl[0]) / np.diff(xl), (py - yl[0]) / np.diff(yl)
        ax.annotate(nm, xy=(px, py), xytext=(u + du, vv + dv), xycoords="data",
                    textcoords="axes fraction", ha=ha, va="center", fontsize=6.0,
                    fontweight="bold", color=col, bbox=WB,
                    arrowprops=dict(arrowstyle="->", color=col, lw=0.8, shrinkA=1, shrinkB=1.5))
    cb = fig.colorbar(hb, cax=axc)
    cb.outline.set_linewidth(0.5)
    cb.set_ticks([round(-v, 2), 0, round(v, 2)])
    cb.ax.tick_params(labelsize=6.0, length=2, width=0.5, pad=1.5)
    cb.set_label("norm_score", fontsize=6, labelpad=1.5)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.02); plt.close(fig)


# ---------------------------------------------------------------- panel D
DT = pd.read_csv(os.path.join(OUT, "fig6D_delta1_74traits_per_tech.csv"))
rk = DT.pivot(index="trait", columns="tech", values="rank_rna_adj")[list(TECH.values())]
order = rk.mean(axis=1).sort_values().index.tolist()
Pd = {t: DT[DT.tech == t].set_index("trait").loc[order] for t in TECH.values()}
NT = len(order)
LIM = 1.40; RMAX = float(NT)
CD = mpl.colormaps["coolwarm"].copy(); CD.set_bad("none")
CR = mpl.colormaps["cividis_r"].copy(); CR.set_bad("none")


def z(a):
    a = np.asarray(a, float); return (a - a.mean()) / a.std(ddof=0)


def render_D(w, h, path):
    TOP, BOT, LM, GAPC, RMC, WG, FS = 0.42, 0.30, 1.42, 0.30, 1.02, 0.13, 8.0   # 8.0 pt = the delivered D's tier
    colw = (w - LM - GAPC - RMC - 2 * WG) / 12.0
    RH = (h - TOP - BOT) / NT
    xs = []
    for i in range(3):
        xs += [colw, colw]
        if i < 2:
            xs.append(WG)
    xedge = np.concatenate([[0.0], np.cumsum(xs)])
    YEDGE = np.arange(NT + 1) * RH
    POS = [(3 * i, 3 * i + 1) for i in range(3)]
    TW, TH = w, NT * RH + TOP + BOT
    fig = plt.figure(figsize=(TW / CM, TH / CM))

    def axpos(x0):
        return [x0 / TW, BOT / TH, xedge[-1] / TW, (TH - TOP - BOT) / TH]
    for ch, ttl, ax in (("rna", "Gene expression", fig.add_axes(axpos(LM))),
                        ("atac", "Gene score", fig.add_axes(axpos(LM + xedge[-1] + GAPC)))):
        Cd = np.full((NT, len(xs)), np.nan); Ck = np.full((NT, len(xs)), np.nan)
        key, rc = ("d1_rna", "rank_rna_adj") if ch == "rna" else ("d1_atac", "rank_atac_adj")
        for i, t in enumerate(TECH.values()):
            Cd[:, POS[i][0]] = z(Pd[t][key].to_numpy())
            Ck[:, POS[i][1]] = 1.0 - (Pd[t][rc].to_numpy() - 1) / (RMAX - 1)
            Xe, Ye = np.meshgrid(xedge, YEDGE)
            ax.pcolormesh(Xe, Ye, Cd, cmap=CD, vmin=-LIM, vmax=LIM, shading="flat", zorder=2)
            ax.pcolormesh(Xe, Ye, Ck, cmap=CR, shading="flat", zorder=2)
        ax.set_xlim(0, xedge[-1]); ax.set_ylim(YEDGE[-1], 0.0)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        for i in range(3):
            x0, x1 = xedge[POS[i][0]], xedge[POS[i][1] + 1]
            ax.add_patch(Rectangle((x0, 0), x1 - x0, YEDGE[-1], fill=False, ec="#666666",
                                   lw=0.5, clip_on=False, zorder=4))
        ax.add_patch(Rectangle((0, 0), xedge[-1], YEDGE[-1], fill=False, ec="#333333", lw=0.7,
                               clip_on=False))
        ax.set_title(ttl, fontsize=FS, color="black", pad=2)
        for i in range(3):
            cx = (xedge[POS[i][0]] + xedge[POS[i][1] + 1]) / 2
            ax.text(cx, YEDGE[-1] + 0.06, NAME[i], ha="center", va="top", fontsize=FS, color="black")
        for tr, nm in MARK.items():
            if tr not in order or ch != "rna":
                continue
            r = order.index(tr); tc = TRAIT_COLOR[nm]; yc = YEDGE[r] + RH / 2
            ax.plot([-0.32, -0.10], [yc, yc], color=tc, lw=1.0, clip_on=False, zorder=7)
            ax.text(-0.38, yc, nm, ha="right", va="center", fontsize=FS, color=tc)
    bx = (LM + 2 * xedge[-1] + GAPC + 0.30) / TW; bw = 0.11 / TW
    for y0, cmap, lab, ticks, fmt in ((0.575, CD, "Δ1 (z)", [-LIM, 0, LIM], lambda v, p: "%+.1f" % v),
                                      (0.185, CR, "Rank", [0, 0.5, 1], None)):
        cax = fig.add_axes([bx, y0, bw, 0.245])
        sm_ = mpl.cm.ScalarMappable(cmap=cmap, norm=(mpl.colors.Normalize(0, 1) if cmap is CR
                                                     else mpl.colors.Normalize(-LIM, LIM)))
        sm_.set_array([-LIM, LIM] if cmap is CD else [0, 1])
        mpl.colorbar.ColorbarBase(cax, mappable=sm_, orientation="vertical", ticks=ticks)
        cax.yaxis.set_ticks_position("right"); cax.yaxis.set_label_position("right")
        if fmt:
            cax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(fmt))
        else:
            cax.set_yticklabels([str(NT), str(NT // 2), "1"])
        cax.tick_params(labelsize=FS, length=1.6, width=0.5, pad=1.4)
        for s in cax.spines.values():
            s.set_linewidth(0.5)
        cax.text(0.5, 1.16, lab, transform=cax.transAxes, ha="center", va="bottom",
                 fontsize=FS, color="black")
    fig.savefig(path, bbox_inches="tight", pad_inches=0.02); plt.close(fig)


# ---------------------------------------------------------------- iterate to the target frames
import fitz
from arial_style import patch_fitz_arial
patch_fitz_arial(fitz)


def size_of(p):
    d = fitz.open(p); r = d[0].rect; d.close(); return r.width * CM / 72, r.height * CM / 72


def fit(fn, wtarget, htarget, name):
    path = os.path.join(DST, name)
    w, h = wtarget, htarget
    for it in range(4):
        fn(w, h, path)
        aw, ah = size_of(path)
        if abs(ah - htarget) < 0.03 and abs(aw - wtarget) < 0.05:
            break
        w *= wtarget / aw; h *= htarget / ah
    aw, ah = size_of(path)
    print("  %-28s %.2f x %.2f cm (target %.2f x %.2f)" % (name, aw, ah, wtarget, htarget))
    return path, aw, ah


H = 5.40
GAP = 0.24
PANELS = os.environ.get("PANELS", "ABD")
TOTAL = float(os.environ.get("TOTAL", 17.86))
WA = 5.90 if STATS == "full" else 4.90
WD = 6.90
if PANELS == "AD":                                  # UMAP dropped: A/D both return to comfortable widths
    WA, WD = 8.40, TOTAL - 8.40 - GAP
    WB_ = 0.0
else:
    WB_ = TOTAL - WA - WD - (2 if "B" in PANELS else 1) * GAP
print("PANELS=%s STATS=%s ASPECT=%s  targets: A %.2f + B %.2f + D %.2f + gaps = %.2f cm high %.2f"
      % (PANELS, STATS, ASPECT, WA, WB_, WD, WA + WB_ + WD + (2 if "B" in PANELS else 1) * GAP, H))
os.makedirs(DST, exist_ok=True)
sfx = ("_%s" % PANELS) + ("_full" if STATS == "full" else "_slim") + ("" if ASPECT == "auto" else "_eq")
rots = {}
pa, wa_, ha_ = fit(lambda w, h, p: render_A(w, h, p, ROT_A=(0 if w > 7.0 else 90)),
                   WA, H, "fig6v36_A_box_row%s.pdf" % sfx)
if "B" in PANELS:
    pb, wb2, hb_ = fit(render_B, WB_, H, "fig6v36_B_umap_row%s.pdf" % sfx)
pd_, wd_, hd_ = fit(render_D, WD, H, "fig6v36_D_ladder_row%s.pdf" % sfx)

items = [(pa, wa_, ha_), (pd_, wd_, hd_)]
if "B" in PANELS:
    items = [(pa, wa_, ha_), (pb, wb2, hb_), (pd_, wd_, hd_)]
PT = 72 / 2.54                                   # cm -> points (fitz works in points)
doc = fitz.open()
tw = sum(i[1] for i in items) + (len(items) - 1) * GAP
HP = max(i[2] for i in items)                    # page height = the tallest panel; NO scaling, or the type would drift from 8.0 to 7.97
page = doc.new_page(width=tw * PT, height=HP * PT)
x = 0.0
for k, (p, pw, ph) in zip(PANELS, items):
    src = fitz.open(p)
    sc = 1.0
    page.show_pdf_page(fitz.Rect(x * PT, (HP - ph) * PT, (x + pw) * PT, HP * PT), src, 0)
    src.close()
    page.insert_text((x * PT + 1, 9.5), k, fontname="Arial-Regular", fontfile=FONTS_PDF["normal"],
                     fontsize=8, color=(0, 0, 0))
    x += pw * sc + GAP
OUTPDF = os.path.join(DST, "Figure_6_v36_row%s.pdf" % sfx)
doc.save(OUTPDF)
d = fitz.open(OUTPDF)
print("composed %s: %.2f x %.2f cm" % (os.path.basename(OUTPDF),
                                   d[0].rect.width * 2.54 / 72, d[0].rect.height * 2.54 / 72))
d[0].get_pixmap(dpi=300).save(OUTPDF.replace(".pdf", ".png"))
