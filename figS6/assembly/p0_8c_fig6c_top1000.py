# -*- coding: utf-8 -*-
"""p0_8c_fig6c_top1000.py - Fig6C v32: gene-level two-channel landing heatmap of the MAGMA top1000 (= the Fig 6A/B scoring set)

The author noted that v31 went off domain: A/B use the MS GWAS MAGMA top1000 gene set, while 38/331 are FigS6A/B/E *cis*-window
locus calls - a different domain. Hence the main-figure C domain must match A/B = top1000.

Statistics domain (no circularity: scDRS scores are never touched, so no self-inclusion)
  rows = the 223 top-1000 genes measurable in this dataset (875 mapped - 652 near-zero, meta.json convention)
  display = the tested quantity: Δ[g] = 0.2810·L[MOL3/6] + 0.7190·L[MOL6] − 0, L[g,s] = mean(g,s) − mean(g,COP/NFOL)
  test    = set-mean Δ vs a detection-matched non-MS background (decile-stratified permutation, 10,000 draws) -> p08c_dual_setlevel_matched.csv
"""
import os, sys
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
import fitz

CM = 1 / 2.54
HERE = os.path.dirname(os.path.abspath(__file__))
# Default = three levels above the repo GitHub/figS6/assembly (the analysis-tree root); override with DSCI_WORK_ROOT.
WORK = os.environ.get("DSCI_WORK_ROOT", os.path.abspath(os.path.join(HERE, "..", "..", "..")))
BASE = WORK + "/figures/routeA_out/P0_8C"
# In the repo the notebook points P08C_TABLES/S6_OUT at tables_v16 and fig6_out; the defaults keep the analysis-tree paths.
TAB = os.environ.get("P08C_TABLES", BASE + "/tables")
PAN = os.environ.get("S6_OUT", BASE + "/panels"); os.makedirs(PAN, exist_ok=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arial_style import use_arial, use_two_tier
import fitz as _fz          # column labels need measured text widths; fitz is only imported for real in the legend section, so alias it here


def tw(t, fs=6.0):
    return _fz.get_text_length(t, fontname="helv", fontsize=fs) * 2.54 / 72.0
use_arial(); use_two_tier()

SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
W = np.array([203, 490, 1016, 688, 1006, 2574], float)
fA = W[4] / (W[4] + W[5])
CONTRAST = os.environ.get("PANEL_CONTRAST", "pool")   # mol6 = the main-figure convention: Δ = MOL6 − COP/NFOL
JC = 5 if CONTRAST == "mol6" else None                 # None = cell-count weighting of the two enriched columns (old convention)
HAS_D = os.environ.get("PANEL_DCOLUMN", "1") == "1"    # the rightmost Δ column; the author decided to drop it once the single-column convention was set
SHOW_NAMES = os.environ.get("PANEL_NAMES", "1") == "1" # D (the null panel) stays unnamed
NSUB = 6
DETECT = 10.0
# ---- z display mode (author 2026-11-17, for the supplement) ----
# In the old Δ-vs-COP display the COP/NFOL column was identically 0 (the reference subtracted from itself) - a cell with no information;
# after per-gene standardization across the six subtypes the reference column itself enters the comparison (risk genes lower in COP -> blue),
# the six columns share one scale, so the RNA/ATAC blocks share ONE colour bar and D/E can share one legend.
ZMODE = os.environ.get("PANEL_ZMODE", "0") == "1"
DRAW_LEG = os.environ.get("PANEL_LEGEND", "1") == "1"
ZCLIP = float(os.environ.get("PANEL_ZCLIP", 2.0))
SC = {"rna": (-1.2, 1.2), "atac": (-0.35, 0.35)}      # different units: RNA log2(CP10k+1), ATAC gene score
GOLD, BLACK = "#b8860b", "#000000"
RNAH, ATAC_H = os.path.join(os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')), 'oligo/3k_dsci_oligo_subset.h5ad'), \
               os.path.join(os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')), 'oligo/fig6_rerun_20261019/dualomics_probe/genescore_oligo_5977.h5ad')

# ---------------------------------------------------------------- cell-level 223-gene sub-matrix (cached)
G = pd.read_csv(TAB + "/p08c_genelevel_all6_pool.csv")[["gene_human", "gene", "mean_pct"]]
cache = TAB + "/p08c_dual_celllevel223.npz"
if os.path.exists(cache):
    z = np.load(cache, allow_pickle=True)
    XG, AG, key = z["XG"], z["AG"], z["key"]
    print("using cache p08c_dual_celllevel223.npz", flush=True)
else:
    def grab(path, tag):
        a = ad.read_h5ad(path, backed="r")
        vn = np.array(list(map(str, a.var_names)))
        ci = np.array([int(np.flatnonzero(vn == g)[0]) if (vn == g).any() else -1 for g in G.gene])
        have = np.flatnonzero(ci >= 0)
        key = pd.Series(a.obs["cell_type"].astype(str)).str.upper().values   # keep '/', matching SUB
        acc = []
        for i in range(0, a.shape[0], 1024):
            X = a.X[i:i + 1024][:, ci[have]]
            X = X.toarray() if sp.issparse(X) else np.asarray(X)
            acc.append(X.astype(np.float32))
        full = np.full((a.shape[0], len(ci)), np.nan, np.float32); full[:, have] = np.vstack(acc)
        print("%s built %dx%d (ATAC missing %d genes)" % (tag, a.shape[0], len(ci), (ci < 0).sum()), flush=True)
        return full, key
    XG, key = grab(RNAH, "RNA")
    AG, _ = grab(ATAC_H, "ATAC")
    np.savez_compressed(cache, XG=XG, AG=AG, key=key, genes=G.gene.values)

rows = {s: np.flatnonzero(key == s) for s in SUB}


def summarize(M):
    with np.errstate(invalid="ignore"):
        Mn = np.column_stack([np.nanmean(M[rows[s]], 0) for s in SUB])
        Pc = np.column_stack([(np.nan_to_num(M[rows[s]]) > 0).mean(0) * 100 for s in SUB])
    return Mn, Pc


rng = np.random.default_rng(20261116)
res = {}
for tag, M in (("rna", XG), ("atac", AG)):
    Mn, Pc = summarize(M)
    L = Mn - Mn[:, :1]
    d = L[:, JC] if JC is not None else fA * L[:, 4] + (1 - fA) * L[:, 5]
    bs = []
    for _ in range(500):
        idx = {s: rng.choice(rows[s], len(rows[s]), replace=True) for s in SUB}
        with np.errstate(invalid="ignore"):
            p = np.column_stack([np.nanmean(M[idx[s]], 0) for s in SUB])
        bs.append((p[:, JC] - p[:, 0]) if JC is not None else
                  fA * (p[:, 4] - p[:, 0]) + (1 - fA) * (p[:, 5] - p[:, 0]))
    se = np.nanstd(np.array(bs), 0)
    res[tag] = dict(L=L, d=d, pct=Pc, se=se, ok=np.abs(d) > 2 * se)
    q = np.nanquantile(np.abs(L), [.5, .9, .95, 1.0])
    print("%s: Δ>0 %d/%d | sign-resolvable %d | |Δ| p50/p90 %.3f/%.3f | <10%%-detected cells %d/%d | landing |L| p50/90/95/max %.3f/%.3f/%.3f/%.3f"
          % (tag, (d > 0).sum(), d.size, res[tag]["ok"].sum(), np.nanmedian(np.abs(d)),
             np.nanquantile(np.abs(d), .9), (Pc < DETECT).sum(), Pc.size, *q), flush=True)

order = np.argsort(-res["rna"]["d"])
G = G.iloc[order].reset_index(drop=True)
for t in res:
    for k in ("L", "d", "pct", "se", "ok"):
        res[t][k] = res[t][k][order] if res[t][k].ndim == 2 else res[t][k][order]
NR = len(G)
rho = pd.DataFrame({"r": res["rna"]["d"], "a": res["atac"]["d"]}).corr(method="spearman").iloc[0, 1]
print("RNA x ATAC Δ Spearman = %.3f | sign-agreeing genes %d/%d" % (rho, ((np.sign(res["rna"]["d"]) == np.sign(res["atac"]["d"])) & res["rna"]["ok"] & res["atac"]["ok"]).sum(), NR))

# ---------------------------------------------------------------- layout (the v31 compliant geometry)
PW = float(os.environ.get("PANEL_W", 10.72)); PH = float(os.environ.get("PANEL_H", 9.60))
# In the narrow version (FigS6, 8.38 cm) the title must wrap: one 62-character line @8pt is ~9.2 cm > the 6.84 cm available,
# 'across' actually reached 8.43 cm and ran off the page (matplotlib text defaults to clip_on=False, so it was invisible and gates could not catch it).
WIDE = PW >= 10.0
TITLE = ("MS GWAS top-1000 gene set, per-gene landing across OL subtypes" if WIDE
         else "MS GWAS top-1000 gene set:\nper-gene landing across OL subtypes")
# Title externalized (author 2026-11-17): with PANEL_TITLE=0 no title is drawn on the face; the title text goes into the caption;
# the saved top space goes to the body, not to page height (row pitch is the only readability bottleneck of a 223-row heatmap).
TITLE_ON = os.environ.get("PANEL_TITLE", "1") == "1"
YA = PH - 0.35 if WIDE else PH - 0.20
YG = PH - 0.63 if WIDE else PH - 1.00
BODY_T = (PH - 0.70 if WIDE else PH - 1.08) if TITLE_ON else PH - 0.34
NARROW = PW < 10.0          # FigS6 narrow version: subtype labels hang vertically (45° would cut diagonally into the gene-name column)
BODY_B = float(os.environ.get("PANEL_BODY_B", (BODY_T - NR * 0.030) if not NARROW else 2.60))
RH = (BODY_T - BODY_B) / NR
# With the title externalized the body moves up, so the channel names "RNA"/"ATAC" must rise above the body - otherwise they sit on
# the first row's cells (verify_text_on_data would report DOT-UNDER-TEXT 'RNA').
if not TITLE_ON:
    # With va="baseline" the box BOTTOM edge sits ~0.045 cm below the baseline (the PDF emits boxes by the font's ascender/descender,
    # taller than matplotlib's extent, see SKILL §17.2); +0.03 would be pressed by the first row; top edge = BODY_T+0.25,
    # so the body needs at least 0.34 cm headroom to the page edge, else show_pdf_page trims at the CropBox (§17.6).
    YG = BODY_T + 0.06
LBLW, x_lbl, x_body = 1.14, 0.30, 1.54
CW = float(os.environ.get("PANEL_CW", 0.70))
# Baseline column width CW0: only the identically-zero COP/NFOL was narrowed on the raw scale; under z it carries information,
# and with the Δ column gone the six columns are equal, so it defaults to the same width as the rest (author 2026-11-17).
CW0 = float(os.environ.get("PANEL_CW0", CW if (ZMODE and not HAS_D) else 0.26))
NC = NSUB + (1 if HAS_D else 0)
XOFF = [0.0] + [CW0 + i * CW for i in range(NC - 1)]
BW = XOFF[-1] + CW
x_blk = [x_body, x_body + BW]
LEG_Y = 0.52
# With its own two-line legend the body needed >= 1.9 at the bottom; once the legend moved to the shared band, only the vertical subtype labels
# (longest "COP/NFOL" @6pt ~ 0.93 cm + 0.12 spacing) need to hang, so 1.20 suffices; the saved height all goes to the 223 rows.
assert BODY_B > (1.9 if DRAW_LEG else 1.20), BODY_B
print("body %.2f-%.2f (RH %.4f) | right edge %.2f" % (BODY_B, BODY_T, RH, x_blk[1] + BW))


def pick_labels(gap_cm=0.26, minpct=30.0, cap=14):
    need = max(1, int(np.ceil(gap_cm / RH)))
    cand = np.flatnonzero((res["rna"]["pct"].max(1) >= minpct) & res["rna"]["ok"])
    keep = []
    for i in sorted(cand, key=lambda x: -abs(res["rna"]["d"][x])):
        if all(abs(i - j) >= need for j in keep):
            keep.append(i)
        if len(keep) >= cap:
            break
    return sorted(keep), need


lab_rows, need = pick_labels()
if SHOW_NAMES:
    print("naming %d (spacing >= %d rows): %s" % (len(lab_rows), need, ", ".join(G.gene_human.iloc[lab_rows])))
else:
    print("PANEL_NAMES=0: this panel is unnamed (null panel); %d candidates recorded but not drawn: %s"
          % (len(lab_rows), ", ".join(G.gene_human.iloc[lab_rows])))

def rowz_disp(Ld):
    """Per gene: center on the six-subtype mean, then standardize; the denominator adds a shrinkage term s0 (= the median of per-gene sd).
    Without shrinkage, "nearly flat across the six columns" genes divide by a tiny sd and amplify noise into arbitrary strong colour - the classic z-heatmap trap.
    With shrinkage, flat genes sit at z ~ 0 (grey-white) and only structured genes open up into red/blue.
    The Δ column applies the SAME cell weights to the display → Δz = 0.2810·z[MOL3/6] + 0.7190·z[MOL6],
    still the weighted average of the two displayed columns, so a reader can verify it by eye against the colour bar; the invariant survives."""
    # A gene undetected everywhere is NaN for the whole row (greyed below by pct<DETECT), so centering/spread/shrink
    # must all use the nan* versions: one NaN would drag np.median of s0 to NaN and whiten the whole heatmap.
    mu = np.nanmean(Ld, 1, keepdims=True); sd = np.nanstd(Ld, 1, ddof=0, keepdims=True)
    s0 = float(np.nanmedian(sd))
    assert np.isfinite(s0) and s0 > 0, ("shrinkage s0 not finite", s0)
    Z = (Ld - mu) / np.sqrt(sd ** 2 + s0 ** 2)
    fz = W / W.sum()
    return Z, (Z[:, JC] if JC is not None else fz[4] * Z[:, 4] + fz[5] * Z[:, 5]), s0


if ZMODE:
    for _t in ("rna", "atac"):
        _Z, _DZ, _s0 = rowz_disp(res[_t]["L"])
        res[_t]["L"], res[_t]["d"] = _Z, _DZ
        print("[zmode] %s shrinkage s0 = %.4g | display switched to per-gene across-subtype z" % (_t, _s0), flush=True)
    SC = {"rna": (-ZCLIP, ZCLIP), "atac": (-ZCLIP, ZCLIP)}

fig = plt.figure(figsize=(PW / 2.54, PH / 2.54))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, PW); ax.set_ylim(0, PH); ax.axis("off")
CMAP = LinearSegmentedColormap.from_list("by", ["#2166ac", "#f7f7f7", "#b2182b"])
if TITLE_ON:
    ax.text(x_body, YA, TITLE, fontsize=8.0, ha="left",
            va="baseline" if WIDE else "top", linespacing=1.3, color=BLACK)


def cell_y(i):
    return BODY_T - (i + 1) * RH


for nm, x0 in (("RNA", x_blk[0]), ("ATAC", x_blk[1])):
    LO, HI = SC[nm.lower()]
    R = res[nm.lower()]; L, d, pct, ok = R["L"], R["d"], R["pct"], R["ok"]
    for i in range(NR):
        y, h = cell_y(i), RH - 0.003
        for j in range(NSUB):
            v, wx, ww = L[i, j], x0 + XOFF[j], (CW0 if j == 0 else CW)
            if pct[i, j] < DETECT:
                ax.add_patch(Rectangle((wx, y), ww, h, fc="#d9d9d9", ec="none"))
            else:
                t = max(min(v, HI), LO)
                ax.add_patch(Rectangle((wx, y), ww, h, fc=CMAP((t - LO) / (HI - LO)), ec="none"))
                if abs(v) > HI:   # the reference column (j=0) is marked too: under z it is no longer identically 0 and all off-scale cells sit there (the old j>0 was a leftover that swallowed every mark)
                    ax.plot(wx + (ww - 0.05 if v > 0 else 0.05), y + h / 2, marker=">" if v > 0 else "<",
                            ms=1.2, color="#ffffff", mew=0.3, zorder=6)
        if not HAS_D:
            continue
        wx = x0 + XOFF[NC - 1]
        if ok[i]:
            t = max(min(d[i], HI), LO)
            ax.add_patch(Rectangle((wx, y), CW, h, fc=CMAP((t - LO) / (HI - LO)), ec="none"))
            if abs(d[i]) > HI:
                ax.plot(wx + (CW - 0.05 if d[i] > 0 else 0.05), y + h / 2, marker=">" if d[i] > 0 else "<",
                        ms=1.2, color="#ffffff", mew=0.3, zorder=6)
        else:
            ax.add_patch(Rectangle((wx, y), CW, h, fc="#d9d9d9", ec="none"))
    _n = int(((np.abs(L) > HI) & (pct >= DETECT)).sum())
    print("  %s block off-scale triangles %d (of which reference column COP/NFOL %d) | |z|>%.1f and not greyed"
          % (nm, _n, int(((np.abs(L[:, 0]) > HI) & (pct[:, 0] >= DETECT)).sum()), HI), flush=True)
    ax.add_patch(Rectangle((x0, BODY_B), BW, BODY_T - BODY_B, fc="none", ec="#333333", lw=0.4))
    if HAS_D:
        ax.add_patch(Rectangle((x0 + XOFF[NC - 1], BODY_B), CW, BODY_T - BODY_B, fc="none", ec=BLACK, lw=0.7))
    ax.text(x0 + BW / 2, YG, nm, fontsize=6.0, ha="center", va="baseline", color=BLACK)
    for j, s in enumerate(SUB + (["\u0394"] if HAS_D else [])):
        ww = CW0 if j == 0 else CW
        if NARROW:
            # The MEASURED semantics of rotation=90 + ha=right + va=center: the anchor lands at the text's vertical centre,
            # horizontal right edge. To hang it top-aligned under the body, drop the anchor by L/2 per label's measured length.
            ax.text(x0 + XOFF[j] + ww / 2 + 0.115, BODY_B - 0.12 - tw(s) / 2.0, s,
                    fontsize=6.0, rotation=90, ha="right", va="center", color=BLACK)
        else:
            ax.text(x0 + XOFF[j] + ww / 2, BODY_B - 0.07, s, fontsize=6.0, rotation=45,
                    ha="right", va="top", color=BLACK)
    for k in ((JC,) if JC is not None else (4, 5)):
        ax.add_patch(Rectangle((x0 + XOFF[k] + 0.06, BODY_B - 0.055), CW - 0.12, 0.035, fc=GOLD, ec="none"))
    if not HAS_D:   # with no Δ column, "which is the numerator and which the denominator" rests entirely on these two identity bars
        ax.add_patch(Rectangle((x0 + XOFF[0] + 0.06, BODY_B - 0.055), CW0 - 0.12, 0.035, fc=BLACK, ec="none"))

for i in (lab_rows if SHOW_NAMES else []):
    ax.text(x_lbl + LBLW - 0.04, cell_y(i) + (RH - 0.003) / 2, G.gene_human[i], fontsize=6.0,
            ha="right", va="center", color=BLACK)

# ---------------------------------------------------------------- legend
def tw(t, fs=6.0):
    return fitz.get_text_length(t, fontname="helv", fontsize=fs) * 2.54 / 72.0


def legend(x0, x1, y, items, gap=0.42):
    MKW = {"grad": 1.45, "sw": 0.26, "car": 0.30, "bar": 0.42}
    wid = [MKW[k] + (0.0 if k == "grad" else 0.06 + tw(t)) for k, t in items]
    tot = sum(wid) + gap * (len(items) - 1)
    if x0 + tot > x1:
        raise SystemExit("legend too wide %.2f > %.2f" % (x0 + tot, x1))
    x = x0 + (x1 - x0 - tot) / 2.0
    for (kind, t), w in zip(items, wid):
        if kind == "grad":
            ax.imshow(np.linspace(HI, LO, 64)[None, :], extent=(x, x + MKW[kind], y - 0.075, y + 0.075),
                      cmap=CMAP, aspect="auto", zorder=3)
            ax.add_patch(Rectangle((x, y - 0.075), MKW[kind], 0.15, fc="none", ec="#333333", lw=0.3))
            lo_, hi_ = SC[{"RNA log2": "rna", "ATAC score": "atac", "z": "rna"}[t]]
            for px, lab in ((x, "%+.2f" % hi_ if hi_ < 1 else "%+.1f" % hi_),
                            (x + MKW[kind] / 2, "0"),
                            (x + MKW[kind], ("\u2212%.2f" % abs(lo_)) if abs(lo_) < 1 else ("\u2212%.1f" % abs(lo_)))):
                ax.text(px, y - 0.115, lab, fontsize=6.0, ha="center", va="top", color=BLACK)
            ax.text(x + MKW[kind] / 2, y + 0.115, t, fontsize=6.0, ha="center", va="bottom", color=BLACK)
        elif kind == "sw":
            ax.add_patch(Rectangle((x, y - 0.075), MKW[kind], 0.15, fc="#d9d9d9", ec="#8c8c8c", lw=0.3))
            ax.text(x + MKW[kind] + 0.06, y, t, fontsize=6.0, ha="left", va="center", color=BLACK)
        elif kind == "car":
            ax.plot(x + 0.07, y, marker="<", ms=2.2, color="#2166ac", mew=0.3)
            ax.plot(x + 0.23, y, marker=">", ms=2.2, color="#b2182b", mew=0.3)
            ax.text(x + MKW[kind] + 0.06, y, t, fontsize=6.0, ha="left", va="center", color=BLACK)
        else:
            ax.add_patch(Rectangle((x, y - 0.055), MKW[kind], 0.075, fc=GOLD, ec="none"))
            ax.text(x + MKW[kind] + 0.06, y, t, fontsize=6.0, ha="left", va="center", color=BLACK)
        x += w + gap
    return x - gap


if not DRAW_LEG:
    print("legend handed to the shared band (PANEL_LEGEND=0)", flush=True)
    xe = x_body
    _it = []
else:
  _it = ([("grad", "z"), ("sw", "not detected"), ("car", "off scale"), ("bar", "tested")]
         if ZMODE else
         [("grad", "RNA log2"), ("grad", "ATAC score"), ("sw", "not detected"),
          ("car", "off scale"), ("bar", "tested")])
if _it and PW < 10.0:      # narrow two-line legend (as above)
    xe = legend(x_body, x_blk[1] + BW, LEG_Y, _it[2:], gap=0.90)
    # the two gradient bars' end ticks centre on the bar ends; too small a gap reads as "−1.2+0.35" (gates silent, eyes misread) -> separate them
    xe = legend(x_body, x_blk[1] + BW, LEG_Y + 0.44, _it[:2], gap=1.10)
elif _it:
    xe = legend(x_body, x_blk[1] + BW, LEG_Y, _it)
if _it:
    print("legend %.2f-%.2f" % (x_body, xe))

out = PAN + "/fig6v32_C_top1000_landing%s.pdf" % os.environ.get("PANEL_TAG", "")
fig.savefig(out, dpi=600)
d = fitz.open(out); p = d[0]
Wp, Hp = p.rect.width * 2.54 / 72, p.rect.height * 2.54 / 72
ws = p.get_text("words")
worst = max((w[2] * 2.54 / 72 for w in ws), default=0); bot = max((w[3] * 2.54 / 72 for w in ws), default=0)
sp_ = [s for b in p.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]]
sizes = sorted({round(s["size"], 2) for s in sp_})
grey = sorted({"%06x" % s["color"] for s in sp_ if s["color"] not in (0, 0xffffff)})
print("size %.2f x %.2f | text right %.2f (limit %.2f) | bottom %.2f | sizes %s | non-black text %s"
      % (Wp, Hp, worst, Wp - 0.20, bot, sizes, grey))
d.close()
