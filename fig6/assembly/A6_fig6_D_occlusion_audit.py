# -*- coding: utf-8 -*-
"""Fig 6 · D heatmap: pixel-level audit of whether the digits cover the heatmap (author 11-03: a backing is legible but hides the heatmap)

No eyeballing: the delivered table + colour map give the colour each rank cell SHOULD have; the rendered bitmap is then measured for how many
pixels went white (>248 in all three channels). White inside a cell that should carry colour (expected-channel sum < 600, excluding the lightest colour-bar end)
means the digit/white pad covered data. It also measures the named rows' white-band height and digit-ink height to judge whether the glyph sits fully inside its own row.

Usage: python3 A6_fig6_D_occlusion_audit.py <D panel PDF> [--tall 0.22]
      with --tall set, row edges follow the "3 named rows grow, the rest shrink proportionally" geometry (the NUMMODE=cell delivered version).
Output: one stdout verdict line + per-row detail (only the polluted rows)
"""
import os, sys

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # arial_style.py, check_text_overlap.py
if os.environ.get("DSCI_PYLIBS"): sys.path.insert(0, os.environ["DSCI_PYLIBS"])
import numpy as np, pandas as pd
import fitz
from PIL import Image
from matplotlib.colors import LinearSegmentedColormap

CM = 2.54
PT = 72 / 2.54
DPI = 1200
TABD = os.environ.get("TABD", os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "fig6_out"))      # fig6D_delta1_74traits_per_tech.csv written by the notebook
TECH = ["dsciCAT-seq (own)", "10x Multiome", "HT-scCAT-seq"]
MARKED_TRAITS = ["PASS_Multiple_sclerosis", "UKB_460K.body_HEIGHTz", "PASS_Schizophrenia_Pardinas2018"]
# Geometry must stay in sync with make_fig6_v21_D_heatmap.py: compact judged by file name, then the same default values
argv = sys.argv[1:]
TALL = float(argv[argv.index("--tall") + 1]) if "--tall" in argv else 0.22   # default = the cell delivery geometry; the halo version passes --tall 0
PDF = next((a for a in argv if not a.startswith("--") and a != "--tall"),
           os.path.join(os.environ.get('DSCI_FY_ROOT', WORK + "/FY"), 'oligo/fig6_rework_v2/fig6D_heatmap_pair_compact.pdf'))
COMPACT = "compact" in os.path.basename(PDF)
WD = float(os.environ.get("WD", 0.55 if COMPACT else 0.92))
WK = float(os.environ.get("WK", 0.55 if COMPACT else 0.92))
RH = float(os.environ.get("RH", 0.058 if COMPACT else 0.112))

# ---- row order and expected colours (audit only the left heatmap's first-chemistry rank column, the easiest to locate; the verdict holds for all 18 columns) ----
D = pd.read_csv(TABD + "/fig6D_delta1_74traits_per_tech.csv")
rk = D.pivot(index="trait", columns="tech", values="rank_rna_adj")[TECH]
order = rk.mean(axis=1).sort_values().index.tolist()
NT = len(order)
CR = LinearSegmentedColormap.from_list("rk", ["#F4F6F8", "#8FA0B8", "#1B2A45"])
vals = 1.0 - (rk.loc[order, TECH[0]].to_numpy() - 1) / (NT - 1)
EXP = (CR(vals)[:, :3] * 255).round(0)
marked = [order.index(t) for t in MARKED_TRAITS if t in order]

# ---- row geometry (cm) ----
if TALL:
    base = (NT * RH - len(marked) * TALL) / (NT - len(marked))
    HROW = np.where(np.isin(np.arange(NT), marked), TALL, base)
else:
    base = RH
    HROW = np.full(NT, RH)
YE = np.concatenate([[0.0], np.cumsum(HROW)])

doc = fitz.open(PDF)
pg = doc[0]
PAGE_H_CM = pg.rect.height / 72 * 2.54   # pt -> cm
PAGE_W_CM = pg.rect.width / 72 * 2.54
pix = pg.get_pixmap(dpi=DPI)
pix.save("/tmp/_occ.png")
doc.close()
im = np.array(Image.open("/tmp/_occ.png").convert("RGB")).astype(int)
H, W, _ = im.shape
pxcm = H / PAGE_H_CM                      # px per cm (whole-page scale; the heatmap and the page share the ratio)
gray = (np.abs(im[:, :, 0] - im[:, :, 1]) < 12) & (np.abs(im[:, :, 1] - im[:, :, 2]) < 12) & (im.sum(axis=2) < 420)
rows = np.where(gray.sum(axis=1) > W * 0.30)[0]
cols = np.where(gray.sum(axis=0) > H * 0.30)[0]
TOP, BOT = float(rows[0]), float(rows[-1])          # heatmap frame top/bottom (px)
LEFT = float(cols[0])                                # frame left = the Δ column's left edge
heat_cm = (BOT - TOP) / pxcm
white = (im[:, :, 0] > 248) & (im[:, :, 1] > 248) & (im[:, :, 2] > 248)
red = np.abs(im - np.array([213, 31, 38])).sum(axis=2) < 90

# left heatmap's first-chemistry rank column: x = LEFT + WD .. LEFT + WD + WK (inset 5 px each side to dodge the column frame lines)
xa, xb = int(LEFT + WD * pxcm) + 5, int(LEFT + (WD + WK) * pxcm) - 5
bad, ink = [], []
GUARD = 6          # row edges carry a 0.5-0.7 pt grey frame stroke (8-12 px at 1200 dpi); step back on both sides so the stroke is not counted into the neighbouring row
for r in range(NT):
    y0, y1 = int(TOP + YE[r] * pxcm) + GUARD, int(TOP + YE[r + 1] * pxcm) - GUARD
    if y1 <= y0:
        continue
    fw = white[y0:y1, xa:xb].mean()
    dark_expected = EXP[r].sum() < 600               # this cell should carry visible colour
    if dark_expected and fw > 0.02:
        bad.append((r, order[r], round(100 * fw, 1)))
    if r in marked:
        w0, w1 = y0 - 20, y1 + 20
        seg = red[w0:w1, xa:xb]
        ys = np.where(seg.any(axis=1))[0]
        if len(ys):
            lo, hi = w0 + ys.min(), w0 + ys.max()
            ink.append((r, "in-row %.0f px  glyph ink %.0f px (%.2f pt)  overflow up %d px / down %d px  %s"
                        % (y1 - y0, hi - lo, (hi - lo) / pxcm * PT,
                           max(0, y0 - lo), max(0, hi - y1),
                           "glyph fully inside its own row" if lo >= y0 and hi <= y1 else "glyph presses the neighbouring row")))

print("PDF %s | heatmap %.3f cm | row geometry: regular %.4f cm (%.1f px) named %.3f cm | rank column x %d-%d"
      % (os.path.basename(PDF), heat_cm, base, base * pxcm, TALL, xa, xb))
print("rows that should carry colour but are covered by white (digit/pad): %d%s" % (len(bad),
      "  -> " + "; ".join("row%d %s %.1f%%" % b for b in bad[:8]) if bad else ""))
for t in ink:
    print("  named row %d: %s" % t)
