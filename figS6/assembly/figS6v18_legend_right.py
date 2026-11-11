#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""figS6v18_legend_right.py - the VERTICAL legend shared by D/E (author 2026-11-17: DE shift left overall, the legend squeezed into the rightmost column).

Why vertical: the horizontal shared legend band (figS6v17, 17.00 x 1.05) took a whole row, page height 20.11 cm.
The vertical version joins row 3 as its third column, saving 1.05 row height + 0.15 attachment gap = **1.20 cm**.
The price: D/E give up this column's width, so the subtype-column width CW shrinks from 0.55 to ~0.48 (the body still has 6 columns x 2 blocks,
cells only narrow 13%, row height and type size unchanged - house spec ranks added whitespace above shrinking type).

The width is not guessed: every text is measured with fitz at its real width, the longest item fixes LW, then CW is solved back,
so that D + E + 2 x GAP + LW is **exactly 17.00** (the same idea as spreading Fig6's column gaps, SKILL §17.9).

Two font sizes 8.0/6.0; face text pure black (the colour-bar end numbers too); gradient direction conventional (blue− left/bottom, red+ right/top).
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import fitz

HERE = os.path.dirname(os.path.abspath(__file__))
# Default = three levels above the repo GitHub/figS6/assembly (the analysis-tree root); override with DSCI_WORK_ROOT.
WORK = os.environ.get("DSCI_WORK_ROOT", os.path.abspath(os.path.join(HERE, "..", "..", "..")))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arial_style import use_arial, use_two_tier

PAN = os.environ.get("S6_OUT", WORK + "/figures/routeA_out/P0_8C/panels"); os.makedirs(PAN, exist_ok=True)
CM = 2.54
BLACK, GOLD, GREY, EDGE = "#000000", "#b8860b", "#d9d9d9", "#333333"
ZCLIP = float(os.environ.get("PANEL_ZCLIP", 2.0))
LW = float(os.environ.get("LEG_W", 2.00))          # canvas width (>= what is needed; tight trims to the ink width)
LH = float(os.environ.get("LEG_H", 5.75))          # matches the D/E row height so they share a top-aligned row
YCENTER = float(os.environ.get("LEG_YC", 3.38))    # the stack's vertical centre: aligns with the body mid-line (1.35+5.41)/2 of both cells
CMAP = LinearSegmentedColormap.from_list("by", ["#2166ac", "#f7f7f7", "#b2182b"])
use_arial(); use_two_tier()


def tw(t, fs=6.0):
    return fitz.get_text_length(t, fontname="helv", fontsize=fs) * CM / 72.0


# At STEP 0.36 a two-line 6pt key (0.42 together) would touch itself; raised to 0.44; stack height 3.97,
# still centred on the body mid-line 3.38 (1.40-5.37, inside the 5.75 page).
BARW, SWH, GAPX, STEP = 0.22, 0.15, 0.07, 0.44
# The off-scale key's swatch size and triangle: MS_OFF is the same 1.2 as the body triangles in the panel scripts,
# not "a bit bigger in the legend" - enlarging would lie: a reader measuring the legend against the face could not match the cells.
CHW, CHH, MS_OFF = 0.16, 0.11, 1.2
TITLE, TBAR = "z (per gene)", 1.60
KEYS = [("sw", "not detected"), ("car", "off scale"),
        ("bar", "tested\n(MOL6)"), ("bar", "reference\n(COP/NFOL)")]
# per-item width = swatch + gap + the widest text line
WID = {"sw": BARW + GAPX + tw("not detected"), "car": 0.34 + GAPX + tw("off scale"),
       "bar": None}
wbar = [BARW + GAPX + max(tw(l) for l in t.split("\n")) for k, t in KEYS if k == "bar"]
wcol = BARW + GAPX + tw("+%.1f" % ZCLIP)
tot = max(wcol, wbar[0], wbar[1], tw(TITLE), BARW + GAPX + tw("not detected"),
          0.34 + GAPX + tw("off scale"))
print("legend needs width %.2f cm (longest item = %s)" % (tot, ["colour scale", "not detected", "off scale",
      "tested(MOL6)", "reference(COP/NFOL)", "title"][
          [wcol, BARW + GAPX + tw("not detected"), 0.34 + GAPX + tw("off scale"),
           wbar[0], wbar[1], tw(TITLE)].index(tot)]))
if LW < tot + 0.10:
    raise SystemExit("LEG_W=%.2f does not fit the longest item %.2f (needs >= %.2f)" % (LW, tot, tot + 0.10))

fig = plt.figure(figsize=(LW / CM, LH / CM))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, LW); ax.set_ylim(0, LH); ax.axis("off")
X0 = 0.10
stack = 0.21 + 0.10 + TBAR + 0.30 + STEP * len(KEYS)      # title + colour scale + gap + four keys
y = YCENTER + stack / 2.0                                   # stack from the top down
ax.text(X0, y - 0.21, TITLE, fontsize=6.0, ha="left", va="bottom", color=BLACK)
y -= 0.21 + 0.10
ax.imshow(np.linspace(-ZCLIP, ZCLIP, 64)[:, None],
          extent=(X0, X0 + BARW, y - TBAR, y), cmap=CMAP, aspect="auto", zorder=3)
ax.add_patch(Rectangle((X0, y - TBAR), BARW, TBAR, fc="none", ec=EDGE, lw=0.3))
for dy, lab in ((0.0, "+%.1f" % ZCLIP), (TBAR / 2.0, "0"), (TBAR, "−%.1f" % ZCLIP)):
    ax.text(X0 + BARW + GAPX, y - dy, lab, fontsize=6.0, ha="left", va="center", color=BLACK)
y -= TBAR + 0.30
for kind, t in KEYS:
    yy = y - STEP / 2.0
    if kind == "sw":
        ax.add_patch(Rectangle((X0, yy - SWH / 2), BARW, SWH, fc=GREY, ec="#8c8c8c", lw=0.3))
    elif kind == "car":
        # The key must look exactly like the face: a WHITE triangle ms=1.2 stuck to the off-scale cell's left edge (negative) / right edge (positive),
        # on the saturated blue/red clipped at ±ZCLIP. This used to draw a 2.2 pt blue/red triangle -
        # matching the face in neither colour nor size (a reader hunting for blue triangles finds none).
        # Now a "swatch + same-size white triangle": swatch width 0.16, the two together still 0.34, so the legend width is unchanged and CW needs no re-solving.
        for k, (cx, fc, mk, mx) in enumerate((
                (X0, CMAP(0.0), "<", 0.030),                 # negative off-scale: blue swatch + white < on the left edge
                (X0 + 0.18, CMAP(1.0), ">", 0.130))):        # positive off-scale: red swatch + white > on the right edge
            ax.add_patch(Rectangle((cx, yy - CHH / 2), CHW, CHH, fc=fc, ec="none"))
            ax.plot(cx + mx, yy, marker=mk, ms=MS_OFF, color="#ffffff", mew=0.3, zorder=6)
    else:
        ax.add_patch(Rectangle((X0, yy - 0.038), BARW, 0.075,
                               fc=GOLD if "tested" in t else BLACK, ec="none"))
    lines = t.split("\n")
    for k, ln in enumerate(lines):        # two-line keys: centred on y as a block
        ax.text(X0 + (BARW if kind != "car" else 0.34) + GAPX,
                yy + (0.09 if len(lines) > 1 and k == 0 else (-0.09 if len(lines) > 1 else 0.0))
                - 0.0, ln, fontsize=6.0, ha="left", va="center", color=BLACK)
    y -= STEP
print("  stack %.2f-%.2f, centre %.2f (target %.2f)" % (y, y + stack, y + stack / 2.0, YCENTER))
out = PAN + "/figS6v18_DE_legend_v.pdf"
fig.savefig(out, dpi=600)
d = fitz.open(out); p = d[0]; ws = p.get_text("words")
sp = [s for b in p.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]]
print("  %-30s page %.2f x %.2f cm | text right %.2f bottom %.2f top %.2f | sizes %s | non-bw %s"
      % (os.path.basename(out), p.rect.width * CM / 72, p.rect.height * CM / 72,
         max(w[2] * CM / 72 for w in ws), max(w[3] * CM / 72 for w in ws),
         min(w[1] * CM / 72 for w in ws), sorted({round(s["size"], 1) for s in sp}),
         [{"#%06x" % s["color"]} for s in sp if s["color"] not in (0, 16777215)]))
d.close()
