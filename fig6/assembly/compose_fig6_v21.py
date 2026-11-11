import os
# -*- coding: utf-8 -*-
"""Fig 6 v21 compose: D becomes the cross-chemistry heatmap pair (Gene expression | Gene score, each 6 columns = 3 chemistries x [Δ1 | rank]).
C switches to CLEAN=2 (the in-figure five-line statistics and per-trait ranks move to the caption; axis labels use formal terms).
A/B unchanged. The four cells align on a 2x2 equal-width/equal-height grid; this script stamps the panel letters, so the panels must not carry letters.
Output fig6_composed_v21.pdf (a PNG is only written to /tmp-like directories when PNGOUT is set; the delivery directory keeps PDFs only).
"""
import sys, os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # arial_style.py, check_text_overlap.py
if os.environ.get("DSCI_PYLIBS"): sys.path.insert(0, os.environ["DSCI_PYLIBS"])
import fitz
from arial_style import patch_fitz_arial, FONTS_PDF
patch_fitz_arial(fitz)
CM = 72 / 2.54
FIG = os.environ.get("FIG6_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "fig6_out"))   # the panel PDFs the notebook wrote
GAP = float(os.environ.get("GAP", 0.24))
# DPDF/OUT exist only for "candidate-version comparison" (NUMMODE=cell/none outputs must not touch the delivery); unset = the current delivery path
P = {"A": FIG + "/fig6v21_A_box.pdf", "B": FIG + "/fig6v21_B_umap.pdf",
     "C": FIG + "/fig6v21_C_contrast.pdf",
     "D": os.environ.get("DPDF", FIG + "/fig6D_heatmap_pair_compact.pdf")}
sz = {k: (lambda d: (d[0].rect.width, d[0].rect.height, d.close()))(fitz.open(v)) for k, v in P.items()}
w = {k: s[0] / CM for k, s in sz.items()}
h = {k: s[1] / CM for k, s in sz.items()}
row1_y, row2_y = 0.0, max(h["A"], h["B"]) + GAP
total_w = max(w["A"] + GAP + w["B"], w["C"] + GAP + w["D"])
total_h = row2_y + max(h["C"], h["D"])
place = dict(A=(0, row1_y), B=(w["A"] + GAP, row1_y), C=(0, row2_y), D=(w["C"] + GAP, row2_y))
doc = fitz.open()
page = doc.new_page(width=total_w * CM, height=total_h * CM)
for k, (x, y) in place.items():
    src = fitz.open(P[k])
    page.show_pdf_page(fitz.Rect(x * CM, y * CM, x * CM + sz[k][0], y * CM + sz[k][1]), src, 0)
    src.close()
    page.insert_text((x * CM + 1, y * CM + 9.5), k, fontname="Arial-Regular",
                     fontfile=FONTS_PDF["normal"], fontsize=8, color=(0, 0, 0))
OUTPDF = os.environ.get("OUTPDF", FIG + "/fig6_composed_v21.pdf")
doc.save(OUTPDF)
d = fitz.open(OUTPDF)
if os.environ.get("PNGOUT"):
    d[0].get_pixmap(matrix=fitz.Matrix(300 / 72, 300 / 72)).save(os.environ["PNGOUT"] + "/fig6_composed_v21.png")
print("Fig6 v21: %.2f x %.2f cm | row1 %.2f/%.2f (h %.2f) row2 %.2f/%.2f (h %.2f) | row-width diff %.2f row-height diff %.2f"
      % (total_w, total_h, w["A"], w["B"], max(h["A"], h["B"]), w["C"], w["D"], max(h["C"], h["D"]),
         abs((w["A"] + GAP + w["B"]) - (w["C"] + GAP + w["D"])),
         abs(max(h["A"], h["B"]) - max(h["C"], h["D"]))))
