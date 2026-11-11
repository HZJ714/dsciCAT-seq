# -*- coding: utf-8 -*-
"""FigS6 v16 page compose - FOUR rows (the original S6 row-1/row-2 layout + the new layout's row-3/row-4)
=========================================================================
row 1  A locus tracks 10.78 | B window-permutation bar 5.38  <- the original S6a/S6b wide-narrow pair; a histogram never takes half a page
row 2  C bin-size UMAP 17.00 full width  <- the Fig 6B supplement (the original S6c)
row 3  D top-1000 per-gene landing | E 38-gene landing (8.38 cells each, column heads aligned with row 4)
row 4  F RNA per-nucleus boxes | G ATAC per-nucleus boxes (8.38 cells each)
Rules: panels at their original size, ZERO scaling; GAP=0.24; inter-row V_GAP=0.34; top-aligned within a row;
      the leftmost column's letters share x=0.04 (§14.10), the rest sit at their own cell heads (§14.3).
Usage: DSCI_WORK_ROOT=<work> python3 figS6v16_compose.py
"""
import os
import fitz
from arial_style import FONTS_PDF, patch_fitz_arial
patch_fitz_arial(fitz)

HERE = os.path.dirname(os.path.abspath(__file__))
# Default = three levels above the repo GitHub/figS6/assembly (the analysis-tree root); override with DSCI_WORK_ROOT.
WORK = os.environ.get("DSCI_WORK_ROOT", os.path.abspath(os.path.join(HERE, "..", "..", "..")))
PAN = os.environ.get("S6_OUT", WORK + "/figures/routeA_out/P0_8C/panels")
OUTD = os.environ.get("S6_SHEET", WORK + "/figures/routeA_out/S6_v16_20261117")
os.makedirs(OUTD, exist_ok=True)
PT = 72 / 2.54
CM = 1.0 / PT
GAP, VGAP, LETX = 0.24, 0.34, 0.04
MAXW = 17.00
# Page height is a SOFT gate: house spec only hard-constrains width (two-column <= 17.00); the S-family precedent S7 v15 = 17.00 x 17.90.
# The author (2026-11-17) laid S6 out in four rows (A|B / C / D|E / F|G); row 4 must hold "two-channel per-nucleus boxes + the 74-trait
# dual heatmap", naturally ~5.8 cm, the whole page ~19.8 cm. Still one page, but beyond 18.00 it warns explicitly instead of passing silently.
SOFT_H, MAXH = 18.00, 20.50

# (file, letter, cell width None=use the ink width); a row may prefix a "gap from the previous row", default VGAP.
# Row 3 switched to z display (PANEL_ZMODE=1) + externalized legend (PANEL_LEGEND=0),
# the two cells share the 17.00 x 1.05 legend band below: attached closer than the regular row pitch (0.15),
# so the eye reads "this legend belongs to the two cells above", not as a standalone panel.
LEG_GAP = 0.15
ROWS = [
    [("figS6v16_A_loci.pdf", "A", None), ("figS6v16_B_perm.pdf", "B", None)],
    [("figS6v16_C_umapbins.pdf", "C", None)],
    # Round 5 (author 2026-11-17): the horizontal legend band became a vertical column at this row's right edge; D/E shift left overall.
    # The 1.20 cm saved = the legend row height 1.05 + attachment gap 0.15. The three widths are solved by figS6v18_row3.sh,
    # so 2 x PW + 2 x GAP + LW = 17.00 holds exactly (not a shrink, a re-render: SKILL §3.3 zero vector scaling).
    [("fig6v32_C_top1000_landing_s6v18D.pdf", "D", None),
     ("fig6v31_C_landing_s6v18E.pdf", "E", None),
     ("figS6v18_DE_legend_v.pdf", None, None)],
    # F = RNA+ATAC merged into one cell; G = the delivered Fig6D 74-trait dual heatmap (expression Δ / score ranks).
    # Narrative: D/E ask "which subtype do risk genes land on", F asks "MS risk nucleus by nucleus", G asks "where in the 74 traits".
    [("figS6v16_F_cellboxes.pdf", "F", None), ("fig6D_heatmap_pair_compact_s6v16G.pdf", "G", None)],
]
# Row pitch: with the legend inside the D/E row there is no "attachment row" any more; all four rows use the regular VGAP (LEG_GAP kept only as a historical note)
VGAPS = [VGAP, VGAP, VGAP]


def size(p):
    d = fitz.open(p); r = d[0].rect; w, h = r.width * CM, r.height * CM; d.close(); return w, h


cells, missing = [], []
for r in ROWS:
    line = []
    for name, k, cell in r:
        p = os.path.join(PAN, name)
        if not os.path.exists(p):
            missing.append(name); continue
        w, h = size(p)
        cw = cell if cell else w
        assert w <= cw + 1e-6, "%s ink width %.2f > cell width %.2f" % (name, w, cw)
        line.append((p, k, cw, h))
    if line:
        cells.append(line)

widths = [sum(c[2] for c in r) + GAP * (len(r) - 1) for r in cells]
heights = [max(c[3] for c in r) for r in cells]
W, H = max(widths), sum(heights) + sum(VGAPS[:len(cells) - 1])
assert W <= MAXW + 0.005, "page width %.2f > %.2f" % (W, MAXW)
assert H <= MAXH, "page height %.2f > %.2f (four rows do not fit; split the figure or compress the row-3 landing height)" % (H, MAXH)
# Letter strip: S6's F has a two-line vertical ylabel ("per-cell norm_score / RNA") whose top edge would hit the letter
# (collide_gate measured 'norm_score' x 'F' 0.10 x 0.27 cm). Same treatment as Fig6: the panels move down as a block
# by LSTRIP so the letters own a top strip and all page letters share one baseline. The only cost is +0.34 cm of page height.
LSTRIP = 0.34
doc = fitz.open(); pg = doc.new_page(width=W * PT, height=(H + LSTRIP) * PT)
y = LSTRIP
# Self-calibrating left margin: some panels are tight-cropped, so their leftmost ink is their own vertical ylabel (x=0 in the source page),
# and placed at cell head x=0 the PDF text box pokes 0.002 cm past the paper edge (F's "per-cell norm_score", measured).
# No hard-coded offset: it is measured from the SOURCE PAGE's text boxes, then the row is re-checked <= 17.00 (SKILL §17.6: clipping is invisible).
EDGE = 0.035
for ri, r in enumerate(cells):
    d0 = fitz.open(r[0][0]); ws0 = d0[0].get_text("words")
    x0min = min(w[0] * CM for w in ws0) if ws0 else 0.0; d0.close()
    inset = max(0.0, EDGE - x0min)
    assert sum(c[2] for c in r) + GAP * (len(r) - 1) + inset <= MAXW + 0.005, \
        "row %d exceeds the page width after adding left margin %.3f" % (ri + 1, inset)
    x = inset
    for p, k, cw, h in r:
        src = fitz.open(p)
        pg.show_pdf_page(fitz.Rect(x * PT, y * PT, (x + cw) * PT, (y + h) * PT), src, 0)
        src.close()
        lx = LETX if abs(x) < 1e-6 else x
        if k is not None:      # the shared legend band gets no panel letter
            pg.insert_text((lx * PT + 1, (y - LSTRIP) * PT + 9.5), k, fontname="Arial-Regular",
                           fontfile=FONTS_PDF["normal"], fontsize=8, color=(0, 0, 0))
        x += cw + GAP
    y += max(c[3] for c in r)                 # gap after a row; the last row needs none
    if ri < len(cells) - 1:
        y += VGAPS[ri]
out = os.path.join(OUTD, "Figure_S6_v16.pdf")
doc.save(out, garbage=3, deflate=True); doc.close()
d = fitz.open(out); p = d[0]
print("FigS6 v16  %.2f x %.2f cm | page-width gate %.2f -> %s | page-height soft gate %.2f -> %s"
      % (p.rect.width * CM, p.rect.height * CM, MAXW,
         "PASS" if p.rect.width * CM <= MAXW + 0.005 else "FAIL", SOFT_H,
         "PASS" if p.rect.height * CM <= SOFT_H + 0.005 else
         ("past soft gate (still one page; hard cap %.2f)" % MAXH)))
print("  row heights %s | row widths %s" % (["%.2f" % h for h in heights], ["%.2f" % w for w in widths]))
ws = p.get_text("words")
print("  text rightmost %.2f/%.2f, bottom %.2f/%.2f" % (max(w[2] * CM for w in ws), p.rect.width * CM,
                                                max(w[3] * CM for w in ws), p.rect.height * CM))
p.get_pixmap(dpi=200).save(out.replace(".pdf", ".png"))
if missing:
    print("  missing panels:", missing)
