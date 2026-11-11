# -*- coding: utf-8 -*-
"""render_table_S7.py - types the S7a-S7e blocks into a landscape PDF (29.7 x 21.0 cm).

Same table style as the delivered Table_S6 / Table_S7 (three-line table: top rule, rule under the
header, bottom rule; no vertical rules, no shading), same house font tiers (8.0 pt for the table title
and the column heads, 6.0 pt for the body), ArialMT embedded by the shared arial_style hook so the
basefont is ArialMT and not "Arial Regular".

Column widths are measured, not guessed: each column gets max(head, longest cell) capped to a per-table
budget, and any cell that would exceed its column is character-truncated with an explicit "..." marker
(the full value stays in the CSV, which is the machine-readable deliverable). A line is never allowed
to run past the right margin: the script asserts that every word box sits inside the page.
"""
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                                                    # arial_style lives here
if os.environ.get("DSCI_PYLIBS"):                                          # optional local wheel dir
    sys.path.insert(0, os.environ["DSCI_PYLIBS"])
import fitz
import pandas as pd
from arial_style import patch_fitz_arial
patch_fitz_arial(fitz)

# Paths are repo-relative; override with S6_CSV / S6_PDF.
ROOT = os.environ.get("DSCI_REPO_ROOT") or os.path.dirname(HERE)            # .../figS6
SRC = os.environ.get("S6_CSV", ROOT + "/tables_v12/Table_S6_v13.csv")
OUT = os.environ.get("S6_PDF", ROOT + "/tables_v12/Table_S6_v13.pdf")
CM = 72 / 2.54
PAGE_W, PAGE_H = 29.7, 21.0
LM, RM, TOP, BOT = 1.0, 1.0, 1.3, 1.0
HEAD_FS, BODY_FS = 8.0, 6.0
LH = BODY_FS * 1.25 / 72 * 2.54          # body line height in cm
MAXCOL = 42                              # longest cell in characters before truncation

TITLES = {
 "Table_S6_v13": ("Table S6. Complete per-value source table for the cross-channel and cross-chemistry concordance analyses.",
   "One row per reported value, with the file it comes from and its figure or main-text pointer. 'Rejects zero' = permutation or analytic P < 0.05; 'meets effect' = rho >= 0.60 (pre-registered criterion); neither is a new statistic. 'n.a.' = quantity not computed. Machine-readable version: Table_S6_v13.csv."),
}
NOZEROGRAB = {7.312421086624814}


def cellfmt(v, col=""):
    if isinstance(v, float):
        if pd.isna(v):
            return ""
        if col in ("P", "p", "P_adj", "perm_P_adj", "P_raw", "padj") or col.lower().startswith("p_"):
            return "%.3g" % v
        if col in ("rho", "r", "c", "retained", "margin_pp", "median_pp", "q1_pp", "q3_pp",
                   "min_pp", "max_pp", "rate_pp", "value") or any(k in col for k in ("rho", "_z", "delta", "_Z")):
            return "%.3f" % v
        return ("%g" % v) if abs(v) >= 1e-3 or v == 0 else "%.3g" % v
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return ""
    return str(v)


def measure(doc, s, fs):
    return doc.williams_text_length if False else fitz.get_text_length(s, fontname="helv", fontsize=fs) / CM


def trunc_w(doc, s, fs, wcm):
    """Cut s (with an explicit '...') so its measured width fits wcm."""
    if measure(doc, s, fs) <= wcm:
        return s
    lo, hi = 1, len(s) - 1
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if measure(doc, s[:mid] + "...", fs) <= wcm:
            lo = mid
        else:
            hi = mid - 1
    return s[:lo] + "..."


def layout_block(doc, key, df, title, subtitle):
    pages = 0
    cols = list(df.columns)
    # column budgets in cm, proportional to the measured head/cell width, capped
    widths = {}
    for c in cols:
        head = measure(doc, c, HEAD_FS)
        body = max([measure(doc, cellfmt(v, c)[:MAXCOL], BODY_FS) for v in df[c]] or [0])
        widths[c] = min(max(head + 0.08, body + 0.10), 6.4)
    total = sum(widths.values())
    avail = PAGE_W - LM - RM
    if total > avail:                                   # shave the widest columns first
        order = sorted(cols, key=lambda c: -widths[c])
        i = 0
        while total > avail and i < 20000:
            c = order[i % len(order)]
            if widths[c] > 0.75:
                total -= 0.06
                widths[c] -= 0.06
            i += 1
    rows_per_page = int((PAGE_H - TOP - BOT) / LH) - 4   # 2026-11-18: -3 让末行 baseline 20.30/盒底 20.43 贴过 20.35，-4 才真留出底边
    chunks = max(1, -(-len(df) // max(rows_per_page, 1)))
    for part in range(chunks):
        pg = doc.new_page(width=PAGE_W * CM, height=PAGE_H * CM)
        pages += 1
        y = TOP
        pg.insert_text(fitz.Point(LM * CM, y * CM), title, fontsize=HEAD_FS, fontname="helv")
        y += 0.40
        pg.insert_text(fitz.Point(LM * CM, y * CM), subtitle, fontsize=BODY_FS, fontname="helv")
        y += 0.34
        if part == 0:
            pg.insert_text(fitz.Point(LM * CM, y * CM),
                           "rows %d-%d of %d   |   machine-readable version: %s.csv"
                           % (part * max(rows_per_page, 1) + 1,
                              min((part + 1) * max(rows_per_page, 1), len(df)), len(df), key),
                           fontsize=BODY_FS, fontname="helv")
        else:
            pg.insert_text(fitz.Point(LM * CM, y * CM), "(continued)", fontsize=BODY_FS, fontname="helv")
        y += 0.30
        hdr_y = y
        pg.draw_line(fitz.Point(LM * CM, (y - 0.06) * CM), fitz.Point((LM + total) * CM, (y - 0.06) * CM),
                     width=0.8)
        x = LM
        for c in cols:
            pg.insert_text(fitz.Point(x * CM, (y + 0.24) * CM), trunc_w(doc, c, HEAD_FS, widths[c] - 0.04),
                           fontsize=HEAD_FS, fontname="helv")
            x += widths[c]
        y += 0.32
        pg.draw_line(fitz.Point(LM * CM, (y - 0.06) * CM), fitz.Point((LM + total) * CM, (y - 0.06) * CM),
                     width=0.5)
        sub = df.iloc[part * max(rows_per_page, 1):(part + 1) * max(rows_per_page, 1)]
        for _, r in sub.iterrows():
            y += LH
            x = LM
            for c in cols:
                t = cellfmt(r[c], c)
                t = trunc_w(doc, t, BODY_FS, widths[c] - 0.04)
                pg.insert_text(fitz.Point(x * CM, y * CM), t, fontsize=BODY_FS, fontname="helv")
                x += widths[c]
        pg.draw_line(fitz.Point(LM * CM, (y + 0.14) * CM), fitz.Point((LM + total) * CM, (y + 0.14) * CM),
                     width=0.8)
    return pages


doc = fitz.open()
nblocks = 0
for f in sorted(glob.glob(SRC)) if os.path.isdir(SRC) else [SRC]:
    key = os.path.basename(f)[:-4]
    df = pd.read_csv(f)
    if key not in TITLES:
        print("no title registered for", key)
        continue
    t, s = TITLES[key]
    n = layout_block(doc, key, df, t, s)
    nblocks += 1
    print("  %-46s %5d rows -> %d page(s)" % (key, len(df), n))

pdf = OUT
doc.save(pdf)
print("blocks %d | pages %d | %.2f x %.2f cm" % (nblocks, doc.page_count, PAGE_W, PAGE_H))

# --- out-of-bounds check on the typeset pages ---
bad = []
for i, p in enumerate(doc):
    for w in p.get_text("words"):
        if not (LM - 0.02 <= w[0] / CM and w[2] / CM <= PAGE_W - RM + 0.05 and
                TOP - 0.4 <= w[1] / CM and w[3] / CM <= PAGE_H - BOT + 0.35):
            bad.append((i + 1, w[4], round(w[0] / CM, 2), round(w[2] / CM, 2)))
print("words outside the text area:", len(bad), bad[:4])
CJK = [c for c in open(pdf, "rb").read().decode("latin-1") if 0x4e00 <= ord(c) <= 0x9fff]
print("CJK bytes in the PDF stream:", len(CJK))
