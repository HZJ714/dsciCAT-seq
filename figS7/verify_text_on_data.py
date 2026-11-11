#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""verify_text_on_data.py — extend the self-audit from "text vs text" to "text vs data ink".

Why this gate exists: a text-vs-text gate can report zero collisions while seven text blocks
sit on top of data marks, because a text-to-text check cannot see points. So the marks
themselves are tested, with two rules run separately to avoid false positives:

  R1 filled small marks (dots, squares, bars, dumbbell ends) with area 0.0004-0.35 cm2:
     report DOT-UNDER-TEXT if text covers more than 50 % of the mark's area.
  R2 thin segments (links, error bars, reference lines; stroke only, length < 3.0 cm,
     width >= 0.5): report LINE-ACROSS-TEXT if a segment crosses more than 6 % of a text
     cell (a diagonal that only clips a corner is not reported).

Not treated as data marks: anything below 0.0004 cm2 (grid dots, background stipple),
anything above 0.35 cm2 (axvspan fills, colourbar strips), and frame lines that hug the
page edge (spines, tick marks).
"""
import sys, os, fitz
CM = 2.54 / 72.0
DOT_COVER, LINE_COVER = 0.50, 0.06

def classify(page):
    dots, lines = [], []
    for it in page.get_drawings():
        r = fitz.Rect(it["rect"])
        if r.is_empty:
            continue
        a = r.get_area() * CM * CM
        if a < 0.0004 or a > 0.35:
            continue
        if (r.y0 < 2 and r.x0 < 2) or (r.y1 > page.rect.height - 2) or (r.x1 > page.rect.width - 2):
            continue                                    # not data marks: < 0.0004 cm2 (grids / background dots), > 0.35 cm2 (axvspan fills /
        f = it.get("fill")
        if f is not None:
            if min(f) > 0.90:
                continue                            #   colourbar strips), and frame lines hugging the page edge (spines, tick marks)
            dots.append(r)
            continue
        elif it.get("color") is not None:
            ln = (r.width ** 2 + r.height ** 2) ** 0.5 * CM
            if ln < 3.0 and float(it.get("width") or 0.5) >= 0.5:
                lines.append(r)
    return dots, lines

def spans(page):
    out = []
    for b in page.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                if s["text"].strip():
                    out.append((fitz.Rect(s["bbox"]), s["text"].strip()))
    return out

def cover(target, others):
    # near-white fills are a label's own backing box, not data ink
    tot = 0.0
    for o in others:
        ix = fitz.Rect(target); ix.intersect(o)
        if not ix.is_empty:
            tot += ix.get_area()
    return tot / max(target.get_area(), 1e-9)

def audit(pdf):
    d = fitz.open(pdf); bad = 0
    for pi, page in enumerate(d):
        dots, lines = classify(page)
        sp = spans(page)
        for r, t in sp:
            for k in dots:
                if cover(k, [r]) > DOT_COVER:
                    print("DOT-UNDER-TEXT (%d, %.1f cm, %r)" % (pi, r.y0 * CM, t[:34])); bad += 1; break
            if cover(r, lines) > LINE_COVER:
                print("LINE-ACROSS-TEXT (%d, %.1f cm, %r, %.0f%%)"
                      % (pi, r.y0 * CM, t[:34], 100 * cover(r, lines))); bad += 1
    d.close(); return bad

if __name__ == "__main__":
    files = sys.argv[1:]
    tot = 0
    for f in files:
        n = audit(f)
        print("%-42s data marks under text: %d  %s"
              % (os.path.basename(f), n, "FAIL" if n else "OK"))
        tot += (n > 0)
    print("\n%d/%d files pass the text-vs-data gate" % (len(files) - tot, len(files)))
    sys.exit(1 if tot else 0)
