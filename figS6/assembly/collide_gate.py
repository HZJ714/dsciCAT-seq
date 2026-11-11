import os
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""collide_gate.py - the generic text-collision gate (house spec §4's collision gate, independent of any per-figure verifier).

Criterion: two DIFFERENT words on the same page whose bounding boxes overlap with area > 0 and an overlap >= 0.02 cm in both directions
count as one collision. End-to-end touching of neighbouring words on a line (overlap 0) does not - but §9 warns:
**intersections only, not touch**, so borderline cases still need an eyeball pass.

Usage: python3 collide_gate.py a.pdf b.pdf ...   (no arguments: run the four DRAFT_v22 drafts + six new panels)
Exit code: 0 = no collisions; 1 = collisions found.
"""
import sys, os, glob, fitz

CM = 28.3465
MIN_OVERLAP = 0.02          # cm; both directions must exceed it to count as real overlap


def words_of(page):
    H = page.rect.height
    return [(w[0] / CM, (H - w[3]) / CM, (w[2] - w[0]) / CM, (H - w[1]) / CM, w[4])
            for w in page.get_text("words") if w[4].strip()]


def audit(path):
    d = fitz.open(path)
    hits = []
    for pi, pg in enumerate(d):
        ws = words_of(pg)
        ws.sort(key=lambda t: t[0])
        for i, (x0, y0, x1, y1, t0) in enumerate(ws):
            for j in range(i + 1, len(ws)):
                a0, b0, a1, b1, t1 = ws[j]
                if a0 >= x1:
                    break                       # x is sorted; no further overlap is possible
                ox = min(x1, a1) - max(x0, a0)
                oy = min(y1, b1) - max(y0, b0)
                if ox > MIN_OVERLAP and oy > MIN_OVERLAP:
                    hits.append((pi + 1, round(ox, 2), round(oy, 2), t0[:22], t1[:22]))
    d.close()
    return hits


DEFAULT = sorted(glob.glob(os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'routeA_out/DRAFT_v22/pdf/*.pdf'))) + \
          sorted(glob.glob(os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'routeA_out/DRAFT_v22/panels/fig*.pdf')))

if __name__ == "__main__":
    files = [a for a in sys.argv[1:] if a.endswith(".pdf")] or DEFAULT
    bad = 0
    for f in files:
        h = audit(f)
        bad += len(h)
        print("%-58s collisions %d%s" % (os.path.basename(f)[:58], len(h),
                                      "" if not h else "  top 3: %s" % "; ".join(
                                          "p%d '%s'×'%s' (%.2f×%.2f cm)" % (x[0], x[3], x[4], x[1], x[2])
                                          for x in h[:3])))
    print("\ntotal collisions %d across %d files" % (bad, len(files)))
    print("reminder: this gate catches intersections only, not end-to-end touch; pass the gate and still eyeball the borders between neighbours (§9).")
    sys.exit(1 if bad else 0)
