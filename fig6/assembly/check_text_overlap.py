# -*- coding: utf-8 -*-
"""Hard gate for on-figure text overlap: reads the PDF's span bboxes; intersecting horizontal intervals on the same baseline (y diff < 0.6 pt) are collisions.
Table / multi-column text panels (Fig 6 v18's D) must pass this, or eyeballed column widths always leak.
Usage: python3 check_text_overlap.py <pdf> [<pdf2> ...]
"""
import os, sys
if os.environ.get("DSCI_PYLIBS"): sys.path.insert(0, os.environ["DSCI_PYLIBS"])
import fitz


def overlaps(pdf, tol=0.6, gap=0.35):
    d = fitz.open(pdf); bad = []
    for pno, pg in enumerate(d):
        sp = []
        for b in pg.get_text("dict")["blocks"]:
            for l in b.get("lines", []):
                for s in l["spans"]:
                    if s["text"].strip():
                        sp.append(s)
        for i in range(len(sp)):
            for j in range(i + 1, len(sp)):
                a, c = sp[i], sp[j]
                rot = [q for q in (a, c) if (q["bbox"][3] - q["bbox"][1]) > 2.0 * (q["bbox"][2] - q["bbox"][0])]
                if rot:                                          # vertical text: bbox judgement does not apply, skip (avoid false positives)
                    continue
                # a same-glyph ghost is not a collision: identical text, nearly identical vertical centre, one layer pure white - this is
                # the rank digits' "white underlay outline" (path_effects would turn text into vector outlines, invisible to the text layer,
                # so two layers of real text are used instead). A true collision never shows the "same glyph + pure-white backing" combination.
                if (a["text"].strip() == c["text"].strip()
                        and abs((a["bbox"][1] + a["bbox"][3]) / 2 - (c["bbox"][1] + c["bbox"][3]) / 2) <= 1.0
                        and 0xffffff in (a["color"], c["color"])):
                    continue
                vy = min(a["bbox"][3], c["bbox"][3]) - max(a["bbox"][1], c["bbox"][1])
                hh = min(a["bbox"][3] - a["bbox"][1], c["bbox"][3] - c["bbox"][1])
                if hh <= 0 or vy / hh < 0.30:                    # not on the same line
                    continue
                x0 = max(a["bbox"][0], c["bbox"][0]); x1 = min(a["bbox"][2], c["bbox"][2])
                if x1 - x0 > gap:                                # horizontal overlap exceeds the tolerance
                    bad.append((pno, round(x1 - x0, 1), a["text"][:34], c["text"][:34]))
        d2 = d
    d.close()
    return bad


if __name__ == "__main__":
    rc = 0
    for f in sys.argv[1:]:
        bad = overlaps(f)
        print("%-52s collisions %d" % (f.split("/")[-1], len(bad)))
        for b in bad[:12]:
            print("   overlap %.1fpt | %r  x  %r" % (b[1], b[2], b[3])); rc = 1
    sys.exit(rc)
