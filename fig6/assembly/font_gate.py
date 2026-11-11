#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""font_gate.py — the two submission gates every delivered page of this paper passes.

1. every visible span is at least FMIN pt (5.5 pt here, the journal floor used throughout);
2. the page carries Arial only (Liberation Sans is accepted as metric-compatible), embedded as
   live text (`pdf.fonttype = 42`), and the page is not wider than WMAX_CM.

Usage:  python3 font_gate.py <page.pdf> [...] [--fmin 5.5] [--wmax 18.0]
Exit code 1 if any page fails. This is the repo-side subset of the checks in
A6_fig6_v21_verify.py; that script additionally re-reads every printed number against the
delivered tables and needs the full working tree.
"""
import os
import sys

if os.environ.get("DSCI_PYLIBS"): sys.path.insert(0, os.environ["DSCI_PYLIBS"])
import fitz  # noqa: E402

OK_FONTS = ("Arial", "LiberationSans", "Liberation")


def gate(path, fmin, wmax):
    doc = fitz.open(path)
    page = doc[0]
    w_cm = page.rect.width / 72 * 2.54
    h_cm = page.rect.height / 72 * 2.54
    spans = [s for b in page.get_text("dict")["blocks"] for l in b.get("lines", [])
             for s in l["spans"] if s["text"].strip()]
    sizes = sorted(round(s["size"], 2) for s in spans)
    fonts = sorted({s["font"] for s in spans})
    doc.close()
    small = [s for s in spans if round(s["size"], 2) < fmin]
    badfont = [f for f in fonts if not any(k in f for k in OK_FONTS)]
    print("%-34s %6.2f x %5.2f cm | %3d spans | font %.2f-%.2f pt | faces %s"
          % (os.path.basename(path), w_cm, h_cm, len(spans),
             sizes[0] if sizes else float("nan"), sizes[-1] if sizes else float("nan"),
             ",".join(fonts)))
    fail = []
    if small:
        fail.append("%d span(s) below %.1f pt: %s" % (len(small), fmin,
                  " | ".join("%r@%.2f" % (s["text"].strip()[:16], s["size"]) for s in small[:5])))
    if badfont:
        fail.append("non-Arial face(s): %s" % ", ".join(badfont))
    if w_cm > wmax + 1e-6:
        fail.append("page %.2f cm wider than %.2f cm" % (w_cm, wmax))
    for f in fail:
        print("   FAIL  " + f)
    return not fail


if __name__ == "__main__":
    a, opts, args = sys.argv[1:], {}, []
    i = 0
    while i < len(a):                       # accepts both --opt=value and --opt value
        if a[i].startswith("--"):
            k, v = (a[i][2:].split("=", 1) + [""])[:2] if "=" in a[i] else (a[i][2:], a[i + 1])
            opts[k] = v
            i += 2 if "=" not in a[i] else 1
        else:
            args.append(a[i]); i += 1
    fmin = float(opts.get("fmin", 5.5)); wmax = float(opts.get("wmax", 18.0))
    if not args:
        print(__doc__)
        sys.exit(2)
    ok = all(gate(p, fmin, wmax) for p in args)
    print("font_gate:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
