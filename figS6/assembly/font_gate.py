#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""font_gate.py - turns house spec §9's "two font sizes only" into a runnable generic hard gate.

Background: the two-size policy (8.0 = panel letters + row-head titles; 6.0 = everything else) previously lived only in figS7's verifier,
so figS6/fig6/fig7 slipped through in droves. This script audits ANY PDF at word level,
serving both as a pre-delivery gate and as an inventory tool for "which figures are still non-compliant".

Why not simply change arial_style's BASE: historical scripts must be able to re-render the delivered layouts (§14.13),
and moving the baseline defaults would make every old figure unrereaderable. Compliance rides on new figures explicitly opting into two tiers (use_two_tier()),
while old figures first get an inventory here and pass gate-by-gate when revised.

Usage:
  python3 font_gate.py                      # inventory every delivered PDF of this project
  python3 font_gate.py a.pdf b.pdf --strict # audit only these two; exit 1 on any violation
Exit code: 0 = compliant; 1 = violations found (with --strict).
"""
import os, sys, glob, fitz

TWO_TIER = (8.0, 6.0)      # the only allowed two sizes
TOL = 0.05                 # measured PDF sizes jitter by 0.01-0.03
MIN_OK = 5.5               # the §1 floor (should not appear under the two-size policy; flagged here as a second alarm)
CM = 28.3465

# §1 font whitelist (journal preflight reads the PDF's /BaseFont, so "looks like Arial" does not pass)
WHITELIST = {"ArialMT", "Arial-BoldMT", "Arial-ItalicMT", "Arial-BoldItalicMT"}

# delivered-file locations (edit as needed; glob only, no hard-coded single files)
DEFAULT_PATTERNS = [
    os.path.join(os.environ.get('DSCI_REPO_ROOT', '.'), 'fig6/fig6_out/*.pdf'),
    os.path.join(os.environ.get('DSCI_REPO_ROOT', '.'), 'figS6/fig6_out/*.pdf'),
    os.path.join(os.environ.get('DSCI_REPO_ROOT', '.'), 'fig7/fig7_out/*.pdf'),
    os.path.join(os.environ.get('DSCI_REPO_ROOT', '.'), 'figS7/fig7_out/*.pdf'),
    os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'routeA_out/SUPPORT_v1/pdf/*.pdf'),
    os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'routeA_out/SUPPORT_v2/pdf/*.pdf'),
    os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'routeA_out/SUPPORT_v2/pdf_v12/*.pdf'),
    os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'routeA_out/FIG*_*/pdf/*.pdf'),
]


def audit(path, verbose=True):
    """Returns a dict: size tier -> {char count, compliant?, examples}; plus the page size."""
    try:
        d = fitz.open(path)
    except Exception as e:
        return dict(path=path, error=str(e))
    hist, ex, sizes = {}, {}, []
    for pi, pg in enumerate(d):
        for b in pg.get_text("dict")["blocks"]:
            if b["type"]:
                continue
            for l in b["lines"]:
                for s in l["spans"]:
                    t = s["text"].strip()
                    if not t:
                        continue
                    sz = round(s["size"], 2)
                    sizes.append(sz)
                    k = sz
                    hist[k] = hist.get(k, 0) + len(t)
                    if k not in ex and not any(abs(k - a) < TOL for a in TWO_TIER):
                        ex[k] = (pi + 1, t[:34])
    ok = all(any(abs(k - a) < TOL for a in TWO_TIER) for k in hist)
    pgs = d[0]
    # ---- §1 gate: font names + Type3
    fonts, type3 = set(), 0
    for pg in d:
        for fo in pg.get_fonts(full=True):
            fonts.add(str(fo[3]))
    for x in range(1, d.xref_length()):
        try:
            o = d.xref_object(x, compressed=True) or ""
        except Exception:
            continue
        if "/Font" in o and "/Type3" in o:
            type3 += 1
    basefont_off = sorted(n for n in fonts if n.split("+")[-1] not in WHITELIST)
    out = dict(path=path, npages=len(d), fonts=sorted(fonts), basefont_off=basefont_off, type3=type3,
               page_cm=(round(pgs.rect.width / CM, 2), round(pgs.rect.height / CM, 2)),
               tiers=sorted(hist), bad_chars=sum(v for k, v in hist.items()
                                                 if not any(abs(k - a) < TOL for a in TWO_TIER)),
               total_chars=sum(hist.values()), min_size=min(sizes) if sizes else None,
               compliant=ok and not basefont_off and not type3,
               size_ok=ok, font_ok=not basefont_off and not type3,
               examples=ex, hist=hist)
    d.close()
    return out


def kind(path):
    """Figures vs tables. Tables follow the reportlab three-line book-table style (its own tiers: titles 8.4 / body 6.2),
    the two-size policy applies to FIGURES only; auditing them together would let Table_S6's 38k characters drown the real signal."""
    b = os.path.basename(path)
    return "table" if b.startswith(("Table_", "table_")) else "figure"


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    strict = "--strict" in sys.argv
    files = args or sorted({f for pat in DEFAULT_PATTERNS for f in glob.glob(pat) if os.path.isfile(f)})
    rows = [audit(f) for f in files]
    bad = [r for r in rows if not r.get("error") and not r["compliant"]]
    print("fitz %s | auditing %d PDFs | two tiers = %s (tolerance %.2f)" % (fitz.VersionBind, len(rows), TWO_TIER, TOL))
    print("%-58s %6s %10s %9s  %s" % ("file", "pages", "bad chars", "total", "bad tier (examples)"))
    for r in sorted(rows, key=lambda x: -x.get("bad_chars", 10**9)):
        if r.get("error"):
            print("%-58s open failed: %s" % (os.path.basename(r["path"])[:58], r["error"][:60])); continue
        tag = ""
        if r.get("basefont_off"):
            tag += " basefont violations=%s" % (r["basefont_off"][:2],)
        if r.get("type3"):
            tag += " Type3=%d" % r["type3"]
        ev = "; ".join("%s→p%s '%s'" % (k, v[0], v[1]) for k, v in sorted(r["examples"].items())[:3])
        print("%-58s %6d %10d %9d  %s%s" % (os.path.basename(r["path"])[:58], r["npages"],
                                            r["bad_chars"], r["total_chars"], (ev[:60] or "—"), tag))
    fig = [r for r in rows if kind(r["path"]) == "figure" and not r.get("error")]
    tab = [r for r in rows if kind(r["path"]) == "table" and not r.get("error")]
    fbad = [r for r in fig if not r.get("size_ok", True)]
    nfbad = [r for r in fig if not r.get("font_ok", True)]
    tbf = [r for r in tab if not r.get("font_ok", True)]
    print("\n[figures] %d files | size violations %d files (%d chars) | basefont/Type3 violations %d files"
          % (len(fig), len(fbad), sum(r["bad_chars"] for r in fbad), len(nfbad)))
    if nfbad:
        print("  basefont gate failed (§1 whitelist; journal preflight hard fail):")
        for r in sorted(nfbad, key=lambda x: os.path.basename(x["path"]))[:14]:
            print("    %-46s off-whitelist=%s Type3=%d" % (os.path.basename(r["path"])[:46], r["basefont_off"][:3], r["type3"]))
    if tbf:
        print("  table basefont gate failed: %s" % [os.path.basename(r["path"]) for r in tbf][:6])
    print("[tables] %d files (three-line table style, not subject to the two-tier gate), chars %d" % (len(tab), sum(r["total_chars"] for r in tab)))
    # floor gate: no PDF may contain text below 5.5 pt (§1)
    tiny = [(os.path.basename(r["path"]), r["min_size"]) for r in rows if r.get("min_size") and r["min_size"] < MIN_OK - TOL]
    if tiny:
        print("\n!! files below the 5.5 pt floor %d: %s" % (len(tiny), "; ".join("%s(%.2f)" % t for t in tiny[:10])))
    if not fbad and not nfbad:
        print("PASS: all audited figures use the two-tier policy (8.0 / 6.0) and whitelisted basefonts")
    else:
        print("figures to fix (by bad-char count, top 15):")
        for r in sorted(fbad, key=lambda x: -x["bad_chars"])[:15]:
            print("  %-52s bad %5d / %5d chars  tier=%s" % (os.path.basename(r["path"])[:52], r["bad_chars"],
                                                           r["total_chars"], r["tiers"]))
    if strict and (fbad or nfbad):
        sys.exit(1)


if __name__ == "__main__":
    main()
