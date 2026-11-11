# -*- coding: utf-8 -*-
"""Shared settings for the S6/S7 supplement figures (2026-11-03, this session = the parallel-figure side).

Discipline (per the FIG6_V11 review receipt §5 / ROUTEA_D receipt):
  - every statistic on a figure is READ ONLY from the delivered tables, never recomputed in a drawing script; derived quantities are first written as
    SUPPORT_v1/tables/*.csv (with a source column), then read back by the panel.
  - all Arial, visible glyphs >= 5.5 pt, no Type3 (pdf.fonttype=42), open spines.
  - each panel is saved as native PDF; composition uses fitz show_pdf_page (vector, never rasterized).
"""
import os, sys, textwrap

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
sys.path.insert(0, WORK + "/scripts/MS_GWAS/fig6_repro/deps")
# Third-party deps (PyMuPDF et al.) come from the active environment; DSCI_PYLIBS is the
# documented escape hatch for a wheel-only prefix.  No personal site-packages paths.
if os.environ.get("DSCI_PYLIBS"): sys.path.insert(0, os.environ["DSCI_PYLIBS"])
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from arial_style import use_arial

ARIAL = WORK + "/fonts/ARIAL.TTF"
mpl.font_manager.fontManager.addfont(ARIAL)
mpl.rcParams["font.family"] = mpl.font_manager.FontProperties(fname=ARIAL).get_name()
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
use_arial({"xtick.labelsize": 5.6, "ytick.labelsize": 5.6})
mpl.rcParams.update({"font.size": 6.0, "axes.linewidth": 0.55,
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                     "xtick.major.size": 2.2, "ytick.major.size": 2.2,
                     "axes.spines.top": False, "axes.spines.right": False})

R = WORK + "/figures/routeA_out"
SUP = R + "/SUPPORT_v1"
PDFDIR = SUP + "/pdf"
TABDIR = SUP + "/tables"
LOGDIR = SUP + "/logs"
for _d in (PDFDIR, TABDIR, LOGDIR):
    os.makedirs(_d, exist_ok=True)

# manuscript-wide canonical colours (fig6_repro/make_fig6_panels.py 2026-10-28 note; red = disease/condition)
RED, BLU, ORG, BLU2, GRY, MUT, GRID = "#D51F26", "#272E6A", "#FF7F0E", "#1F77B4", "0.42", "0.30", "0.86"
OLIGO, OPC = "#1F77B4", "#FF7F0E"          # shared with Fig 3/4/S4/GA: Oligo blue / OPC orange
WB = dict(fc="white", ec="none", alpha=0.92, pad=0.3, zorder=8)


def save(fig, name, note="", tight=False):
    """Saves a PDF (vector). Default is a FIXED canvas (no tight bbox), so overflowing text exposes itself instead of stretching the canvas."""
    out = os.path.join(PDFDIR, name + ".pdf")
    if tight:
        fig.savefig(out, bbox_inches="tight", pad_inches=0.03)
    else:
        fig.savefig(out)
    d = plt.get_current_fig_manager()
    w, h = fig.get_size_inches()
    plt.close(fig)
    import fitz
    doc = fitz.open(out); r = doc[0].rect; doc.close()
    if os.environ.get("PNGOUT"):
        import fitz as _f
        _dd = _f.open(out)
        _dd[0].get_pixmap(matrix=_f.Matrix(300 / 72, 300 / 72)).save(
            os.path.join(os.environ["PNGOUT"], name + ".png"))
        _dd.close()
    PT2CM = 2.54 / 72.0
    print("[panel] %-34s %.2f x %.2f cm  %s" % (name, r.width * PT2CM, r.height * PT2CM, note), flush=True)
    return name, r.width * PT2CM, r.height * PT2CM


CM = 2.54


def fmt_p(p):
    """Canonical P printing (three decimals at >= 1e-3; a.bde-x below)."""
    if pd.isna(p):
        return "NA"
    if p >= 0.001:
        return "%.3f" % p
    for dec in (1, 2, 3):
        e = int(np.floor(np.log10(p)))
        m = p / 10 ** e
        if round(m, dec) < 10:
            return ("%%.%dfe%%d" % dec) % (m, e)
    return "%.3f" % p


def read(src_rel, **kw):
    """Reads a delivered table from a routeA_out-relative path."""
    return pd.read_csv(os.path.join(R, src_rel), **kw)


def panel(fig, gs, row, col, title, sub=None):
    ax = fig.add_subplot(gs[row, col])
    ax.set_title(title, fontsize=6.4, loc="left", pad=(12 if sub else 3))
    if sub:
        ax.text(0.0, 1.005, sub, transform=ax.transAxes, fontsize=5.5, color=MUT, va="bottom")
    return ax


def note(ax, x, y, s, fs=5.5, ha="left", va="top", col="0.25"):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=fs, ha=ha, va=va, color=col, bbox=WB, zorder=9)


def wrapblock(ax, lines, fs=5.8, width=74, x=0.0, y=1.0, ls=1.5):
    """ls is tunable: once a text block's line count grows it falls off the figure (measured: S7c's last paragraph vanished, absent from the PDF),
    and 1.25 line spacing stays legible at 5.5 pt - cheaper than deleting content."""
    body = "\n".join(("\n".join(textwrap.wrap(b, width)) if b.strip() else "") for b in lines)
    ax.text(x, y, body, transform=ax.transAxes, va="top", ha="left", fontsize=fs,
            linespacing=ls, color="0.15")


def headroom(ax, top=0.0, bot=0.0):
    """Reserves a BLANK BAND for annotations: they may never sit on data points (receipt S6F7_PANEL_BEAUTY §6, the DOT-UNDER-TEXT gate).
    note() already carries a white box; dropping it on a point erases the point - worse than a cosmetic blemish."""
    lo, hi = ax.get_ylim(); sp = hi - lo
    ax.set_ylim(lo - bot * sp, hi + top * sp)


def gridspec(fig, nrow, ncol, w_cm=18.0, h_cm=17.4):
    fig.set_size_inches(w_cm / CM, h_cm / CM)
    return fig.add_gridspec(nrow, ncol, left=0.135, right=0.985, top=0.945, bottom=0.055,
                            wspace=0.36, hspace=0.52)


def _text_w_pt(s, fs, weight="normal", rend=None):
    """Measures a text line's width (pixels). Table column widths used to be guessed in characters; tall thin columns heavy in capitals/parentheses overran,
    and whole tables collided (measured: 401 in S7a) - so the renderer now measures."""
    global _REND
    fig = plt.gcf()
    if rend is None:
        if _REND is None:
            fig.canvas.draw()
            _REND = fig.canvas.get_renderer()
        rend = _REND
    t = fig.text(0.0, 0.0, s, fontsize=fs, weight=weight)
    w = t.get_window_extent(rend).width
    t.remove()
    return w


_REND = None


def draw_table(ax, rows, cols, fs=5.5, lh=0.030, pad=0.014, top=1.0, auto=True, gap_pt=6.0):
    """Wraps by column width (characters) and lays out by the actual line count, eliminating table-text overlap.
    cols = [(header, value key, width in chars, x-start fraction)], rows = list[dict]. Returns the total height used (axes fraction).
    with auto=True the x starts are not hand-written fractions but accumulated per column from the renderer-measured widest line + gap_pt,
    (the previous hand-written fractions disagreed with the real glyph widths, 401 collisions across the table)."""
    import math
    ax.axis("off")
    xs = [c[3] for c in cols]
    if auto:
        # first compute each row's wrapping (same rule as the drawing below), then measure; iteratively narrow wch if it does not fit
        ax.figure.canvas.draw()                      # take this figure's own renderer, not the previous figure's
        rend = ax.figure.canvas.get_renderer()
        axw = ax.get_window_extent(rend).width
        gap = gap_pt * rend.dpi / 72.0          # get_window_extent returns pixels, not points
        wchs = [c[2] for c in cols]
        cells_all, widths, need, _it = [], [], 0.0, 0
        while True:
            _it += 1
            cells_all, widths = [], []
            for ck, (title, key, _w, _x) in enumerate(cols):
                mx = _text_w_pt(title, fs, "bold", rend)
                col_lines = []
                for r in rows:
                    lines = []
                    for para in str(r.get(key, "")).split("\n"):
                        lines += (textwrap.wrap(para, wchs[ck]) or [""]) if para.strip() else [""]
                    col_lines.append(lines)
                    for ln in lines:
                        mx = max(mx, _text_w_pt(ln, fs, "normal", rend))
                cells_all.append(col_lines)
                widths.append(mx)
            need = sum(widths) + gap * (len(cols) - 1)
            if need <= axw or _it >= 6:
                break
            wchs = [max(6, int(w * axw / need * 0.97)) for w in wchs]
        if _it > 1:
            print("  [table auto-shrink] after round %d: %.1f/%.1f px" % (_it, need, axw), flush=True)
        x = 0.0
        xs = []
        for w in widths:
            xs.append(x / axw)
            x += w + gap
    y = top
    for k, (title, key, wch, _x) in enumerate(cols):
        ax.text(xs[k], y, title, fontsize=fs, weight="bold", color="0.25", va="top", ha="left",
                transform=ax.transAxes)
    y -= lh + pad
    for ri, r in enumerate(rows):
        nlines = 1
        cells = []
        for k, (title, key, wch, _x) in enumerate(cols):
            if auto:
                lines = cells_all[k][ri]
            else:
                txt = str(r.get(key, ""))
                lines = []
                for para in txt.split("\n"):
                    lines += (textwrap.wrap(para, wch) or [""]) if para.strip() else [""]
            nlines = max(nlines, len(lines))
            cells.append(lines)
        h = nlines * lh
        for k, lines in enumerate(cells):
            for i, ln in enumerate(lines):
                ax.text(xs[k], y - i * lh, ln, fontsize=fs, color="0.12", va="top", ha="left",
                        transform=ax.transAxes)
        y -= h + pad * 2
        ax.plot([0, 1], [y + pad, y + pad], lw=0.35, color="0.86", transform=ax.transAxes, zorder=0)
    return top - y


def say(*a):
    print("[%s] %s" % (pd.Timestamp.now().strftime("%H:%M:%S"), " ".join(str(x) for x in a)), flush=True)
