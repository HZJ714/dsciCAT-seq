# -*- coding: utf-8 -*-
"""figstyle.py — shared figure setup for the Figure S7 pipeline.

Three things are fixed here so that every panel script agrees:
  * fonts: Arial only, ``pdf.fonttype = 42`` (no Type 3 in the delivered PDF);
  * paths: three environment variables, all optional
      ``FIG7_ROOT``  directory that holds the delivered input tables,
                   laid out as ``<cohort>/tables/*.csv`` (default ``<this dir>/data``)
      ``FIG7_OUT``   where the PDF / PNG / number register are written
                   (default ``<this dir>/figS7_out``)
      ``FIG7_FONT``  path to an Arial TrueType file (Arial is **not**
                   redistributed by this repository; see README "Fonts")
  * one discipline: panels never recompute a statistic that a delivered table
    already contains — they read the table. Derived quantities are written to
    ``data/`` by the scripts in ``upstream/`` first, then read here.

Only the pieces used by ``Figure_S7.ipynb`` / ``verify_figS7.py`` are kept.
"""
import os
import sys

import matplotlib as mpl

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("FIG7_ROOT", os.path.join(HERE, "data"))
OUTDIR = os.environ.get("FIG7_OUT", os.path.join(HERE, "figS7_out"))
os.makedirs(OUTDIR, exist_ok=True)

# ---------------------------------------------------------------- fonts
FONT_CANDIDATES = [
    os.environ.get("FIG7_FONT", ""),
    os.path.join(HERE, "fonts", "ARIAL.TTF"),
    WORK + "/fonts/ARIAL.TTF",
    "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf",
]
FONT = next((f for f in FONT_CANDIDATES if f and os.path.exists(f)), None)
if FONT:
    mpl.font_manager.fontManager.addfont(FONT)
    FAMILY = mpl.font_manager.FontProperties(fname=FONT).get_name()
else:                                    # fall back to whatever the system calls Arial
    FAMILY = "Arial"
    sys.stderr.write("[figstyle] Arial TTF not found; relying on the system font "
                     "named 'Arial'. Set FIG7_FONT=/path/to/ARIAL.TTF.\n")
# House style baseline — the values the author's pipeline actually ran with
# (`arial_style.BASE`, then `support_common` overrides, then the per-figure rcParams in
# the drawing script). These are NOT cosmetic: `xtick/ytick.major.pad` and
# `xtick/ytick.major.size` decide where every tick label sits, and `savefig.dpi` decides
# the resolution at which the heatmap rasters are embedded. Reproducing the signed page
# point-for-point needs this block, not just the font name.
# The TrueType subsetter used by the PDF backend logs a benign
# "meta NOT subset; don't know how to subset; dropped" note once per embedded font per
# output file. Not a font problem: gate G2a still reads the embedded face from the PDF.
import logging
logging.getLogger("fontTools.subset").setLevel(logging.ERROR)

mpl.rcParams.update({
    "font.family": FAMILY,
    "mathtext.fontset": "custom",        # math glyphs follow the body font, never DejaVu
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.size": 6.0,                    # BASE 7.0, overridden to the house 6 pt body text
    "axes.titlesize": 7.0, "axes.labelsize": 7.0,
    "xtick.labelsize": 5.6, "ytick.labelsize": 5.6,
    "legend.fontsize": 6.0, "axes.titleweight": "normal",
    "axes.linewidth": 0.55,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2,
    "xtick.minor.size": 1.0, "ytick.minor.size": 1.0,
    "xtick.major.pad": 1.6, "ytick.major.pad": 1.6,
    "lines.linewidth": 0.8, "patch.linewidth": 0.5,
    "savefig.bbox": None, "savefig.dpi": 600, "figure.dpi": 120,
    "axes.spines.top": False, "axes.spines.right": False,
})
# house fallback rule: if the TTF is missing, at least keep the face name auditable
if FONT is None:
    mpl.rcParams["mathtext.fontset"] = "dejavusans"

CM = 2.54                                # cm per inch used by the canvas helpers
RED, BLU, ORG, BLU2, GRY, MUT, GRID = ("#D51F26", "#272E6A", "#FF7F0E", "#1F77B4",
                                       "0.42", "0.30", "0.86")
OLIGO, OPC = "#1F77B4", "#FF7F0E"        # lineage identity colours, shared with Figs 3/4


def read(rel, **kw):
    """Read a delivered table by its ``<cohort>/tables/...`` relative path."""
    return pd.read_csv(os.path.join(ROOT, rel), **kw)


def exists(rel):
    return os.path.exists(os.path.join(ROOT, rel))


def note(ax, x, y, s, fs=5.5, ha="left", va="top", col="0.25"):
    """Small grey annotation inside an axes (used for provenance notes)."""
    ax.text(x, y, s, transform=ax.transAxes, fontsize=fs, ha=ha, va=va, color=col,
            bbox=dict(fc="white", ec="none", alpha=0.92, pad=0.3, zorder=8), zorder=9)


def say(*a):
    print(" ".join(str(x) for x in a), flush=True)


def save(fig, name, note_txt="", tight=False):
    """Write ``<OUTDIR>/<name>.pdf`` on a fixed canvas (overflow stays visible),
    optionally a 300-dpi PNG when ``PNGOUT`` is set, and report the size in cm."""
    out = os.path.join(OUTDIR, name + ".pdf")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.03) if tight else fig.savefig(out)
    plt.close(fig)
    w_cm, h_cm = pdf_size_cm(out)                     # needs PyMuPDF, else MediaBox parse
    if os.environ.get("PNGOUT"):
        try:
            import fitz
            dd = fitz.open(out)
            dd[0].get_pixmap(matrix=fitz.Matrix(300 / 72, 300 / 72)).save(
                os.path.join(os.environ["PNGOUT"], name + ".png"))
            dd.close()
        except ImportError:
            say("[save] PyMuPDF missing - PNG preview skipped (pip install PyMuPDF)")
    say("[panel] %-34s %.2f x %.2f cm  %s" % (name, w_cm, h_cm, note_txt))
    return name, w_cm, h_cm


def pdf_size_cm(path):
    """Page size in cm. PyMuPDF when available, otherwise a MediaBox regex parse."""
    try:
        import fitz
        doc = fitz.open(path)
        r = doc[0].rect
        doc.close()
        return r.width * 2.54 / 72.0, r.height * 2.54 / 72.0
    except ImportError:
        import re as _re
        box = _re.search(rb"/MediaBox\s*\[\s*([\d.\s-]+)\]", open(path, "rb").read())
        v = [float(x) for x in box.group(1).split()]
        return (v[2] - v[0]) * 2.54 / 72.0, (v[3] - v[1]) * 2.54 / 72.0
