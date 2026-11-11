#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build ``Figure_6.ipynb`` — Figure 6, v36 final (single row A | B | C, 2026-11-17).

The delivered figure is one row: A the scDRS MS score across the six oligodendroglial subtypes
(box plot + the one-line MWU *P* table), B the same score on the UMAP of the same nuclei
(aspect-locked hexbin, horizontal colour bar), C the 74-trait cross-channel (expression vs
chromatin) Δ1 ladder with the MS / SCZ / HEIGHT labels. The gene-level landing panel that used
to be main-figure C was pulled in v36 (no defensible positive); the row is composed at exactly
17.00 cm with two-tier type (8.0/6.0) and zero vector scaling.

How the notebook is built: every code cell is the *verbatim text* of the delivered pipeline
script (`assembly/fig6v36_abd_final.py`), preceded by a banner that injects ``__file__`` so the
script's own HERE/WORK resolution behaves exactly as when it runs as a file. Retyping the panel
code would create a second implementation that can drift from the delivered page, so the
notebook *is* the delivered pipeline. The generator asserts the embedded file exists.

The notebook ends at the three panels (`fig6_out/fig6v36_{A_score,B_umap,C_ladder}.pdf`); the
row composition is the delivered driver's job (`FIG6_PANELS_ONLY=1` stops the embedded script
before the fitz sheet assembly). Acceptance gates and audit scripts are meta, so they live in
`assembly/` and are not part of the notebook.

Run:  cd fig6 && python3 build_fig6_notebook.py
Then: cd fig6 && jupyter nbconvert --to notebook --execute --inplace Figure_6.ipynb
"""
import os
import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = HERE
ASSY = os.path.join(REPO, "assembly")

FINAL = os.path.join(ASSY, "fig6v36_abd_final.py")   # A | B | C + compose, the delivered script

SHEET = "17.00 x 4.89 cm"
SIZES = {"A": "5.80 x 4.89", "B": "5.75 x 4.89", "C": "4.85 x 4.89"}

cells = []


def md(s):
    cells.append(nbf.v4.new_markdown_cell(s.strip()))


def code(s):
    cells.append(nbf.v4.new_code_cell(s.strip("\n")))


def src(path):
    assert os.path.exists(path), "missing delivered script: %s" % path
    return open(path, encoding="utf-8").read()


def banner(relpath):
    return ('# ------------------------------------------------------------ %s (verbatim)\n'
            '# The cell body below IS the delivered Figure 6 script. `__file__` is injected so\n'
            "# the script's own HERE/WORK resolution is identical to running it as a file.\n"
            '__file__ = os.path.join(REPO, "%s")\n' % (relpath, relpath)).replace(
                "The cell body below IS", "# The cell body below IS").replace(
                "the script's own", "# the script's own")


md(r"""# Figure 6 — Cross-channel concordance of the Oligo–OPC contrast (v36 FINAL, single row)

**What the figure claims.**  The MS GWAS scDRS score separates mature oligodendrocytes from
OPC/COP-NFOL precursors in the own dsciCAT-seq data (A, B). The same *contrast* — not the
absolute score — can be measured independently in the expression channel and in the chromatin
(gene-score) channel of the *same nuclei*: panel C plots the 74-trait per-trait Δ1 pair with the
three traits that matter (MS, SCZ, HEIGHT) named. Every statistic on the face is a Spearman
correlation over the 74 traits or a number read from a delivered table; nothing is recomputed by
the drawing code.

| Panel | File | Size (cm) | Content |
|---|---|---|---|
| A | `fig6v36_A_score.pdf` | %(A)s | scDRS MS `norm_score` per nucleus across the six Fig 2 subtypes; boxes = median/IQR, whiskers 5–95 %%, jittered dots = nuclei, MWU *P* vs COP/NFOL above |
| B | `fig6v36_B_umap.pdf` | %(B)s | the same score on the UMAP of the same nuclei (hexbin of the binned mean, aspect locked to the data's own 1.63 ratio) |
| C | `fig6v36_C_ladder.pdf` | %(C)s | 74 GWAS traits: Δ1 in the chromatin channel (x) vs the expression channel (y), own same-nucleus data; MS / SCZ / HEIGHT named with their ranks |

Delivered sheet: `fig6_out/Figure_6_v36_ABC_row.pdf`, **%(SHEET)s** — three panels, one row,
letters stamped by the compose layer, two-tier type 8.0/6.0, zero vector scaling, page-width
gate 17.00 cm. The row is composed by the delivered driver outside this notebook
(`assembly/fig6v36_abd_final.py` run by `run_all.sh`): per the author instruction the
**composition code is not part of the notebook** — the notebook reproduces the three panels
at the delivered sizes, and the driver keeps producing the sheet.
""" % dict(SIZES, SHEET=SHEET))

md(r"""---

## 0. Environment, paths and versions

The notebook reads **delivered tables and one small cache**; the analysis-tree objects are needed
only to *rebuild* the cache.

| Variable | Default (in-repo) | Points at |
|---|---|---|
| `DSCI_WORK_ROOT` | three levels above `assembly/` | the analysis tree; needed only when the cache is missing |
| `FIG6_OUT` | `fig6_out/` | delivered tables (`fig6A_stats.csv`, `fig6C_*`) and everything this notebook writes |
| `FIG6_PANELS` | `fig6_out/` | where the panel PDFs and the composed row are written |
| `FIG6_CACHE` | `fig6_out/caches/` | the per-nucleus score table (`fig6_cellscore.csv`) |

With the cache present the notebook runs without the sequencing tree; `upstream/` produced every
input. The script itself validates A's statistics against the delivered `fig6A_stats.csv` and the
ladder's ρ against the delivered `fig6C_stats.csv`.""")

md("""### 0.1 Where the inputs come from

| File | Produced by | Used by |
|---|---|---|
| `caches/fig6_cellscore.csv` | this script's first run (from the h5ad + MS score) | A, B |
| `oligo/3k_dsci_oligo_subset.h5ad` (5,977 nuclei; `obs.cell_type` = Fig 2 subtype, `obs.UMAP_1/2`) | **Fig 2 lineage curation: converted from the author's own Seurat RDS** | A, B (only when the cache is missing) |
| `oligo/scdrs_ms_out/PASS_Multiple_sclerosis.score.gz` + `.scdrs_group.cell_type` | `scdrs compute-score` (see `upstream/README.md`) | A (score + group P; score only when the cache is missing) |
| `fig6A_stats.csv` | the delivered v21 chain (same-nucleus statistics) | A (the assert that the face matches the ledger) |
| `fig6C_delta1_own74.csv`, `fig6C_stats.csv` | the delivered v21 chain (same-nucleus 74-trait scores, both channels) | C |

**RDS → h5ad chain (the two own-data files).**  `upstream/README.md` Step 0 documents both
conversions, each re-run from the RDS and verified bit-for-bit against the delivered file:
`FY/oligo/3k_dsci_oligo_subset.RDS` (794 MB, the Fig 2 lineage Seurat object) →
`3k_dsci_oligo_subset.h5ad` via `rds_to_dump.R` + `build_h5ad.py` (5,977 nuclei × 31,263 genes;
`cell_type` = the six subtypes, UMAP inherited, never recomputed); and
`FY/3k_dsci_ht_10x_gs_subset.RDS` (13.2 GB, Signac multiome) → `genescore_oligo_5977.h5ad` via
`export_genescore.R` + `build_genescore_h5ad.py` (same nuclei in the same order × 24,333 genes).
The MS scDRS score comes from `scdrs compute-score` on the first h5ad (`upstream/README.md`,
step 1); no third-party data enters Figure 6 — the three library chemistries of the deposit
(CNGBdb **CNP0010263**) appear in the supplementary S6m, not here.""")

md("""### 0.2 Software versions""")

code("""import os, sys, time
import numpy as np, pandas as pd
import matplotlib, scipy, fitz

REPO = os.path.abspath(os.getcwd())              # run the notebook from fig6/
ASSY = os.path.join(REPO, "assembly")
OUTD = os.path.join(REPO, "fig6_out")
CACH = os.path.join(OUTD, "caches")
for _d in (OUTD, CACH):
    os.makedirs(_d, exist_ok=True)

os.environ.setdefault("DSCI_WORK_ROOT", os.path.abspath(os.path.join(REPO, "..", "..")))
os.environ["FIG6_OUT"], os.environ["FIG6_PANELS"] = OUTD, OUTD
os.environ["FIG6_CACHE"] = CACH
os.environ["FIG6_PANELS_ONLY"] = "1"      # stop at the three panels; the row stays the driver's product
sys.path.insert(0, ASSY)

print("python   ", sys.version.split()[0])
for m in (np, pd, scipy, matplotlib, fitz):
    print("%-10s %s" % (m.__name__, m.__version__))
# Hard kernel gate (2026-11-18). The delivered panel geometry is a matplotlib OUTCOME, not a request:
# the tight bounding box of the colour-bar tick row moves with the renderer, and under matplotlib 3.6.2
# (python 3.8) the three panels land 0.05 cm SHORT of the delivered sizes - identical words, identical
# drawing counts, wrong page. The delivered v36 row was rendered with matplotlib 3.9.x, so the run
# refuses to continue on anything older instead of quietly emitting a figure that then fails the row gate.
_mv = tuple(int(x) for x in matplotlib.__version__.split("+")[0].split(".")[:2])
assert _mv >= (3, 9), (
    "matplotlib %s is too old to reproduce the delivered geometry (panels land ~0.05 cm short); run this "
    "notebook in the python 3.9 kernel (matplotlib >= 3.9), see README 'How to run'" % matplotlib.__version__)
print("kernel gate OK: matplotlib %s reproduces the delivered sizes" % matplotlib.__version__)
print("PyMuPDF is load-bearing: every size claim below is measured from the produced PDF,")
print("never taken from the nominal canvas (tight bounding box != figsize).")""")

md("""### 0.3 Frozen criteria (fixed before the run, not chosen afterwards)

| | |
|---|---|
| contrast | Δ1 = `median(norm_score \| MOL6) − median(norm_score \| COP/NFOL)` per trait, per channel; labels are the Fig 2 lineage calls of the own h5ad (the two nucleus pools of A/B and C are the same object here) |
| main statistic | Spearman ρ over the 74 traits between the two channels' Δ1 — **0.819** as delivered, asserted against `fig6C_stats.csv` in the cell below |
| named traits | MS, SCZ, HEIGHT (the same three as the old C panel, same colours from `fig6_palette`); labels placed by the delivered leader-search so no label sits on the y = x line or on a point |
| typography | two sizes only, 8.0 pt (letters, axis titles, row heads) and 6.0 pt (ticks, annotations, statistics); Arial; face text black except the three named-trait labels |
| geometry | row height solved from B's aspect lock (B rendered first, its ink height measured, A/C re-rendered to match); three panels then spread to exactly **17.00 cm** by widening the gaps — no panel is ever scaled |
| jitter / RNG | panel A jitter: `RandomState(7)`, Gaussian σ = 0.07; hexbin `gridsize=25`, `mincnt=5`, mean reduce — the same parameters as the delivered page |""")

md("""---

## 1. The figure — one cell, the delivered script

Everything above set the environment. The single cell below is the delivered
`assembly/fig6v36_abd_final.py` verbatim: it loads the same data the delivered page was made
from, renders A, B and C, asserts A's statistics and the ladder's ρ against the delivered
tables, and solves the row height from B's aspect lock so each panel lands at the delivered
size. The notebook exports `FIG6_PANELS_ONLY=1`, so the cell ends at the three panel PDFs —
the row composition, the letters and the 17.00 cm page-width gate stay in the delivered driver
(per the author instruction the composition code is not part of the notebook).""")

code(banner("assembly/fig6v36_abd_final.py") + src(FINAL))

md("""---

## 2. What this notebook produces

Three panel PDFs, written to `fig6_out/` by the cell above. Sizes are measured from the produced PDFs
(the row height is solved from B's aspect lock, so these are outcomes, not requests):

| output | measured size (cm) | content |
|---|---|---|
| `fig6_out/fig6v36_A_score.pdf` | 5.82 x 4.54 | scDRS MS `norm_score` per nucleus across the six Fig 2 subtypes; boxes median/IQR, whiskers 5-95 %, jittered dots, MWU *P* vs COP/NFOL |
| `fig6_out/fig6v36_B_umap.pdf` | 5.76 x 4.53 | the same score on the UMAP of the same nuclei (hexbin of binned mean, aspect locked to the data's 1.63 ratio) |
| `fig6_out/fig6v36_C_ladder.pdf` | 4.87 x 4.55 | cross-chemistry ladder: per-trait rho between channels, size-adjusted, with the MS / SCZ / HEIGHT points named |

**The kernel is part of the specification.** Those sizes are outcomes of matplotlib's tight-bbox. Under
matplotlib 3.6.2 (python 3.8) all three panels land 0.05 cm short with otherwise identical ink, so the
version cell hard-asserts matplotlib >= 3.9 and the kernelspec names the python 3.9 kernel: running this
notebook under the base python 3.8 kernel is expected to stop on that assert, not to produce a page that
then fails the 17.00 cm row gate.

Composing the row, drawing the panel letters and enforcing the 17.00 cm page-width gate stay in the
delivered driver (`assembly/fig6v36_abd_final.py` run by `run_all.sh`), because composition is journal
production and not content (house spec SS1, author instruction). The delivered page is

```
figures/routeA_out/P0_8C/panels/Figure_6_v36_ABC_row.pdf      17.00 x 4.89 cm
```

**Delivery condition.** The cell above asserts A's per-subtype statistics against
`tables_v36/fig6A_stats.csv` and the ladder's rho against `tables_v36/fig6C_stats.csv` before it
exports, and it refuses to write a row whose panels do not land at the measured sizes above. A notebook
that runs with any assertion failing is not the delivered figure. The gene-level landing panel that used
to be main-figure C was pulled in v36 (no defensible positive); its scripts and tables remain in
`upstream/` and `tables_v36/` for provenance.""")

# =============================================================================
nb = nbf.v4.new_notebook()
nb.cells = cells
nb.metadata = {
    "kernelspec": {"display_name": "py39 (matplotlib >= 3.9)", "language": "python", "name": "py39"},
    "language_info": {"name": "python", "version": "3.9"},
    "dsci": {"figure": "6", "version": "v36", "sheet": SHEET,
             "generator": "build_fig6_notebook.py",
             "embedded_scripts": ["assembly/fig6v36_abd_final.py"]},
}
out = os.path.join(REPO, "Figure_6.ipynb")
with open(out, "w") as f:
    nbf.write(nb, f)
print("wrote %s (%d cells)" % (out, len(cells)))
