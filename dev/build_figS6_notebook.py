#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build ../figS6/Figure_S6.ipynb — Supplementary Figure S6, v16 (final 2026-11-18).

Layout of the delivered sheet (four rows, seven panels + one shared legend):

    row 1  A four MS loci (cis-window tracks)   |  B window permutation of the 331 platform genes
    row 2  C bin-size robustness of the Fig 6B RNA UMAP (17.00 cm wide)
    row 3  D MS top-1000 gene-level landing     |  E lead-SNP-window (cis38) landing  | legend
    row 4  F per-nucleus MS score, 2 channels x 3 chemistries  |  G 74-trait cross-channel pair

v16 replaced v11 (2026-11-06): panels D/E of the old notebook (subtype profiles, permutation
ladder) were dropped, the two gene-level "landing" heatmaps took their place, and the panel that
used to be Figure 6 D (the 74-trait expression/chromatin heatmap pair) moved here as panel G.
The superseded v11 notebook, its panels and its generator are archived in
`figS6/_superseded_20261118/`.

HOW THIS NOTEBOOK IS BUILT — and why it is not hand-written
-----------------------------------------------------------
Every code cell of the notebook is the *verbatim text* of the delivered pipeline script that
produces that panel (see `SRC` / `embed` below), preceded only by a two-line banner that injects
`__file__` so the script's own HERE/WORK resolution behaves exactly as when it is run as a file.

That is a deliberate change from the v11 generator, which retyped the panel code into cell
strings. Retyping created a second implementation: the notebook and the delivered page could
drift, and when they did, nobody knew which one was the figure. Here the notebook *is* the
delivered pipeline, so:

  * `python3 assembly/figS6v16_run_all.sh` and running this notebook give the same PDFs;
  * the generator asserts that every embedded file still exists and still carries the section
    banners it splits on, so a renamed or restructured script fails the build instead of
    silently generating a stale notebook.

The notebook ends at the panels and the legend: environment -> one section per panel (inputs,
core parameters, the delivered script). The four-row page composition is the delivered
driver's job (`figS6v16_run_all.sh`), not the notebook's. Acceptance gates and the
word-level comparison against a reference page are meta, so they live in `run_all.sh` and
dev, not in the notebook.

Run:  cd dev && python3 build_figS6_notebook.py        ->  ../figS6/Figure_S6.ipynb
Then: cd ../figS6 && jupyter nbconvert --to notebook --execute --inplace Figure_S6.ipynb
"""
import os
import re
import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(os.path.dirname(HERE), "figS6")
ASSY = os.path.join(REPO, "assembly")

PANELS = os.path.join(ASSY, "figS6v16_panels.py")        # A, B, C, F
LAND1000 = os.path.join(ASSY, "p0_8c_fig6c_top1000.py")  # D
LAND38 = os.path.join(ASSY, "p0_8c_fig6c_landing.py")    # E
LEGEND = os.path.join(ASSY, "figS6v18_legend_right.py")  # row-3 legend column
GHEAT = os.path.join(ASSY, "make_figS6v16_G_heatmap.py")  # G (ex Fig 6 D)

# ---------------------------------------------------------------- the numbers of the delivery
SHEET = "17.00 x 18.91 cm"
SIZES = {"A": "10.78 x 2.79", "B": "5.38 x 3.05", "C": "17.00 x 3.00",
         "D": "7.45 x 5.75", "E": "7.45 x 5.75", "LEG": "1.62 x 5.75",
         "F": "8.08 x 5.75", "G": "8.64 x 5.39"}
OFFSCALE = {"D_rna": 3, "D_atac": 11, "E_rna": 2, "E_atac": 1, "legend": 2}

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
            'The cell body below IS the delivered pipeline script. `__file__` is injected so\n'
            "the script's own HERE/WORK resolution is identical to running it as a file.\n"
            'Everything after the CLI entry point (`if __name__ == "__main__"`) is dropped:\n'
            'one kernel runs many cells, so a script that re-launches its sibling panels would\n'
            'redraw them under leaked rcParams. The panel functions are called explicitly.\n'
            '__file__ = os.path.join(REPO, "%s")\n' % (relpath, relpath)).replace(
                "The cell body below IS", "# The cell body below IS").replace(
                "one kernel runs", "# one kernel runs").replace(
                "the script's own", "# the script's own").replace(
                "Everything after", "# Everything after").replace(
                "redraw them under", "# redraw them under")


def embed(path, relpath):
    """Whole delivered script, verbatim, with the __file__ banner."""
    return banner(relpath) + src(path)


def split_panels():
    """figS6v16_panels.py is one file with five banner-delimited sections; hand them to the
    notebook as separate cells so each panel gets its own section (asserted, not assumed)."""
    t = src(PANELS)
    parts = re.split(r"^# =+ ", t, flags=re.M)
    assert len(parts) == 6, "figS6v16_panels.py lost its section banners (%d parts)" % len(parts)
    head, a, b, c, boxes, fcomb = parts
    pick = lambda s, n: "# ---- section %s of the delivered file ----\n" % n + s.split("\n", 1)[1]
    # the last chunk carries the script's CLI entry point; drop it (see the cell banner)
    fcomb = fcomb.split('if __name__ ==')[0].rstrip() + "\n"
    return (head, pick(a, "A"), pick(b, "B"), pick(c, "C"), pick(boxes, "per-channel boxes"),
            pick(fcomb, "F combined"))


HEAD, PA, PB, PC, PBOX, PF = split_panels()

# =============================================================================
md("""# Figure S6 — MS risk genes, from GWAS loci to the mature oligodendrocyte subtype (v16 FINAL)

**What the figure claims.**  MS GWAS signal reaches the oligodendrocyte lineage in a way that
depends on how the risk set is *defined at the gene level*.  Positional evidence (A: four loci
resolved to single genes by their lead SNP; B: how many platform genes fall in those windows by
chance) is followed by two gene-level landing heatmaps drawn on the same scale: the MS
MAGMA **top-1000** set restricted to genes measurable on the platform (D, 223 rows) and the
**cis-window** set of 38 genes (E).  D shows nothing that a reader can pick out — that is the
point — while E's set lands on MOL6 in both channels.  F and G then place the same contrast in
per-nucleus scores across three library chemistries and in the 74-trait cross-channel ranking.

| Panel | File | Size (cm) | Content |
|---|---|---|---|
| A | `figS6v16_A_loci.pdf` | %(A)s | four MS loci: gene TSS ticks within lead-SNP +/- 500 kb, hit gene in red, lead-SNP *P*, signed TSS distance |
| B | `figS6v16_B_perm.pdf` | %(B)s | window permutation of the hit count: 38/331 = 11.5 %%, *P* = 0.277 (10,000 draws) |
| C | `figS6v16_C_umapbins.pdf` | %(C)s | bin-size robustness of the Fig 6B RNA UMAP (single cells / 15 / 25 / 40), one shared colour scale, aspect locked |
| D | `fig6v32_C_top1000_landing_s6v18D.pdf` | %(D)s | **MS top-1000, 223 measurable genes** x 6 subtypes, RNA + ATAC, per-gene z; **deliberately unnamed** |
| E | `fig6v31_C_landing_s6v18E.pdf` | %(E)s | **cis-window 38 genes** x 6 subtypes, same scale; the four genes resolved in A are named |
| – | `figS6v18_DE_legend_v.pdf` | %(LEG)s | the one legend shared by D and E, as a vertical column at the right edge of row 3 |
| F | `figS6v16_F_cellboxes.pdf` | %(F)s | per-nucleus MS scDRS score, OPC vs Oligo, 2 channels x 3 chemistries, MWU *P* per panel |
| G | `fig6D_heatmap_pair_compact_s6v16G.pdf` | %(G)s | the 74-trait expression/chromatin pair (this was Figure 6 D before the move) |

Delivered sheet: `figS6_out/Figure_S6_v16.pdf`, **%(SHEET)s**, four rows (heights 3.05 / 3.00 / 5.75 /
5.75 cm), page-width gate 17.00 cm, page-height soft gate 18.00 cm (hard cap 20.50). The sheet
is composed by the delivered driver outside this notebook (`assembly/figS6v16_run_all.sh`):
per the author instruction the **composition code is not part of the notebook** — the notebook
reproduces the seven panels and the shared legend at the delivered sizes, and the driver
keeps producing the page.

**Display vs tested quantity — read this before the numbers.**  The heatmaps are coloured by a
*per-gene z across the six subtypes*; the statistic that is tested, and that orders the rows, is
`%(DELTA)s`.  The z is an affine, per-gene transform of Δ, so colour depth is comparable **within** a
row and never **between** rows.  The tested contrast therefore still has to be stated in the
legend even though its column was deleted from the face.
""" % dict(SIZES, SHEET=SHEET, DELTA="Δ = mean(MOL6) − mean(COP/NFOL)"))

md("""---

## 0. Environment, paths and versions

The notebook reads **only delivered tables and small caches**; nothing is re-scored, re-tested or
re-permuted by the drawing code (house rule: every statistic on the face is read from a CSV that
`upstream/` produced, so the figure cannot disagree with the ledger).

| Variable | Default (in-repo) | Points at |
|---|---|---|
| `DSCI_WORK_ROOT` | three levels above `assembly/` | the analysis tree; needed only to *rebuild* the caches |
| `S6_OUT` | `figS6_out/` | where the panel PDFs are written |
| `S6_V8` | `tables_v16/` | delivered statistic tables (window permutation, MWU ledger) |
| `S6_SCMS` | `tables_v16/` | the locus hit table (panel A) |
| `S6_CACHE` | `figS6_out/caches/` | column projections of the big objects (GTF, h5ad, 6 GB score sweeps) |
| `P08C_TABLES` | `tables_v16/` | the P0_8C gene-level tables behind D and E |

With the caches present the notebook runs on a laptop with no access to the sequencing tree.
`upstream/S6v16_export_caches.py` regenerates them from `DSCI_WORK_ROOT`; the caches are
projections only — no statistic is computed in that script.""" % ())

code("""import os, sys, glob, time
import numpy as np, pandas as pd

REPO = os.path.abspath(os.getcwd())            # run the notebook from figS6/
ASSY = os.path.join(REPO, "assembly")
TABD = os.path.join(REPO, "tables_v16")
OUTD = os.path.join(REPO, "figS6_out")
CACH = os.path.join(OUTD, "caches")
for _d in (OUTD, CACH):
    os.makedirs(_d, exist_ok=True)

os.environ.setdefault("DSCI_WORK_ROOT", os.path.abspath(os.path.join(REPO, "..", "..")))
os.environ["S6_OUT"] = OUTD
os.environ["S6_V8"], os.environ["S6_SCMS"] = TABD, TABD
os.environ["S6_CACHE"] = CACH
os.environ["P08C_TABLES"] = TABD
sys.path.insert(0, ASSY); sys.path.insert(0, os.path.join(REPO, "upstream"))

# --- module pinning (a shared-kernel hazard that per-script processes never show) ---------
# The repo carries its own arial_style / two_tier_check / support_common / fig6_palette; the
# analysis tree has same-named files too (<work>/scripts/arial_style.py and
# .../fig6_repro/deps/arial_style.py; neither has use_two_tier). The delivered scripts each
# run in their own process with a fixed resolution order; the notebook shares one kernel, and
# support_common inserts .../fig6_repro/deps at sys.path[0], so later cells would pick the
# wrong copy. Pinning the repo copies into sys.modules here makes every later
# `from arial_style import ...` hit the cache, immune to path order.
import importlib
for _n in ("arial_style", "two_tier_check", "support_common", "fig6_palette"):
    if os.path.exists(os.path.join(ASSY, _n + ".py")):
        sys.modules.pop(_n, None)
        importlib.import_module(_n)
print("pinned   :", ", ".join("%s=%s" % (k, os.path.basename(sys.modules[k].__file__))
                              for k in ("arial_style", "two_tier_check", "support_common")
                              if k in sys.modules))
print("repo     :", REPO)
print("analysis :", os.environ["DSCI_WORK_ROOT"], "(only used if a cache is missing)")
print("caches   :", ", ".join(sorted(os.listdir(CACH))) or "NONE — panels fall back to the analysis tree")""")

md("""### 0.1 Where the inputs come from

| File (in `tables_v16/` unless noted) | Produced by | Used by |
|---|---|---|
| `ms_locus_overlap_table.csv` | `upstream/ms_locus_overlap.py` (lead SNP +/- 500 kb against GENCODE v45) | A |
| `figS6v8_window_permutation{,_null}.csv` | `upstream/S6v8_window_permutation.py` (10,000 draws, seed fixed in the script) | B |
| `p08c_genelevel_all6_pool.csv`, `p08c_dual_celllevel223.npz` | `assembly/p0_8c_gene_enrichment.py` / `p0_8c_fig6c_full223.py` (scDRS gene-level scores per subtype, both channels) | D |
| `p08c_c1_rna_atac_by_subtype.csv`, `p08c_c2_celllevel38.npz` | `assembly/p0_8c_c1_cluster_probe.py` / `p0_8c_c2_stability.py` | E |
| `p08c_dual_setlevel_matched2d_mol6.csv` | `assembly/p0_8c_match2d.py` — **the single code path for every set-level P on D/E** (size-matched 2D permutation, 5x5 strata, 10,000 draws) | D, E, legend text |
| `figS6v10_G_opc_oligo_mwu.csv` | `upstream/S6v10_G_mwu_p.py` | F |
| `fig6D_delta1_74traits_per_tech.csv`, `fig6D_crosschem_ladder.csv` | `upstream/D6_crosschem_analysis.py` (same-nucleus scores of the three chemistries) | G |
| `caches/s6a_gene_coords.csv` | `upstream/S6v16_export_caches.py` | A |
| `caches/s6c_umap_score.csv` | same (UMAP_1/2 + cell_type + MS RNA score, 5,962 nuclei) | C |
| `caches/s6f_<job>.csv` | same (per-nucleus MS `norm_score` + label, 6 arms) | F |

The three chemistries are the deposited ones (CNGBdb **CNP0010263**): dsciCAT-seq (own),
10x Multiome, HT-scCAT-seq. The chromatin arms of F and G come from the **log1p-corrected**
re-run of 2026-11-04 (`upstream/README.md`, "Input transform"); the RNA arms are the delivered tree.

**Data origins — from which RDS/h5ad every number here comes.**  The own multiome traces to the
authors' RDS objects converted under `fig6/upstream/README.md` Step 0 (each verified bit-for-bit):
`FY/oligo/3k_dsci_oligo_subset.RDS` → `3k_dsci_oligo_subset.h5ad` (`rds_to_dump.R` + `build_h5ad.py`;
5,977 nuclei, the six Fig 2 subtypes, UMAP inherited — the panels' 5,962 = the paired subset that has
both-channel scores) and `FY/3k_dsci_ht_10x_gs_subset.RDS` → `genescore_oligo_5977.h5ad`
(`export_genescore.R` + `build_genescore_h5ad.py`; 24,333 ATAC gene scores, same nuclei, same order).
The per-nucleus MS `norm_score` trees behind C/F/G are `*.score.gz` / `PASS_*.full_score.gz` sweeps
run by `upstream/D5_run_scdrs.sh` (scdrs 1.0.2): the RNA arms are the delivered tree
(`D_crosschemistry/scdrs/`), the ATAC arms the **log1p-corrected** rerun (`scdrs_log1p_20261104/`).
The only external dataset is the mouse EAE **GSE250589** (panel F mouse rung, `M2b_formal_authoranno.py`).
The full production table is `figS6/upstream/README.md`.""")

md("""### 0.2 Software versions""")

code("""import matplotlib, scipy, fitz
print("python   ", sys.version.split()[0])
for m in (np, pd, scipy, matplotlib, fitz):
    print("%-10s %s" % (m.__name__, m.__version__))
print("PyMuPDF is load-bearing: every size claim below is measured from the produced PDF,")
print("never taken from the nominal canvas (tight bounding box != figsize).")""")

md("""### 0.3 Frozen criteria (fixed before the run, not chosen afterwards)

| | |
|---|---|
| tested contrast | `Δ = mean(MOL6) − mean(COP/NFOL)` per gene, per channel — the same poles as Fig 6 A/B, so the supplement and the main figure test one thing |
| display | per-gene z across the six subtypes: `z = (L − mean_s L) / sqrt(sd_s^2 + s0^2)`, `s0` = the row-scale floor (D: RNA 0.07755 / ATAC 0.02619; E: 0.1476 / 0.07267), clipped at ±2.0, one colour scale shared by RNA and ATAC and by D and E |
| detection floor | a cell is grey if the gene is detected in <10 %% of that subtype's nuclei: D 255/1338 (RNA) and 81/1338 (ATAC) cells; E 49/228 and 7/228 |
| off scale | a white triangle at the clipped edge of the cell: **D %d RNA + %d ATAC = %d cells, E %d + %d = %d cells**, all of them in the COP/NFOL reference column (see §4) |
| set-level test | size-matched 2D permutation, 5x5 strata (n_hit x mean detection), 10,000 draws, one seed, one code path (`p08c_dual_setlevel_matched2d_mol6.csv`): top-1000 P = 0.096 (RNA) / 0.745 (ATAC); cis38 P = 0.0017 / 0.0007; myelin anchor 0.0107 / 0.0080; 38 vs the 293 non-window platform genes 0.0127 / 0.0066 |
| row order | by \|Δ\| descending — display only, never a test |
| naming | D is **not** named (a null panel with named rows invites readers to find a pattern; and a name list produced by \|Δ\| plus a row-spacing rule changes when the page height changes). E names exactly the four genes that panel A resolves to single loci — positional evidence that is not on the face of E |
| jitter / RNG | panel F dots: `default_rng(20261117)`; panel C hexbins: deterministic |
| typography | two sizes only, 8.0 pt (letters, axis titles) and 6.0 pt (everything else); Arial; text on the face is black except subtype labels and arrows, which keep the identity colours |

The page geometry is **solved, not guessed**: `PW = 1.78 + 12·CW` and `2·PW + 2·0.24 + LW = 17.00`,
with `LW` measured from the legend's own ink — the same arithmetic as `assembly/figS6v16_run_all.sh`."""
    % (OFFSCALE["D_rna"], OFFSCALE["D_atac"], OFFSCALE["D_rna"] + OFFSCALE["D_atac"],
       OFFSCALE["E_rna"], OFFSCALE["E_atac"], OFFSCALE["E_rna"] + OFFSCALE["E_atac"]))

md("""---

## 1. Row 1 — panels A and B: positional evidence

**A — four loci.**  For each of PHLDB1 (chr11), PHGDH (chr1), BCAS1 (chr20) and PLLP (chr16) the
lead SNP is plotted with every GENCODE gene TSS inside ±500 kb; the hit gene is the red tick, and
the signed distance from lead SNP to gene TSS is printed in the corner. These four are the genes
that the window analysis resolves to a *single* locus, which is why they — and only they — are
named in panel E.

**B — is 38/331 more than chance?**  The number of platform genes inside the lead-SNP windows of
the MS hit list, permuted 10,000 times against size-matched random windows: observed 38 of 331
(11.5 %), permutation *P* = 0.277. The positional evidence is therefore *locus-level*, not
*count-level* — which is exactly why the supplement then asks the gene-level question in D and E.

Both panels come from `assembly/figS6v16_panels.py`; the file also draws C and F, so its header
(the two-size style, the ink-measuring `savefig`, the colour constants) is executed here.""")

code(banner("assembly/figS6v16_panels.py") + HEAD + "\n" + PA + "\npanel_A()\n")
code(banner("assembly/figS6v16_panels.py") + PB + "\npanel_B()\n")

md("""---

## 2. Row 2 — panel C: the UMAP behind Fig 6B, at four bin sizes

**Data chain**: `obs.UMAP_1/2` + `cell_type` of the lineage object joined to the per-nucleus MS
RNA score (`caches/s6c_umap_score.csv`, 5,962 nuclei). Panel C exists because a reader should not
have to trust one rendering choice: the same nuclei are hex-binned at single-cell resolution and
at gridsize 15, 25 (the Fig 6B choice) and 40, on one shared colour scale, with the aspect locked
to the data's own length-to-width ratio (1.63) so no tile is stretched.

Nothing about the embedding is re-derived here, and no conclusion depends on its shape.""")

code(banner("assembly/figS6v16_panels.py") + PC + "\npanel_C()\n")

md("""---

## 3. Row 3 — geometry is solved before D and E are drawn

Row 3 carries three columns: D, E and the shared vertical legend. The legend used to be its own
row (a horizontal band, 1.05 cm plus a 0.15 cm attachment gap); making it a column of the same row
recovered **1.20 cm of page height** (20.11 → 18.91 cm) at the cost of 14 % of the subtype-column
width. Row height, cell pitch and both font sizes did not move — house spec ranks added whitespace
above shrinking type, and shrinking type is the one thing that is not allowed.

The width is not a rounded guess. `PW = 1.78 + 12·CW` (gene-name column 1.54 + right pad 0.24,
twelve subtype columns) and `2·PW + 2·0.24 + LW = 17.00`, where `LW` is the legend's own width
estimate + 0.11 — the same number the delivered `figS6v16_run_all.sh` reads off the legend
script's stdout. The cell below measures it, solves the equation, and writes the result into the
environment that D, E and the legend read.""")

code('# measure pass: no LEG_W in the environment yet -> the script uses its own default canvas\n'
    'os.environ.pop("LEG_W", None)\n'
    + embed(LEGEND, "assembly/figS6v18_legend_right.py"))

code("""# The delivered recipe reads the legend's OWN width estimate off its stdout ("legend needs width ... cm")
# and takes LW = estimate + 0.11 (figS6v16_run_all.sh step 2). The legend cell just ran in this
# kernel, so its module-level `tot` is in the namespace. Re-measuring the PDF ink instead would
# add the X0=0.10 left margin (ink = X0 + tot) and give LW=1.72 - a geometry 0.05 cm wider per
# D/E column than the delivered page.
ink = tot
LW = round(ink + 0.11, 2)
CW = (17.00 - 2 * 0.24 - LW - 3.56) / 24.0
PW = 1.78 + 12 * CW
assert abs(2 * PW + 2 * 0.24 + LW - 17.00) < 0.02, (PW, LW)
for k, v in dict(PANEL_W="%.4f" % PW, PANEL_CW="%.4f" % CW, LEG_W="%.2f" % LW,
                 PANEL_H="5.75", PANEL_BODY_B="1.35", PANEL_ZMODE="1", PANEL_LEGEND="0",
                 PANEL_TITLE="0", PANEL_DCOLUMN="0", PANEL_CONTRAST="mol6").items():
    os.environ[k] = v
print("legend estimate %.2f -> LW %.2f | CW %.4f -> PW %.2f | row width %.2f (= 17.00)"
      % (ink, LW, CW, PW, 2 * PW + 2 * 0.24 + LW))""")

md("""## 4. Panel D — the MS top-1000 lands nowhere a reader can see

223 genes: the MS MAGMA top-1000, restricted to genes measurable on this platform (mapped minus
near-zero), each drawn across the six subtypes in both channels. The tested result is
*P* = 0.096 (RNA) and 0.745 (ATAC) against size-matched background sets — a null. The panel is
deliberately **not** labelled: a null panel with ten gene names on it invites the reader to hunt
for a pattern that the test says is not there, and a name list produced by |Δ| plus a minimum
row spacing changes whenever the page height changes (three geometries of the same data gave
8/10/10 names with only two in common).

**Why every off-scale triangle sits in the reference column.**  Under the per-gene z the
COP/NFOL column is the row-offset column: `L[:,0] ≡ 0` while the row mean is set by the other five
columns, so its |z| is systematically the largest (mean |z| 0.85/1.01 in the reference column vs
0.49/0.47 elsewhere). That is a property of the transform, not of the data — and it is *not* a
systematic colour bias: only 18 % (RNA) / 34 % (ATAC) of reference-column cells are negative, so
both directions occur. An alternative centring (five columns, reference excluded) was measured and
rejected: the reference column's mean |z| rises to 1.42/2.16 and the off-scale cells go from
3/11 to 54/86, i.e. the column becomes a stripe of triangles.

Historical note, kept because it is the kind of bug a gate cannot see: the triangles were once
drawn only when `abs(v) > HI and j > 0`, a leftover from the era when the reference column was a
constant-zero placeholder. Under the z display every off-scale cell *is* in column 0, so the
condition silently swallowed all of them and the legend key became a promise nothing fulfilled.
The guard is gone; the panel scripts now count what they draw and print it.""")

code('# run recipe for this panel (the same env the delivered driver exports)\n'
    'os.environ["PANEL_TAG"] = "_s6v18D"\n'
    'os.environ["PANEL_NAMES"] = "0"        # D stays unnamed (a null panel must not hand readers a pattern to hunt for)\n'
    + embed(LAND1000, "assembly/p0_8c_fig6c_top1000.py"))

md("""## 5. Panel E — the cis-window set lands on MOL6

38 genes whose platform gene lies inside a lead-SNP ±500 kb window, same scale, same contrast.
Here the tested result is positive in both channels: mean Δ = +0.145 (RNA) / +0.097 (ATAC)
against detection-breadth × signal-matched backgrounds, *z* = +2.87 / +3.07,
*P* = 0.0017 / 0.0007; the myelin/oligodendrocyte anchor set gives 0.0107 / 0.0080 and the same
panel's 293 non-window platform genes give 0.0127 / 0.0066 (leave-one-gene-out worst 0.0127
dropping PRR5L, 0.0038 dropping KCTD13).

**E is not a subset of D.**  Only 8 of the 38 window genes appear in the MS top-1000 at all, and
6 of those are measurable here, so the two panels share 6 rows. Two independent gene-level
definitions of the same GWAS, disagreeing — which is the finding the two panels are there to show.

The four named genes are the panel-A loci (PHLDB1, PHGDH, PLLP, BCAS1). Naming is independent of
the test: row order and the permutation both use Δ, and the names come from positional evidence
that is not on this panel's face. BCAS1 is the one that a reader may question: it is the
*progenitor*-end gene of the four (top subtype COP/NFOL, Δ_MOL6 = −1.71, ranked 38/38 within the
set, and −1.14 in the independent human mature-vs-precursor contrast), it is **not** in the MS
top-1000, so its negative Δ is not a counter-example to Fig 6A/B — it is the reason the manuscript
calls the BCAS1 coincidence a coincidence of gene identity.""")

code('# run recipe for this panel\n'
    'os.environ["PANEL_TAG"] = "_s6v18E"\n'
    'os.environ["PANEL_CONTRAST"] = "mol6"\n'
    'os.environ["PANEL_ZMODE"] = "1"\n'
    'os.environ["PANEL_DCOLUMN"] = "0"\n'
    'os.environ["PANEL_NAMES"] = "1"        # E names its genes; the D cell leaves the shared env at 0, so this restores 1 explicitly\n'
    'os.environ["PANEL_PRESPEC"] = "PHLDB1,PHGDH,BCAS1,PLLP"   # same env as figS6v16_run_all.sh\n'
    + embed(LAND38, "assembly/p0_8c_fig6c_landing.py"))

md("""## 6. Row 3 legend, row 4 — panels F and G

The legend is re-rendered at the solved column width so its canvas equals its ink and it can share
row 3's top edge. Its `off scale` key is drawn the way the mark actually looks on the body — a
1.2 pt **white** triangle at the clipped edge of a saturated cell — because a key that shows a
2.2 pt blue/red triangle is a different glyph from the one the reader is being asked to find.

**F** puts the same OPC-vs-Oligo contrast that Fig 6 A/B states for the own cohort into the three
chemistries, per nucleus, both channels, in one panel: two rows are the channels, three columns
are the chemistries, and the MWU *P* per panel is read from the ledger and *asserted* to equal the
label drawn ("P = " + p_fig), so a rewritten table cannot leave a stale P on the face.

**G** is the panel that used to be Figure 6 D: the 74-trait expression/chromatin Δ1 pair with the
rank columns, re-cut for the supplement by narrowing columns through environment variables only —
zero vector scaling, no font change. It sits in row 4 next to F so that "one nucleus at a time"
and "all 74 traits at once" are read as the same claim at two resolutions.""")

code(banner("assembly/figS6v16_panels.py") + PBOX + "\n" + PF + "\npanel_F_combined()\n")

code('# final pass at the solved column width (LEG_W was written by the geometry cell)\n'
    + embed(LEGEND, "assembly/figS6v18_legend_right.py"))

code("""os.environ.update(COMPACT="1", SUF="_s6v16G", WD="0.40", WK="0.40", WG="0.12",
                  LM="1.55", GAPC="0.42", RMC="1.10", TABD=TABD, OUTDIR=OUTD)
""" + embed(GHEAT, "assembly/make_figS6v16_G_heatmap.py"))

md("""---

## 7. What this notebook produces

Seven panel PDFs + the shared legend, all at the delivered sizes:

| output | size (cm) |
|---|---|
| `figS6_out/figS6v16_A_loci.pdf` | %(A)s |
| `figS6_out/figS6v16_B_perm.pdf` | %(B)s |
| `figS6_out/figS6v16_C_umapbins.pdf` | %(C)s |
| `figS6_out/fig6v32_C_top1000_landing_s6v18D.pdf` | %(D)s |
| `figS6_out/fig6v31_C_landing_s6v18E.pdf` | %(E)s |
| `figS6_out/figS6v18_DE_legend_v.pdf` | %(LEG)s |
| `figS6_out/figS6v16_F_cellboxes.pdf` | %(F)s |
| `figS6_out/fig6D_heatmap_pair_compact_s6v16G.pdf` | %(G)s |

The four-row page `Figure_S6_v16.pdf` (17.00 x 18.91 cm) is composed from these eight files by
the delivered driver `assembly/figS6v16_run_all.sh` (step 4/5, `figS6v16_compose.py`) — the
composition code is deliberately **not** part of this notebook (author instruction). The
word-level gate against the delivered page lives in `dev/`, not here.
""" % SIZES)

# =============================================================================
nb = nbf.v4.new_notebook()
nb.cells = cells
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.9"},
    "dsci": {"figure": "S6", "version": "v16", "sheet": SHEET,
             "generator": "dev/build_figS6_notebook.py",
             "embedded_scripts": [os.path.relpath(x, REPO) for x in
                                  (PANELS, LAND1000, LAND38, LEGEND, GHEAT)]},
}
out = os.path.join(REPO, "Figure_S6.ipynb")
with open(out, "w") as f:
    nbf.write(nb, f)
print("wrote %s (%d cells)" % (out, len(cells)))
