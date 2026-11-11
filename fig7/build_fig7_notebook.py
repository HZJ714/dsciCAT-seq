#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""build_fig7_notebook.py -- generate Figure_7.ipynb from assembly/make_fig7_v17_panels.py.

Why generated: the notebook and the panel script must never disagree.  The script holds the code
(one "# %% PANEL X" block per panel); this builder splits it and wraps it in the markdown below.
Edit the code in assembly/make_fig7_v17_panels.py, edit the prose here, then re-run:

    python3 build_fig7_notebook.py            # writes Figure_7.ipynb
    jupyter nbconvert --to notebook --execute --inplace Figure_7.ipynb   # run it (see README)

Format follows Figure_6.ipynb / Figure_S7.ipynb: one section per panel, single-panel PDFs only,
no composing, no QA gates, no caption text on the face, English throughout.
"""
import os
import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "assembly", "make_fig7_v17_panels.py")
OUT = os.path.join(HERE, "Figure_7.ipynb")

src = open(SCRIPT).read()
parts = src.split("# %% ")
setup = parts[0]
# the script's shebang and module docstring are file-level furniture; in a notebook they are noise,
# so the setup cell keeps everything from the first import onward
i = setup.index("import os")
setup = ("# Shared setup: paths, the house two-tier font policy, palette constants and the geometry "
         "helpers every panel uses.\n" + setup[i:])
# notebook context: there is no __file__, the kernel's cwd is the fig7/ directory
setup = setup.replace(
    'sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # arial_style.py, the font policy\n'
    'from arial_style import use_arial, use_two_tier, FAMILY, open_spines   # noqa: E402\n\n'
    'REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))',
    'REPO = os.path.abspath(os.getcwd())           # run this notebook from the fig7/ directory\n'
    'sys.path.insert(0, os.path.join(REPO, "assembly"))    # arial_style.py, the font policy\n'
    'from arial_style import use_arial, use_two_tier, FAMILY, open_spines\n\n'
    '# Every savefig emits one line: "meta NOT subset; do not know how to subset; dropped".\n'
    '# That is fontTools speaking while matplotlib subsets the embedded Arial/Liberation face for\n'
    '# pdf.fonttype=42; the table it drops is the font metadata table "meta", not a glyph.  The PDFs\n'
    '# still embed only the Arial family (ArialMT / Arial-BoldMT), which is what the font gate checks.')
assert "__file__" not in setup, "setup cell still references __file__"
blocks = {}
for p in parts[1:]:
    name, _, body = p.partition("\n")
    blocks[name.strip().upper()] = body

MD = {}

MD["title"] = """# Figure 7 — External replication of the mature-OL vs OPC scDRS contrast (v17 FINAL, two-column sheet)

**What the figure claims.**  The mature-OL-minus-OPC contrast of the MS GWAS scDRS score, measured in
four independent human post-mortem MS cohorts (GSE277430, GSE118257, GSE301908, GSE231585) and in mouse
EAE spinal cord (GSE250589).  Panel A gives the contrast per donor / library, Panel B ranks the MS set
against the other 73 GWAS sets on that contrast, Panel C shows that the MS peak is specific to the
OL-maturation axis by re-projecting the same nuclei onto two other identity axes, and Panel D asks which
of the 18 gold-standard drug targets are actually reachable in the mature-OL compartment.

| Panel | File | Size (cm) | Content |
|---|---|---|---|
| A | `fig7v17_A_donors.pdf` | 6.58 x 5.90 | unit-level Δ1 = median(mature OL) − median(OPC) of the MS scDRS score: four human cohorts (donor- or library-level dots) + mouse GSE250589 libraries split Early/Peak/Late; Mann-Whitney *P* above each column (bare number, the "P =" prefix lives in the caption) |
| B | `fig7v17_B_adjranks.pdf` | 3.78 x 5.85 | 74 GWAS gene sets x 4 cohorts, portrait: colour = rank of cohort-level Δ1 **after residualising Δ1 on ln(regulatory-domain size)**, rows ordered by the mean of those four ranks; MS / AD / SCZ named. The raw ranking and the ±15 % matched-window control are in Fig S7 |
| C | `fig7v17_C_threeaxes.pdf` | 5.84 x 5.90 | the same nuclei, three identity axes (OL maturation / astro−immune / neuron−glia) as robust *z* within cohort x axis, 2 x 2 tiles; MS (red) and SCZ (blue) dots carry `#rank` on the axis they sit on; the single `*` hangs **directly above the neuron−glia violin of GSE231585** |
| D | `fig7v17_D_targets.pdf` | 17.00 x 7.90 | 18 audited drug targets: two within-cohort log2-OR blocks (mature OL vs OPC; OL lineage vs non-lineage) over three cohorts, a drug-target flag block (label / stated MOA / used in MS), and the % of nuclei detecting each gene in GSE301908 by population. The red box marks the one row that passes the detection rule |

The journal sheet (A | B | C on row 1, D on row 2, 17.00 x 14.34 cm, panel letters, final crop) is
assembled outside this notebook; nothing here draws letters or composes panels.

**About the `*`.**  It is not a significance marker.  GSE231585 carries 62,911 nuclei of which only 717
are annotated neurons (1.1 %; 4-133 per sample), so under the frozen ">= 10 nuclei per group" rule two of
its 15 samples drop out of the neuron−glia axis (nB = 13 while n1 = nA = 14).  The caveat is specific to
that axis, which is why the mark sits over that violin and not on the tile title.  The other cohorts'
neuron pools are 7-100x larger (GSE277430 76,436; GSE118257 5,031; GSE301908 capped at <= 12,000 per
lineage by the stratified draw), so no other tile is marked."""

MD["env"] = """---

## 0. Environment, paths and versions

Three roots, all overridable; nothing in the tables below is produced by this notebook, `upstream/`
contains the scripts that produced it (see `upstream/README.md`).

| Variable | Default | Points at |
|---|---|---|
| `DSCI_WORK_ROOT` | two levels above `fig7/` | the analysis tree; needed here only for the Arial TTF. Every table the notebook reads is vendored (2026-11-18, provenance audit G1) |
| `FIG7_DATA` | `data/` | the vendored A/B/C inputs, laid out like `routeA_out/` (see `data/MANIFEST.csv`) |
| `FIG7_OUT` | `fig7_out/` | where the four panel PDFs are written |
| `FIG7_TABLES` | `tables_v17/` | the small vendored tables that Panel D reads |

Geometry switches (defaults reproduce the delivered sheet): `W_A/H_A`, `W_B/H_B`, `W_C/H_C` request each
panel's canvas in cm, `FIG7D_FW` sets D's canvas width, `FIG7A_PVALS=bare|P=|off` chooses how A prints
its *P* values, `STAR_POS=violin|label|title` moves the GSE231585 marker (only `violin` is delivered;
the two rejected positions are kept switchable so the archived PDFs stay reproducible),
`FIG7D_LETTER=1` opts into printing the letter "D"."""

MD["inputs"] = """### 0.1 Where the inputs come from

Every input this notebook reads is **inside the bundle**: the A/B/C inputs are vendored under `data/`
(mirroring the `routeA_out/` layout) and Panel D's under `tables_v17/`. `data/MANIFEST.csv` records, per
file, the bytes/rows/md5, the producing script, which panel reads it and what the columns mean; the
vendored copies are byte-identical to the analysis-tree originals, and re-running this notebook after
vendoring reproduced all four delivered panels (page size, every word, coordinates to 0.02 cm, drawing
counts, fonts). Nothing here needs the analysis tree except the Arial TTF.

| Used by | File | Produced by | Content |
|---|---|---|---|
| A | `data/FIG7_v11/tables/fig7A_donors_delta1.csv` | `upstream/D9_fig7_assembly.py` | per-unit \u03941 for the four human cohorts + mouse GSE250589 libraries |
| A | `data/FIG7_v11/tables/fig7A_between_donor_P.csv` | `upstream/D9_fig7_assembly.py` | between-unit Mann-Whitney MS vs Ctrl per cohort (recomputed in-cell and asserted) |
| B | `data/FIG7_v1/tables/seven4_confound_ALL_cohorts.csv` | `upstream/A8c_confound_unified.py` | \u03941 matrix + raw ranks + ranks after residualising on ln n_peak; **the same file is a Fig S7 input** (`figS7/data/FIG7_v1/tables/`), which is why the two figures cannot quote different ranks |
| C | `data/A2_GSE277430_v1/tables/trait_sweep_dual_axis.csv` | `upstream/A2_GSE277430_delta2_mirror.py` | the \u03941 sweep merged with the two mirror axes (n_ctrl=100, seed=0) |
| C | `data/A3_GSE118257_v1/tables/trait_sweep_delta2_GSE118257_74.csv` | `upstream/A3_GSE118257_delta2_full74.py` | three axes, MIN_CELL=20 nuclei per donor per group |
| C | `data/A4_GSE301908_v1/tables/trait_sweep_delta2_GSE301908_74.csv` | `upstream/A4_GSE301908_delta2_full74.py` | three axes from one n_ctrl=200 pass, \u2264 12,000 nuclei per lineage |
| C | `data/A7_GSE231585_rna_v1/tables/trait_sweep_delta2_GSE231585_74.csv` | `upstream/A7_GSE231585_delta2_full74.py` | three axes; its dB axis rests on 717 annotated NEU nuclei and is the one Panel C marks `*` |
| D | `tables_v17/S13_two_contrasts.csv` | `upstream/s13_two_contrasts.py` | within-cohort log2 OR, mature OL vs OPC and lineage vs non-lineage, BH within cohort and direction |
| D | `tables_v17/S23_new_targets_stats.csv` | `upstream/s23_expand_new_targets.py` + `data/EXPANSION_33/expansion33_rows.csv` | the expansion atoms scored against the released populations |
| D | `tables_v17/S26_rowset_C12.csv`, `S26_target_ledger.csv` | `upstream/s27_target_ledger.py` | the 83-atom ledger (80 scored drug-target atoms) and the 18 rows the panel draws |
| D | `tables_v17/S29_flags_v2.csv` | `upstream/s32_flags_v2.py` | the 18 drawn rows x 3 evidence flags, author-side openFDA + Europe PMC re-audit |
| D | `tables_v17/S15_det_GSE301908.csv`, `S15_n_GSE301908.csv` | `upstream/s15b_allgene_pop_counts.py` | detected / total nuclei per gene x population, pooled over the whole released matrix (293,970 nuclei). The same script also wrote GSE277430 / GSE118257 / GSE231585 files; only GSE301908 is drawn, because it is the one cohort where every population reaches >= 2,000 nuclei (2,123 in the smallest, LYM) - GSE118257 bottoms out at 352 and GSE231585 at 299, while the GSE277430 OL subset carries no non-lineage populations at all |
| D, B | the 74 MAGMA gene sets behind every sweep | the shared scDRS resource `$DSCI_FY_ROOT/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs` | **the same gene-set file Fig 6 scores**, so a trait number means the same thing in both figures |

Two inputs are deliberately **not** vendored because this notebook never opens them:

* `A6_GSE231584_atac_v1/tables/trait_chromatin_sizeMatched.csv` (`upstream/A6_stepC2_sizematched.py`)
  feeds `A8c`, whose *output* is vendored above - Panel B reads ln n_peak through `seven4_confound`;
* `M1b_authoranno_v1/tables/per_lib_delta1_author_vs_self.csv` (`upstream/M1b_authoranno.py`) is the
  mouse \u03941 source behind `fig7A_donors_delta1.csv`, which is vendored.

`upstream/s15b_allgene_pop_counts.py` is the only step in the chain that reads the cohort matrices
(~1.1 GB); everything else is table-only, so re-deriving Panel D's inputs from a clone costs one heavy
run plus minutes. `s23`/`s27` read `data/EXPANSION_33/expansion33_rows.csv`, the ASCII extract of the
author's expansion memo: the memo is an internal Chinese-language document and is not pushed, and the
English status strings in the CSV are exactly the ones the delivered ledger carries (before 2026-11-18
that translation had only ever happened ad hoc, so the delivered `S26_target_ledger.csv` was not
re-derivable by any committed script).

The cohort matrices and annotations behind the sweep tables are listed in `upstream/README.md`; the
score tables are read as delivered, so this notebook re-derives ranks and robust *z* only, never scores.
The donor-level re-test of Panel D's detection rule is **not** part of Figure 7 and is handed to
Fig S7 (`figS7/upstream/`, see that bundle)."""

MD["versions"] = """### 0.2 Software versions

Two environments are involved and both are reported, because the numbers on the face were produced in
the first one and only rendered in the second:

| | scoring / analysis environment (all sweep tables, Fig S7, the donor test) | this notebook (rendering) |
|---|---|---|
| python | 3.9.23 (`$DSCI_PY`) | whatever kernel runs the notebook |
| scdrs | **1.0.2** (`score_cell`, `n_ctrl=200`, `pp.preprocess`) | not imported - scores are read as delivered |
| anndata / numpy / pandas / scipy / statsmodels | 0.10.9 / 1.26.4 / 2.3.1 / 1.13.1 / 0.14.5 | printed below |
| matplotlib | 3.9.4 | printed below |

`upstream/README.md` pins the same versions for every script that produced an input table."""

MD["frozen"] = """### 0.3 Frozen criteria (fixed before the analysis was run, not chosen afterwards)

| | |
|---|---|
| contrast | Δ1 = median(MS scDRS `norm_score` of mature OL) − median(OPC), per unit (donor / library / sample); cohort-level Δ1 = median over units |
| scoring | `scdrs` 1.0.2 `score_cell(n_ctrl=200)` on the same 74 MAGMA 10-kb gene sets (1000 genes each) in every cohort, after `scdrs.pp.preprocess`; seeds 0 (A2/A4), 20261029 (A3), 20261030 (A7) |
| unit floor | a group enters a per-unit median only with >= 10 nuclei (MIN_CELL); this is what removes 2 of GSE231585's 15 samples from Δ2B |
| rank | descending rank over the 74 traits, `method="min"`, rank 1 = strongest OL contrast |
| B colour | rank of Δ1 **after residualising on ln n_peak_raw** (OLS on the traits present in that cohort). Rationale: Δ1 tracks regulatory-domain size in all four cohorts, and the size-adjusted ordering is what the main figure reports; the four Spearman rho values (−0.838 / −0.691 / −0.778 / −0.873) are recomputed in the shared setup from the same table and asserted there, so the caption cannot drift from the figure. MS ranks 1 / 1 / 3 / 1; the raw ordering and the ±15 % matched-window control are in Fig S7 |
| C transform | robust *z* = (Δ − median) / (1.4826 x MAD) within cohort x axis, so the four tiles share one y-axis; violins are KDEs evaluated only inside the visible window (±4.7 / 5.5), ranks and IQR boxes use all 74 traits |
| C marker | exactly one `*`, over GSE231585's neuron−glia violin: 717 annotated neurons, 1.1 % of that cohort, nB = 13 vs n1 = nA = 14. It is a data-thinness caveat, not a p-value |
| D rule (pooled, as drawn) | red box iff mature-OL detection % > every non-lineage population's by >= 1.0 pp AND >= 1.5x, on the pooled GSE301908 rates; only S1PR5 passes |
| D rule (pooled, as drawn) | red box iff mature-OL detection % > every non-lineage population's by >= 1.0 pp AND >= 1.5x, on the pooled GSE301908 rates; only S1PR5 passes. Resolving the same rule per donor is a Fig S7 question, not a Fig 7 one |
| mouse | GSE250589 spinal cord, FACS GFP+ (Sox10-cre); its 74-trait sweep is not in Panel C (78-neuron neuronal layer) and is reported in S7 |"""

MD["house"] = """### 0.4 House style baseline (geometry, not cosmetics)

* Two font tiers only: 8.0 pt (panel titles, axis titles, row heads) and 6.0 pt (ticks, annotations,
  statistics). Arial family; face text black except the named traits and the axis colours.
* Panels are drawn at their final width. No delivered PDF is rescaled to fit (a rescale would move
  every 6 pt glyph off the tier); the sheet is composed from the four files above.
* Rotated labels are anchored by measured text width (`rot_anchor`), so a 45/60 degree label sits
  centred under its own column instead of drifting into the neighbour's.
* `#rank` labels are placed by a bounded ring search (`place_rank_labels`) that tests the measured box
  against every other label, every dot and the figure edge - collisions are resolved by moving labels,
  never by shrinking fonts or dropping the numbers."""

MD["panelA"] = """---

## 1. Panel A — unit-level Δ1 across four human cohorts and mouse EAE stages

**Data chain**: `fig7A_donors_delta1.csv` (one row per donor / library, built by `upstream/D9_...` from
the per-cohort score tables) + the between-unit *P* table. The cell below recomputes the human *P*
values from the unit table and asserts them against the delivered file, so the number on the face is
tied to the data rather than to a spreadsheet.

Switches: `W_A`/`H_A` canvas, `FIG7A_PVALS` (`bare` prints `0.091`, `P=` prints `P=0.091`, `off` prints
nothing and moves the value to the caption)."""

MD["panelB"] = """---

## 2. Panel B — 74 GWAS gene sets x 4 cohorts, portrait rank heatmap

**Data chain**: the four sweep tables are re-read here to rebuild the cohort-level Δ1 matrix; the
size-adjusted ranks are computed in-cell (OLS residual of Δ1 on ln n_peak, then re-ranked) and asserted
against the delivered `seven4_confound_ALL_cohorts.csv`. Rows are ordered by the mean of the four
adjusted ranks, so the strongest OL contrasts sit at the top; MS / AD / SCZ are named in the margin.

The regulatory-domain size `n_peak_raw` is counted on the public human snATAC reference (GSE231584,
GRCh38, 121,505 deduplicated peaks) inside TSS ± 100 kb of each gene set's genes."""

MD["panelC"] = """---

## 3. Panel C — the same nuclei, three cell-identity axes

**Data chain**: the same four sweep tables, which also carry Δ on two further axes computed from the same
per-nucleus MS scores: astro−immune (median AST − median IMM) and neuron−glia (median NEU − median GLIA,
GLIA = AST + IMM + OLIGO + OPC + VAS). Within cohort x axis the 74 Δ values become robust *z*; the dots
carry their rank on the axis they sit on.

The marker logic is in `place_axis_marks`: `TILES` stores the **axis index** whose pool is thin (2 for
GSE231585), and the `*` is placed over that column, outside the top spine, after a measured search that
rejects any position touching a label, the tile title, a dot or the figure edge."""

MD["tail"] = """---

## 1. What this notebook produces

Four panel PDFs, one per section, written to `FIG7_OUT` (default `fig7_out/`). Sizes are measured from
the produced PDFs, not from the nominal canvas (`savefig.bbox="tight"` means the two differ):

| file | panel | measured size | content |
|---|---|---|---|
| `fig7v17_A_donors.pdf` | A | 6.58 x 5.90 cm | donor/library-level \u03941 per cohort, between-unit MWU *P* |
| `fig7v17_B_adjranks.pdf` | B | 3.78 x 5.85 cm | 74 traits x 4 cohorts, size-adjusted rank colour |
| `fig7v17_C_threeaxes.pdf` | C | 5.84 x 5.98 cm | three identity axes per cohort, MS/SCZ/AD ranks, `*` on the thin NEU axis |
| `fig7v17_D_targets.pdf` | D | 16.20 x 7.90 cm | 18 audited drug-target rows x detection tiers |

Assembling these into the two-column sheet, adding the panel letters and running the QA gates is
deliberately **not** in this notebook (house spec \u00a71: assembly and verification are journal
production, not content). The delivered composed sheet is

```
figures/SUPPL_ANALYSES_20261112/FIG7_2COL_20261117/pdf/fig7_composed_v17c_2col.pdf   17.00 x 14.42 cm
```

with A | B | C on row 1 and D on row 2; the composing script is `s38_compose_fig7_2col.py` in that
directory, and the gates are `font_gate.py` / `collide_gate.py` / `verify_text_on_data.py` plus the
face-number register `s39_face_register.py`.

**Delivery condition.** Every `assert` in the four panel sections must pass: Panel A recomputes the
between-unit *P* values from the vendored per-donor table, Panel B recomputes the size-adjusted ranks,
Panel C re-derives each printed rank from the plotted column and places the single `*` only where a free
position exists, Panel D checks the 18 drawn rows against the ledger. A notebook that runs with any
assertion failing, or that writes a PDF whose measured size differs from the table above, is not the
delivered figure."""

MD["panelD"] = """---

## 4. Panel D — which targets are reachable in the mature-OL compartment

**Data chain**: 18 audited targets (`S26_rowset_C12.csv`), two within-cohort log2 odds ratios per cohort
(`S13_two_contrasts.csv`, `S23_new_targets_stats.csv`), drug-target flags (`S29_flags_v2.csv`: does a
label / a stated MOA name this target, and is it used in MS through this target - the flags describe the
drug-target relation, not localisation), and the % of nuclei detecting each gene by population in
GSE301908 (`S15_det_*.csv` / `S15_n_*.csv`, pooled over 293,970 nuclei).

The red box is a display rule on those pooled rates (>= 1.0 pp and >= 1.5x over every non-lineage
population); S1PR5 is the only row that passes it. The ledger is asserted here so that "18 of 80 scored"
in the caption is a counted statement, not a remembered one: 83 rows, 3 of them non-drug atoms."""

nb = nbf.v4.new_notebook()
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.9"},
}
cells = [nbf.v4.new_markdown_cell(MD["title"]),
         nbf.v4.new_markdown_cell(MD["env"]),
         nbf.v4.new_code_cell(setup.split('def savefig')[0].rstrip() + "\n"),
         nbf.v4.new_markdown_cell("The block below is the shared setup: paths, the two-tier font policy "
                                 "and the geometry helpers every panel uses."),
         nbf.v4.new_code_cell("def savefig" + setup.split('def savefig')[1]),
         nbf.v4.new_markdown_cell(MD["inputs"]),
         nbf.v4.new_markdown_cell(MD["versions"]),
         nbf.v4.new_code_cell("""import sys, scipy, statsmodels, anndata, matplotlib
print("python   ", sys.version.split()[0])
for m in (np, pd, scipy, statsmodels, anndata, matplotlib):
    print("%-12s %s" % (m.__name__, m.__version__))
try:                              # scoring runs upstream; this notebook only reads the score tables
    import scdrs
    print("scdrs      ", scdrs.__version__, "(this kernel)")
except ImportError:
    print("scdrs       1.0.2 - not in this kernel; the upstream sweeps ran it in")
    print("              $DSCI_PY (python 3.9, scdrs 1.0.2), see upstream/README.md")
print("PyMuPDF is load-bearing: every size printed below is measured from the produced PDF,")
import fitz
print("              never taken from the nominal canvas (tight bbox != figsize).")
# Kernel gate (2026-11-18, learned from Fig 6): panel sizes are matplotlib tight-bbox OUTCOMES and the
# outcome moves with the renderer - Fig 6's delivered row is a matplotlib 3.9 product and lands 0.05 cm
# short under 3.6. This sheet is the opposite case: v17 was delivered under python 3.8 / matplotlib 3.6.2
# and the vendored-input re-run reproduces it exactly there. So the version is pinned, not assumed: a
# reader on a newer matplotlib gets a clear message instead of a quietly drifted page.
_mv = tuple(int(x) for x in matplotlib.__version__.split("+")[0].split(".")[:2])
assert _mv[:2] == (3, 6), (
    "this sheet was delivered under matplotlib 3.6.x (python 3.8 kernel); on %s the tight-bbox geometry "
    "drifts - re-measure the panels and re-run the gates before accepting any output"
    % matplotlib.__version__)
print("kernel gate OK: matplotlib %s is the delivered renderer for this sheet" % matplotlib.__version__)"""),
         nbf.v4.new_markdown_cell(MD["frozen"]),
         nbf.v4.new_markdown_cell(MD["house"]),
         nbf.v4.new_markdown_cell(MD["panelA"]),
         nbf.v4.new_code_cell(blocks["PANEL A"].strip()),
         nbf.v4.new_markdown_cell(MD["panelB"]),
         nbf.v4.new_code_cell(blocks["PANEL B"].strip()),
         nbf.v4.new_markdown_cell(MD["panelC"]),
         nbf.v4.new_code_cell(blocks["PANEL C"].strip()),
         nbf.v4.new_markdown_cell(MD["panelD"]),
         nbf.v4.new_code_cell(blocks["PANEL D"].strip()),
         nbf.v4.new_markdown_cell(MD["tail"])]
nb.cells = cells
nbf.write(nb, OUT)
print("wrote %s | %d cells | %d markdown, %d code"
      % (OUT, len(cells), sum(c.cell_type == "markdown" for c in cells),
         sum(c.cell_type == "code" for c in cells)))
