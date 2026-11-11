#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""build_figS7_notebook.py -- generate Figure_S7.ipynb from assembly/make_figS7_v30_panels.py.

Why generated: the notebook and the panel script must never disagree.  The script holds the code
(one "# %% PANEL X" block per panel); this builder splits it and wraps it in the markdown below.
Edit the code in assembly/make_figS7_v30_panels.py, edit the prose here, then re-run:

    python3 build_figS7_notebook.py
    jupyter nbconvert --to notebook --execute --inplace Figure_S7.ipynb

Format follows Figure_7.ipynb / Figure_6.ipynb: one section per panel, single-panel PDFs only,
no composing, no QA gates, no caption text on the face, English throughout.
"""
import os

import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "assembly", "make_figS7_v30_panels.py")
OUT = os.path.join(HERE, "Figure_S7.ipynb")

src = open(SCRIPT).read()
parts = src.split("# %% ")
setup = parts[0]
i = setup.index("import json")
setup = setup[i:]
# the notebook has no __file__: the kernel's cwd is the figS7/ directory
setup = setup.replace(
    'REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))',
    'REPO = os.path.abspath(os.getcwd())           # run this notebook from the figS7/ directory')
assert "__file__" not in setup, "setup cell still references __file__"

blocks = {}
for p in parts[1:]:
    name, _, body = p.partition("\n")
    blocks[name.strip()] = body
# the runner at the end of the script is compose-time furniture; the notebook calls each panel itself
for k in blocks:
    blocks[k] = blocks[k].split("# ------------------------------------------------------------------ run everything")[0]

MD = {}

MD["title"] = """# Figure S7 — the rank basis, the cross-channel controls, and the donor-level detection test (v30)

**What the figure claims.**  Figure 7 reports the MS risk signature as a mature-oligodendrocyte (OL)
contrast on the *size-adjusted* scDRS rank.  Figure S7 exists to show what that adjustment does and
what remains without it, in every channel that can address the same question: four human RNA cohorts
(GSE277430, GSE118257, GSE301908, GSE231585), a human snATAC cohort (GSE231584, same donors as
GSE231585, profiled separately), and mouse EAE spinal cord multiome (GSE250589, RNA and ATAC measured
in the *same nuclei*).  It also carries the donor-level re-test of Figure 7D's S1PR5 detection rule,
which is a per-person statement and therefore does not belong in the main figure.

| Panel | File | Size (cm) | Content |
|---|---|---|---|
| K | `S7v30_K_keyband.pdf` | 14.13 x 0.31 | legend artwork: the seven channels, one row, each named by accession |
| A | `S7v30_A_ladder.pdf` | 4.38 x 4.03 | 74 traits x 4 human RNA cohorts, ranked twice: as computed, and after the peak-domain-size trend is removed; MS ranks under each half |
| B | `S7v30_B_tightpool.pdf` | 3.27 x 3.24 | the same four cohorts inside a tight pool (MS + its 15 nearest neighbours on domain size, 16 traits), where MS is rank 1 in every cohort |
| C | `S7v30_C_human_samedonor.pdf` | 3.32 x 3.57 | RNA d1 (GSE231585) vs ATAC d1 (GSE231584) for the 74 traits: rho 0.722 raw, 0.137 after size residualisation (P 0.24); MS 2 -> 36 on the plotted axes |
| D | `S7v30_D_mouse_samenucleus.pdf` | 2.91 x 3.68 | mouse EAE multiome, RNA d1 vs ATAC d1 measured in the same 14,928 nuclei: rho 0.248, P 0.033; MS 11 -> 36 |
| E | `S7v30_E_ledger.pdf` | 15.17 x 3.00 | the 74 x 21 ledger: four human RNA cohorts + human ATAC + mouse RNA + mouse ATAC, three delta axes each, every cell a rank out of 74 |
| RB | `S7v30_RB_rankbar.pdf` | 4.06 x 0.56 | the rank colour bar shared by A, B and E |
| L2 | `S7v30_L2_axiskey.pdf` | 8.35 x 0.31 | what d1 / d2A / d2B subtract, in one row, in the same symbols Fig 7C uses |
| F | `S7v30_F_detection.pdf` | 8.35 x 3.49 | five Fig 7D targets x 15 populations in GSE118257: area = % detecting, fill = mean log2(CP10k+1) |
| G | `S7v30_G_donor.pdf` | 6.16 x 3.36 | the Fig 7D detection rule re-tested per donor: 27 evaluable donors in three cohorts, 22 pass |

The journal sheet (K above row 1; A | B | C | D on row 1; the rank bar and delta key on the band above
E; E full width; F | G on row 3; 17.00 x 14.36 cm; panel letters drawn by the composer) is assembled
outside this notebook: nothing here draws letters, composes panels, or runs QA gates."""

MD["env"] = """---

## 0. Environment, paths and versions

Three roots, all overridable.  Every table this notebook reads is **vendored** under `data/` or
`tables_v29/`, so the figure reproduces without the analysis tree; `DSCI_WORK_ROOT` is needed only to
find the Arial TTF.

| Variable | Default | Points at |
|---|---|---|
| `DSCI_WORK_ROOT` | two levels above `figS7/` | the analysis tree; here only for `fonts/ARIAL.TTF` |
| `FIGS7_OUT` | `figS7_out/panels_v30` | where the ten panel PDFs are written |
| `FIGS7_DATA` | `data/` | the vendored inputs (see `data/MANIFEST.csv`) |
| `FIGS7_TABLES` | `tables_v29/` | the two donor tables, produced by the handoff pipeline |

`export()` crops every panel to its own ink and **refuses to write a panel whose ink touches the
canvas edge** - that assert is what catches a clipped label, and it is the reason the scatter panels
carry 0.25 cm of head room."""

MD["inputs"] = """### 0.1 Where the inputs come from

| Used by | File | Produced by | Content |
|---|---|---|---|
| A, B, E, C | `data/FIG7_v1/tables/seven4_confound_ALL_cohorts.csv` | `upstream/confound_all_cohorts.py` | 74 traits x 4 cohorts: d1, ln(peak-domain size), raw and size-adjusted ranks. **The same file is a Fig 7 input**, which is why the two figures cannot quote different ranks |
| A, B | `data/SUPPORT_v1/tables/fig7Sb_matched_window.csv` | `upstream/matched_window_table.py` | the +-15 % domain-size window around MS (40 traits) and rho(ln n_peak, d1) per cohort |
| E | `data/SUPPORT_v1/tables/fig7Sc_rank_matrix_4cohorts.csv` | `upstream/axis_rho_and_rank_matrix.py` | the raw rank matrix (c0_1 ... c3_B); its `_1` columns are the as-computed d1 ranks, and their mean is the ledger's row order |
| E | `data/SUPPORT_v1/tables/S7_atac_crossaxis_human_full.csv` | `upstream/human_atac_threeaxis.py` | GSE231584 ATAC, 74 traits x three axes + split-half |
| E, D | `data/SUPPORT_v1/tables/S7_atac_crossaxis_mouse_full.csv`, `S7_rna_crossaxis_mouse_full.csv` | `upstream/mouse_atac_threeaxis.py` (+ composition fix), `upstream/mouse_rna_threeaxis.py` | mouse EAE GSE250589, ATAC and RNA channels on identical nuclei, labels, gene sets and seed (20261030, 14,928 nuclei, 25 libraries) |
| C | `data/A6_GSE231584_atac_v1/tables/delta1_vs_peakdomain_allchannels.csv` | `upstream/delta1_vs_peakdomain.py` | the size-matched ATAC d1 (`median_groupdiff`) plotted on C's y axis |
| C | `data/SUPPORT_v1/tables/figS7c_crosschannel_pairs.csv` | `upstream/crosschannel_pairs_table.py` | the rho / P that C is allowed to print, for the pair GSE231585 <-> GSE231584 |
| C | `data/A6_GSE231584_atac_v1/tables/same_donor_RNA_vs_CHROMATIN_delta1_74traits.csv` | delivered input (donor-paired ATAC chain) | the donor-paired ranks; kept for the record and for Table S7 - **not** the series C plots, which is why C's label changed in v29 (see section 4) |
| D | `data/SUPPORT_v1/tables/S7_mouse_same_nucleus_delta1.csv` + `.json` | `upstream/export_mouse_same_nucleus_delta1.py` | per-trait d1 for both channels of the same nuclei, plus rho / P / MS ranks / n_nuclei / n_libs recomputed from the score files |
| F | `data/A3_GSE118257_v1/tables/A3_GSE118257_target_stats.csv` | `upstream/target_stats_GSE118257_v7.py` | 30 targets x 15 populations: % detecting and mean depth-corrected expression, over 17,799 nuclei (9 donors, 20 libraries) |
| G | `tables_v29/S41_S1PR5_per_donor.csv`, `S41_summary_by_cohort.csv` | `handoff/donor_s1pr5_20261117/s41_donor_all_human_cohorts.py` | one row per donor: OL % , the strongest other population with >= 30 nuclei in that donor, the margin, the fold, the rule verdict |

The scoring itself is never re-run here: every number this notebook draws is read from a delivered
score table, and the scripts that made them are under `upstream/` and `handoff/`."""

MD["frozen"] = """### 0.2 Frozen criteria (fixed before the analysis, not chosen afterwards)

| | |
|---|---|
| scores | scDRS 1.0.2 `score_cell(n_ctrl=200)`, the same 74 MAGMA 10-kb gene sets (1000 genes each) as Fig 6; seeds 0 (A2/A4), 20261029 (A3), 20261030 (A7 and both mouse channels) |
| contrast | d1 = median(mature OL) - median(OPC); d2A = median(AST) - median(immune pool); d2B = median(neuron) - median(all non-neuronal nuclei) |
| group floor | a population counts in a unit only with >= 10 nuclei (the sweeps); the donor panel uses >= 30 nuclei per comparator population |
| size adjustment | rank after residualising d1 on ln(regulatory-domain size) = the *trend removed* half of A and Fig 7B |
| tight pool (B) | MS + its 15 nearest neighbours on ln(peak-domain size), ranked within those 16 |
| detection rule (Fig 7D, tested per donor in G) | margin >= 1.0 pp AND (OL >= 1.5x the comparator OR comparator = 0); margin = OL %% - the strongest other population's %% in the same donor |
| evaluable donor (G) | OL >= 30 nuclei and at least one comparator population >= 30 nuclei |

Statistical unit: donor where the cohort is donor-resolved, library where it is not.  GSE301908's 17
released libraries are 11 donors and are pooled to donors before G is drawn; GSE118257 is 9 donors
from 20 libraries (17,799 nuclei); GSE231585 is 8 donors from 15 libraries."""

MD["house"] = """### 0.3 House rules this figure implements

* Arial only, `pdf.fonttype = 42`, two size tiers (8.0 / 6.0 pt); the 8.0 pt tier is used by the panel
  letters, which the composer - not this notebook - draws.
* No title beside a panel letter: the information that used to sit there went into an axis title
  (G: `S1PR5: OL - strongest other population (pp)`), into a centred column-group head (F:
  `GSE118257`), or into the legend.
* No informal `(i)/(ii)` codes: seven consecutive letters.
* A legend sits where the ink it decodes is: the rank bar and the delta key form one band directly
  above the ledger and are left-aligned to the ledger's *measured* matrix edge, not to the panel edge.
* No dead padding: the ledger's head margin is 0.24 cm because the rotated cohort heads it was
  reserved for no longer exist; the detection map's size guide sits under its own `% detecting`
  heading instead of 2.1 cm below it.
* **One cell-type vocabulary across the figure family (v30).** The delta key prints
  `d1 = OL/MOL - OPC`, `d2A = AST - (MG+LYM)`, `d2B = NEU - non-NEU`, character for character the
  strings Fig 7C's column labels use, and panel F's column heads now use the same symbols for the four
  populations that appear in those definitions (AST / MG / LYM / NEU). Before v30 the key said
  `AST - immune` and `neuron - non-neuron` while Fig 7D already printed AST/MG/LYM/NEU: two spellings
  of one contrast in one paper. The oligodendrocyte subtypes (`OPCs`, `COPs`, `ImOlGs`, `Oligo1-6`) and
  the vascular labels keep the annotation table's own names - they are a source vocabulary, not part of
  any contrast. This costs two rows in the manuscript abbreviation table (`NEU`, `LYM`); see
  `NEEDS_ABBREVIATIONS.md`.
* Colour is asserted, not remembered: `check_palette()` recomputes the smallest CIE76 distance
  between adjacent strips, between any pair, and against the six semantic colours, under normal,
  deuteranopic and protanopic vision, and the render stops if a floor is broken."""

MD["K"] = """---

## 1. Panel K - identity key band

Legend artwork, not a panel: every colour the figure spends, in one row, each entry an accession.
The channel is stated once per block (RNA / ATAC) rather than once per entry, which is what buys the
room, and it is also what distinguishes GSE250589's two columns - the same accession legitimately
appears in both blocks because RNA and ATAC were measured in the same nuclei."""

MD["A"] = """---

## 2. Panel A - the rank, as computed and with the trend removed

**The claim.**  MS is rank 2 / 4 / 7 / 2 out of 74 on the raw d1 rank, and 1 / 1 / 3 / 1 once the
peak-domain-size trend is removed.  Both halves are shown side by side because Fig 7B prints only the
adjusted one; the reader needs the un-adjusted basis to see what the adjustment does.

**Row order is not a rank.**  Rows are ordered by peak-domain size (`n_peak_raw` descending), the
delivered sheet's key, which puts MS on row 37 and PLATELET on row 56 - both asserted.  The ledger E
uses a third order (mean of the four raw d1 ranks, MS on row 2) and Fig 7B a fourth (mean of the four
*size-adjusted* ranks, MS on row 1).  The legend names all of them; the face cannot."""

MD["B"] = """---

## 3. Panel B - the tight pool

**The claim.**  Inside a pool of 16 traits whose regulatory-domain size is within MS's neighbourhood
(the 15 nearest neighbours on ln n_peak), MS is rank 1 in all four cohorts.  This is the size-matched
statement in its own right, complementary to the residualisation in A's right half; the +-15 % window
itself (40 traits) and its rho are asserted here and tabulated in Table S7."""

MD["C"] = """---

## 4. Panel C - human RNA versus human ATAC in the same donors

**The claim.**  The two channels order the 74 traits similarly (rho 0.722) but the similarity is
mostly the shared size trend: after residualising on ln n_peak, rho falls to 0.137 (P 0.24).

**v29 correction.**  The delivered face said `MS 2 -> 21`.  2 is MS's rank in the plotted x-series,
but 21 is its rank in a *different* ATAC series - the donor-paired column of
`same_donor_RNA_vs_CHROMATIN_delta1_74traits.csv`, which is not what is plotted.  The plotted y is the
size-matched ATAC d1 (`median_groupdiff`), the series the rho is computed on, and MS is rank 36 there.
Both numbers are now computed from the plotted columns and asserted, so the arrow describes the axes it
sits on.  The donor-paired ranks stay in Table S7.

The dot colour is the delivered one (blue for the human pair, green for the mouse pair): a cross-channel
cloud does not belong to a cohort and must not borrow a cohort colour from the key band."""

MD["D"] = """---

## 5. Panel D - mouse EAE multiome, both channels in the same nucleus

**The claim.**  Where RNA and ATAC are measured in one nucleus there is no pairing to argue about, and
the two channels still disagree: rho 0.248 (P 0.033) over 74 traits, 14,928 nuclei, 25 libraries.  MS
sits at rank 11 in RNA and rank 36 in ATAC - the same direction as C, on a different species and a
different assay.  rho, P and both ranks are read from the JSON the export script wrote and re-asserted
against the CSV."""

MD["E"] = """---

## 6. Panel E - the 74 x 21 ledger, with the band that decodes it

**The claim.**  One table holds every rank the figure talks about: four human RNA cohorts, human ATAC,
mouse RNA, mouse ATAC, each on d1 / d2A / d2B.  MS is the red rule; AD and SCZ are the green and blue
rules.  Reading across the MS row gives the profile the earlier prototype carried as two separate
panels - those panels are withdrawn, and the numbers live here and in Table S7.

The rank bar (RB) and the delta key (L2) are drawn as their own PDFs and placed by the composer in one
band directly above the matrix, left-aligned to the matrix edge measured out of the produced PDF
(`matrix_x0`).  That measurement is not decoration: cropping to ink means the panel's left edge is set
by whatever text is furthest left - here the right-aligned red `MS` label - so the matrix does **not**
begin where the source constants say it does."""

MD["F"] = """---

## 7. Panel F - which cells carry the Fig 7D targets

**The claim.**  In GSE118257 (17,799 nuclei, 9 donors from 20 libraries) S1PR5 is detected in a small
minority of nuclei and is OL-biased, while NFE2L2 and CD52 are broad and ITGA4 is vascular; the
per-population denominators are printed under each column because the percentages mean nothing without
them.  Area is the detection rate, fill is mean log2(CP10k+1) (depth-corrected), both asserted against
the table by `verify_detection_ink.py` in the QA bundle."""

MD["G"] = """---

## 8. Panel G - the detection rule re-tested one donor at a time

**The claim.**  Figure 7D's S1PR5 rule is a pooled statement; here it is re-tested per donor.  Each dot
is a donor's own within-donor contrast: the % of that donor's OL nuclei detecting S1PR5 minus the % in
that donor's strongest *other* population.  22 of 27 evaluable donors meet the rule (8/11 GSE301908,
7/8 GSE118257, 7/8 GSE231585).

| Encoding | Meaning |
|---|---|
| colour | cohort (the key band's colour) |
| filled / open | rule met / not met; all five failures are thin-comparator donors (31-410 nuclei in the population that beat or nearly beat OL) |
| dot area | the comparator population's nucleus count in that donor |
| red dashed line | the +1.0 pp margin arm of the rule; the tinted band above it is where that arm holds |
| grey dotted line | 0, i.e. OL exactly ties the strongest other population |
| thick bar | cohort median margin |

GSE277430 is absent for a stated reason: its released object is an OL-enriched subset with no
non-lineage comparator, and its scores were built on normalised counts.  Diagnosis (16 MS, 11 Ctrl) is
deliberately **not** encoded: the between-group test on this within-donor statistic is null (pooled
Mann-Whitney P = 0.336) and the two cohorts trend in opposite directions, so it belongs in Table S7 and
the legend, not on the face."""

nb = nbf.v4.new_notebook()
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.9"},
}

# split the setup into the four blocks a reader needs in order: paths/style, colour, inputs, canvas
j = setup.index("# ------------------------------------------------------------------ the seven channel colours")
k = setup.index("# ------------------------------------------------------------------ inputs, with the assertions")
m = setup.index("# ------------------------------------------------------------------ canvas mechanics")
s_style, s_colour, s_inputs, s_canvas = setup[:j], setup[j:k], setup[k:m], setup[m:]

cells = [
    nbf.v4.new_markdown_cell(MD["title"]),
    nbf.v4.new_markdown_cell(MD["env"]),
    nbf.v4.new_code_cell(s_style.rstrip() + "\n"),
    nbf.v4.new_markdown_cell(MD["inputs"]),
    nbf.v4.new_code_cell("""import sys, matplotlib, numpy as _np, pandas as _pd
print("python    ", sys.version.split()[0])
for _m in (_np, _pd, matplotlib):
    print("%-12s %s" % (_m.__name__, _m.__version__))
# Kernel gate (2026-11-18, learned from Fig 6, where the delivered geometry turned out to be a
# matplotlib 3.9 product).  Panel sizes here are tight-bbox outcomes and the outcome moves with the
# renderer; v30 was delivered under python 3.8 / matplotlib 3.6.2 and the notebook reproduces the
# delivered ink exactly there.  Pin it, so a newer matplotlib fails loudly instead of quietly
# shifting every canvas by a few hundredths of a centimetre and then tripping the row gate.
_mv = tuple(int(x) for x in matplotlib.__version__.split("+")[0].split(".")[:2])
assert _mv[:2] == (3, 6), (
    "Fig S7 v30 was delivered under matplotlib 3.6.x (python 3.8 kernel); on %s the tight-bbox "
    "geometry drifts - re-measure the panels and re-run gates_v30 before accepting any output"
    % matplotlib.__version__)
print("kernel gate OK: matplotlib %s is the delivered renderer for this figure" % matplotlib.__version__)"""),
    nbf.v4.new_code_cell(s_colour.rstrip() + "\ncheck_palette()\n"),
    nbf.v4.new_markdown_cell(MD["frozen"]),
    nbf.v4.new_code_cell(s_inputs.rstrip() + "\n"),
    nbf.v4.new_markdown_cell("The helpers below are the only canvas machinery the panels use.  "
                             "`tw()` measures a string instead of estimating it, and `export()` crops "
                             "to ink while refusing any panel whose ink reaches the canvas edge."),
    nbf.v4.new_code_cell(s_canvas.rstrip() + "\n"),
    nbf.v4.new_markdown_cell(MD["house"]),
]
SECTION = [
    ("K", MD["K"]), ("A", MD["A"]), ("B", MD["B"]), ("C", MD["C"]), ("D", MD["D"]), ("E", MD["E"]),
    ("RB", "### 6.1 The band's left half (RB)\n\nThe rank bar is drawn as its own PDF so it can be "
           "cropped and placed independently; the composer puts it, with the delta key, in one band "
           "directly above the ledger."),
    ("L2", "### 6.2 The band's right half (L2)\n\nThe delta key takes whatever width the band has "
           "left after the rank bar, which is why it is drawn after E: the band's width comes from "
           "E's produced PDF."),
    ("F", MD["F"]), ("G", MD["G"]),
]
CALL = {
    "K": "key_band()", "A": "ladder_A()", "B": "tight_pool_B()", "C": "same_donor_C()",
    "D": "same_nucleus_D()", "E": "ledger_E(AVAIL)", "RB": "rank_bar()",
    "L2": 'axis_key(AVAIL - matrix_x0("E") - FR["RB"]["w"] - 2 * GAP + 0.30)',
    "F": "detection_map()", "G": "donor_s1pr5(6.40)",
}
NOTE = {
    "E": "\n# ledger_E prints the measured matrix edge that RB and L2 are then aligned to",
    "L2": "\n# the delta key takes whatever width the band has left after the rank bar",
}
for key, md in SECTION:
    cells.append(nbf.v4.new_markdown_cell(md))
    cells.append(nbf.v4.new_code_cell(blocks["PANEL " + key].strip() + "\n\n" + CALL[key]
                                      + NOTE.get(key, "") + "\n"))

nb.cells = cells
nbf.write(nb, OUT)
print("wrote %s | %d cells | %d markdown, %d code"
      % (OUT, len(cells), sum(c.cell_type == "markdown" for c in cells),
         sum(c.cell_type == "code" for c in cells)))
