# Figure S7 — legend and face-number contract (FINAL = v30, 2026-11-18)

Delivered sheet: **17.00 x 14.28 cm** (composed from the ten panels; the sheet PDF is working-tree
only, not versioned), panels **A–G** plus the key band
(K) and the legend band (RB + L2) that the composer places.  Every panel PDF in `figS7_out/panels_v30/`
is produced by `assembly/make_figS7_v30_panels.py`, which is also the source of `Figure_S7.ipynb`
(see `build_figS7_notebook.py`).  Gate battery: `FAILURES: none` (section 4).

This file has five parts: (1) the journal legend, ready to paste; (2) the numbers that are deliberately
*not* on the face and must therefore be quoted from here or from Table S7 (house spec §14.7); (3) the
mapping from the old manuscript references (S7A…S7I of v22) to the final letters; (4) typesetting and
compliance record; (5) what changed in v29 and v30 and why.

---

## 1. Journal legend (ready to paste)

**Supplementary Figure S7 | The oligodendrocyte-lineage risk axis is size-robust, channel-specific, and
the S1PR5 detection call survives being resolved one donor at a time.**

**(A)** Rank ladder of the 74 GWAS gene sets on the mature-oligodendrocyte-minus-precursor axis (Δ1) in
the four human post-mortem cohorts (GSE277430, GSE118257, GSE301908, GSE231585), shown twice: the rank
as computed, and the rank after the regulatory-domain-size trend is removed (colour = rank of 74; the
shared rank bar sits above **E**). Rows are ordered by regulatory-domain size (*n*~peak~ descending, the
ordering of the delivered sheet), which places MS on row 37 and PLATELET on row 56; this is a size
ordering, **not** a rank ordering. The MS rank under each column is **2 / 4 / 7 / 2** as computed and
**1 / 1 / 3 / 1** with the size trend removed.

**(B)** Size-matched tight pool: for each cohort, MS plus the 15 gene sets whose peak-domain size is
nearest to MS's, ranked **within that pool of 16**. MS ranks first in all four cohorts (1/16 in each).
The ±15 % window used for the sensitivity analysis (40 gene sets) and its ρ(ln *n*~peak~, Δ1) are in
Table S7.

**(C)** Δ1 concordance between the two channels profiled from the **same donors** (RNA: GSE231585;
chromatin: GSE231584 snATAC, ArchR gene scores; the two assays were built separately, so this is a
trait-level correspondence, not a paired measurement), 74 gene sets. Spearman ρ falls from **0.722**
(raw) to **0.137** (*P* = 0.24) once both channels are residualised on regulatory-domain size, so the
unadjusted cross-channel correlation is dominated by gene-set size. On the plotted axes MS is rank
**2** in RNA and rank **36** in the size-matched chromatin series. (In the donor-paired table the same
two ranks are **2 → 21**; that series is not plotted here and the paired ranks are in Table S7.)

**(D)** The same comparison in the mouse EAE spinal-cord multiome (**GSE250589**), where both channels
were measured **in the same nuclei** (14,928 nuclei, 25 libraries, identical cell-type labels, gene sets
and seed): ρ = **0.248**, *P* = **0.033**; MS ranks **11** in RNA and **36** in chromatin.

**(E)** Full axis ledger: 74 gene sets × 21 columns — seven channels (the four human RNA cohorts, human
ATAC GSE231584, mouse RNA GSE250589, mouse ATAC GSE250589) × three axes (Δ1, Δ2A, Δ2B); colour = rank of
74 on the shared rank scale; the seven channel strips above the matrix and the three axis strips below
it carry the colours keyed in the top band and in the band above **E**. Rows are ordered by the mean of
the four human cohorts' **as-computed** Δ1 rank (MS on row 2). Note that this is a third row key: Fig 7B
orders its rows by the mean **size-adjusted** rank (MS on row 1), and **A** orders by domain size (MS on
row 37). The MS row is marked in red; AD and SCZ rows in green and blue.

**(F)** Drug-target × cell-type detection map for **GSE118257** (17,799 nuclei; 9 donors from 20
libraries; per-population *n* printed under each column: OPCs 352, COPs 242, ImOlGs 207, Oligo1 1,129,
Oligo2 1,839, Oligo3 775, Oligo4 1,579, Oligo5 1,167, Oligo6 1,484, AST 1,242, MG 796, LYM 423,
Endo 836, Peri 697, NEU 5,031 - v30 prints the four populations that appear in the axis definitions
with the same symbols the delta key and Fig 7C use; the OL subtypes and the vascular labels keep the
annotation table's own names, because they are a source vocabulary rather than a contrast). Marker diameter = **percentage of nuclei in which the target is
detected** (≥ 1 read; guide circles 1 % / 10 % / 50 %); colour = **mean log2(CP10k+1) among detecting
nuclei** (2.2–4.4 scale, depth-corrected as log2(count × 10,000 / detected genes + 1)); an open circle
marks a population with no detecting nucleus. S1PR5 is detected in 6.8–14.9 % of differentiated
oligodendrocyte states and in COPs (14.9 %) but in only 0.85 % of OPCs and ≤ 3.6 % of any non-lineage
population; S1PR1 is endothelial/pericyte-biased (24.8 % / 19.9 %) and ≤ 1.2 % in the lineage; the
lymphocyte targets CD52 (0 % in all nine lineage columns, 15.8 % in immune) and ITGA4 (≤ 1.7 % in the
lineage, 9.5 % in immune) are the specificity controls; NFE2L2 is broadly detected (6.4–38.8 %) and is
the non-specific reference. MS4A1 is absent from the GSE118257 count matrix, so the fifth row is ITGA4.

**(G)** The **F** call resolved per donor. For each of the 27 evaluable donors of the three human
cohorts that carry non-lineage populations, the plotted value is that donor's own within-donor contrast:
S1PR5 detection percentage in mature oligodendrocytes minus its percentage in the **strongest other
population of the same donor** (percentage points). A population enters a donor's comparison only if
that donor contributes ≥ 30 nuclei of it, and a donor is plotted only if its OL pool is ≥ 30 nuclei and
at least one comparator survives. The red dashed line is the **+1.0 pp** margin arm of the rule Fig 7D
applies to pooled data, and the tinted band above it is where that arm holds; the grey dotted line is
zero, i.e. OL exactly ties the strongest other population. Filled dots meet the full rule
(margin ≥ 1.0 pp **and** OL ≥ 1.5× the comparator, or the comparator is 0): **22 of 27** donors
(8/11 GSE301908, 7/8 GSE118257, 7/8 GSE231585). Dot area is the comparator population's nucleus count;
the thick bar is the cohort median. All five failures are thin-comparator donors (31–410 nuclei in the
population that beat or nearly beat OL), i.e. they are low-information donors rather than counter-
examples. GSE277430 is not represented: its released object is an OL-enriched subset with no
non-lineage comparator, and its scores were built on normalised counts. Diagnosis is not encoded — see
section 2.

---

## 2. Numbers carried by the legend / Table S7 instead of the face

| Item | Value | Where it comes from |
|---|---|---|
| MS profile across the 21 ledger columns (Δ1 / Δ2A / Δ2B per channel) | human RNA 2/65/44, 4/71/59, 7/66/60, 2/67/74; mean over the four human RNA cohorts **3.8 / 67.2 / 59.2**; human ATAC **31 / 65 / 43**; mouse ATAC **36 / 64 / 67**; mouse RNA **11 / 74 / 59** | the MS row of `data/SUPPORT_v1/tables/fig7Sc_rank_matrix_4cohorts.csv` and the ATAC/mouse rank columns — i.e. read straight off **E**; these are the numbers the withdrawn v22 profile panels (F, G) used to print |
| three-basis ranks of the four identity sets | MS 2/1/1, 4/1/2, 7/3/4, 2/1/1; PLATELET 1/11/2, 2/5/1, 1/8/3, 1/7/2; AD 25/18, 41/49, 26/27, 9/12; SCZ 67/68, 73/74, 63/42, 67/71 (raw / size-stratified / matched-window) | v22 panel C, withdrawn; the matched-window ranks are defined for MS and PLATELET only |
| size trend on the face of the ranking | ρ(ln *n*~peak~, Δ1) per cohort, ±15 % window = 40 traits | `fig7Sb_matched_window.csv` |
| axis orthogonality | Δ1 vs Δ2A ρ −0.84 / −0.69 / −0.78 / −0.87 | `fig7Sc_axis_pair_rho.csv` (the withdrawn v22 panel's content) |
| donor rule, no floor | GSE301908 **6/11** rule met, OL highest 11/11 (sign test P = 4.9 × 10⁻⁴) | `handoff/.../tables/S40_summary.csv` — the 2026-11-17 delivered figure; kept for the record because the panel now uses the ≥ 30-nucleus comparator rule |
| floor sweep (GSE301908, OL-subset convention) | rule met 4/7 at ≥ 30 nuclei, 3/5 at ≥ 50, 2/3 at ≥ 100 | `handoff/.../tables/S40b_summary_floors.csv` |
| diagnosis, per donor | 16 MS, 11 Ctrl among the 27 evaluable; median margin MS 5.08 pp vs Ctrl 6.22 pp, pooled Mann-Whitney **P = 0.336**; per cohort P = 0.571 (GSE118257), 1.000 (GSE231585), 0.376 (GSE301908); rule met 12/16 MS, 10/11 Ctrl | `handoff/donor_s1pr5_20261117/tables/S41_S1PR5_per_donor_with_dx.csv` and that directory's README. The two cohorts trend in opposite directions, which is why diagnosis is a covariate in the table and not an encoding on the face |
| donor / library counts | GSE301908 11 donors from 17 released libraries (293,970 nuclei); GSE231585 8 donors from 15 libraries (62,911 nuclei); GSE118257 9 donors from 20 libraries (17,799 nuclei) | the cohort READMEs; GSE301908's `ms1…ms16` labels are subject identifiers, **not** MS status (8 patients, 3 controls) |
| GSE231585 neuron pool | 717 annotated neurons of 62,911 nuclei (1.1 %), 4–133 per sample; 2 of 15 samples drop from the Δ2B axis under the ≥ 10-nuclei rule (nB = 13 vs n1 = nA = 14) | the Fig 7C asterisk; the same caveat belongs to this figure's Δ2B columns |

---

## 3. Panel letters: v22 reference → v29 letter

| v22 (delivered, A–I) | v29 | note |
|---|---|---|
| A size-stratified rank ladder | **A** | now shows both bases side by side |
| B tight-pool ranks | **B** | unchanged content |
| C three-basis rank table | withdrawn | → Table S7 (section 2) |
| D human same-donor scatter | **C** | MS label corrected (section 5) |
| E mouse same-nucleus scatter | **D** | unchanged |
| F cross-channel rank profile | withdrawn | → the MS row of **E** + Table S7 |
| G MS rank table | withdrawn | → Table S7 |
| H 74 × 21 ledger | **E** | full text width; head margin 0.86 → 0.24 cm |
| I detection map | **F** | legend block moved under its own heading |
| — (new) | **G** | donor-level S1PR5, moved here out of Fig 7 |
| LEG axis key | RB + L2 | one band directly above **E**, aligned to the measured matrix edge |
| (top of sheet) | **K** | identity key band, seven accessions, one row |

Manuscript pointers that must be repointed when this version is pushed: the text still cites S7F and
S7G (both withdrawn) → those claims go to Table S7; S7C/S7D → **C**/**D**; S7H → **E**; S7I → **F**;
the donor-level sentence in the Fig 7D legend → **G** plus Table S7.

---

## 4. Typesetting and compliance record

* Sheet 17.00 cm wide (double-column cap) × 14.36 cm; Arial only (`pdf.fonttype = 42`, no Type 3).
* Two face tiers only: 8.0 pt (panel letters, drawn by the composer) and 6.0 pt (everything else);
  `font_gate.py --strict` PASS.
* Panel letters A–G, ArialMT 8.0 pt, alone, in a reserved 0.34 cm strip above their row; the leftmost
  letter of each row shares one x (0.835 cm), which is also, measured, the x of the ledger's own
  leftmost text. No title beside any letter (`no title beside letter` PASS).
* Gates: collisions **0**; text-on-data **2** (the MS callouts on their own highlighted dots, identical
  in the signed v22 sheet); ink outside page **0**; ink on page border **0**; run-on word pairs **0**;
  word ledger **sheet 167 = panels 160 + letters 7**, missing {}, phantom is exactly the 7 letters;
  clipped ink **0 of 10 panels**; detection-ink audit **CONSISTENT** (every bubble re-derived from
  `A3_GSE118257_target_stats.csv`).
* Palette asserted at render time: adjacent strips min CIE76 **68.6**, any pair **40.7**, vs the six
  semantic colours **18.6**, worst of normal / deuteranopia / protanopia; the render stops if a floor
  (50 / 35 / 17) is broken.
* Geometry measured on the composed page: band ink 1.213–13.61 × 5.58–6.08 cm; ledger matrix
  1.213–15.94 × 6.84–9.40 cm (band left x equals matrix left x to 0.000 cm; band bottom to matrix top
  0.77 cm); row 1 matrices aligned to a common axes top; bottom page margin 0.34 cm.
* Reproducibility: `Figure_S7.ipynb` is generated from `assembly/make_figS7_v30_panels.py` and was
  executed end-to-end with no error cells; its ten panel PDFs carry the same words and the same drawing
  geometry as the delivered ones, and panels A, B, E, F, G, K, RB, L2 are ink-identical to the approved
  v28 prototype (raster pixels identical, internal geometry ≤ 0.02 cm; the rank bar's gradient is
  0.012 cm wider, same ramp).

---

## 5. What changed in v29 and v30, and why

1. **Panel C's MS annotation is corrected: `2 → 21` becomes `2 → 36`.** The plotted y is the
   size-matched ATAC Δ1 series — the very series the printed ρ 0.722 is computed on — and MS is rank 36
   in it. The delivered number 21 is MS's rank in the *donor-paired* ATAC column of
   `same_donor_RNA_vs_CHROMATIN_delta1_74traits.csv`, which is not plotted. Both ranks are now computed
   from the plotted columns and asserted, and the donor-paired values move to Table S7. This was found
   while checking the audit receipt of 2026-11-18; the audit's own §4 had accepted "2 → 21 ✓" against
   the table rather than against the plotted axes.
2. **Panels C and D are drawn from data instead of being carried over as v22 PDFs.** A notebook must be
   able to produce every panel it shows. `verify_v30_ink.py` registers each re-drawn panel on its axes
   frame and compares all 74 dots, the highlighted MS dot and the tick marks against the delivered v22
   PDFs: maximum deviation **0.03 pt** (C) and **0.02 pt** (D), with the only word differences being the
   corrected rank and the two baked letters that the composer now draws.
3. **Panel letters moved out of the panels entirely.** The prototype chain baked a letter into each
   panel and redacted it afterwards; the repository script draws none, so the composer places them in
   the row strip. Row 1 is therefore aligned by measured content top (matrix / axes) rather than by
   canvas top, which also fixes a latent misalignment: with the letters baked, A's and B's matrices were
   sharing a line only by accident of the crop.
4. **GSE231585's Δ2B wording.** The audit's H3 asked for "neuron − non-neuron" on the face. v29
   spelled it out (`Δ2B = neuron − non-neuron`, `Δ2A = AST − immune`).
   What still says "neuron − glia" is **Fig 7C's axis label** and the manuscript abbreviation table
   ("neuron minus glia control axis"); the computation in all four cohorts is median(neuron) −
   median(non-neuronal) (`A2_GSE277430_delta2_mirror.py:48`, `A3_GSE118257_delta2_full74.py:54`,
   `A4_GSE301908_delta2_full74.py:135`, `A7_GSE231585_delta2_full74.py:96`). That fix belongs to the
   Fig 7 bundle and the manuscript, not here.

5. **v30 - one vocabulary across the figure family.** The band above **E** now reads
   `Δ1 = OL/MOL − OPC`, `Δ2A = AST − (MG+LYM)`, `Δ2B = NEU − non-NEU`, character for character
   the strings Fig 7C's three violin columns carry, and panel F's column heads use AST / MG / LYM / NEU.
   The page is 0.08 cm shorter than v29 (17.00 x 14.28 vs 14.36) because F's rotated labels got shorter;
   `gates_v30.py` is clean, and the only raster difference against v29 is inside the key band and F's
   label row (1.0% of pixels, x 1.22-15.93, y 5.90-13.97). `OL/MOL` rather than plain `OL` because the
   key decodes the mouse and own-multiome arms too, whose mature side is MOL. Cost: the manuscript
   abbreviation table must gain `NEU` and `LYM` (see NEEDS_ABBREVIATIONS.md).

## 6. Push prerequisites

* Table S7 must be written with the section-2 rows (three-basis table, MS profile, per-donor detail with
  the diagnosis column, the floor sweep, the no-floor 6/11) — the withdrawn panels have nowhere to live
  until it exists.
* Manuscript pointers S7F/S7G → Table S7; S7C/S7D → S7C/S7D unchanged in letter but re-worded; S7H →
  S7E; S7I → S7F; add the S7G pointer for the donor-level statement.
* Fig 7D legend: state that GSE301908's 11 donors are 8 patients + 3 controls.
* Face-number register: **done** - `assembly/register_figS7_faces.py` writes `figS7_out/figS7v30_face_numbers.csv`
  (38 entries, every value computed from the vendored table it cites) and gates both directions: registered strings
  all present on the page, and every numeric token on the page either registered or an explicitly listed tick.
* Abbreviation table additions/consistency: `pp`, `ImOlGs`, `Peri`, `Myeloid`, `PLATELET`, channel
  heads `RNA` / `ATAC` are on the face but not in the table; `Astro` / `Endo` / `CP10k` on the face vs
  `AST` / `EC` / `CP10K` in the table (see `NEEDS_ABBREVIATIONS.md`).
