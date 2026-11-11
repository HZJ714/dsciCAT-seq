# Figure 7 — journal legend (FINAL, v17 two-column sheet, 2026-11-17)

Paste-ready. Produced by `Figure_7.ipynb` (panels) and assembled outside it; the sheet is
**17.00 × 14.34 cm**, A | B | C on row 1 and D on row 2. Face typography is Arial at two sizes only
(8.0 / 6.0 pt); panel letters are 8.0 pt regular.

---

**Figure 7 | The mature-oligodendrocyte-versus-precursor risk contrast replicates in independent human
MS cohorts and in mouse EAE, and is specific to the oligodendrocyte-maturation axis.**
**(A)** Unit-level Δ1 = median scDRS score of mature oligodendrocytes − median score of OPC/COP
precursors, computed within each donor (GSE277430, GSE118257, GSE231585), within each library
(GSE301908) and within each mouse library (GSE250589, spinal cord, FACS GFP+ from Sox10-cre animals,
lesion phase as labelled). Dots are units, bars are group medians, colours are MS (red), EAE (orange)
and control (dark blue). The number above each column is the Mann-Whitney rank-sum *P* for that column
(MS or EAE versus control); the "P =" prefix is omitted on the face to keep the 0.75 cm column pitch
clear. Values: GSE277430 0.500, GSE118257 1.000, GSE301908 0.091, GSE231585 0.393; mouse Early 0.100,
Peak 0.333, Late 0.857.
**(B)** The same contrast ranked across 74 GWAS gene sets (MAGMA, 10 kb windows, 1,000 genes each) in
the four human cohorts. Rows are the 74 gene sets ordered top-down by the mean of their four
within-cohort ranks; columns are cohorts; colour is the rank of cohort-level Δ1 within that cohort
(1 = strongest mature-OL contrast, 74 = weakest). Because Δ1 tracks the size of each gene set's
regulatory domain in every cohort (Spearman ρ of Δ1 with ln n_peak = −0.84 / −0.69 / −0.78 / −0.87;
n_peak counted on public human snATAC GSE231584, GRCh38, 121,505 deduplicated peaks, TSS ± 100 kb),
the ranks shown are those obtained after residualising Δ1 on ln n_peak, so that a gene set is not
ranked ahead of another merely for covering more chromatin. On this basis MS ranks 1 / 1 / 3 / 1,
Alzheimer's disease 18 / 49 / 27 / 12 and schizophrenia 68 / 74 / 42 / 71 (GSE277430 / GSE118257 /
GSE301908 / GSE231585). The unadjusted ranking (MS 2 / 4 / 7 / 2), the ±15 % matched-window control
(MS 2 / 3 / 5 / 2) and the ρ strip are reported in Fig S7. Ranks are always computed within a cohort
and never from the transformed values in (C).
**(C)** The same nuclei re-projected onto two further cell-identity axes, printed with the cell-type
symbols the figure's own panel D uses: `AST−(MG+LYM)` (median astrocyte − median myeloid +
lymphocyte) and `NEU−non-NEU` (median neuron − median of every non-neuronal nucleus: astrocyte, myeloid,
lymphocyte, oligodendrocyte, OPC and vascular), beside the OL-maturation axis `OL−OPC` of (A) and (B).
The old face text "neuron−glia" was wrong on the comparator side, and "astro−immune" spelled out
what panel D had already abbreviated. Values are
robust *z* of Δ, standardised within cohort × axis by median and MAD, so the four tiles share one
y-axis; violins are kernel densities of the 74 gene-set values in that cohort × axis, evaluated only
inside the visible window (truncated at −4.7 / +5.5), while ranks and IQR boxes use all 74. Red dots
are MS and dark-blue dots schizophrenia, each labelled `#n` with its rank on the axis it sits on.
The single `*` marks **that axis, not that cohort**: GSE231585 contains 62,911 nuclei of which 717
(1.1 %) are annotated neurons, 4-133 per sample, so under the frozen rule that a group must contribute
at least 10 nuclei to a per-unit median, two of its 15 samples drop out of the neuron−glia axis
(13 units there, 14 on the other two axes). Nothing is subsampled: 717 is the annotated neuron pool of
the cohort. The other cohorts' neuron pools are 7-100× larger (GSE277430 76,436; GSE118257 5,031;
GSE301908 ≤ 12,000 per lineage by the stratified draw), so no other tile carries a mark.
**(D)** The 18 gold-standard drug targets that survived audit out of 80 scored (83 audited atoms, 3 of
them not drugs). Left block: within-cohort log2 odds ratio of detection in mature OL versus OPC, over all
four cohorts; middle block: OL lineage versus all non-lineage cells, over the three full atlases only —
GSE277430 is an OL-subset release and carries no non-lineage cells, so that contrast is undefined there
and no column is drawn for it. A hatched cell is the key's `not in matrix`: the gene is absent from that
cohort's tested universe; a pale cell is one whose within-cohort contrast did not reach q < 0.05. Right
block: the percentage of nuclei in which the gene is detected in GSE301908, shown separately for each cell
population — a within-population detection rate, not a composition. The red box marks the single row that
passes the pre-set rule (mature OL exceeds every non-lineage population by
>= 1.0 percentage point and by >= 1.5-fold): S1PR5. The same rule re-applied inside individual donors is
a robustness check and is reported in Fig S7, not here. The three drug-target columns (label / stated
mechanism of action / used in MS through this target) describe the drug-target relation itself and are
not localisation evidence.

Statistical unit is the donor (library where a cohort is library-resolved), never the nucleus;
n per unit, cohort sources, platforms and annotation provenance are in Table S6 and Fig S7A.
scDRS 1.0.2 (`score_cell`, n_ctrl = 200, fixed seeds per cohort); MATLAB/Python environment details in
the repository README of this figure.

---

## Notes for the production editor

* No panel is rescaled: A 6.58 × 5.90 cm, B 3.78 × 5.85 cm, C 5.84 × 5.90 cm, D 17.00 × 7.90 cm,
  placed as A | B | C (row 1) + D (row 2) inside 17.00 × 14.34 cm.
* Fonts: Arial (delivered as metrically identical Liberation Sans; swap ArialMT → Arial in one step).
  Two sizes only, 8.0 and 6.0 pt; letters A-D are 8.0 pt regular, not bold.
* Nothing in the figure uses colour as the only carrier of meaning: MS/EAE/control are also separated
  by column, the three axes in (C) are also labelled, and (B) names its three traits with leaders.

## Changes against the previously finalised single-column version (v16)

| item | v16 (single column, 7.93 cm) | v17 (two column, 17.00 cm) |
|---|---|---|
| layout | A / B / C stacked | A \| B \| C on row 1, D added on row 2 |
| B | landscape band, 74 columns × 4 rows, unadjusted ranks | portrait heatmap, 74 rows × 4 columns, ranks after size residualisation (MS 1 / 1 / 3 / 1) |
| A column heads | "GSE / 277430" on two lines | full accession at 45° |
| A *P* values | "P=0.091" | "0.091" (prefix moved to this legend) |
| sizes | 5.5-8 pt mixed | 8.0 / 6.0 pt only |
| D | absent | 18 targets × three blocks |
| the `*` | printed after the tile title, worded as "thinned to 717" | printed above the neuron−glia violin it qualifies, worded as the cohort's annotated neuron pool |
| per-donor re-test of D's detection rule | quoted in this legend | removed from Figure 7 entirely and handed to Fig S7 |
