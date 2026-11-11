# Handoff to Fig S7 — donor-level re-test of the Fig 7D detection rule (S1PR5), 2026-11-17

**Why it is here and not in Figure 7.** The author's decision (2026-11-17): a donor-level robustness
test of a Panel D claim belongs to the supplement whether it comes out positive or negative, so it was
removed from `GitHub/fig7/` (notebook section 5 deleted, scripts and tables moved here). Fig 7D keeps
the pooled detection heatmap and its red box exactly as delivered; nothing on the Fig 7 face changed.

## What was computed

| script | what it does | tables it writes (into `tables/`) |
|---|---|---|
| `s1pr5_per_donor_detection.py` | re-counts detected / total nuclei per **gene x population x donor** straight off the released GSE301908 matrix (`routeA_probe/GSE301908_counts/full_csr.bin`; 293,970 nuclei, 555,975,985 nonzero, 14 chunks, 22 genes = the 18 audited targets + PLP1/PDGFRA/GAPDH/ACTB). Read-only, ~5 min. | `S40_det_donor_pop.csv`, `S40_rate_donor_pop.csv`, `S40_summary.csv`, `S40_rule_per_donor.csv` |
| `s1pr5_per_donor_summary.py` | applies the Panel D rule (`ol > max_nonlin + 1.0` AND `ol > 1.5 x max_nonlin`, verbatim from the Fig 7D panel code) **inside each donor**, with a pool-size floor sweep (30/50/100 nuclei per donor x population) and a per-donor one-sided Fisher test against the pooled non-lineage count | `S40b_summary_floors.csv`, `S40b_S1PR5_per_donor.csv` |

Run with `DSCI_WORK_ROOT=<work>` (or `FIG7_PROBE` pointing at `GSE301908_counts/`); outputs go to
`tables/` unless `S7_DONOR_OUT` says otherwise.

## Result (GSE301908, 17 libraries = 11 donors)

| statement | value |
|---|---|
| OL has the highest S1PR5 detection rate among the five non-lineage populations | **11 / 11 donors**, one-sided sign test p = 4.9e-4 |
| full red-box margin (>= 1.0 pp AND >= 1.5x) | **6 / 11 donors**; 4 / 7 among donors with >= 30 nuclei in every comparator, 3 / 5 at >= 50, 2 / 3 at >= 100 |
| per-donor median rates | OL 8.32 % vs strongest non-lineage 3.73 %, median fold 2.5x |
| the five misses | ms10 (LYM n=9), ms7 (NEU n=15), ms8 (LYM n=31) = tiny comparator pools; ms11 (OL 3.17 % vs VASC 2.68 %, n=410) and ms15 (OL 5.20 % vs NEU 4.22 %, n=237) = genuine low-margin donors |
| the other 17 targets | none passes the rule in any evaluable donor; PLP1 is "OL highest" in 6 / 7 evaluable donors but with fold 1.07, i.e. it is near-saturated in every population |
| pooled rule, as drawn in Fig 7D | only S1PR5 passes — unchanged |

Two unit corrections that matter for any text written from this:

1. **11 donors, not 17 libraries.** `library_group_map.csv` shows ms1/ms4/ms6/ms12/ms15/ms16 contribute
   two libraries each. The three libraries GEO labels "Ctrl" (ms6r1, ms12r1, ms16r1) are second
   libraries of donors that also appear as MS: "Ctrl" there is a second site of the same subject.
   The earlier wording in the S8 bubble draft ("17/17 libraries, sign test p = 7.6e-06") counts technical
   replicates as independent units; quote 11/11 and p = 4.9e-4 instead.
2. **717 neurons in GSE231585 is not a thinning.** Unrelated to this test but noticed while auditing:
   the Fig 7 `*` caveat is the cohort's annotated neuron count; `A7_GSE231585_delta2_full74.py` applies
   no cap (unlike A4's <= 12,000 per lineage). Fig 7's README/caption were corrected accordingly.

## What this bundle still owes before the panel can be drawn

The relayout receipt (`figures/FIGS7_RELAYOUT_FOR_FIG7V17_20261117.md`, section 4.3) asks for a
**four-cohort** per-donor view (x = OL minus non-OL detection rate, one column per cohort).

* A4 / GSE301908: done, here (293,970 nuclei, 7 populations, smallest population 2,123 nuclei).
* A3 / GSE118257: **feasible with the identical statistic** - `s15b` already reads its raw matrix
  (`routeA_probe/GSE118257_matrix.txt.gz`) and the released annotation table carries a per-nucleus
  `Sample` column, so the donor split is a re-aggregation of the same pass. Pool sizes per donor are
  small though: the cohort has 17,799 OL-subset nuclei over 8 populations, with OPC n = 352 and
  LYM n = 423 overall, so a per-donor OL-vs-non-lineage margin will be noisy by construction.
* A7 / GSE231585: feasible (raw counts in `derived/`, per-nucleus `sample` column, 15 donor x lesion
  units) but the same population-depth problem: VASC n = 299, LYM n = 395.
* A2 / GSE277430: **cannot carry this rule.** The released OL subset contains only OL and OPC labels
  (35,836 / 8,931 nuclei), i.e. there is no non-lineage comparator, and its stored values are
  `normCount`, not raw counts, so "detected = stored value > 0" is a different statistic there. Leave
  A2 out of the panel and say so in one line.

So the four-cohort version of the receipt's section 4.3 is not achievable with one consistent statistic:
the honest panel is A3 + A4 + A7, or A4 alone at donor resolution (11 donors) with the other two cohorts
shown at their pooled depth. Suggest deciding between those two before drawing anything.

Disciplines from the receipt that this data satisfies: only detection rates are used (no absolute
expression across cohorts), every donor row carries its pool sizes, and the wording stays
"OL-biased detection replicates across donors within cohorts" - no cross-species label transfer, no
"pharmacological engagement".

## Diagnosis per donor (added 2026-11-18, asked by the author)

`tables/S41_S1PR5_per_donor_with_dx.csv` = the same 28 rows plus `dx` (MS / Ctrl) and `dx_source`.
Nothing else was touched; the original `S41_S1PR5_per_donor.csv` is unchanged.

**The donor ids in GSE301908 are not MS status.** `ms1 … ms16` are subject labels; the case/control
field is `patient_info` in `routeA_probe/GSE301908_counts/full_meta.tsv`, and it says **8 patient /
3 ctrl** (ms6, ms12, ms16 are controls). Anyone reading "11 ms donors" as "11 MS patients" would be
wrong, and Fig 7D's pooled GSE301908 block mixes both. GSE118257 uses `Condition` from the released
annotation table (4 MS / 5 Ctrl, `SD48/16` is Ctrl); GSE231585 uses the sample prefix (5 MS / 3 Ctrl).
Of the 27 evaluable donors: **16 MS, 11 Ctrl**.

Does the marker's margin differ by diagnosis? It is not tested anywhere in the manuscript, and on this
data it does not separate - which is exactly why panel G does not colour by it:

| cohort | MS median (pp) | Ctrl median (pp) | Mann-Whitney P | rule met |
|---|---|---|---|---|
| GSE118257 | 4.69 (n=3) | 7.48 (n=5) | 0.571 | 3/3 vs 4/5 |
| GSE231585 | 6.87 (n=5) | 6.22 (n=3) | 1.000 | 4/5 vs 3/3 |
| GSE301908 | 4.31 (n=8) | 6.03 (n=3) | 0.376 | 5/8 vs 3/3 |
| pooled (cohort-confounded) | 5.08 (n=16) | 6.22 (n=11) | 0.336 | 12/16 vs 10/11 |

Two of three cohorts point the *other* way (controls with the larger OL-vs-rest margin), none reaches
significance, and the smallest n is 3. The panel's y value is a within-donor contrast, so diagnosis is
a covariate of it, not the quantity plotted; drawing it would put a null (with a direction) at the
centre of a panel whose only job is "the pooled Fig 7D rule holds donor by donor".
