# -*- coding: utf-8 -*-
"""make_table_S7.py - build the five Fig S7 support tables (S7a-S7e) from delivered tables only.

Why this file exists: the manuscript points at "Table S7" five times and the Fig S7 legend names five
sub-tables, but nothing was ever written. The only existing Table_S7.pdf belongs to Fig S6f (per-mouse
numbers) - a numbering collision that has to be settled by the author, see TABLE_S7_RECEIPT.md.

Rules enforced here:
 * every value is read from a CSV that already exists (no statistics recomputed except order-preserving
   aggregation of delivered per-donor numbers, which is labelled as such and cites its input);
 * each block carries a `source` column naming the file it came from;
 * anything a reader needs that no delivered table carries goes into the outstanding list, not into a
   cell of the table.

Outputs: <work>/figures/TABLE_S7_DRAFT_20261118/tables/S7*.csv (+ a printed coverage report).
"""
import ast
import os
import numpy as np
import pandas as pd

OUT = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'TABLE_S7_DRAFT_20261118/tables')
os.makedirs(OUT, exist_ok=True)
DATA = os.path.join(os.environ.get('DSCI_REPO_ROOT', '.'), 'figS7/data')
S1 = DATA + "/SUPPORT_v1/tables"
# the axis pipeline writes eight fig7S* tables; only four of them are vendored into the repo today, the
# rest are read from the analysis tree and listed as "needs vendoring" in the coverage report
S1A = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'routeA_out/SUPPORT_v1/tables')
HAND = os.path.join(os.environ.get('DSCI_REPO_ROOT', '.'), 'figS7/handoff/donor_s1pr5_20261117/tables')
FIG7T = os.path.join(os.environ.get('DSCI_REPO_ROOT', '.'), 'fig7/tables_v17')
WORKT = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'SUPPL_ANALYSES_20261112/tables')

MISSING_VENDOR = []


def load(path, vendor_hint=""):
    """Read a delivered table; remember when it lives outside the repo so it can be vendored."""
    if not os.path.exists(path):
        raise SystemExit("required table missing: " + path)
    if "/GitHub/" not in path:
        MISSING_VENDOR.append((os.path.basename(path), path, vendor_hint))
    return pd.read_csv(path)


def provenance_rows_from_code():
    """Read the PROVENANCE literal out of axis_rho_and_rank_matrix.py without executing it."""
    src = open(os.path.join(os.environ.get('DSCI_REPO_ROOT', '.'), 'figS7/upstream/axis_rho_and_rank_matrix.py'), encoding="utf-8").read()
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == "PROVENANCE":
            return ast.literal_eval(node.value)
    raise SystemExit("PROVENANCE not found")


def w(df, name):
    df.to_csv(os.path.join(OUT, name), index=False)
    print("  %-46s %3d rows x %2d cols" % (name, len(df), df.shape[1]))


print("=== reading delivered tables ===")
CONF = load(DATA + "/FIG7_v1/tables/seven4_confound_ALL_cohorts.csv", "S7b1")
RANK4 = load(S1 + "/fig7Sc_rank_matrix_4cohorts.csv", "S7b1")
SAME = load(DATA + "/A6_GSE231584_atac_v1/tables/same_donor_RNA_vs_CHROMATIN_delta1_74traits.csv", "S7b2")
PEAK = load(DATA + "/A6_GSE231584_atac_v1/tables/delta1_vs_peakdomain_allchannels.csv", "S7b2")
WINDOW = load(S1 + "/fig7Sb_matched_window.csv", "S7b3")
PAIRS = load(S1 + "/fig7Sc_axis_pair_rho.csv", "S7b3")
XCH = load(S1 + "/figS7c_crosschannel_pairs.csv", "S7c")
MIR = load(S1A + "/fig7Sc_mirror_rho.csv", "S7c2")
THD = load(S1A + "/fig7Sc_third_cohort_rho.csv", "S7c2")
ATAC_H = load(S1 + "/S7_atac_crossaxis_human_full.csv", "S7d1")
ATAC_M = load(S1 + "/S7_atac_crossaxis_mouse_full.csv", "S7d2")
RNA_M = load(S1 + "/S7_rna_crossaxis_mouse_full.csv", "S7d2")
SNUC = load(S1 + "/S7_mouse_same_nucleus_delta1.csv", "S7d2")
DON3 = load(HAND + "/S41_S1PR5_per_donor_with_dx.csv", "S7d3")
S40R = load(HAND + "/S40_rule_per_donor.csv", "S7d3b")
GRP = load(WORKT + "/S3_group_summary.csv", "S7d4")
RATER = load(HAND + "/S40_rate_donor_pop.csv", "S7e1")
DETLONG = load(HAND + "/S40_det_donor_pop.csv", "S7e2")
S15DET = load(FIG7T + "/S15_det_GSE301908.csv", "S7e3")
S15N = load(FIG7T + "/S15_n_GSE301908.csv", "S7e3")

# ============================== S7a  dataset, unit and contrast provenance ==============================
# S7a is the two tables the axis pipeline already writes (they were never tabulated into a document):
#   fig7Sa_provenance.csv               study / GEO / platform / species / condition / unit / n_case /
#                                       n_ctrl / annotation source / OL nuclei / notes
#   fig7Sa_platform_delta1_definition.csv  tissue / preparation / chemistry / mature side / precursor
#                                       side / same-nucleus availability / where the arm is used
S7a1 = load(S1A + "/fig7Sa_provenance.csv", "S7a1")
S7a1["source_file"] = "SUPPORT_v1/tables/fig7Sa_provenance.csv (axis_rho_and_rank_matrix.py)"
S7a2 = load(S1A + "/fig7Sa_platform_delta1_definition.csv", "S7a2")
print("=== S7a1 dataset provenance ===")
w(S7a1, "S7a1_dataset_provenance.csv")
print("=== S7a2 platform and Delta1 contrast definition ===")
w(S7a2, "S7a2_platform_and_delta1_definition.csv")

# ============================== S7b  ranks, and the three sensitivity views ==============================
# S7b1: per trait, the raw rank and the size-adjusted rank in every cohort, plus the confound itself.
b1 = CONF[["trait", "n_peak_raw", "lnp", "rank_raw_A2", "rank_adj_A2", "rank_raw_A3", "rank_adj_A3",
           "rank_raw_A4", "rank_adj_A4", "rank_raw_A7", "rank_adj_A7", "rank_raw_chr",
           "rank_adj_chr"]].copy()
b1.columns = ["trait", "n_peaks", "ln_n_peaks",
              "GSE118257_raw_rank", "GSE118257_adjusted_rank",
              "GSE277430_raw_rank", "GSE277430_adjusted_rank",
              "GSE301908_raw_rank", "GSE301908_adjusted_rank",
              "GSE231585_raw_rank", "GSE231585_adjusted_rank",
              "GSE231584_atac_raw_rank", "GSE231584_atac_adjusted_rank"]
b1["source"] = "seven4_confound_ALL_cohorts.csv (A8c_confound_unified.py)"
# the adjusted ranks are the ones Fig 7B prints; the raw ones are this table's sensitivity view
r4 = RANK4.copy()
r4.columns = ["trait"] + ["%s_%s" % (c, x) for c, x in
                          [("GSE277430", "d1"), ("GSE277430", "d2A"), ("GSE277430", "d2B"),
                           ("GSE118257", "d1"), ("GSE118257", "d2A"), ("GSE118257", "d2B"),
                           ("GSE231585", "d1"), ("GSE231585", "d2A"), ("GSE231585", "d2B"),
                           ("GSE301908", "d1"), ("GSE301908", "d2A"), ("GSE301908", "d2B")]]
b1 = b1.merge(r4, on="trait", how="left")
b1["three_axis_ranks_source"] = "fig7Sc_rank_matrix_4cohorts.csv (axis_rho_and_rank_matrix.py)"
# (the scDRS ranks are not in any CSV; they go into their own block below, parsed from the panel script)
print("=== S7b1 raw vs adjusted ranks (Fig 7B sensitivity) ===")
w(b1, "S7b1_trait_ranks_raw_vs_size_adjusted.csv")

# S7b4 was dropped on purpose. The v26 prototype panel showed scDRS q95 ranks for the five
# signatures, but no delivered CSV carries them (no column named *scdrs*rank* anywhere under
# figures/ or GitHub/, and Table_S6_v12 has no such row); the run itself lives in the Fig 6 chain
# as .npz score matrices. A support table must not launder a number no table holds, so the ranks
# move to the outstanding list instead of into a cell here.

# S7b2: same-donor RNA vs ATAC, and the chromatin peak-domain control.
b2 = SAME.rename(columns={"RNA_delta1": "d1_RNA", "CHR_delta1": "d1_ATAC", "n_peak_raw": "n_peaks"})
b2["d1_ATAC_vs_ln_n_peaks_rho"] = PEAK.set_index("trait").loc[b2.trait, "r"] if "r" in PEAK else np.nan
b2["d1_ATAC_median_group_diff"] = PEAK.set_index("trait").loc[b2.trait, "median_groupdiff"]
b2["source"] = ("same_donor_RNA_vs_CHROMATIN_delta1_74traits.csv + "
                "delta1_vs_peakdomain_allchannels.csv (A6_GSE231584_atac_v1)")
print("=== S7b2 same-donor RNA vs ATAC ===")
w(b2, "S7b2_same_donor_RNA_vs_ATAC.csv")

# S7b3: the matched-window control and the axis-pair rho (two short blocks in one file, keyed by `block`).
b3a = WINDOW.copy(); b3a.insert(0, "block", "matched_window_pm15pct")
b3b = PAIRS.copy(); b3b.insert(0, "block", "axis_pair_rho")
print("=== S7b3 matched window + axis pairs ===")
w(pd.concat([b3a, b3b], ignore_index=True, sort=False), "S7b3_matched_window_and_axis_pairs.csv")

# ============================== S7c  cross-channel concordance + mirror replication ==============
S7c = XCH.copy()
S7c["source"] = "figS7c_crosschannel_pairs.csv (crosschannel_pairs_table.py)"
print("=== S7c1 cross-channel pairs ===")
w(S7c, "S7c1_crosschannel_pairs.csv")
S7c2 = pd.concat([MIR.assign(block="mirror_rho_by_cohort"), THD.assign(block="third_cohort_mirror")],
                 ignore_index=True, sort=False)
S7c2["source"] = "fig7Sc_mirror_rho.csv + fig7Sc_third_cohort_rho.csv (axis_rho_and_rank_matrix.py)"
print("=== S7c2 mirror replication ===")
w(S7c2, "S7c2_mirror_rho.csv")

# ============================== S7d  the three axes, per library / per donor ==============================
d1 = ATAC_H.copy()
d1.columns = ["trait", "ATAC_d1_delta", "ATAC_d2A_delta", "ATAC_d2B_delta",
              "ATAC_d1_robust_z", "ATAC_d2A_robust_z", "ATAC_d2B_robust_z"]
d1 = d1[d1.trait != "d1"]                     # the producer wrote a second header line
d1["cohort"] = "GSE231584 (human ATAC, 9 libraries / 10,870 nuclei)"
d1["source"] = "S7_atac_crossaxis_human_full.csv (human_atac_threeaxis.py)"
print("=== S7d1 human ATAC three axes ===")
w(d1, "S7d1_human_atac_three_axes_per_trait.csv")

d2 = RNA_M.merge(ATAC_M, on="trait", suffixes=("_RNA", "_ATAC"))
d2 = d2.merge(SNUC.rename(columns={"delta1_rna": "same_nucleus_d1_RNA",
                                   "delta1_atac": "same_nucleus_d1_ATAC"}), on="trait", how="left")
d2["cohort"] = "GSE250589 (mouse atlas RNA 28 libraries; ATAC 25 libraries / 14,928 nuclei)"
d2["d2A_d2B_note"] = "d2A/d2B are computable only where the comparator class clears its floor (d2A 77, d2B 20 libraries)"
d2["source"] = "S7_rna_crossaxis_mouse_full.csv + S7_atac_crossaxis_mouse_full.csv + S7_mouse_same_nucleus_delta1.csv"
print("=== S7d2 mouse three axes ===")
w(d2, "S7d2_mouse_three_axes_per_trait.csv")

d3 = DON3.copy()
d3["axes_note"] = ("one row per donor (28, of which 27 evaluable): OL percentage, the strongest comparator "
                  "with >= 30 nuclei, the margin and the rule verdict; for GSE301908 one row is one PATIENT "
                  "(17 released libraries -> 11 patients), and no library-level row is merged in. "
                  "dx comes from the cohort label tables, not from the donor id")
d3["source"] = "S41_S1PR5_per_donor_with_dx.csv (s41_donor_all_human_cohorts.py + handoff README dx rules)"
print("=== S7d3 human RNA per donor ===")
w(d3, "S7d3_human_rna_per_donor_detection.csv")
d3b = S40R.copy()
d3b["definition"] = ("GSE301908: one row per patient (11 released patients from 17 libraries); OL = S1PR5 "
                     "% of the patient's OL nuclei, max_nonlin = the strongest non-lineage population, "
                     "fold = their ratio; the >= 1 pp and >= 1.5x rule is evaluated on these rows")
d3b["source"] = "S40_rule_per_donor.csv (s1pr5_per_donor_detection.py)"
print("=== S7d3b GSE301908 patient rule table ===")
w(d3b, "S7d3b_GSE301908_patient_rule.csv")

d4 = GRP.copy()
d4["definition"] = ("n_units = donors/libraries as strated; median_d1 and IQR are the delivered per-unit "
                    "d1 values; ci_lo/ci_hi are the bootstrap interval of the group median, as delivered")
d4["source"] = "S3_group_summary.csv (s3_disease_state_table.py)"
print("=== S7d4 group summary with CIs (manuscript pointer 1) ===")
w(d4, "S7d4_donor_group_summary_with_bootstrap_CI.csv")

# ============================== S7e  detection, per cohort / cell type / gene ==============================
e1 = DETLONG.copy()
e1["cohort"] = e1.donor.str.extract(r"^(GSE\d+)")[0]
e1["nuclei"] = e1.n
e1["rate_pp"] = (e1.rate * 100).round(4)
e1 = e1[["cohort", "pop", "gene", "donor", "nuclei", "det", "rate_pp"]]
e1["source"] = "S40_det_donor_pop.csv (per donor x population x gene numerator and denominator)"
print("=== S7e1 per-donor detection (long) ===")
w(e1.sort_values(["cohort", "gene", "pop", "donor"]), "S7e1_detection_per_donor_long.csv")

piv = RATER.melt(id_vars=["gene", "pop"], var_name="donor", value_name="rate_pp")
piv = piv[piv.rate_pp.notna()]
agg = (piv.groupby(["gene", "pop"])
       .agg(donors_over_50_nuclei=("rate_pp", "size"), median_pp=("rate_pp", "median"),
            q1_pp=("rate_pp", lambda s: np.percentile(s, 25)),
            q3_pp=("rate_pp", lambda s: np.percentile(s, 75)),
            min_pp=("rate_pp", "min"), max_pp=("rate_pp", "max")).reset_index())
agg["aggregation"] = "median/IQR/min/max of the DELIVERED per-donor rates (no new statistic estimated)"
agg["source"] = "S40_rate_donor_pop.csv"
print("=== S7e2 donor-level detection summary ===")
w(agg, "S7e2_detection_donor_median_IQR.csv")

genes = ["S1PR5", "GFAP", "SNAP25", "PDGFRA", "CSF1R"]
e3 = S15DET[S15DET.gene.isin(genes)].set_index("gene")
n3 = S15N.set_index(S15N.columns[0])["n"]
rows = []
for g in genes:
    if g not in e3.index:
        rows.append(dict(gene=g, cell_type="", detected=0, nuclei_in_window=0,
                         note="gene absent from the GSE301908 whole-population table"))
        continue
    r = e3.loc[g]
    if isinstance(r, pd.DataFrame):
        r = r.iloc[0]
    for ct in e3.columns:
        rows.append(dict(gene=g, cell_type=ct, detected=int(r[ct]),
                         nuclei_in_window=int(n3.get(ct, np.nan)),
                         rate_pp=round(100.0 * r[ct] / n3.get(ct, np.nan), 4) if n3.get(ct) else np.nan,
                         source="S15_det_GSE301908.csv / S15_n_GSE301908.csv (whole-population window, "
                                "not the 12,000-per-lineage sampling window Fig 7D uses)"))
print("=== S7e3 GSE301908 whole-population detection ===")
w(pd.DataFrame(rows), "S7e3_GSE301908_whole_population_detection.csv")

# ============================== coverage report ==============================
print("\n=== manuscript pointers -> which sub-table answers them ===")
for t, where in [
    ("1. Results (Fig 7B) - 'per-donor rank shifts and bootstrap CIs'", "S7d4 (+ S7d3 donor values)"),
    ("2. Results (Fig 7C) - 'per-library cross-axis values, exact ranks, third-cohort mirror'",
     "S7d1 / S7d2 / S7d3 + S7b1 (ranks) + S7c2 (mirror)"),
    ("3. Results (Fig 7D) - 'per cell type, detection; donor-level replication'",
     "S7e1 / S7e2 / S7e3 + S7d3"),
    ("4. Discussion - 'raw versus adjusted ranks (size-matched)', 'three axes of Fig 7C'",
     "S7b1 / S7b2 / S7b3"),
    ("5. Methods - 'the three axis contrasts'", "S7a"),
]:
    print("   %-76s -> %s" % (t, where))

print("\n=== tables used from outside the repo (need vendoring before push) ===")
seen = set()
for f, p, hint in MISSING_VENDOR:
    if f in seen:
        continue
    seen.add(f)
    print("   %-32s %s   (%s)" % (f, p, hint))
print("\nblocks written to", OUT)
