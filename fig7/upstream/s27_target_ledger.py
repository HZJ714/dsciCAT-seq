#!/usr/bin/env python
"""s27_target_ledger.py  --  build the merged MS drug-target ledger (Table 3 candidate).

Implements the rules of
  online/figures/FIG7D_C10_REVIEW_EXPANSION33_RECEIPT_20261117.md  sections 2.1, 2.4, 3.1-3.7, 8.5
mechanically, so that the figure row set is a FUNCTION OF A LEDGER instead of a hand-written list.

inputs
  tables/S12_curation_summary.csv          47 drug genes + 18 anchors, tier as drawn
  tables/S12_target_curation_audit_strict.csv  EuropePMC verdicts (local verification)
  tables/S17_target_evidence_flags.csv     label/MOA/MS flags for rows already audited
  tables/S23_gene_universe_status.csv      class (S13 rule) for all 102 scored atoms
  tables/S15_det_/S15_n_ GSE301908/GSE231585   detection rates for the floor rule
  <work>/figures/MS_DRUG_TARGET_EXPANSION_33.md   the 33-row expansion (status per row)

outputs
  tables/S26_target_ledger.csv   one row per scored atom: source, tier_as_drawn, tier_audited,
                                 class, flags, floor, drawn_in_main, reason
  tables/S26_rowset_C12.csv      the rows the main figure will draw, in draw order
"""
import os
import re, sys, pandas as pd, numpy as np
# Root of the author's analysis tree (figures/, routeA_out/, routeA_probe/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
T = WORK + "/figures/SUPPL_ANALYSES_20261112/tables/"
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
# Vendored ASCII extract of the expansion memo (see data/EXPANSION_33/): the memo is an
# internal Chinese-language document, and fda_status_en holds the exact English strings the
# delivered ledger carries, so this chain reproduces the delivered table instead of translating it.
EXP_CSV = os.environ.get("FIG7_EXPANSION_CSV", os.path.join(REPO, "data", "EXPANSION_33", "expansion33_rows.csv"))

# ---------------------------------------------------------------- 1. the 47-side
S12 = pd.read_csv(T + "S12_curation_summary.csv")
OLD = S12[~S12.tier_as_drawn.str.startswith("T6")].copy()          # 47 drug genes (anchors excluded)
OLD["tier47"] = OLD.tier_as_drawn.str[:2]
OLD["source"] = "curation47"
F0 = pd.read_csv(T + "S17_target_evidence_flags.csv").set_index("gene")
DRUGFIX = {"BTK": "tolebrutinib / fenebrutinib", "DHODH": "teriflunomide", "TOP2A": "mitoxantrone",
           "S1PR3": "fingolimod (off-target)", "S1PR4": "fingolimod (off-target)", "DCK": "cladribine",
           "S1PR2": "(no MS drug engages S1PR2)", "ADA": "pegademase (not used in MS)"}
OLD["drug"] = OLD.gene.map(lambda g: F0.drug.get(g, "") or DRUGFIX.get(g, ""))

# ---------------------------------------------------------------- 2. the 33-row expansion -> single genes
X = pd.read_csv(EXP_CSV)
assert len(X) == 33, len(X)
SPLIT = {"IMPDH1/2": ["IMPDH1", "IMPDH2"], "CD80/CD86": ["CD80", "CD86"], "SLC6A4/2": ["SLC6A4", "SLC6A2"],
         "SCN1A/8A": ["SCN1A", "SCN8A"], "ESR1/2": ["ESR1", "ESR2"]}                       # 3.1/3.2
exp = []
for r in X.itertuples():
    for g in SPLIT.get(r.claim_cell, [r.claim_cell]):           # this ledger's own atom split (no HTR2B)
        exp.append(dict(gene=g, exp_no=int(r.exp_no), drug=r.drug, status=r.fda_status_en))
E = pd.DataFrame(exp)
E["source"] = "expansion33"
assert len(E) == 38, len(E)                                     # V4
E["dup_of_47"] = E.gene.isin(set(OLD.gene))                     # V5  CD38, TNFSF13B
NONDRUG = {"TNFSF13": "ligand, not a drug target", "MTNR1A": "melatonin = food supplement, not a licensed drug",
           "IL7R": "genetic risk locus / pathway, no MS drug"}  # 3.3
E["nondrug"] = E.gene.map(NONDRUG)

# ---------------------------------------------------------------- 3. status -> tier, mechanical rule (3.7)
def tier_from_status(g, status):
    s = status
    if g in NONDRUG:
        return "X"
    if g in ("NR3C1", "MC2R", "IFNAR2", "RRM2B"):               # 3.7 explicit
        return "T1"
    if g == "CD40LG":                                           # 3.7 explicit
        return "T2"
    if "approved: MS" in s:                                      # incl. "approved: MS exacerbation"
        return "T1"
    if "worsened relapse" in s:                                 # 3.5 IFNGR1
        return "T5"
    if s.startswith("not approved") or "in trials" in s:
        return "T4"                                             # 3.7: "the rest" never enter T2
    if "not a drug" in s or "dietary supplement" in s or "no FDA-approved product" in s:
        return "T4"
    return "T4"                                                 # approved elsewhere / preclinical
E["tier_exp"] = [tier_from_status(g, s) for g, s in zip(E.gene, E.status)]

# ---------------------------------------------------------------- 4. merge
L = pd.concat([
    OLD[["gene", "tier47", "source", "status_claimed"]].rename(columns={"tier47": "tier_raw", "status_claimed": "status"}),
    E.loc[~E.dup_of_47, ["gene", "tier_exp", "source", "drug", "status"]].rename(columns={"tier_exp": "tier_raw"}),
], ignore_index=True)
L["drug"] = L.drug.fillna("")
assert len(L) == 83, len(L)                                    # V6

# ---------------------------------------------------------------- 5. audited tier corrections (3.8 + S12 must-fix)
FIX = {
 "BTK":   ("T2", "no BTK inhibitor is approved in MS (tolebrutinib/fenebrutinib phase III); S12 must-fix, "
                 "EuropePMC screen could not verify an approval"),
 "IL2RA": ("T5", "daclizumab withdrawn 2018; S12 strict verdict NOT VERIFIED / pass2 FAIL -> negative control, not an approved target"),
 "DCK":   ("T1", "cladribine is an approved MS drug and its label names DCK as the activating kinase "
                 "(expansion receipt 7.3, Mavenclad clinical pharmacology) -> target of an approved MS drug"),
 "TOP2A": ("T1", "kept T1 but flagged: mitoxantrone approved for SPMS only, cardiotoxic, largely abandoned"),
}
L["tier_audited"] = [FIX.get(g, (t, ""))[0] if t != "X" else "X" for g, t in zip(L.gene, L.tier_raw)]
L["tier_fix_reason"] = [FIX.get(g, ("", ""))[1] for g in L.gene]
L["tier_fix"] = L.tier_raw != L.tier_audited

# ---------------------------------------------------------------- 6. class (S13 rule) + floor rule
CLS = pd.read_csv(T + "S23_gene_universe_status.csv").set_index("gene")
L["class"] = L.gene.map(CLS.class_new)
assert L["class"].notna().all(), L[L["class"].isna()].gene.tolist()

def rates(acc):
    d = pd.read_csv(T + "S15_det_%s.csv" % acc, index_col=0).groupby(level=0).max()
    n = pd.read_csv(T + "S15_n_%s.csv" % acc, index_col=0)["n"]
    return d.div(n.reindex(d.columns).values.astype(float), axis=1) * 100
R4, R7 = rates("GSE301908"), rates("GSE231585")
L["in_A4"] = L.gene.isin(R4.index); L["in_A7"] = L.gene.isin(R7.index)
L["max_rate_A4"] = L.gene.map(R4.max(axis=1)); L["max_rate_A7"] = L.gene.map(R7.max(axis=1))
L["floor"] = [(np.nanmax([a, b]) < 0.10) for a, b in zip(L.max_rate_A4.fillna(0), L.max_rate_A7.fillna(0))]
L.loc[~(L.in_A4 | L.in_A7), "floor"] = False                   # absent everywhere != below floor

# ---------------------------------------------------------------- 7. flags for the rows that have none
F = pd.read_csv(T + "S17_target_evidence_flags.csv").set_index("gene")
NEWFLAG = {  # state: 1 documented here, 2 asserted by the parallel env's label query (receipt 7.2/7.3), 0 not documented
 "BTK":    (0, 1, 0, "tolebrutinib/fenebrutinib MOA names BTK (EuropePMC VERIFIED); no approved MS drug -> label/MS blank"),
 "DHODH":  (2, 1, 1, "Aubagio label MOA 'inhibits dihydroorotate dehydrogenase' (receipt 7.3, not re-read here); MS approval VERIFIED locally"),
 "TOP2A":  (0, 1, 2, "mitoxantrone = topo-II/DNA poison (S12 pass2 note); Novantrone label does not name TOP2A; SPMS approval not verified by the screen"),
 "S1PR3":  (0, 0, 0, "off-target of fingolimod by affinity tables (S12 VERIFIED); no label, no stated MOA, no MS use through S1PR3"),
 "S1PR4":  (0, 0, 0, "not engaged by approved modulators (S12 VERIFIED, allergy literature); kept in the space as an off-target row"),
 "DCK":    (2, 2, 1, "Mavenclad label clinical pharmacology names deoxycytidine kinase (receipt 7.3); cladribine approved in MS"),
 "RRM2B":  (2, 2, 1, "same drug as DCK; the label also names deoxyguanosine kinase in mitochondria -> ledger row, not a main-figure row"),
 "NR3C1":  (0, 2, 2, "GR agonism is the stated corticosteroid MOA and Solu-Medrol carries MS-relapse labels (receipt 7.2, 18 labels); neither re-read here"),
 "IFNAR2": (0, 2, 2, "IFNbeta signals through type-I IFN receptor; AVONEX/REBIF/PLEGRIDY MS labels (receipt 7.2); not re-read here"),
 "MC2R":   (0, 2, 2, "ACTH acts through MC2R; ACTHAR GEL MS label (receipt 7.2, 1 label); not re-read here"),
 "CD40LG": (0, 2, 0, "frexalimab = anti-CD40L, phase II/III in MS (receipt 7.2: any-labels = 0) -> no label, no approved MS use"),
}
for g, (a, b, c, note) in NEWFLAG.items():
    if g not in F.index:
        F.loc[g] = dict(drug="", f_label=a, f_moa=b, f_msuse=c, note=note, n_met=a + b + c,
                       grade_from_flags="", grade_audited="", grade_as_drawn="", flags_vs_drawn="")
L["label"] = L.gene.map(F.f_label); L["MOA"] = L.gene.map(F.f_moa); L["MS"] = L.gene.map(F.f_msuse)
L["flag_note"] = L.gene.map(F.note)
L["flag_provenance"] = np.where(L.label.notna(),
    np.where(L.gene.isin(NEWFLAG), "audited here + parallel-env label query", "audited here (S17)"), "")
# rows nobody has signed yet: mechanical proxy, clearly marked, never drawn (they are all T3/T4/T5)
STRICT = pd.read_csv(T + "S12_target_curation_audit_strict.csv").set_index("gene")
miss = L.label.isna()
L.loc[miss, "MOA"] = [int(bool(str(STRICT.pass2_verdict.get(g, "")).startswith("PASS")) or
                    "VERIFIED" in str(STRICT.verdict.get(g, ""))) for g in L.gene[miss]]
L.loc[miss, "label"] = 0
L.loc[miss, "MS"] = [int(t in ("T1",)) for t in L.tier_audited[miss]]
L.loc[miss, "flag_provenance"] = "mechanical proxy, unsigned (Table 3 only)"
assert L[["label", "MOA", "MS"]].notna().all().all(), L[L[["label", "MOA", "MS"]].isna().any(axis=1)].gene.tolist()

# ---------------------------------------------------------------- 8. the row rule (2.1) + one-drug-one-row (3.4)
L["space"] = L.tier_audited != "X"                              # 80 scored drug-target rows
L["drug_key"] = L.drug.str.lower().str.replace(r"[^a-z0-9]+", "|", regex=True).str.split("|").map(
    lambda p: sorted(set(x for x in p if x))[:2])
L["dup_drug"] = [g in ("RRM2B",) for g in L.gene]               # 3.4: cladribine's main row is DCK
L["drawn_as_drawn"] = (L.tier_raw.isin(["T1", "T2"])) & ~L.dup_drug
L["drawn_audited"] = (L.tier_audited.isin(["T1", "T2"])) & ~L.dup_drug
L["not_drawn_reason"] = np.where(L.drawn_audited, "", np.where(L.dup_drug, "same drug as DCK (3.4)",
                         np.where(L.tier_audited == "X", "non-drug atom (3.3)", "tier T3/T4/T5 -> Table 3 only")))

ORDER = {"1": 0, "2": 1, "3": 2, "6": 3, "4": 4, "5": 5}
L["cls_num"] = L["class"].str[0]
L = L.sort_values(["drawn_audited", "cls_num", "gene"], ascending=[False, True, True])
L.to_csv(T + "S26_target_ledger.csv", index=False)

R = L[L.drawn_audited].copy()
R["ord"] = R.cls_num.map(ORDER)
R = R.sort_values(["ord", "gene"])
R[["gene", "drug", "class", "tier_audited", "label", "MOA", "MS", "floor", "max_rate_A4", "max_rate_A7"]].to_csv(T + "S26_rowset_C12.csv", index=False)

print("ledger atoms %d | scored drug-target space %d | non-drug %d" % (len(L), L.space.sum(), (~L.space).sum()))
print("tier_audited counts:\n", L[L.space].tier_audited.value_counts().sort_index().to_string())
print("\nrow rule T1+T2 (audited) = %d rows   [as-drawn tiers would give %d]" % (len(R), L.drawn_as_drawn.sum()))
print("tier corrections:", [(g, a, b) for g, a, b in L[L.tier_fix][["gene", "tier_raw", "tier_audited"]].values])
print("\n", R[["gene", "class", "tier_audited", "label", "MOA", "MS", "floor", "max_rate_A4", "max_rate_A7"]].to_string(index=False))
print("\nclass tally of the %d drawn rows:" % len(R), R.cls_num.value_counts().sort_index().to_dict())
print("floor rows:", R[R.floor].gene.tolist())
fr = []
for g in R.gene:
    if g in R4.index:
        s = R4.loc[g]; best = s[["AST", "MG", "LYM", "VASC", "NEU"]].max()
        if s["OL"] > best + 1 and s["OL"] > 1.5 * best:
            fr.append(g)
print("framed rows (A4 rule OL>best non-lineage +1pp and >1.5x):", fr)
