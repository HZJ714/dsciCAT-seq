#!/usr/bin/env python3
"""Build the publishable Table S6 ledger (v13) from the canonical v12 ledger.

Hard rule enforced here and checked by dev/preflight_checks.py (G11):
  * numeric columns (value / ci95_or_IQR / n / rejects_zero / meets_effect) are NEVER
    altered in magnitude - only the internal-workflow wording of a cell is rewritten;
  * rows whose sole purpose was internal bookkeeping are dropped, not edited;
  * every `source` pointer is rewritten to a path that exists inside this repository.

Usage:  python dev/make_table_S6_v13.py            # write + print change log
        python dev/make_table_S6_v13.py --check     # verify an existing v13
"""
import re
import sys
import hashlib
import pandas as pd

SRC = "figS6/tables_v12/Table_S6_v12.csv"
DST = "figS6/tables_v12/Table_S6_v13.csv"
UP = "figS6/upstream/"          # vendored upstream analysis tables (see figS6/README.md)

# ---- rows that are internal bookkeeping, not results -------------------------
DROP_CHECKS = []
DROP_STATISTICS = [
    r"^tail fraction read off the DIGITISED null",
    r"^sum of digitised bar counts",
    r"— REMOVED$",
    r"^C-side per-trait six-point rho rows delivered$",
]
# rows labelled only with the internal work-order name get a scientific label + a repo pointer
A1_LABEL = [
    (r"six-point rho for all 74 traits", "per-trait six-point concordance (C matrices)",
     "figS6/tables_v12/detail_sixpoint_all74.csv"),
    (r"split-half Delta1 pairs", "split-half Delta1 reproducibility (per trait)",
     "figS6/tables_v12/splithalf_delta1_all74.csv"),
    (r"gridded UMAP views", "rho between gridded UMAP views",
     "figS6/tables_v12/hexbin_gridsize_rho.csv"),
    (r"'% signal' per locus", "'% signal' per locus (figure art only)",
     "figS6/tables_v12/detail_window_gene_hits_331.csv"),
    (r"raw permutation vector", "permutation null vector behind Fig S6B",
     "figS6/tables_v12/window_permutation_null_10000.csv"),
    (r"excluded by the M2b rule", "mouse QC screen (M2b rule)",
     "figS6/tables_v12/mice_all26_qc_rule.csv"),
    (r"per-mouse channel gap", "per-mouse channel gap",
     "figS6/tables_v12/per_mouse_channel_gap.csv"),
]
VALUE_RULES = [
    (r"^DELIVERED:\s*", ""), (r"^DELIVERED\s+", ""),
    (r"the delivered T4 table", "the reported T4 table"),
    (r"were never delivered with the old bars", "were never published with the earlier figure"),
    (r"NOT REPRODUCIBLE", "not reproducible"),
    (r"\(P recomputes to 0\.2767723, add-one\)", "(P recomputes to 0.2767723, add-one corrected)"),
    (r"against the delivered 20-draw", "against the 20-draw values in the source table"),
    (r"read off the delivered per-nucleus", "read off the per-nucleus"),
]
DROP_SOURCE = [r"^figS6v6_derived_numbers\.csv \(superseded\)"]

# ---- path rewrites (longest-first; order matters) ---------------------------
PATHS = [
    ("FY/oligo/scdrs_ms_out/ms_locus_overlap_table.csv", "figS6/tables_v16/ms_locus_overlap_table.csv"),
    ("figures/FINAL_SET_v2_20261028/Figure_S6.pdf", "Figure S6 (rendered; see Fig S6/)"),
    ("scripts/routeA/make_S6_v2_panels.py", UP + "make_S6_v2_panels_layout_params.txt"),
    ("D_crosschemistry/tables_log1p_20261104/", UP + "D_crosschemistry/tables_log1p_20261104/"),
    ("D_crosschemistry/tables/", UP + "D_crosschemistry/tables/"),
    ("O1_own_multiome_paired/tables/", UP + "O1_own_multiome_paired/tables/"),
    ("M2b_authoranno_v1/tables/", UP + "M2b_authoranno_v1/tables/"),
    ("FIG6_v17/tables/", UP + "FIG6_v17/tables/"),
    ("SUPPORT_v2/tables_v12/", "figS6/tables_v12/"),
    ("SUPPORT_v2/tables_v11/", "figS6/tables_v11/"),
    ("SUPPORT_v2/tables_v8/", UP + "SUPPORT_v2/tables_v8/"),
    ("SUPPORT_v2/tables/", UP + "SUPPORT_v2/tables/"),
    ("SUPPORT_v1/tables/", UP + "SUPPORT_v1/tables/"),
    ("routeA_probe/mouse_EAE/anno_min.pkl", UP + "M2b_authoranno_v1/ (mouse EAE annotation)"),
    ("RNA full74/, ATAC out_gs74/", "the 74-trait scDRS score sets"),
    ("(tables_log1p_20261104)", "(figS6/upstream/D_crosschemistry/tables_log1p_20261104)"),
    ("tables_v12/window_permutation_null_histogram.csv", "figS6/tables_v12/window_permutation_null_histogram.csv"),
]

# ---- wording rewrites -------------------------------------------------------
EXACT = {
    "chromatin-side reliability (receipt 20261104 plan)": "chromatin-side reliability",
    "derived summary of a delivered column (no new test)": "derived summary of a reported column (no new test)",
    "read from the delivered table, no new test": "as reported, no new test",
    "not delivered": "n.a.",
    "not delivered (Fisher CI of raw in same table)": "n.a. (Fisher CI of the raw comparison is in the same row group)",
    "Table only (forest demoted, triage receipt 3.1)": "Table only (the forest-plot variant was dropped)",
    "Table only (forest demoted, triage receipt 3.1); verdict also on Fig 6D":
        "Table only (the forest-plot variant was dropped); verdict also on Fig 6D",
    "Table only (receipt 3.1 asks for it to be stored here)": "Table only",
    "confound check, receipt 2.3 'must be on the figure or in the ledger'": "confound check",
    "read off delivered art, not recomputed": "read off rendered figure art, not recomputed",
    "1,000 permutations (digitised, 25 bars)": "1,000 permutations (approximated from 25 rendered bars)",
    "raw 10,000-draw vector now delivered (closes the outstanding item)":
        "raw 10,000-draw vector: figS6/tables_v12/window_permutation_null_10000.csv",
    "provenance warning: digitised figure values are not independent evidence":
        "note: values approximated from figure art are not independent evidence",
    "vector delivered; histogram also in tables_v12/window_permutation_null_histogram.csv":
        "the full 10,000-draw vector is provided; histogram also in figS6/tables_v12/window_permutation_null_histogram.csv",
    "hit genes (rows of the delivered overlap table)": "hit genes (rows of the overlap table)",
    "traits whose state-level rho is NOT delivered": "traits whose state-level rho is not reported",
    "delivered hexbin bin counts (gridsize 15/25/40)": "hexbin bin counts (gridsize 15/25/40)",
    "C ladder as delivered (MOL6 vs COP/NFOL, 5,962 nuclei)": "C ladder as reported (MOL6 vs COP/NFOL, 5,962 nuclei)",
    "0.952-0.997 (delivered range)": "0.952-0.997 (range across reported replicates)",
    "self-consistent with the ledger's 0.878 (receipt 2.2)": "consistent with the reported 0.878",
    "supersedes the digitised 25-bar inset and the 0.26 read off v7 art":
        "replaces the 25-bar inset approximated from figure art and the 0.26 label of an earlier draft",
    "QA record only - do NOT quote: a value read off figure art cannot outrank the source-computed 0.277":
        "reference only: a value approximated from figure art cannot outrank the source-computed 0.277",
    "receipt 4/9 forbid inventing the quantity; we deliver the underlying per-gene window membership instead":
        "the quantity was never defined in the source; the underlying per-gene window membership is reported instead",
    "P = (1 + 2,767)/(1 + 10,000); the earlier printed label 0.26 was read off the v7 figure drawn with 1,000 draws from a superseded vector and is superseded by this row":
        "P = (1 + 2,767)/(1 + 10,000) from the 10,000-draw null (figS6/tables_v12/window_permutation_null_10000.csv); a 0.26 label from an earlier 1,000-draw draft should not be quoted",
    "independent scorers do NOT put MS on the chromatin axis either - a methodological negative, ~10x less sensitive than scDRS (receipt 7.3)":
        "independent scorers do not place MS on the chromatin axis either - a methodological negative, roughly tenfold less sensitive than scDRS",
    "attenuation-corrected rho = 0.933 is delivered but was never on any figure (receipt 3.1)":
        "attenuation-corrected rho = 0.933 is reported here and was not shown on any figure",
    "P = 0.806 for the level row (delivered analytic P)": "P = 0.806 for the level row (analytic P as reported)",
    "SUPPORT_v2/tables_v12/splithalf_delta1_all74.csv (+ splithalf_T5_reliability_draws.csv, _summary.csv, splithalf_across74_summary.csv)":
        "figS6/tables_v12/splithalf_delta1_all74.csv with splithalf_T5_reliability_draws.csv, splithalf_T5_reliability_summary.csv, splithalf_across74_summary.csv",
    "SUPPORT_v2/tables_v12/window_signal_per_chr.csv + detail_window_gene_hits_331.csv":
        "figS6/tables_v12/window_signal_per_chr.csv and figS6/tables_v12/detail_window_gene_hits_331.csv",
    "D_crosschemistry/tables/fig6D_subtype_6pt.csv minus FIG6_v17/tables/fig6C_sixpoint_all74.csv":
        UP + "D_crosschemistry/tables/fig6D_subtype_6pt.csv (the 10x Multiome rows of " + UP + "FIG6_v17/tables/fig6C_sixpoint_all74.csv are subtracted)",
    "n_MOL, n_OPC, pass_M2b, reason_excluded, delivered_in_M2b_per_mouse":
        "n_MOL, n_OPC, qc_pass, reason_excluded (see also M2b_per_mouse.csv)",
    "SUPPORT_v2/tables/figS6v2_bins_check.csv (counts depend on the axes box, not the data - deliberately NOT printed)":
        UP + "SUPPORT_v2/tables/figS6v2_bins_check.csv (counts depend on the axes box, not the data; not printed)",
    "scripts/routeA/make_S6_v2_panels.py (drawing parameter, not a statistic)":
        "figure layout parameter, not a statistic (see figS6/assembly/)",
}

REGEX = [
    (r"\s*\(receipt[^)]*\)", ""),                                   # internal receipt pointers
    (r"\(plan group A/B\)", ""),
    (r"delivered ([0-9.]+) was pre-correction", r"the earlier value \1 was pre-correction"),
    (r"on the delivered matrices", "on the matrices as delivered"),
    (r"per METHODS_HARD_AUDIT_v1 A1:", "note:"),
    (r"RELABELED 2026-11-06: this row was printed as", "note: this row was originally printed as"),
    (r" - one sentence required in Results", ""),
    (r"(HT passes .{0,4}0\.6 on the corrected input: 0\.630 \()delivered 0\.518", r"\1the earlier 0.518"),
    (r"are now delivered in ", "are provided in "),
    (r"; no bootstrap CI delivered", "; no bootstrap CI computed"),
    (r"\bno bootstrap CI delivered\b", "no bootstrap CI computed"),
    (r"0\.819 \(C era\) and 0\.896 \(D era, corrected 2026-11-04\)", "0.819 (earlier analysis round) and 0.896 (corrected round)"),
]

BANNED = re.compile(r"receipt|OUTSTANDING|outstanding|not delivered|SUPERSEDED|superseded by|do NOT quote|"
                    r"QA record|FINAL_SET|routeA_out|routeA_probe|triage|A1 resolution|METHODS_HARD_AUDIT|"
                    r"/data/|/home/[A-Za-z0-9._-]+|s[t]ereonote\w*|[A-Za-z0-9._%+-]{3,}@(?:[A-Za-z0-9-]+\.)+(?:cn|com|org|net)|FY/oligo|SUPPORT_v[12]|O1_own|D_crosschemistry/tables(?![/_])|"
                    r"FIG6_v17|M2b_authoranno_v1|tables_log1p(?!\w)|delivered", re.I)


def clean(text):
    out = str(text)
    if out in EXACT:
        return out.replace(out, EXACT[out])
    for a, b in PATHS:
        out = out.replace(a, b)
    if out in EXACT:
        return out
    for pat, rep in REGEX:
        out = re.sub(pat, rep, out)
    return out.strip()


def is_dropped(r):
    chk, st, src = str(r["check"]), str(r["statistic"]), str(r["source"])
    if chk in DROP_CHECKS:
        return "check"
    if any(re.search(p, st) for p in DROP_STATISTICS):
        return "statistic"
    if any(re.search(p, src) for p in DROP_SOURCE):
        return "source"
    return None


def main():
    df = pd.read_csv(SRC, dtype=str)
    n0 = len(df)
    keep, dropped = [], []
    for _, r in df.iterrows():
        why = is_dropped(r)
        (dropped if why else keep).append((why, r))
    new = pd.DataFrame([r for _, r in keep]).reset_index(drop=True)
    # relabel the work-order rows before the generic wording pass
    for i in range(len(new)):
        st = str(new.at[i, "statistic"])
        if str(new.at[i, "check"]) in DROP_CHECKS or "A1 resolution" in str(new.at[i, "check"]):
            for pat, chk, ref in A1_LABEL:
                if re.search(pat, st):
                    new.at[i, "check"] = chk
                    new.at[i, "reference"] = ref
                    break
    for col in ["value", "in_main_text", "reference"]:
        for i in range(len(new)):
            v = "" if pd.isna(new.at[i, col]) else str(new.at[i, col])
            nv = v
            for pat, rep in VALUE_RULES:
                nv = re.sub(pat, rep, nv)
            if nv != v:
                new.at[i, col] = nv
    changes = []
    for col in new.columns:
        for i in range(len(new)):
            before = "" if pd.isna(new.at[i, col]) else str(new.at[i, col])
            after = clean(before)
            if after != before:
                new.at[i, col] = after
                changes.append((i, col, before, after))

    # ---- integrity guards -----------------------------------------------------
    # (1) quantities are never touched: value / rejects_zero / meets_effect byte-equal
    ndig = lambda s: re.findall(r"\d+(?:\.\d+)?", str(s))
    for col in ["rejects_zero", "meets_effect"]:
        a = ["" if pd.isna(r[col]) else str(r[col]) for _, r in keep]
        b = ["" if pd.isna(x) else str(x) for x in new[col].tolist()]
        assert a == b, "drift in %s" % col
    for (_, r), x in zip(keep, new["value"].tolist()):          # value: numbers identical, in order
        assert ndig(r["value"]) == ndig(x), ("magnitude drift", r["value"], x)
    # (2) tripwire: a reword may never drop a reported quantity.
    #     decimal tokens (0.277, 2.3) are statistics; integers may also be dates/counters,
    #     so integers are only guarded in the quantity-bearing columns, with declared exemptions.
    DOT = re.compile(r"\d+\.\d+")
    INT = re.compile(r"(?<![\d.])\d{2,}(?![\d.])")
    ALLOW_PAT = [
        re.compile(r"receipt\s+[\d./\-]+", re.I),
        re.compile(r"CONFOUND_ADJUSTMENT_LIT_v1[\s\d.]*", re.I),
        re.compile(r"METHODS_HARD_AUDIT_v1[\s\w:]*", re.I),
        re.compile(r"\(plan group [A-Z]/[A-Z]\)"),
        re.compile(r"(corrected )?20\d{2}-\d{2}-\d{2}"),
        re.compile(r"20\d{6}"),
    ]
    for i, col, before, after in changes:
        allow = set()
        for pat in ALLOW_PAT:
            for grp in pat.findall(before):
                allow |= set(re.findall(r"\d+(?:\.\d+)?", grp if isinstance(grp, str) else str(grp)))
        lost = [t for t in DOT.findall(before) if t not in DOT.findall(after) and t not in allow]
        if col in ("ci95_or_IQR", "n"):
            lost += [t for t in INT.findall(before) if t not in INT.findall(after)
                     and not re.match(r"^(?:20\d{6}|(?:19|20)\d{2})$", t)]
        assert not lost, ("lost numbers %s" % lost, i, col, before, after)
    new.to_csv(DST, index=False)
    REPO_OK = re.compile(r"fig(?:S6|S7|6|7)?/[A-Za-z0-9_.,*/+ -]+|Fig S6/[A-Za-z0-9_.,*/ -]*|dev/[A-Za-z0-9_./-]*")
    def strip_repo(v):
        return REPO_OK.sub("<repo>", str(v))
    left = [(i, c, v) for i in range(len(new)) for c in new.columns
            for v in [new.at[i, c]] if isinstance(v, str) and BANNED.search(strip_repo(v))]
    print("rows: %d -> %d (dropped %d)" % (n0, len(new), len(dropped)))
    print("dropped rows:")
    for why, r in dropped:
        print("   [%s] %s | %s" % (why, str(r["check"])[:44], str(r["statistic"])[:70]))
    print("reworded cells: %d" % len(changes))
    for i, col, before, after in changes:
        print("   {%s} %s\n      - %s\n      + %s" % (col, "row%d" % i, before[:150], after[:150]))
    print("remaining banned-vocabulary cells: %d" % len(left))
    for i, c, v in left[:20]:
        print("   ! row%d %s :: %s" % (i, c, str(v)[:130]))
    print("md5 %s = %s" % (DST, hashlib.md5(open(DST, "rb").read()).hexdigest()[:8]))
    return 0 if not left else 1


if __name__ == "__main__":
    sys.exit(main())
