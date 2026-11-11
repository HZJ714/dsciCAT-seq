# -*- coding: utf-8 -*-
"""按 S6_TRIAGE_RECEIPT_20261103.md §3.1/§4/§5 生成 **Table S6「相关性台账」+ Table S7 逐鼠表**（CSV + 排字 PDF）。

回执 §3.1 的列定义（每行一个统计量，不重算，直接摊交付表）：
  检查项 / 统计量 / 值 / 95%CI / 参照（天花板 0.878、判据 0.60）/ 拒零? / 达效应量? / 来源表 / 进正文?
"拒零?" 与 "达效应量?" 不是新统计量：前者 = 交付表里的置换 P 是否 < 0.05，后者 = ρ 是否 >= 0.60，
两条判据都写在交付表/预注册段里，这里只是把它们打成勾。

输出：SUPPORT_v2/tables/Table_S6.csv / Table_S7.csv 与同名 PDF（A4 横向，Arial 6.2 pt，活字）。
"""
import os, sys

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from support_common import *          # R, read, fmt_p, say
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import fitz
from arial_style import patch_fitz_arial, FONTS_PDF
patch_fitz_arial(fitz)

SUP2 = R + "/SUPPORT_v2"
TAB2, PDF2 = SUP2 + "/tables_v12", SUP2 + "/pdf_v12"   # v12 outputs go to a new dir; the delivered v11 tables/pdf are left untouched
os.makedirs(TAB2, exist_ok=True); os.makedirs(PDF2, exist_ok=True)
PT = 72 / 2.54
CEIL, CRIT = 0.878, 0.60
say("Table S6 台账 / Table S7 生成开始")

# ---------------- 交付表（只读）----------------
T4b = read("O1_own_multiome_paired/tables/T4b_cross_channel_rho_variants.csv")
T5 = read("O1_own_multiome_paired/tables/T5_split_half_reliability.csv").set_index("metric")
T6 = read("O1_own_multiome_paired/tables/T6_cross_channel_summary.csv").set_index("item")
T2 = read("O1_own_multiome_paired/tables/T2_MS_per_nucleus_rna_vs_atac.csv")
T3 = read("O1_own_multiome_paired/tables/T3_state_profile_concordance.csv")
T3b = T3[T3.test.astype(str).str.startswith("T3b")].iloc[0]
T4c = read("O1_own_multiome_paired/tables/T4c_level_effect_control.csv").iloc[0]
D3 = read("D_crosschemistry/tables/D3_crossconcordance_domainadjust.csv")
D6TAB = os.environ.get("D6TAB", "D_crosschemistry/tables_log1p_20261104")  # 2026-11-04 修正重跑
LAD = read(D6TAB + "/fig6D_crosschem_ladder.csv")
V2 = read(D6TAB + "/V2_bridge_ladder_dsciCAT.csv")
MSR = read(D6TAB + "/fig6D_MS_row.csv")
MG = read("SUPPORT_v1/tables/figS6h_matched_genes_summary.csv")
SIXC = read("FIG6_v17/tables/fig6C_sixpoint_all74.csv")
SIXD = read(D6TAB + "/fig6D_subtype_6pt.csv")
AUC = read(D6TAB + "/fig6S_aucell_addmodule_crosscheck.csv")
M2B = read("M2b_authoranno_v1/tables/M2b_per_mouse.csv")
PCT = read("SUPPORT_v1/tables/figS6g_pctcells_fdr02.csv")
DER = read("SUPPORT_v2/tables/figS6v2_derived_numbers.csv")
BIN = read("SUPPORT_v2/tables/figS6v2_bins_check.csv")
PDIG = read("SUPPORT_v2/tables/figS6_perm_hist_digitized.csv")
HIT = pd.read_csv(WORK + "/FY/oligo/scdrs_ms_out/ms_locus_overlap_table.csv")
row_o1 = D3[D3["pair"].astype(str).str.startswith("O1 自家")].iloc[0]
ORDER = ["dsciCAT-seq (own)", "10x Multiome", "HT-scCAT-seq"]
Li = LAD.set_index("tech").loc[ORDER].reset_index()
SRC4B = "O1_own_multiome_paired/tables/T4b_cross_channel_rho_variants.csv"

ROWS = []
def add(check, stat, value, n="", source="", col="", ci="", ref="", rz="", me="", main="", fig=""):
    ROWS.append(dict(check=check, statistic=stat, value=value, ci95_or_IQR=ci, reference=ref,
                     rejects_zero=rz, meets_effect=me, n=n, source=source, source_column=col,
                     in_figure=fig, in_main_text=main))


def verdict(v, p=None):
    """把交付判据打成两列勾：拒零 = 交付 P < 0.05；达效应量 = rho >= 0.60。"""
    me = "yes" if np.isfinite(v) and v >= CRIT else "no"
    if p is None or (isinstance(p, float) and np.isnan(p)):
        rz = "n.d."
    else:
        rz = "yes" if float(p) < 0.05 else "no"
    return rz, me


# ============ 台账带 A：对比定义敏感性（回执 §7.1 的标签更正）============
for _, r in T4b.iterrows():
    v = float(r.rho_cross_channel)
    rz, me = verdict(v, r.P_perm)
    is_sens = str(r.variant).startswith("sens")
    add("A contrast-definition sensitivity", str(r.variant), "%.3f" % v,
        n="74 traits, 5,962 nuclei", source=SRC4B,
        col="rho_cross_channel, P, P_perm, null_975",
        ci="null 97.5%% = %.3f; no bootstrap CI delivered" % float(r.null_975),
        ref="ceiling 0.878 / criterion 0.60", rz=rz, me=me,
        main=("0.424 = ALTERNATIVE POLE DEFINITION, not a domain adjustment; P_perm = 5.0e-4 still "
              "rejects zero, analytic P = 1.6e-4, only the effect-size criterion is missed (receipt 7.1/7.2)"
              if is_sens else ""),
        fig="Table only (forest demoted, triage receipt 3.1)")
v = float(row_o1.rho_adj)
rz, me = verdict(v, row_o1.perm_P_adj)
add("A contrast-definition sensitivity", "gene-set-size adjusted (ln n_hit_atac, D3 O1 row)",
    "%.3f" % v, n="74 traits", source="D_crosschemistry/tables/D3_crossconcordance_domainadjust.csv",
    col="rho_adj, ci_lo_adj, ci_hi_adj, perm_P_adj",
    ci="[%.3f, %.3f]" % (row_o1.ci_lo_adj, row_o1.ci_hi_adj),
    ref="ceiling 0.878 / criterion 0.60", rz=rz, me=me,
    main="RELABELED 2026-11-06: this row was printed as 'domain-size adjusted' but the D3 O1 row "
         "residualises on ln n_hit_atac (sensitivity on ln n_hit_rna gives 0.823); the regulatory-domain "
         "(ln n_peak) adjustment is the row below (receipt 7.1 superseded by CONFOUND_ADJUSTMENT_LIT_v1 8.7)",
    fig="Table only (forest demoted, triage receipt 3.1)")
_npk = pd.read_csv(SUP2 + "/tables_v11/figS6v11_unit_summary_ext.csv").set_index("unit").loc["trait_O1_npeak"]
add("A contrast-definition sensitivity", "regulatory-domain-size adjusted (ln n_peak, O1 rung)",
    "%.3f" % _npk.rho_point, n="74 traits", source="SUPPORT_v2/tables_v11/figS6v11_unit_summary_ext.csv",
    col="rho_point, ci_lo, ci_hi (B = 1,000, seed 20261104)",
    ci="[%.3f, %.3f]" % (_npk.ci_lo, _npk.ci_hi),
    ref="D10 B = 10,000: 0.550-0.815, perm P 0.0001; residual confound |rho| <= 0.055", rz="yes", me="yes",
    main="the true domain-size adjustment: moves the concordance by -0.111 (0.819 -> 0.708); "
         "shown as the 'trait n_peak adj' rung of Fig S6F",
    fig="Fig S6F rung + Table")

# ============ 台账带 B：测量可靠性 ============
for key, nm in [("rna_within_channel", "RNA split-half (within channel)"),
                ("atac_within_channel", "ATAC split-half (within channel)"),
                ("xchan_half1_vs_half2", "cross-channel, independent nucleus halves"),
                ("xchan_same_half", "cross-channel, same nuclei halves")]:
    r = T5.loc[key]
    v = float(r.median_split_half_rho)
    add("B measurement reliability", nm, "%.3f" % v, n="%d draws" % int(r.n_draw),
        source="O1_own_multiome_paired/tables/T5_split_half_reliability.csv",
        col="median_split_half_rho, q25, q75", ci="IQR [%.3f, %.3f]" % (r.q25, r.q75),
        ref="feeds the 0.878 ceiling", rz="n.a.", me=verdict(v)[1],
        fig="S6c prints the median only; the 1,000 per-draw values behind all four metrics are now "
            "delivered in SUPPORT_v2/tables_v12/splithalf_T5_reliability_draws.csv (20-draw T5 recomputed "
            "to 0.898/0.863/0.775/0.773, delta <= 0.013)")
for item, nm in [("noise ceiling sqrt(relRNA*relATAC)", "noise ceiling sqrt(relRNA x relATAC)"),
                 ("attenuation-corrected rho", "attenuation-corrected rho"),
                 ("cross-channel rho (raw)", "cross-channel rho (raw, = band A raw)"),
                 ("MS off-diagonal residual z", "MS off-diagonal residual z"),
                 ("MS off-diagonal rank of 74", "MS off-diagonal rank (of 74)")]:
    v = float(T6.loc[item, "value"])
    add("B measurement reliability", nm, "%.4f" % v, n="74 traits" if "rho" in item else "74 traits",
        source="O1_own_multiome_paired/tables/T6_cross_channel_summary.csv", col="value",
        ci="not delivered", ref="ceiling 0.878 / criterion 0.60",
        rz="n.a." if "z" in item or "rank" in item else verdict(v)[0],
        me=verdict(v)[1] if "rho" in item else "n.a.",
        main="attenuation-corrected rho = 0.933 is delivered but was never on any figure (receipt 3.1)"
        if item == "attenuation-corrected rho" else "",
        fig=("S6c (printed: ceiling 0.878)" if item == "noise ceiling sqrt(relRNA*relATAC)" else
             "S6c / S6e (printed 0.819)" if item == "cross-channel rho (raw)" else
             "Table only (receipt 3.1 asks for it to be stored here)"))

# ============ 台账带 C：水平 vs 形状 ============
for col, nm, pcol in [("rho_trait_means_cross_channel", "rho on trait absolute levels", "P_trait_means"),
                      ("rho_delta1_after_removing_trait_level", "rho on contrast after removing trait level", None)]:
    v = float(T4c[col])
    p = float(T4c[pcol]) if pcol else np.nan
    rz, me = verdict(v, p)
    if pcol:
        rz = "no"                                  # 交付的解析 P = 0.806，不拒零
    add("C level vs shape", nm, "%.3f" % v, n="74 traits",
        source="O1_own_multiome_paired/tables/T4c_level_effect_control.csv", col="%s%s" % (col, "" if not pcol else ", " + pcol),
        ci="not delivered", ref="criterion 0.60", rz=rz, me=me,
        main="P = 0.806 for the level row (delivered analytic P)" if pcol else "", fig="Table only")

# ============ 台账带 D：逐核共变（阴性）============
for key, nm in [("Multiple_sclerosis", "MS"), ("Schizophrenia", "SCZ"), ("HEIGHT", "HEIGHT")]:
    r = T2[T2.trait.astype(str).str.contains(key, regex=False)].iloc[0]
    v = float(r.rho_raw)
    add("D per-nucleus co-variation", "per-nucleus cross-channel rho, " + nm, "%.3f" % v,
        n="%s nuclei" % format(int(r.n_cell), ","),
        source="O1_own_multiome_paired/tables/T2_MS_per_nucleus_rna_vs_atac.csv",
        col="rho_raw, P_raw, rho_depth_adjusted", ci="not delivered",
        ref="criterion 0.60", rz="yes" if float(r.P_raw) < 0.05 else "no", me=verdict(v)[1],
        main="the three per-nucleus values are the negative set; depth-adjusted 0.079/0.112/0.061 same table"
        if nm == "MS" else "", fig="Table only")

# ============ 台账带 E：跨化学梯子 ============
for _, r in Li.iterrows():
    v = float(r.rho_adj)
    add("E cross-chemistry ladder", "rho_adj, " + str(r.tech), "%.3f" % v,
        n="Oligo %s / OPC %s / %d traits" % (format(int(r.n_oligo), ","), format(int(r.n_opc), ","), int(r.n_traits)),
        source="D_crosschemistry/tables/fig6D_crosschem_ladder.csv",
        col="rho_adj, ci_lo_adj, ci_hi_adj, perm_P_adj, noise_ceiling, ladder_reproduced",
        ci="[%.3f, %.3f]" % (r.ci_lo_adj, r.ci_hi_adj),
        ref="ceiling %.3f / criterion 0.60" % float(r.noise_ceiling),
        rz="yes" if float(r.perm_P_adj) < 0.05 else "no", me=verdict(v)[1],
        main=("ladder value quoted in Results; CI %.2f-%.2f crosses neither bound of the criterion"
              % (r.ci_lo_adj, r.ci_hi_adj)
              if "own" in str(r.tech) else
              ("HT passes \u03c1 \u2265 0.6 on the corrected input: %.3f "
               "(delivered 0.518 was computed from un-logged chromatin matrices), P = 1.0e-4, "
               "%.0f%% of its ceiling" % (r.rho_adj, 100 * r.rho_adj / r.noise_ceiling)
               if "HT" in str(r.tech) else "")),
        fig="Table only (forest demoted, triage receipt 3.1); verdict also on Fig 6D")
    add("E cross-chemistry ladder", "rho_raw, " + str(r.tech), "%.3f" % float(r.rho_raw),
        n="%d traits" % int(r.n_traits), source="D_crosschemistry/tables/fig6D_crosschem_ladder.csv",
        col="rho_raw", ci="not delivered (Fisher CI of raw in same table)",
        ref="criterion 0.60", rz="yes", me=verdict(float(r.rho_raw))[1], fig="Table only")

# ============ 台账带 F：V2 桥接（失败）============
for it, nm in [(0, "C ladder as delivered (MOL6 vs COP/NFOL, 5,962 nuclei)"),
               (1, "V2: same criterion on the D matrices, raw"),
               (2, "V2: same criterion on the D matrices, domain-adjusted"),
               (3, "difference (V2 raw - C)"), (4, "criterion |difference| <= 0.05 ?")]:
    v = V2.value.iloc[it]
    isnum = pd.notna(v)
    add("F bridge QC (FAILED)", nm, ("%.3f" % float(v)) if isnum else "n.d.",
        n=("%s nuclei" % format(int(V2.n.iloc[it]), ",")) if pd.notna(V2.n.iloc[it]) else "n/a",
        source="D_crosschemistry/tables/V2_bridge_ladder_dsciCAT.csv", col="value, n, pass_0.05",
        ci="not delivered", ref="criterion |delta| <= 0.05",
        rz="n.a.", me=("NO" if str(V2["pass_0.05"].iloc[it]) == "NO" else ("yes" if isnum else "n.d.")),
        main=("0.819 (C era) and %.3f (D era, corrected 2026-11-04) are NOT the same quantity; "
               "0.819 -> %.3f/%.3f, delta = %.3f (receipt 7.2) - one sentence required in Results"
               % (float(LAD[LAD.tech == "dsciCAT-seq (own)"].rho_adj.iloc[0]),
                  float(V2.value.iloc[1]), float(V2.value.iloc[2]), float(V2.value.iloc[3]))
               if it == 0 else ""),
        fig="Table only (forest demoted, triage receipt 3.1)")

# ============ 图面上仍然出现的数（v7：S6a/S6b/S6c/S6d/S6e/S6f/S6g）============
add("locus overlap (S6b)", "observed hit rate within +-500 kb of MS lead SNPs", "11.5%",
    n="38 of 331 platform genes", source="FY/oligo/scdrs_ms_out/ms_locus_overlap_table.csv",
    col="n_rows / 331 platform genes", ci="recomputed 11.48% from 38/331", ref="perm P = 0.277 (source-computed, see next row)",
    rz="no", me="n.a.", fig="S6b (printed)")
add("locus overlap (S6b)", "permutation P for the overlap (one-sided, add-one corrected)", "0.2767723",
    n="10,000 permutations, seed=20261104", source="SUPPORT_v2/tables_v8/figS6v8_window_permutation.csv",
    col="stat='permutation P (one-sided, observed vs null)'", ci="null mean 34.33 hits [24, 45] (95%)",
    ref="P = (1 + 2,767)/(1 + 10,000); the earlier printed label 0.26 was read off the v7 figure drawn with "
        "1,000 draws from a superseded vector and is superseded by this row", rz="no", me="n.a.",
    main="one sentence: the overlap is NOT more than expected (P = 0.28, 10,000 permutations)", fig="S6b (printed)")
add("locus overlap (S6b)", "exceedances in the source-computed null (count of draws >= observed 38)", "2767",
    n="10,000 draws", source="SUPPORT_v2/tables_v12/window_permutation_null_10000.csv", col="null_hit_count",
    ci="n.a.", ref="raw 10,000-draw vector now delivered (closes the outstanding item)", rz="no", me="n.a.",
    main="vector delivered; histogram also in tables_v12/window_permutation_null_histogram.csv",
    fig="not printed (kept in the table only)")
add("locus overlap (S6b)", "tail fraction read off the DIGITISED null (SUPERSEDED, kept as QA record)", "0.201",
    n="1,000 permutations (digitised, 25 bars)", source="SUPPORT_v2/tables/figS6_perm_hist_digitized.csv",
    col="sum of bars with bin_left >= observed", ci="n.a.",
    ref="QA record only - do NOT quote: a value read off figure art cannot outrank the source-computed 0.277",
    rz="no", me="n.a.",
    main="provenance warning: digitised figure values are not independent evidence",
    fig="not printed (kept in the table only)")
add("locus overlap (S6b)", "sum of digitised bar counts (self-check of the digitisation)", "1000.0",
    n="25 bars", source="SUPPORT_v2/tables/figS6_perm_hist_digitized.csv", col="n_permutations",
    ci="n.a.", ref="matches the printed n = 1,000", rz="n.a.", me="n.a.", fig="S6b (printed as n = 1,000)")
for _, r in HIT.iterrows():
    add("MS GWAS locus track (S6a)", "%s (chr%s) lead SNP %s" % (r.gene, r["chr"], r.snp),
        "p = %.0e" % r.pvalue, n="1 locus", source="FY/oligo/scdrs_ms_out/ms_locus_overlap_table.csv",
        col="pvalue", ci="lead-SNP distance %+.0f kb" % (r.distance_bp / 1000.0),
        ref="TSS distance, not a test", rz="n.a.", me="n.a.",
        fig="S6a (4 of them drawn)" if r.gene in ["PHLDB1", "PHGDH", "BCAS1", "PLLP"] else "Table only")
for _, r in DER.iterrows():
    v = float(r.value)
    add("derived summary of a delivered column (no new test)", r.quantity, "%.4g" % v,
        n=r.note if isinstance(r.note, str) and r.note else "n/a", source=r.source, col=r.column,
        ci="n.a.", ref="descriptive only", rz="n.a.", me="n.a.", fig=("S6f (printed)" if "AUC" in r.quantity else
         "S6d (printed)" if ("state-level" in r.quantity or "rho_6pt" in r.quantity) else
         "S6b (printed)" if ("hit rate" in r.quantity or "permutation P" in r.quantity
                             or "platform gene" in r.quantity or "hit genes (rows" in r.quantity) else
         "not printed (table only)"))
add("six-point state profile (S6d)", "median rho, C matrices (T3b)", "%.3f" % float(T3b.rho),
    n="74 traits", source="O1_own_multiome_paired/tables/T3_state_profile_concordance.csv", col="rho (T3b row)",
    ci="not delivered", ref="criterion 0.60", rz="n.d.", me=verdict(float(T3b.rho))[1], fig="S6d (printed)")
add("six-point state profile (S6d)", "share rho > 0, C matrices (T3b)", "%.1f%%" % (100 * float(T3b.frac_pos)),
    n="74 traits", source="O1_own_multiome_paired/tables/T3_state_profile_concordance.csv",
    col="frac_pos (T3b row)", ci="n.a.", ref="descriptive", rz="n.a.", me="n.a.", fig="S6d (printed)")
add("six-point state profile (S6d)", "median rho, D matrices",
    "%.3f" % float(SIXD.rho_6pt.median()), n="74 traits",
    source="D_crosschemistry/tables/fig6D_subtype_6pt.csv", col="rho_6pt (column median)",
    ci="not delivered", ref="criterion 0.60", rz="n.d.", me=verdict(float(SIXD.rho_6pt.median()))[1],
    fig="S6d (printed)")
add("six-point state profile (S6d)", "share rho > 0, D matrices",
    "%.1f%%" % (100 * float((SIXD.rho_6pt > 0).mean())), n="74 traits",
    source="D_crosschemistry/tables/fig6D_subtype_6pt.csv", col="rho_6pt > 0 (column share)",
    ci="n.a.", ref="descriptive", rz="n.a.", me="n.a.", fig="S6d (printed)")
for _, r in SIXC.iterrows():
    add("six-point state profile, C matrices (per trait)", str(r.trait),
        "rho %.3f" % float(r.rho), n="%d points" % int(r.n_points),
        source="FIG6_v17/tables/fig6C_sixpoint_all74.csv", col="rho, P",
        ci="P = %.3f" % float(r.P), ref="criterion 0.60", rz="no", me=verdict(float(r.rho))[1],
        fig="S6d (open squares, 3 traits only)")
for _, r in SIXD[SIXD.trait.astype(str).str.contains("Multiple_sclerosis")].iterrows():
    add("six-point state profile, D matrices (per trait)", str(r.trait),
        "rho %.3f" % float(r.rho_6pt), n="%d points" % int(r.n_points),
        source="D_crosschemistry/tables/fig6D_subtype_6pt.csv", col="rho_6pt, P_6pt",
        ci="P = %.3f" % float(r.P_6pt), ref="criterion 0.60", rz="no", me=verdict(float(r.rho_6pt))[1],
        main="MS = -0.143 (C) vs +0.543 (D): the six-point reading is NOT source-invariant (L5 disclosure)",
        fig="S6d (annotated both points)")
add("per-mouse positive control (S6f, GSE250589)", "mice with ATAC split-half AUC >= 0.95", "18 of 18", n="18 mice",
    source="M2b_authoranno_v1/tables/M2b_per_mouse.csv", col="auc_test", ci="0.952-0.997 (delivered range)",
    ref="chance 0.5", rz="yes", me="yes", fig="S6f (printed)")
add("per-mouse positive control (S6f, GSE250589)", "permutation P, every mouse", "%.4f" % float(M2B.p_perm.iloc[0]),
    n="101 permutations", source="M2b_authoranno_v1/tables/M2b_per_mouse.csv", col="p_perm",
    ci="n.a.", ref="P < 0.05", rz="yes", me="n.a.", fig="S6f (printed)")
add("per-mouse positive control (S6f, GSE250589)", "median AUC across the 18 mice", "%.3f" % float(M2B.auc_test.median()),
    n="18 mice", source="M2b_authoranno_v1/tables/M2b_per_mouse.csv", col="auc_test (median)",
    ci="n.a.", ref="descriptive", rz="n.a.", me="n.a.", fig="S6f (printed)")

# ============ 从图面搬下来的标量：%cells / AUCell / 基因数 / 核数 / MS 行 ============
for _, r in PCT[~PCT.subtype.isin(["RANGE_pp", "MAX_OVER_MIN"])].iterrows():
    add("share of nuclei passing scDRS FDR < 0.2", "%s, %s" % (r.subtype, r.channel),
        "%.1f%%" % float(r.pct), n=format(int(r.n_cells), ","),
        source="SUPPORT_v1/tables/figS6g_pctcells_fdr02.csv", col="pct", ci="n.a.",
        ref="FDR < 0.2 within subtype", rz="n.a.", me="n.a.", fig="Table only (was a bar panel)")
for ch, lab in [("RNA scDRS", "expression"), ("ATAC gene-score scDRS", "chromatin")]:
    add("share of nuclei passing scDRS FDR < 0.2", "range across 6 subtypes, " + lab,
        "%.1f pp" % float(PCT[(PCT.subtype == "RANGE_pp") & (PCT.channel == ch)].pct.iloc[0]),
        n=format(int(PCT[PCT.subtype == "RANGE_pp"].n_cells.iloc[0]), ","),
        source="SUPPORT_v1/tables/figS6g_pctcells_fdr02.csv", col="pct (RANGE_pp row)", ci="n.a.",
        ref="descriptive", rz="n.a.", me="n.a.", fig="Table only")
    add("share of nuclei passing scDRS FDR < 0.2", "max/min across 6 subtypes, " + lab,
        "%.2fx" % float(PCT[(PCT.subtype == "MAX_OVER_MIN") & (PCT.channel == ch)].pct.iloc[0]), n="n/a",
        source="SUPPORT_v1/tables/figS6g_pctcells_fdr02.csv", col="pct (MAX_OVER_MIN row)", ci="n.a.",
        ref="descriptive", rz="n.a.", me="n.a.", fig="Table only")
for ch in ["RNA", "gene_score"]:
    for m in ["AUCell", "AddModuleScore", "scDRS(norm_score)"]:
        r = AUC[(AUC.channel == ch) & (AUC.method == m)].iloc[0]
        add("method-independent scorers, MS Delta1", "%s, %s" % (m, ch), "%.4f" % float(r.d1),
            n="n_matched = %d genes" % int(r.n_matched),
            source="D_crosschemistry/tables/fig6S_aucell_addmodule_crosscheck.csv", col="d1",
            ci="not delivered", ref="compare scDRS 0.397 / 0.262", rz="n.a.", me="n.a.",
            main="independent scorers do NOT put MS on the chromatin axis either - a methodological "
                 "negative, ~10x less sensitive than scDRS (receipt 7.3)" if m == "scDRS(norm_score)" else "",
            fig="Table only (was a bar panel)")
for _, r in MG.iterrows():
    add("genes compared across chemistries", "%s: gene set RNA / ATAC" % r.tech,
        "%s / %s" % (format(int(r.n_genes_rna), ","), format(int(r.n_genes_atac), ",")),
        n="%d traits" % int(r.n_traits),
        source="D_crosschemistry/tables/fig6D_matched_genes.csv (flattened into SUPPORT_v1/tables/"
               "figS6h_matched_genes_summary.csv)", col="n_genes_rna, n_genes_atac", ci="n.a.",
        ref="identical in all three arms", rz="n.a.", me="n.a.", fig="Table only (was prose)")
    add("genes compared across chemistries", "%s: hits per trait RNA / ATAC" % r.tech,
        "%d-%d (med %d) / %d-%d (med %d)" % (int(r.hit_rna_min), int(r.hit_rna_max), int(r.hit_rna_med),
                                             int(r.hit_atac_min), int(r.hit_atac_max), int(r.hit_atac_med)),
        n="%d traits" % int(r.n_traits),
        source="D_crosschemistry/tables/fig6D_matched_genes.csv (same summary)", col="hit_*", ci="n.a.",
        ref="descriptive", rz="n.a.", me="n.a.", fig="Table only")
for _, r in Li.iterrows():
    add("nuclei per chemistry", "%s: Oligo / OPC / COP-NFOL" % r.tech,
        "%s / %s / %s" % tuple(format(int(x), ",") for x in [r.n_oligo, r.n_opc, r.n_cop]),
        n="%d traits" % int(r.n_traits), source="D_crosschemistry/tables/fig6D_crosschem_ladder.csv",
        col="n_oligo, n_opc, n_cop", ci="n.a.", ref="HT sparsity is shown on S6g", rz="n.a.", me="n.a.",
        fig="Table only (was a bar panel); sparsity visible in S6g")
for _, r in MSR.iterrows():
    add("MS on the chromatin ladder, per chemistry", "%s: Delta1 RNA / ATAC" % r.tech,
        "%+.3f / %+.3f" % (r.d1_rna, r.d1_atac), n="74 traits",
        source="D_crosschemistry/tables/fig6D_MS_row.csv", col="d1_rna, d1_atac", ci="n.a.",
        ref="sign + magnitude is what S6g shows", rz="n.a.", me="n.a.", fig="S6g (printed on the boxes)")
    add("MS on the chromatin ladder, per chemistry", "%s: off-diagonal z, rank of MS" % r.tech,
        "z = %+.2f, rank %d/74" % (r.offdiag_z, int(r.offdiag_rank)), n="74 traits",
        source="D_crosschemistry/tables/fig6D_MS_row.csv", col="offdiag_z, offdiag_rank", ci="n.a.",
        ref="criterion: rank_atac_adj > rank_rna_adj and z > 0", rz="n.a.",
        me="yes" if (r.offdiag_z > 0 and r.rank_atac_adj > r.rank_rna_adj) else "no", fig="Table only")
    add("MS on the chromatin ladder, per chemistry", "%s: adjusted rank RNA / ATAC" % r.tech,
        "%d / %d" % (int(r.rank_rna_adj), int(r.rank_atac_adj)), n="74 traits",
        source="D_crosschemistry/tables/fig6D_MS_row.csv", col="rank_rna_adj, rank_atac_adj", ci="n.a.",
        ref="both ranks always printed together", rz="n.a.", me="n.a.", fig="Table only")

# ============ 台账带 G：染色质侧可靠性（FIG6C_SUPPLEMENT_RELIABILITY_PLAN_20261104 §8.3/§8.5）============
# 2026-11-04 修正重跑后，DER6 快照里引用 atac 列的行必须从修正表重算（杜绝旧值转印）：
#  - cross-arm rho, ATAC * → 用 tables_log1p_20261104/fig6D_delta1_74traits_per_tech.csv 重算，标 pre-correction 旧值；
#  - cross-channel rho_adj * → 与带 E 重复且带 E 已是修正值，直接跳过；
#  - RNA/ATAC axis-contrast ratio * → METHODS_HARD_AUDIT_v1 A1：跨模态幅度比不可解释，撤除转印。
DL6 = read(D6TAB + "/fig6D_crosschem_ladder.csv")
DT6 = read(D6TAB + "/fig6D_delta1_74traits_per_tech.csv")
TECH6 = {"own": "dsciCAT-seq (own)", "10x": "10x Multiome", "HT": "HT-scCAT-seq"}

def atac_cross_arm(q, old_v):
    """"cross-arm rho, ATAC, a vs b (on ...)" → 从修正表重算"""
    import re as _re
    m = _re.match(r"cross-arm rho, ATAC, (\w+) vs (\w+) \((on (.+))\)", q)
    a, b, on = TECH6[m.group(1)], TECH6[m.group(2)], m.group(3)
    col = "rank_atac_adj" if "n_peak-adjusted" in on else "d1_atac"
    x = DT6[DT6.tech == a][col].rank(); y = DT6[DT6.tech == b][col].rank()
    import numpy as _np
    v = _np.corrcoef(x, y)[0, 1]
    return v, ("recomputed from corrected atac column %s (tables_log1p_20261104); delivered %.4f was "
               "pre-correction (un-logged input)" % (col, old_v))

DER6 = read("SUPPORT_v2/tables/figS6v6_derived_numbers.csv")
SRC6 = "SUPPORT_v2/tables/figS6v6_derived_numbers.csv"
for _, r in DER6.iterrows():
    q, v = str(r.quantity), float(r.value)
    if q.startswith("cross-arm rho, ATAC"):
        v, note = atac_cross_arm(q, v)
        printed = "S6c (printed)"
        src6, col6 = D6TAB + "/fig6D_delta1_74traits_per_tech.csv", "d1_atac/rank_atac_adj by tech"
    elif q.startswith("cross-channel rho_adj"):
        continue            # 带 E 已报修正值，避免同一表两处数字
    elif "ratio" in q and "axis-contrast" in q:
        continue            # A1：跨模态幅度比撤除转印（METHODS_HARD_AUDIT_v1）
    else:
        note, printed, src6, col6 = None, None, r.source, r.column
    printed = printed if printed else (
        "S6c (printed)" if ("cross-arm" in q or "T5" in q or "ceiling" in q or "observed / ceiling" in q
                            or "n_hit" in q) else
        "S6e (printed)" if ("ladder rung" in q or "axis-contrast" in q or "rho(n_hit" in q) else
        "S6d (printed)" if "state-level" in q else
        "S6f (printed)" if "ratio" in q or "AUC" in q else
        "S6g (printed)" if "rho_adj" in q else "not printed (table only)")
    add("chromatin-side reliability (receipt 20261104 plan)", q, "%.4f" % v if note else "%.4g" % v,
        n=(note if note else (r.note if isinstance(r.note, str) and r.note else "74 traits / 20 draws as stated")),
        source=src6, col=col6, ci="n.a.", ref="read from the delivered table, no new test",
        rz="n.a.", me="n.a.", fig=printed)
add("chromatin-side reliability (receipt 20261104 plan)",
    "RNA/ATAC axis-contrast ratio (min/median/max) — REMOVED", "n.a.",
    n="per METHODS_HARD_AUDIT_v1 A1: norm_score scales differ across modalities, cross-modality amplitude "
      "ratios are not interpretable (was 3.34 / 7.00 / 14.77 on the delivered matrices)",
    source="figS6v6_derived_numbers.csv (superseded)", col="—", ci="n.a.", ref="A1", rz="n.a.", me="n.a.",
    fig="not printed (table only)")
for _, r in Li.iterrows():
    ok = bool(r.ladder_reproduced)
    add("per-chemistry cross-channel judgement (S6g)", "%s: rho_adj vs the 0.60 criterion" % r.tech,
        "%.3f" % float(r.rho_adj), n="%s nuclei, 74 traits" % format(int(r.n_oligo) + int(r.n_opc), ","),
        source="D_crosschemistry/tables/fig6D_crosschem_ladder.csv",
        col="rho_adj, ladder_reproduced", ci="95%% CI [%.3f, %.3f]" % (float(r.ci_lo_adj), float(r.ci_hi_adj)),
        ref="criterion 0.60 (cross-channel, within one chemistry)", rz="yes",
        me="yes" if ok else "no", fig="S6g (printed under the box groups as %s)" % ("PASS" if ok else "FAIL"))
add("per-chemistry cross-channel judgement (S6g)", "noise ceiling per chemistry, own arm",
    "%.3f" % float(Li[Li.tech == "dsciCAT-seq (own)"].noise_ceiling.iloc[0]), n="74 traits",
    source="D_crosschemistry/tables/fig6D_crosschem_ladder.csv", col="noise_ceiling", ci="n.a.",
    ref="not compared with 0.60; descriptive", rz="n.a.", me="n.a.", fig="Table only")

# ============ 全局口径 ============
for item, val, src, col in [("n permutations for every rho ladder", "10,000",
                             "D_crosschemistry/tables/fig6D_crosschem_ladder.csv", "n_perm"),
                            ("seed", "20261102", "D_crosschemistry/tables/fig6D_crosschem_ladder.csv", "seed"),
                            ("perm P, C ladder (all T4b variants)", "5.0e-04", SRC4B, "P_perm"),
                            ("hexbin min-count threshold on S6b", "5 nuclei per bin",
                             "scripts/routeA/make_S6_v2_panels.py (drawing parameter, not a statistic)", "-"),
                            ("delivered hexbin bin counts (gridsize 15/25/40)", "147 / 362 / 509",
                             "figures/FINAL_SET_v2_20261028/Figure_S6.pdf panel E", "printed labels"),
                            ("bin counts when the same data are redrawn at our axes aspect", "150 / 645",
                             "SUPPORT_v2/tables/figS6v2_bins_check.csv (counts depend on the axes box, "
                             "not the data - deliberately NOT printed)", "redrawn_*")]:
    add("settings / drawing parameters", item, val, n="n/a", source=src, col=col, ci="n.a.",
        ref="reproducibility", rz="n.a.", me="n.a.", main=val if False else "",
        fig="legend / not printed")

# ============ 回执点名要、过去交付侧没有、现已由 A1 关闭（a1_fill_outstanding_v12.py） ============
V12 = "SUPPORT_v2/tables_v12/"
RES = [("C matrices, per-trait six-point rho for all 74 traits",
        "DELIVERED 74 rows (median rho_6pt = 0.600; BH q and rank per trait)", "74 traits",
        V12 + "detail_sixpoint_all74.csv", "rho_6pt, P_6pt, q_BH_74, rank_rho_desc + 12 subtype means",
        "closes the item: every trait now has a row, so the six-point scatter no longer needs to live on a panel"),
       ("per-trait split-half Delta1 pairs (plan group A/B)",
        "DELIVERED 74 traits x 2 channels = 148 rows (1,000 stratified exact-half draws, seed 20261112); "
        "delta1_full reproduces the delivered T4 table to 0.00000 on every trait and channel; "
        "the same draws also recompute the four T5 reliability metrics (median rho 0.898 / 0.863 / 0.775 / 0.773 "
        "against the delivered 20-draw 0.895 / 0.860 / 0.777 / 0.785)",
        "5,962 same-nucleus cells; 2,574 MOL6 / 203 COP-NFOL, halves of 1,287 / 101",
        V12 + "splithalf_delta1_all74.csv (+ splithalf_T5_reliability_draws.csv, _summary.csv, "
        "splithalf_across74_summary.csv)",
        "delta1_full, delta1_halfA_mean, delta1_halfB_mean, half_diff_mean/se, ci95_lo/hi, sign_concordance, "
        "both_halves_same_sign_as_full; per-draw: rna_within_channel, atac_within_channel, xchan_same_half, "
        "xchan_half1_vs_half2",
        "read off the delivered per-nucleus scDRS norm_score (RNA full74/, ATAC out_gs74/), so no re-scoring was "
        "needed; the halves share the whole-sample empirical normalisation, so this quantifies SAMPLING "
        "variability of Delta1 at n/2 only, not re-normalisation variability. Sign concordance is 1.000 for "
        "74/74 (RNA) and 73/74 (ATAC) traits, which is expected for a median-based statistic on a shared "
        "normalisation and must not be read as independent replication"),
       ("rho between the gridded UMAP views (gridsize 15 / 25 / 40)",
        "DELIVERED 6 pairs (raster view 0.696 / 0.690 / 0.236; nearest-hex view also given)", "5,962 nuclei, 160x160 raster",
        V12 + "hexbin_gridsize_rho.csv", "spearman_rho_raster, spearman_rho_nearesthex",
        "computed from the same H5AD + score.gz the panels are drawn from; both resampling views printed"),
       ("'% signal' per locus (receipt 4 lists 13 values)",
        "NOT REPRODUCIBLE - definition and algorithm were never delivered with the old bars; "
        "replaced by a traceable per-chromosome hit table (38/331 = 11.48%, 23 chromosomes)",
        "331 platform genes", V12 + "window_signal_per_chr.csv + detail_window_gene_hits_331.csv",
        "platform_genes, hits_in_window, pct_of_platform, closest_lead_snp",
        "receipt 4/9 forbid inventing the quantity; we deliver the underlying per-gene window membership instead"),
       ("raw permutation vector behind the S6b",
        "DELIVERED 10,000 draws (P recomputes to 0.2767723, add-one)", "10,000 draws, seed=20261104",
        V12 + "window_permutation_null_10000.csv", "null_hit_count, null_pct",
        "supersedes the digitised 25-bar inset and the 0.26 read off v7 art"),
       ("the 8 mice excluded by the M2b rule (MOL < 20 or OPC < 20 nuclei)",
        "DELIVERED: 26 screened / 18 passed / 8 excluded, all 8 for OPC < 20 "
        "(Early1, EarlyCtrl2, EarlyCtrl3, EarlyFemale3, EarlyMale3, Late1, PeakCtrl2, PeakCtrl3)",
        "26 mice", V12 + "mice_all26_qc_rule.csv",
        "n_MOL, n_OPC, pass_M2b, reason_excluded, delivered_in_M2b_per_mouse",
        "rebuilt from routeA_probe/mouse_EAE/anno_min.pkl; the passing set reproduces M2b_per_mouse.csv exactly "
        "(checked True), so the exclusion is by OPC scarcity, not by MOL"),
       ("mean of the per-mouse channel gap (per-mouse values in Table S7)",
        "DELIVERED: mean -2.108 (all 18), -1.965 (EAE 13), -2.479 (Ctrl 5)", "18 mice",
        V12 + "per_mouse_channel_gap.csv", "gap_atac_minus_rna, abs_gap",
        "no panel prints it; per-mouse values are in Table S7")]
for a1, a2, a3, a4, a5, a6 in RES:
    add("A1 resolution of outstanding items (20261112)", a1, a2, n=a3, source=a4, col=a5, ci="n.a.",
        ref="a1_fill_outstanding_v12.py", rz="n.a.", me="n.a.", main=a6, fig="table only")

# ====== 回执 v16 #1 + 整改单 2026-11-10 §4：正文 §3.10 / 图注 E 把 Fig 7 分析链（Table S6/S7）的数归到 Table S6 ======
# 这些值全部属于 Fig 7 分析链侧（GSE250589 multiome gene-activity channel / three-mouse aggregate /
# EAE composition），不是 Fig S7 的面值；docx v16 正文明确 "These per-library gene-activity
# statistics belong to the Figure 7 analysis chain (Table S6)"。从交付表读取，不硬编码。
_RES = read("SUPPORT_v1/tables/S7_atac_rescore_fullstats.csv")
_ms = _RES[_RES.trait.astype(str).str.contains("Multiple_sclerosis", regex=False)].iloc[0]
_r, _p = spearmanr(_RES.atac_d1_fullstats, _RES.rna_d1_3mouse)
_ratio = _RES.atac_d1_fullstats.std() / _RES.rna_d1_3mouse.std()
_COMP = read("D_crosschemistry/tables/fig6D_composition.csv")
_comp = _COMP[_COMP.quantity.astype(str).str.startswith("composition ratio")].iloc[0]
_eadj = _COMP[(_COMP.quantity.astype(str).str.contains("adjusting log composition")) &
              (_COMP.axis == "RNA")].iloc[0]
_aadj = _COMP[(_COMP.quantity.astype(str).str.contains("adjusting log composition")) &
              (_COMP.axis == "ATAC")].iloc[0]
add("Fig 7 analysis chain (Table S6/S7)",
    "MS gene-activity rank (of 74, chromatin arm, full-data rescore)",
    "%d of 74" % int(_ms.rank_fullstats), n="74 traits",
    source="SUPPORT_v1/tables/S7_atac_rescore_fullstats.csv",
    col="rank_fullstats", ci="n.a.", ref="Fig 7 analysis chain", rz="n.a.", me="n.a.",
    main="docx v16 §3.10: 'gene-activity rank 28 / 74' — a different estimator from Fig S7E "
         "(same-nucleus gene-set, MS 36/74); the text must not equate them",
    fig="Table only (not on Fig S7 face)")
add("Fig 7 analysis chain (Table S6/S7)",
    "cross-channel Spearman rho (three-mouse aggregate RNA vs full-data chromatin), 74 traits",
    "%.3f" % _r, n="74 traits",
    source="SUPPORT_v1/tables/S7_atac_rescore_fullstats.csv",
    col="atac_d1_fullstats vs rna_d1_3mouse",
    ci="P = %.1e" % _p, ref="three-mouse aggregate (M3); P < 0.05", rz="yes", me="no",
    main="0.474 / P = 2e-05 (printed on Fig S7 D-row in earlier drafts; final sheet has no "
         "multiome channel, so this value lives in Table S6)",
    fig="Table only (not on Fig S7 face)")
add("Fig 7 analysis chain (Table S6/S7)",
    "SD ratio of the two channels (chromatin / RNA), three-mouse aggregate",
    "%.2f" % _ratio, n="74 traits",
    source="SUPPORT_v1/tables/S7_atac_rescore_fullstats.csv",
    col="std(atac_d1_fullstats) / std(rna_d1_3mouse)",
    ci="n.a.", ref="docx v16: 'SD ratio 0.26 versus 1.06 in RNA'", rz="n.a.", me="n.a.",
    main="0.26 is this ratio; 1.06 is the same ratio computed on the RNA-side per-library spread "
         "(M3_paired3 summary IQR_y) — both are descriptive, no test",
    fig="Table only (not on Fig S7 face)")
add("Fig 7 analysis chain (Table S6/S7)",
    "domain-size confound rho(Delta1, ln n_peak) in the mouse chromatin arm",
    "+0.050", n="74 traits",
    source="SUPPORT_v1/tables/figS7_v15_face_numbers.csv (D d4, M3_paired3)",
    col="D d4 rho(ln n_peak) mouse", ci="n.a.",
    ref="|rho| < 0.10: no domain-size confound in this arm", rz="no", me="n.a.",
    main="docx v16 §3.10: 'no domain-size confound in this arm, rho = +0.050'",
    fig="Table only (not on Fig S7 face)")
add("Fig 7 analysis chain (Table S6/S7)",
    "EAE composition n_MOL/n_OPC, log scale, EAE vs Ctrl",
    "%.2f vs %.2f" % (np.log(_comp.median_EAE), np.log(_comp.median_Ctrl)), n="18 mice (EAE 13 / Ctrl 5)",
    source="D_crosschemistry/tables/fig6D_composition.csv",
    col="median_EAE / median_Ctrl (log)",
    ci="MWU P = %.3f" % _comp.mwu_P, ref="P < 0.05", rz="yes", me="n.a.",
    main="2.78 vs 4.01 = ln(16.18) vs ln(55.06); the raw ratios are 16.2 vs 55.1 — the log scale "
         "is the one quoted in docx v16",
    fig="Table only (not on Fig S7 face)")
add("Fig 7 analysis chain (Table S6/S7)",
    "RNA-channel EAE effect after composition adjustment (OLS residual)",
    "%.3f" % _eadj.delta, n="18 mice",
    source="D_crosschemistry/tables/fig6D_composition.csv",
    col="delta (axis EAE-Ctrl after adjusting log composition, RNA)",
    ci="P = %.3f" % _eadj.mwu_P, ref="P >= 0.05: no composition-independent RNA effect",
    rz="no", me="n.a.", main="docx v16: 'the RNA-channel EAE effect after composition adjustment "
         "is P = 0.387, n.s.'", fig="Table only (not on Fig S7 face)")
add("Fig 7 analysis chain (Table S6/S7)",
    "ATAC-channel EAE effect after composition adjustment (OLS residual)",
    "+%.3f" % _aadj.delta, n="18 mice",
    source="D_crosschemistry/tables/fig6D_composition.csv",
    col="delta (axis EAE-Ctrl after adjusting log composition, ATAC)",
    ci="P = %.3f" % _aadj.mwu_P, ref="P >= 0.05: no ATAC effect",
    rz="no", me="n.a.", main="docx v16: 'the ATAC channel shows no effect, +0.011, P = 0.775'",
    fig="Table only (not on Fig S7 face)")

T6OUT = pd.DataFrame(ROWS)[["check", "statistic", "value", "ci95_or_IQR", "reference", "rejects_zero",
                            "meets_effect", "n", "source", "source_column", "in_figure", "in_main_text"]]
T6OUT.to_csv(TAB2 + "/Table_S6_v12.csv", index=False)
say("Table S6 台账：%d 行（其中相关性台账带 A-F %d 行）"
    % (len(T6OUT), int(T6OUT.check.str.match(r"^[A-F] ").sum())))

# ---------------- Table S7 ----------------
m = M2B.sort_values(["group", "mouse"]).reset_index(drop=True)
T7 = pd.DataFrame({
    "mouse": m.mouse, "group": m.group,
    "n_MOL": [format(int(x), ",") for x in m.n_MOL], "n_OPC": [format(int(x), ",") for x in m.n_OPC],
    "split_half_AUC": m.auc_test.round(4), "perm_P": m.p_perm.round(4), "qc_pass": m.qc_pass,
    "ATAC_axis_median": m.atac_axis_median.round(3), "RNA_axis_median": m.rna_axis_median.round(3),
    "ATAC_axis_contrast": m.atac_axis_contrast.round(3), "RNA_axis_contrast": m.rna_axis_contrast.round(3),
    "n_panel_genes": m.n_panel, "source": "M2b_authoranno_v1/tables/M2b_per_mouse.csv"})
T7.insert(0, "row", np.arange(1, len(T7) + 1))
T7.to_csv(TAB2 + ("/Table_S7_v11skip.csv" if os.environ.get("SKIP_S7") else "/Table_S7.csv"), index=False)
say("Table S7：%d 行 × %d 列" % (T7.shape))


# ---------------- 排字：A4 横向、Arial 6.2 pt 活字、按真实字宽截断 ----------------
_FONTS = {False: fitz.Font(fontfile=FONTS_PDF["normal"]), True: fitz.Font(fontfile=FONTS_PDF["bold"])}


def tw(txt, fs, bold=False):
    return _FONTS[bool(bold)].text_length(txt, fontsize=fs)


def clip(txt, lim, fs, bold=False):
    if tw(txt, fs, bold) <= lim:
        return txt
    while len(txt) > 1 and tw(txt + "\u2026", fs, bold) > lim:
        txt = txt[:-1]
    return txt.rstrip(" ,;(-/") + "\u2026"


def typeset(df, out_pdf, title, sub, col_w_cm, fs=6.2, lh=7.6):
    left = 0.9 * PT
    PAGE_H = 21.0 * PT
    assert sum(col_w_cm) <= 29.7 - 1.8, "表宽超出 A4 横向：%s cm" % sum(col_w_cm)
    cols = list(df.columns)
    doc = fitz.open()

    def newpage():
        pg = doc.new_page(width=29.7 * PT, height=PAGE_H)
        pg.insert_text((left, 0.7 * PT), title, fontname="Arial-Bold",
                       fontfile=FONTS_PDF["bold"], fontsize=8.4)
        pg.insert_text((left, 0.7 * PT + 10), sub, fontname="Arial-Regular",
                       fontfile=FONTS_PDF["normal"], fontsize=6.2, color=(0.25, 0.25, 0.25))
        yy = 0.7 * PT + 24
        x = left
        for c, w in zip(cols, col_w_cm):
            pg.insert_text((x, yy), clip(str(c), w * PT - 3, fs, bold=True),
                           fontname="Arial-Bold", fontfile=FONTS_PDF["bold"], fontsize=fs)
            x += w * PT
        pg.draw_line(fitz.Point(left, yy + 2.2), fitz.Point(left + sum(col_w_cm) * PT, yy + 2.2),
                     width=0.4, color=(0.3, 0.3, 0.3))
        return pg, yy + lh

    page, yy = newpage()
    for _, r in df.iterrows():
        if yy > PAGE_H - 1.0 * PT:
            page, yy = newpage()
        x = left
        for c, w in zip(cols, col_w_cm):
            txt = str(r[c]).replace("\n", " ")
            lim = w * PT - 3
            if tw(txt, fs) > lim:
                while len(txt) > 1 and tw(txt + "\u2026", fs) > lim:
                    txt = txt[:-1]
                txt = txt.rstrip(" ,;(-/") + "\u2026"
            page.insert_text((x, yy), txt, fontname="Arial-Regular", fontfile=FONTS_PDF["normal"],
                             fontsize=fs, color=(0.1, 0.1, 0.1))
            x += w * PT
        yy += lh
    doc.save(out_pdf)
    d = fitz.open(out_pdf)
    say("[table pdf] %-16s %d 页  %.1f x %.1f cm" % (os.path.basename(out_pdf), len(d),
                                                     d[0].rect.width * 2.54 / 72,
                                                     d[0].rect.height * 2.54 / 72))
    d.close()


T6d = T6OUT.rename(columns={"rejects_zero": "P<0.05?", "meets_effect": "rho>=0.60?",
                            "ci95_or_IQR": "CI / IQR", "source_column": "column",
                            "in_figure": "on figure", "in_main_text": "in main text"})
typeset(T6d, PDF2 + "/Table_S6_v12.pdf",
        "Table S6. Correlation ledger for Fig S6 - one row per value, with its delivered source "
        "(receipt 3.1: the rho forest plot is demoted to this table).",
        "'rejects zero?' = the delivered permutation/analytic P < 0.05; 'meets effect?' = rho >= 0.60 "
        "(pre-registered criterion); neither is a new statistic. Rows marked OUTSTANDING are values the "
        "triage receipt asked for that no delivered table carries. Source: "
        "SUPPORT_v2/tables/Table_S6_v11.csv.",
        [2.7, 3.5, 1.7, 2.4, 2.4, 1.2, 1.4, 1.9, 3.6, 2.0, 2.2, 2.9])
typeset(T7, PDF2 + ("/Table_S7_skip.pdf" if os.environ.get("SKIP_S7") else "/Table_S7.pdf"),
        "Table S7. Per-mouse numbers behind the mouse rung of Fig S6f (external GSE250589 paired multiome; "
        "ATAC-only MOL-vs-OPC split-half classifier).",
        "Source: M2b_authoranno_v1/tables/M2b_per_mouse.csv, sorted by group then name; the 18 rows are "
        "the same order as the 18 bootstrap units of the mouse rung in Fig S6f. The 8 mice excluded by the "
        "M2b rule (MOL < 20 or OPC < 20 nuclei; 26 screened, 18 delivered) are now listed in "
        "SUPPORT_v2/tables_v12/mice_all26_qc_rule.csv - all 8 failed on OPC < 20, none on MOL.",
        [0.9, 3.2, 1.4, 1.6, 1.5, 2.1, 1.5, 1.4, 2.1, 2.1, 2.1, 2.1, 1.7])
say("表格完成")
