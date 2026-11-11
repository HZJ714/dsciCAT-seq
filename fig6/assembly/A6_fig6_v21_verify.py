# -*- coding: utf-8 -*-
"""Fig 6 v21 acceptance hard gate (D = the cross-chemistry heatmap pair; all of C's in-figure text moved to the caption)
 1) every digit on the face must trace to a row in the delivered tables (any gap is a MISS; "roughly right" is not allowed)
 2) fonts/sizes: Arial only, floor 5.5 pt (submission hard gate, final_audit.py:206; v21 once broke it with 5.0 pt in-cell ranks)
 3) text-collision hard gate (the check_text_overlap criteria): the main figure + the S6m supplement
 4) new in v21: the author demanded "all D-face text regular Arial 8pt black, no greys" - audited span by span against D's layout rectangle:
    colour must be pure black, sizes only 8.0 (annotations) and FSR (in-cell ranks; the column is only 0.55 cm wide)
 5) every sentence deleted from the face (the top ladder line, C's five corner lines, C's per-trait ranks, C's grey title)
    must appear verbatim in fig6_v21_legend_draft.md - "moved to the caption" must really move, not vanish
Output FIG6_v21/tables/fig6_v21_audit.txt + fig6_all_numbers_v21.csv
"""
import os, sys, math

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
sys.path.insert(0, WORK + "/scripts/MS_GWAS/fig6_repro/deps")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fitz
import numpy as np, pandas as pd
from check_text_overlap import overlaps

FSR = float(os.environ.get("FSR", "5.5"))     # same source as make_fig6_v21_D_heatmap.py (same default)
FMIN = 5.5                                    # submission hard gate: the smallest allowed visible size
# the rank-digit placement matches the drawing script (halo = white-underlay digits / cell = named rows grow, cells blanked / none = not printed)
NUMMODE = os.environ.get("NUMMODE", "none")    # must match make_fig6_v21_D_heatmap.py's default
FSD = float(os.environ.get("FSD", "6.0"))
FDIG = FSR if NUMMODE == "halo" else FSD      # the legal size for colour rank digits (unused in none mode)
NDIGC = 18 if NUMMODE in ("halo", "cell") else 0   # on-face colour rank digits (none = the author's 11-03 decision: not printed)
NDIGW = 18 if NUMMODE == "halo" else 0        # white-underlay digits exist only in the halo version

CM = 72 / 2.54
FIG = WORK + "/FY/oligo/fig6_rework_v2"
S15 = WORK + "/figures/routeA_out/FIG6_v15/tables"
V17 = WORK + "/figures/routeA_out/FIG6_v17/tables"
V21 = WORK + "/figures/routeA_out/FIG6_v21/tables"
D = os.environ.get("TABD", WORK + "/figures/routeA_out/D_crosschemistry/tables")
PDF = os.environ.get("PDF", FIG + "/fig6_composed_v21.pdf")
S6M = WORK + "/figures/routeA_out/FIG6_v18/tables/fig6S6m_crosschem_scatters.pdf"
LEG = os.path.dirname(V21) + "/fig6_v21_legend_draft.md"
os.makedirs(V21, exist_ok=True)

A = pd.read_csv(S15 + "/fig6A_stats.csv").set_index("cell_type")
CS = pd.read_csv(S15 + "/fig6C_stats.csv").iloc[0]
T = pd.read_csv(WORK + "/figures/routeA_out/O1_own_multiome_paired/tables/T4_delta1_74traits_both_channels.csv")
DIS = pd.read_csv(V17 + "/fig6C_fourline_disclosure.csv")
SIX = pd.read_csv(V17 + "/fig6C_sixpoint_all74.csv")
DL = pd.read_csv(D + "/fig6D_crosschem_ladder.csv")
DM = pd.read_csv(D + "/fig6D_MS_row.csv")
DT = pd.read_csv(D + "/fig6D_delta1_74traits_per_tech.csv")
TECH = ["dsciCAT-seq (own)", "10x Multiome", "HT-scCAT-seq"]
MARKS = [("PASS_Multiple_sclerosis", "MS"), ("PASS_Schizophrenia_Pardinas2018", "SCZ"),
         ("UKB_460K.body_HEIGHTz", "HEIGHT")]


def norm(x):
    """Normalization for comparisons: strip whitespace, strip markdown bold asterisks, unify the minus sign (the caption writes U+2212; the PDF writes U+2212 too)"""
    return "".join(str(x).replace("\u2212", "-").replace("*", "").split())


def fmtp(p):
    if p >= 0.001:
        return "%.3f" % p
    e = int(math.floor(math.log10(p)))
    return "%.1fe%d" % (p / 10 ** e, e)


def disc(prefix, col="value"):
    return float(DIS[DIS.line.str.startswith(prefix)].iloc[0][col])


exp, keys, cap = [], [], []
# ---- A (unchanged from v20) ----
for ct in ["MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]:
    r = A.loc[ct]
    exp += [fmtp(r.mwu_p_vs_COPNFOL), "%.3f" % r.group_P, "%.3f" % r.group_q_BH5, "%+.2f" % r.cliffs_delta]
for ct in A.index:
    exp.append("n=%s" % format(int(A.loc[ct, "n"]), ","))
    keys.append(("A.n." + ct, int(A.loc[ct, "n"]), "fig6A_stats.csv"))
exp += ["scDRS norm_score", "(MS GWAS)", "COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
# ---- B ----
exp += ["MS scDRS score on UMAP", "UMAP_1", "UMAP_2", "norm_score", "MOL6", "COP/NFOL", "MOL2"]
# ---- C: v21 keeps only the named traits + axis labels; everything else moves to the caption ----
exp += ["MS", "SCZ", "HEIGHT", "Gene score (MOL6 \u2212 COP/NFOL)", "Gene expression",
        "(MOL6 \u2212 COP/NFOL)", "\u22122.0", "\u22122"]
for key, nm in [("Multiple_sclerosis", "MS"), ("Schizophrenia", "SCZ"), ("HEIGHTz", "HEIGHT")]:
    r = T[T.trait.str.contains(key, regex=False)].iloc[0]
    keys += [("C.%s_rank_rna_raw" % nm, int(r.rank_rna_raw), "O1/T4_delta1_74traits_both_channels.csv"),
             ("C.%s_rank_atac_raw" % nm, int(r.rank_atac_raw), "O1/T4_delta1_74traits_both_channels.csv"),
             ("C.%s_d1_rna" % nm, float(r.d1_rna), "O1/T4_delta1_74traits_both_channels.csv"),
             ("C.%s_d1_atac" % nm, float(r.d1_atac), "O1/T4_delta1_74traits_both_channels.csv")]
cap += ["z = %+.2f" % disc("L2"), "rank 1/74",       # MS's off-diagonal z is not a rank, so it retreats to the caption (ZLINE=1 restores it)
        "\u03c1 = %.3f" % disc("L1"), "perm P = %.1e" % CS.perm_P_delta1,
        "\u03c1 = %.3f" % disc("L3 six-subtype mean profile, cross-channel (trait=MS)"),
        "\u03c1 = %.3f" % disc("L4"), "ceiling %.3f" % CS.noise_ceiling,
        "level-removed %.3f" % CS.rho_delta1_level_removed, "# = rank within %d" % int(CS.n_traits),
        "trait absolute levels across channels: \u03c1 = %.3f, P = %.2f" % (CS.rho_trait_level, CS.P_trait_level),
        "dsciCAT-seq (own) \u00b7 same nuclei \u00b7 points = %d traits" % int(CS.n_traits)]
for key, nm in [("Multiple_sclerosis", "MS"), ("Schizophrenia", "SCZ"), ("HEIGHTz", "HEIGHT")]:
    r = T[T.trait.str.contains(key, regex=False)].iloc[0]
    exp += ["RNA #%d" % int(r.rank_rna_raw), "ATAC #%d" % int(r.rank_atac_raw)]
keys += [("C.L1_ladder_rho", disc("L1"), "fig6C_fourline_disclosure.csv"),
         ("C.L2_MS_offdiag_z", disc("L2"), "fig6C_fourline_disclosure.csv"),
         ("C.L3_MS_6point_rho", disc("L3 six-subtype mean profile, cross-channel (trait=MS)"), "fig6C_fourline_disclosure.csv"),
         ("C.L4_per_nucleus_rho", disc("L4"), "fig6C_fourline_disclosure.csv"),
         ("C.perm_P", float(CS.perm_P_delta1), "fig6C_stats.csv"),
         ("C.noise_ceiling", float(CS.noise_ceiling), "fig6C_stats.csv"),
         ("C.level_removed", float(CS.rho_delta1_level_removed), "fig6C_stats.csv")]
# ---- D: the heatmap pair. The face shows only ranks, colour-bar endpoints, channel names, chemistry names and the three named traits ----
exp += ["Gene expression", "Gene score", "dsciCAT", "10x", "Multiome", "HT-scCAT", "seq",
        "\u03941 (z)", "Rank", "+1.4", "0.0", "-1.4", "74", "37", "1"]
keys += [("D.clip_LIM", 1.40, "design constant: \u0394 z \u5217\u88c1\u65ad\u4f4d"),
         ("D.n_traits", int(CS.n_traits), "fig6C_stats.csv")]
for t in TECH:
    TT = DT[DT.tech == t]
    for key, nm in MARKS:
        q = TT[TT.trait.str.contains(key, regex=False)].iloc[0]
        # with NUMMODE=none the ranks are off the face and must live in the caption/text (the rank list assembled below),
        # hence they move from exp (must appear on the face) to cap (must appear in the caption) - otherwise the “every face digit traces to a table row” gate misses the 18 ranks.
        if NUMMODE != "none":
            exp += [str(int(q.rank_rna_adj)), str(int(q.rank_atac_adj))]
        else:
            cap += ["%d/%d" % (int(q.rank_rna_adj), int(q.rank_atac_adj))]
        keys += [("D.%s.%s_rank_rna_adj" % (t, nm), int(q.rank_rna_adj), "fig6D_delta1_74traits_per_tech.csv"),
                 ("D.%s.%s_rank_atac_adj" % (t, nm), int(q.rank_atac_adj), "fig6D_delta1_74traits_per_tech.csv"),
                 ("D.%s.%s_d1_rna" % (t, nm), float(q.d1_rna), "fig6D_delta1_74traits_per_tech.csv"),
                 ("D.%s.%s_d1_atac" % (t, nm), float(q.d1_atac), "fig6D_delta1_74traits_per_tech.csv")]
# the ladder verdict line deleted from D's face must appear wholesale in the caption
for i, r in DL.iterrows():
    cap += ["\u03c1_adj = **%.3f**" % r.rho_adj, "%.4f" % r.perm_P_adj,
            "95%% CI %.3f\u2013%.3f" % (r.ci_lo_adj, r.ci_hi_adj)]
    keys += [("D.%s.rho_adj" % r.tech, float(r.rho_adj), "fig6D_crosschem_ladder.csv"),
             ("D.%s.perm_P_adj" % r.tech, float(r.perm_P_adj), "fig6D_crosschem_ladder.csv"),
             ("D.%s.ci_lo_adj" % r.tech, float(r.ci_lo_adj), "fig6D_crosschem_ladder.csv"),
             ("D.%s.ci_hi_adj" % r.tech, float(r.ci_hi_adj), "fig6D_crosschem_ladder.csv"),
             ("D.%s.n_oligo" % r.tech, int(r.n_oligo), "fig6D_crosschem_ladder.csv"),
             ("D.%s.n_opc" % r.tech, int(r.n_opc), "fig6D_crosschem_ladder.csv")]
npass = int(DL.ladder_reproduced.sum())
cap += ["Oligo / OPC): ", "5,813 / 1,488", "540 / 298", "293 / 220", "\u03c1 \u2265 0.6",
        "%d of 3" % npass,
        "30,559-gene intersection", "n_ctrl = 1000", "measured Seurat matrix",
        "row order is the ascending mean of the three domain-adjusted *expression*-channel ranks",
        "within-column z-score", "clipped at \u00b11.4", "domain-adjusted rank within 74",
        "concordance is reported as \u03c1 (text), not inferred from the heatmap",
        "pre-registered effect-size screen", "rather than a confirmatory hypothesis test"]
# the MS cross-chemistry counter-example paragraph derives from the MS-row table, no longer hard-coded (the numbers changed after the 2026-11-04 corrected re-run)
for _, r in DM.iterrows():
    cap += ["%+.3f \u2192 %+.3f" % (r.d1_rna, r.d1_atac),
            "#%d \u2192 #%d" % (r.rank_rna_adj, r.rank_atac_adj),
            "z = %+.2f" % r.offdiag_z]
# the MS detail table must exist and match D's Δ1
SUPP = WORK + "/figures/routeA_out/FIG6_v20/tables/fig6D_MS_supplementary_table.csv"
assert os.path.exists(SUPP), "\u7f3a MS \u660e\u7ec6\u9644\u8868"
_sp = pd.read_csv(SUPP)
assert float(((_sp.d1_rna.values - DM.d1_rna.values).__abs__()).max()) < 5e-4, "\u9644\u8868\u4e0e D \u7684 \u03941 \u4e0d\u4e00\u81f4"

# ---- read the main figure ----
doc = fitz.open(PDF); pg = doc[0]
flat = norm(pg.get_text())
sizes, fonts, spans = [], set(), []
for b in pg.get_text("dict")["blocks"]:
    for l in b.get("lines", []):
        for s in l["spans"]:
            if s["text"].strip():
                sizes.append(round(s["size"], 2)); fonts.add(s["font"]); spans.append(s)
miss = [e for e in exp if norm(e) not in flat]
undersz = [s_ for s_ in spans if s_["text"].strip() and round(s_["size"], 2) < FMIN]

# ---- D-region text-style hard gate: all black + sizes only 8.0 / 6.0 ----
wC = fitz.open(FIG + "/fig6v21_C_contrast.pdf")[0].rect.width / CM
hA = fitz.open(FIG + "/fig6v21_A_box.pdf")[0].rect.height / CM
hB = fitz.open(FIG + "/fig6v21_B_umap.pdf")[0].rect.height / CM
xD = (wC + 0.24) * CM; yD = max(hA, hB) * CM
dsty = [s_ for s_ in spans if s_["bbox"][0] >= xD and s_["bbox"][1] >= yD - 1]
NAMES3 = ("MS", "SCZ", "HEIGHT")
# D-region text falls into four classes, each audited (author 11-03: rank digits take the trait colour, smaller, no white box)
names = [s_ for s_ in dsty if s_["text"].strip() in NAMES3]
dig_c = [s_ for s_ in dsty if s_["text"].strip().isdigit() and s_["color"] not in (0, 0xffffff)]   # colour rank digits
dig_w = [s_ for s_ in dsty if s_["text"].strip().isdigit() and s_["color"] == 0xffffff]            # white-underlay digits
rest = [s_ for s_ in dsty if s_ not in names and s_ not in dig_c and s_ not in dig_w]
bad = [(s_["text"], s_["color"], round(s_["size"], 2)) for s_ in rest
       if s_["color"] != 0 or round(s_["size"], 2) != 8.0]                     # everything else must be pure black 8 pt
bad += [(s_["text"], s_["color"], round(s_["size"], 2)) for s_ in names
        if round(s_["size"], 2) != 8.0]                                         # named traits 8 pt
bad += [(s_["text"], s_["color"], round(s_["size"], 2)) for s_ in dig_c
        if round(s_["size"], 2) != FDIG]                                       # rank digits = FDIG

# ---- cross-panel colour identity: D's named-trait colours must equal C's (comparing the rendered ink, not the source) ----
cspan = {}
for b in pg.get_text("dict")["blocks"]:
    for l in b.get("lines", []):
        for s_ in l["spans"]:
            t = s_["text"].strip()
            if t in NAMES3 and s_["bbox"][0] < xD and s_["bbox"][1] >= yD - 1:
                cspan[t] = s_["color"]
dcol = {s_["text"].strip(): s_["color"] for s_ in names}
mismatch = [t for t in NAMES3 if cspan.get(t) != dcol.get(t)]
# each rank digit's colour must equal its row trait's colour (rows assigned by y centre)
_rows = {t: min((s_["bbox"][1] + s_["bbox"][3]) / 2 for s_ in names if s_["text"].strip() == t) for t in NAMES3}
misnum = [(s_["text"], "#%06x" % s_["color"]) for s_ in dig_c
          if s_["color"] != cspan.get(min(NAMES3, key=lambda k: abs(_rows[k] - (s_["bbox"][1] + s_["bbox"][3]) / 2)))]
# the white box was removed as required: D must contain no white-filled rectangles at glyph size (whole-page white excluded)
doc2 = fitz.open(PDF); dr = doc2[0].get_drawings(); doc2.close()
wbox = [g for g in dr if g["rect"].x0 >= xD and g["rect"].y0 >= yD and g.get("fill")
        and max(g["fill"]) > 0.95 and 3 < g["rect"].width < 14 and 3 < g["rect"].height < 12]
occ = pd.read_csv(V21 + "/fig6C_label_occlusion.csv")
occ_bad = occ[occ.covered_points > 0]

# ---- the S6m supplement ----
exp6m, miss6m = [], []
if os.path.exists(S6M):
    d6 = fitz.open(S6M); f6 = norm(d6[0].get_text()); d6.close()
    for _, r in DL.iterrows():
        exp6m += ["\u03c1 %.3f \u2192 %.3f" % (r.rho_raw, r.rho_adj), "perm P %.4f" % r.perm_P_adj,
                  "95%% CI %.2f\u2013%.2f" % (r.ci_lo_adj, r.ci_hi_adj),
                  "n=%s Oligo" % format(int(r.n_oligo), ","), "%s OPC" % format(int(r.n_opc), ","),
                  "chromatin \u03941 (Oligo \u2212 OPC)", "expression \u03941 (Oligo \u2212 OPC)"]
    for t in TECH:
        TT = DT[DT.tech == t]
        for key, nm in MARKS:
            q = TT[TT.trait.str.contains(key, regex=False)].iloc[0]
            exp6m.append("%s #%d\u2192#%d" % (nm, int(q.rank_rna_adj), int(q.rank_atac_adj)))
    miss6m = [e for e in exp6m if norm(e) not in f6]
    # S6m passes the size gate too (the old three-line provenance at 5.0 pt moved to the caption)
    d6b = fitz.open(S6M); sz6 = [round(sp["size"], 2) for b in d6b[0].get_text("dict")["blocks"]
                                 for l in b.get("lines", []) for sp in l["spans"] if sp["text"].strip()]
    d6b.close()
    under6m = [z for z in sz6 if z < FMIN]
else:
    print("!! \u9644\u56fe S6m \u672a\u751f\u6210")
    under6m = []

# S6m's bottom three-line provenance (5.0 pt, fails the size gate) was deleted from the face -> the caption must carry it
cap += ["HT-scCAT-seq and 10x Multiome data were obtained from the HT-scCAT-seq study",
        "these are independent experiments, from separate animals, and are not the same sample as the dsciCAT-seq libraries",
        "Cell types follow the harmonised 33-class annotation of Fig 2 (WNN space)",
        "per-channel labels were not taken from the source studies",
        "Expression channel = measured Seurat log-normalised data",
        "chromatin channel = ArchR gene score",
        "Ranks shown are domain-adjusted ranks within 74 traits"]

assert os.path.exists(LEG), "\u7f3a\u56fe\u6ce8\u8349\u7a3f %s\uff08\u4e0d\u80fd\u9759\u9ed8\u901a\u8fc7\uff09" % LEG
_lg = norm(open(LEG, encoding="utf-8").read())
misscap = [c for c in cap if norm(c) not in _lg]

coll = {}
for f in [PDF, S6M]:
    if os.path.exists(f):
        coll[os.path.basename(f)] = overlaps(f)

lines = ["Fig 6 v21 acceptance  %s" % pd.Timestamp.now(), "PDF: %s" % PDF,
         "sizes min/max = %.2f / %.2f pt; fonts = %s" % (min(sizes), max(sizes), ", ".join(sorted(fonts))),
         "size hard gate >=%.2f pt: %d spans broken%s" % (FMIN, len(undersz),
             "" if not undersz else " -> " + " | ".join("%r@%.2f" % (u["text"].strip()[:18], u["size"]) for u in undersz[:6])),
         "main figure expects %d items, missing %d" % (len(exp), len(miss))]
if miss:
    lines.append("MISS(main figure): " + " | ".join(repr(m) for m in miss))
lines.append("D-region text %d spans, off-black or off-size %d" % (len(dsty), len(bad)))
lines.append("points covered by C's trait labels = %d (per-point containment in window rects at render time; see fig6C_label_occlusion.csv)" % int(occ.covered_points.sum()))
lines.append("D in-cell ranks: %d colour digits + %d white-underlay digits (NUMMODE=%s, expected %d / %d), colours disagreeing with their row trait %d"
             % (len(dig_c), len(dig_w), NUMMODE, NDIGC, NDIGW, len(misnum)))
if misnum:
    lines.append("BAD(rank colours): " + " | ".join("%s@%s=%s" % m for m in misnum[:6]))
lines.append("white boxes left in D %d (must be 0)" % len(wbox))
lines.append("named-trait colours C=%s  D=%s  mismatch %s" % (
    {k: "#%06x" % v for k, v in sorted(cspan.items())},
    {k: "#%06x" % v for k, v in sorted(dcol.items())}, mismatch or "none"))
if bad:
    lines.append("BAD(D style): " + " | ".join("%r color=%d size=%.2f" % b for b in bad[:8]))
lines.append("caption carries %d sentences (deleted from the face), %d not in the caption" % (len(cap), len(misscap)))
if misscap:
    lines.append("MISS(caption): " + " | ".join(repr(m) for m in misscap))
lines.append("S6m expects %d items, missing %d; spans below %.2f pt: %d" % (len(exp6m), len(miss6m), FMIN, len(under6m)))
if miss6m:
    lines.append("MISS(S6m): " + " | ".join(repr(m) for m in miss6m))
for k, v in coll.items():
    lines.append("text collisions %s: %d%s" % (k, len(v), "" if not v else " -> " + " | ".join("%r\u00d7%r" % (x[2], x[3]) for x in v[:6])))
tab = pd.DataFrame([{"key": k, "value": v, "source": s} for k, v, s in keys])
tab.to_csv(V21 + "/fig6_all_numbers_v21.csv", index=False)
lines.append("fig6_all_numbers_v21.csv: %d rows (every face digit%s -> a table row)" % (len(tab), "" if NUMMODE != "none" else " + the 18 ranks in the caption"))
ok = all([not miss, not misscap, not bad, not mismatch, not misnum, not under6m, not undersz,
          len(dig_c) == NDIGC, len(dig_w) == NDIGW, not wbox, occ_bad.empty, not miss6m,
          all(not v for v in coll.values()), min(sizes) >= FMIN,
          all("Arial" in f for f in fonts)])
lines.append("verdict: %s" % ("PASS" if ok else "FAIL"))
out = "\n".join(lines)
print(out)
open(V21 + "/fig6_v21_audit.txt", "w").write(out + "\n")
sys.exit(0 if ok else 1)
