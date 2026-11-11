# -*- coding: utf-8 -*-
"""Table_S6 台账与 v10 图面（2026-11-05 16:00 重建后）对齐。

背景：v10 在 15:24 那轮按作者意见把 G 的 rho_adj/PASS 行、Δ1 撤下图面，16:00 又把 G
右下角换成精确 P。Table_S6.csv 的 in_figure 列还写着"某格 (printed)"，属于**过期上图声明**。

本脚本做三件事（统计量数值一律不改，只动 in_figure；追加 G 的检验块）：
 1) 核验每条"上图形"声明：声明里那个数字若已不在所声称面板的 PDF 文字层 → 降级
    "Table only (…)"。
 2) **迁移**：声明的面板号是旧血统（"Fig S6 v8 panel X"）或写错面板，但该值确实由
    figS6v10_number_index.csv 登记在某个 v10 面板上 → 改指到那个面板并写明来源可追。
    只认台账（值 + 出处文件都对得上）才算同一条统计量，避免"0.6 / 0.75 撞车"式误迁移。
 3) 追加 G 的 6 行 Mann-Whitney U（右下角图面就是这些 P，写法逐字一致）。

幂等：一切从 tables/Table_S6_v10face_backup.csv（未改动前的 290 行）重放，
所以可以反复跑；改动明细同时写到 audit_20261105/table_s6_face_retags.csv。
"""
import os
import re
import sys

import pandas as pd

import fitz

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

TAB = WORK + "/figures/routeA_out/SUPPORT_v2/tables"
V8 = WORK + "/figures/routeA_out/SUPPORT_v2/tables_v8"
PDF = WORK + "/figures/routeA_out/SUPPORT_v2/pdf"
AUD = WORK + "/figures/routeA_out/SUPPORT_v2/audit_20261105"
BASE = TAB + "/Table_S6_v10face_backup.csv"
OUT = TAB + "/Table_S6.csv"
CELL = {  # 面板字母 -> v10 cell pdf
    "A": "S6a_loci_tracks_v10.pdf", "B": "S6b_permutation_strip_v10.pdf",
    "C": "S6c_umap_bins_v10.pdf", "D": "S6d_subtype_profiles_v10.pdf",
    "E": "S6e_ladder_onepanel_v10.pdf", "F": "S6f_unit_violins_v10.pdf",
    "G": "S6g_boxes_matrix_v10.pdf"}
LETTER = {v: k for k, v in CELL.items()}
TXT = {k: " ".join(" ".join(p.get_text().split()) for p in fitz.open(os.path.join(PDF, v)))
       for k, v in CELL.items()}
NUM = re.compile(r"[-\u2212+]?\d[\d,]*(?:\.\d+)?(?:[eE][-+]?\d+)?")


def fmt_match(val, t):
    """图面写法 token t 是否**就是** val 的合法舍入（按 t 自身的位数格式化后逐字比，
    不用绝对容差——否则 1e-33 这种极小数会集体蒙混）。支持百分号写法。"""
    low = t.lower()
    if "e" in low:
        mant, _, exp = low.partition("e")
        k = len(mant.split(".")[1]) if "." in mant else 0
        s = "%.*e" % (k, float(val))
        m2, _, e2 = s.partition("e")
        try:
            return m2 == mant and int(e2) == int(exp)
        except ValueError:
            return False
    dec = len(low.split(".")[1]) if "." in low else 0
    return ("%.*f" % (dec, float(val)) == low) or ("%.*f" % (dec, 100.0 * float(val)) == low)


def on_face(val, panel):
    """值 val 是否以某个写法出现在 panel 的文字层。"""
    for m in NUM.findall(TXT[panel]):
        t = m.replace("\u2212", "-").replace(",", "").lstrip("+").rstrip("%xX")
        try:
            float(t)
        except ValueError:
            continue
        if fmt_match(val, t):
            return True
    return False


# ---- 台账（唯一可信的"这条统计量印在哪个格"）----
IDX = pd.read_csv(V8 + "/figS6v10_number_index.csv")
IDX["panel"] = IDX.pdf.map(lambda p: LETTER.get(p, "?"))


def ledger_hit(val):
    """在 v10 上图数字台账里找同一个统计量：**相对**容差对齐数值（不用 round(…,4)，
    否则 1.7e-33 和 3.1e-04 都会被归成 0.0），并且它确实印在对应格面上。"""
    import math
    for _, r in IDX.iterrows():
        try:
            if not math.isclose(float(val), float(r.value), rel_tol=1e-4, abs_tol=0.0):
                continue
        except (TypeError, ValueError):
            continue
        if on_face(float(r.value), r.panel):
            return r.panel, r.key, r.source
    return None


T = pd.read_csv(BASE)
print("从备份重放：Table_S6 %d 行（未改动前状态）" % len(T))

PANEL_RE = re.compile(r"panel ([A-G])\b")
S6LET_RE = re.compile(r"\bS6([a-g])\b")
DEMOTE = re.compile(r"^Table only \(.*; was: (.*)\)$")


def claimed_panels(txt):
    ks = set(m.group(1).upper() for m in PANEL_RE.finditer(txt))
    ks |= set(m.group(1).upper() for m in S6LET_RE.finditer(txt))
    return ks


changes = []
for i, r in T.iterrows():
    inf = str(r.in_figure)
    m = DEMOTE.match(inf)                      # 上一轮的降级写法 -> 还原成原始声明再判
    if m:
        inf = m.group(1)
    if inf.strip().lower().startswith(("table only", "no", "nan")):
        continue
    ks = claimed_panels(inf)
    if not ks:
        continue
    v = pd.to_numeric(str(r.value).replace(",", "").replace("%", "").lstrip("+"), errors="coerce")
    if pd.isna(v):
        continue
    v = float(v)
    if any(on_face(v, k) for k in ks) and "v8" not in inf:      # 声明成立，原样保留
        continue
    hit = ledger_hit(r.value)                                    # 声明过期：能迁到别的格吗
    if hit:
        panel, key, src = hit
        new = ("Fig S6 v10 panel %s (%s; value printed there; "
               "was: %s)" % (panel, key, inf))
        T.at[i, "in_figure"] = new
        changes.append((i, r.check, str(r.statistic), r.value, inf, new, "migrate", src))
    else:
        new = ("Table only (value not on the S6 v10 face after the 2026-11-05 detail round; "
               "was: %s)" % inf)
        T.at[i, "in_figure"] = new
        changes.append((i, r.check, str(r.statistic), r.value, inf, new, "demote", ""))

n_mig = sum(1 for c in changes if c[6] == "migrate")
n_dem = sum(1 for c in changes if c[6] == "demote")
print("\n== 改了 %d 条上图声明：迁移 %d、降级 %d ==" % (len(changes), n_mig, n_dem))
for i, chk, st, val, old, new, kind, src in changes:
    print(" [%s] row %-4d %-38s %-46s = %-10s\n        旧 %s\n        新 %s"
          % (kind, i, chk[:38], st[:46], val, old[:96], new[:120]))

# ---- G 的 Mann-Whitney 块（图面右下角印的就是这些 p_label）----
M = pd.read_csv(V8 + "/figS6v10_G_opc_oligo_mwu.csv", dtype={"p_fig": str, "p_label": str})
add = []
for _, m in M.iterrows():
    printed = m.p_label.replace("P = ", "").replace("P < ", "")
    onG = on_face(float(m.p_two_sided), "G") or (" %s" % m.p_label) in TXT["G"]
    add.append(dict(
        check="G OPC vs Oligo per-nucleus contrast (Mann-Whitney U)",
        statistic="%s, %s channel: MWU P (OPC vs Oligo norm_score)" % (m.tech, "RNA" if m.channel == "rna" else "ATAC"),
        value=m.p_label, ci95_or_IQR="median diff (Oligo-OPC) %+0.3f; U %+0.3g" % (m.median_diff, m.mwu_U),
        reference="two-sided Mann-Whitney U, no permutation, deterministic",
        rejects_zero="n.a. (location test)", meets_effect="n.a.",
        n="OPC %d / Oligo %d" % (m.n_opc, m.n_oligo),
        source="tables_v8/figS6v10_G_opc_oligo_mwu.csv", source_column="p_two_sided",
        in_figure=("Fig S6 v10 panel G (printed verbatim '%s' in the facet)" % m.p_label
                   if onG else "Table only (label '%s' not found on the G face)" % m.p_label),
        in_main_text="no"))
T = pd.concat([T, pd.DataFrame(add)], ignore_index=True)
print("\n== 追加 G MWU 块 %d 行；Table_S6 现 %d 行 ==" % (len(add), len(T)))
for a in add:
    print("   %-52s %-12s %s" % (a["statistic"][:52], a["value"], a["in_figure"][:70]))

T.to_csv(OUT, index=False)
os.makedirs(AUD, exist_ok=True)
pd.DataFrame([dict(row=c[0], check=c[1], statistic=c[2], value=c[3], old=c[4], new=c[5],
                   action=c[6], ledger_source=c[7]) for c in changes]).to_csv(
    AUD + "/table_s6_face_retags.csv", index=False)
print("\n已写回 %s\n改动明细 %s/table_s6_face_retags.csv" % (OUT, AUD))
