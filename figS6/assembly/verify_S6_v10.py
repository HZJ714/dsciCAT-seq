# -*- coding: utf-8 -*-
"""S6 **v10** 反向溯源门（严格索引模式，2026-11-05 起）：图面每个数字必须能在
figS6v10_number_index.csv 注册（逐字 printed 或按该 token 有效位对齐 value），
或命中 STRUCT 结构值白名单（坐标刻度/染色体号/GSE 号）；否则判 UNREGISTERED/NO-SOURCE 失败。
历史说明：v8 版用的是"SRC 全值池 ALLOW"，被证明会饱和放行（见 S6v10_gate_allow_probe.py）；
封死"旧血统数字回流"（log1p 修正前的 0.726/0.763/0.518、ATAC 跨臂 0.61/0.53、
数字化 n=1,000 置换、跨模态 1/7 幅度比、旧 G 格 Δ1 0.261/0.430/0.297）。
用法：python verify_S6_v10.py [pdf ...]
"""
import os, re, sys

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from support_common import R, say
import pandas as pd
import fitz

SUP2 = R + "/SUPPORT_v2"
PDF2, TAB2, V8 = SUP2 + "/pdf", SUP2 + "/tables", SUP2 + "/tables_v8"
NUM = re.compile(r"[-\u2212+]?\d[\d,]*(?:\.\d+)?(?:[eE][-+]?\d+)?%?x?")

def norm(t):
    t = t.replace("\u2212", "-").replace(",", "").rstrip("%xX")
    if t.startswith("+"): t = t[1:]
    if t.startswith("-"): return t[1:]
    return t

def toks(txt): return set(norm(m) for m in NUM.findall(txt))

def page_text(path):
    d = fitz.open(path)
    s = " ".join(" ".join(p.get_text().split()) for p in d)
    d.close(); return s

ALLOW = set()
SRC = []
for f in sorted(os.listdir(V8)):
    if f.endswith(".csv"): SRC.append(os.path.join(V8, f))
for f in ["Table_S6.csv", "Table_S7.csv"]:
    SRC.append(os.path.join(TAB2, f))
for f in ["D_crosschemistry/tables_log1p_20261104/fig6D_delta1_74traits_per_tech.csv",
          "D_crosschemistry/tables_log1p_20261104/fig6D_crosschem_ladder.csv",
          "O1_own_multiome_paired/tables/T4b_cross_channel_rho_variants.csv",
          "O1_own_multiome_paired/tables/T5_split_half_reliability.csv",
          "O1_own_multiome_paired/tables/T2_MS_per_nucleus_rna_vs_atac.csv",
          "O1_own_multiome_paired/tables/L1_by_library.csv",
          "O1_own_multiome_paired/tables/cell_level_table.csv",
          "M2b_authoranno_v1/tables/M2b_per_mouse.csv"]:
    SRC.append(os.path.join(R, f))
SRC.append(WORK + "/FY/oligo/scdrs_ms_out/ms_locus_overlap_table.csv")
for p in SRC:
    try:
        df = pd.read_csv(p)
    except Exception as e:
        print("  [warn] 读表失败", p, e); continue
    for c in df.columns:
        for v in df[c].astype(str):
            ALLOW |= toks(v)
STRUCT = {
    "0": "轴原点", "1": "刻度", "2": "刻度", "3": "刻度", "4": "刻度", "5": "刻度", "6": "刻度",
    "7": "刻度", "8": "刻度", "9": "刻度", "10": "刻度", "12": "刻度", "14": "刻度", "16": "刻度",
    "18": "逐鼠 18", "37": "秩轴中点", "74": "性状数", "331": "平台基因数", "6": "状态数",
    "0.0": "刻度", "0.2": "刻度", "0.4": "刻度", "0.6": "预注册判据/刻度", "0.8": "刻度",
    "1.0": "刻度", "1.5": "Delta1 截断刻度", "-1.5": "同上", "0.5": "AUC 机会线",
    "0.75": "刻度", "-0.4": "刻度", "-0.8": "刻度", "-0.2": "刻度", "-1": "刻度",
    "118.4": "位点窗口 Mb 刻度", "119.4": "同上", "119.2": "同上", "120.2": "同上",
    "53.7": "同上", "54.7": "同上", "56.5": "同上", "57.5": "同上",
    "20261104": "seed（角标/表）", "3/3": "三臂达判据", "18/18": "逐鼠同号",
    "500": "窗口半宽 kb", "15": "hexbin gridsize", "25": "hexbin gridsize", "40": "hexbin gridsize", "3000": "零分布直方图刻度",
    "2.5": "bootstrap 百分位", "97.5": "bootstrap 百分位", "2.0": "色标刻度", "1.5": "色标刻度", "271": "距离 kb（A 格）", "76": "同上", "103": "同上", "241": "同上",
    "1000": "group assoc n_ctrl（ttl 文字带）", "5,962": "核数", "900": "抽样显示数",
    # 作者 2026-11-05 严格门上线时补：A 格染色体号、F 格 rho 轴刻度、文库名里的 50/100
    "20": "染色体号（chr20，A 格）", "11": "染色体号（chr11，A 格）",
    "0.00": "F 格 rho 轴刻度", "0.25": "F 格 rho 轴刻度", "0.50": "F 格 rho 轴刻度",
    "1.00": "F 格 rho 轴刻度", "50": "文库名 seq-50k", "100": "文库名 seq-100k",
    "250589": "外部数据 GSE 登录号（mouse EAE，F 格梯级，Fig 7 分析链（Table S6/S7）同源）",
}
ALLOW |= {norm(k) for k in STRUCT}

CELLS = ["S6a_loci_tracks_v10.pdf", "S6b_permutation_strip_v10.pdf", "S6c_umap_bins_v10.pdf",
         "S6d_subtype_profiles_v10.pdf", "S6e_ladder_onepanel_v10.pdf", "S6f_unit_violins_v10.pdf",
         "S6g_boxes_matrix_v10.pdf"]
DRAFT = "figS6_v10_draft.pdf"
BANNED_ON_FIGURE = ["0.726", "0.763", "0.518", "0.261", "0.430", "0.297", "1/7 of", "ratio median",
                    "n = 1,000", "digitis", "P = 0.26", "rho +0.61", "rho +0.53", "rho +0.47",
                    # 作者 2026-11-05：内部叫法与判据文字一律不上图
                    "(own)", "D-own", "FY13", "FY14", "rho_adj", "PASS"]
MUST = {"S6a_loci_tracks_v10.pdf": ["GWAS lead-SNP", "5e-26", "PHLDB1", "PLLP"],
        "S6b_permutation_strip_v10.pdf": ["38/331", "11.5", "10,000", "0.277", "331"],
        "S6c_umap_bins_v10.pdf": ["single cells (no smoothing)", "gridsize 15",
                                  "gridsize 25 (Fig 6B)", "gridsize 40", "MOL6", "COP/NFOL"],
        "S6d_subtype_profiles_v10.pdf": ["MS", "SCZ", "HEIGHT", "RNA channel", "ATAC channel",
                                         "median norm_score", "COP/NFOL", "MOL6"],
        "S6e_ladder_onepanel_v10.pdf": ["+0.819", "0.0001", "0.227", "10,000", "0 / 10,000",
                                        "observed", "null 97.5", "trait pairing permuted"],
        "S6f_unit_violins_v10.pdf": ["0.079", "0.812", "0.729", "0.472", "0.819", "0.424",
                                     "0.920", "descriptive", "1,000", "seq-50k", "seq-100k",
                                     "dsciCAT", "GSE250589"],
        "S6g_boxes_matrix_v10.pdf": ["P = 1.7e-33", "P = 1.8e-18", "P = 6.8e-06", "P = 3.1e-04",
                                     "P = 0.0022", "P = 0.75", "OPC", "Oligo",
                                     "per-cell norm_score", "HT-scCAT-seq", "dsciCAT-seq"]}

ALLOWF = set()
for _v in ALLOW:
    try: ALLOWF.add(float(_v))
    except ValueError: pass

# 严格门（作者 2026-11-05，崩溃会话 01a1aab1/01a1aac1 未竟事项）：
#   旧逻辑把 SRC 全部数值（含 cell_level_table 这种几十万行的原始表）摊成 ALLOW 再按 1e-9 容差匹配，
#   实测任何编造数字几乎都能“找到出处”（S6v10_gate_allow_probe.py：0.817/0.5123/0.123 全命中）。
#   新逻辑：图面数字必须命中 **figS6v10_number_index.csv 的注册条目**（printed 逐字 或 按 token
#   有效位数舍入后等于注册的 value），或 STRUCT 结构值（刻度/坐标轴/布局）。其余：
#     在台账里碰巧有 → UNREGISTERED（严格模式 FAIL；S6_GATE_LENIENT=1 降为警告）
#     台账也没有     → NO-SOURCE（FAIL）
import math
_SIGRE = re.compile(r"(\d+)(?:\.(\d+))?(?:[eE][-+]?\d+)?$")
def _sig(t):
    m = _SIGRE.match(t)
    if not m: return None
    digits = (m.group(1) + (m.group(2) or "")).lstrip("0")
    return len(digits) or 1
def _rr(x, sig):
    if x == 0: return 0.0
    e = math.floor(math.log10(abs(x)))
    return round(x, sig - 1 - e)

IDX = pd.read_csv(V8 + "/figS6v10_number_index.csv")
IDXT = set()
IDXV = []
for _r in IDX.itertuples():
    IDXT |= toks(str(_r.printed))
    try:
        IDXV.append(float(str(_r.value).replace(",", "")))
    except ValueError: pass
STRUCTN = {norm(k) for k in STRUCT}

def classify(t):
    if t in IDXT: return "registered"
    try: f = float(t)
    except ValueError: return "no-source"
    sig = _sig(t)
    if sig is not None:
        for a in IDXV:
            if _rr(a, sig) == f: return "registered"
    if t in STRUCTN or norm(t) in STRUCTN: return "structural"
    if t in ALLOW: return "unregistered"
    try:
        fv = float(t)
        dec = len(t.split(".")[1]) if "." in t and "e" not in t.lower() else 0
        for a in ALLOWF:
            if abs(round(a, dec) - fv) < 1e-9 or abs(round(100 * a, dec) - fv) < 1e-9: return "unregistered"
    except ValueError: pass
    return "no-source"

LENIENT = os.environ.get("S6_GATE_LENIENT") == "1"
def sourced(t):
    return classify(t) in ("registered", "structural") or (LENIENT and classify(t) == "unregistered")

bad = 0
for f in (sys.argv[1:] or [os.path.join(PDF2, c) for c in CELLS + [DRAFT]]):
    name = os.path.basename(f)
    txt = page_text(f)
    ts = toks(txt)
    cls = {t: classify(t) for t in ts if t}
    miss = sorted(t for t, c in cls.items() if c == "no-source")
    unreg = sorted(t for t, c in cls.items() if c == "unregistered")
    banned = [b for b in BANNED_ON_FIGURE if b in txt.replace(",", ".")]
    must = MUST.get(name, [])
    absent = [m for m in must if norm(m) not in ts and m not in txt]
    ok_unreg = unreg and LENIENT
    status = "OK" if not miss and not absent and not banned and (not unreg or ok_unreg) else "FAIL"
    if status != "OK": bad += 1
    print("%-32s 数字 %3d 个，无出处 %d，未注册 %d，缺关键 %d，旧血统回流 %d  %s"
          % (name[:32], len(ts), len(miss), len(unreg), len(absent), len(banned), status))
    for m in miss[:12]: print("      NO-SOURCE     %r" % m)
    for u in unreg[:12]: print("      %s  %r" % ("WARN-UNREG" if ok_unreg else "UNREGISTERED", u))
    for a in absent[:12]: print("      MISSING    %r" % a)
    for b in banned[:24]: print("      BANNED-OLD-VALUE  %s" % b)
tot = (len(CELLS) + 1) if not sys.argv[1:] else (len(sys.argv) - 1)
print("\n%d/%d 通过 v10 反向溯源门" % (tot - bad, tot))
sys.exit(1 if bad else 0)
