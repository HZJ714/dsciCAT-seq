# -*- coding: utf-8 -*-
"""反向溯源门 ALLOW 白名单的"饱和度"诊断（只读）。

verify_S6_v10.py 把 SRC 里每张表的**每一列每一个值**的字符串都塞进 ALLOW，
其中含逐核表 cell_level_table.csv、10,000 次置换抽样 figS6v10_permnull_o1.csv、
bootstrap 明细 figS6v10_bootstrap_units.csv、74x6x2 全表这类**行级/明细表**。
后果：随便编一个 3~4 位小数的 rho/P 都大概率在白名单里撞得到，
"图面每个数字都能反查出处"这条纪律实际是空转的。

三件事：
 1) ALLOW 有多少数字、分别来自哪些表；
 2) 随机编造的 rho/P 通过率（只看汇总表 vs 含明细表）；
 3) 门的 1e-9 绝对容差造成的"极小数自动放行"。
"""
import os
import re
import sys

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from support_common import R

NUM = re.compile(r"[-\u2212+]?\d[\d,]*(?:\.\d+)?(?:[eE][-+]?\d+)?%?x?")


def norm(t):
    t = t.replace("\u2212", "-").replace(",", "").rstrip("%xX")
    return t[1:] if t.startswith(("+", "-")) else t


def toks(txt):
    return set(norm(m) for m in NUM.findall(txt))


SUP2 = R + "/SUPPORT_v2"
V8, TAB2 = SUP2 + "/tables_v8", SUP2 + "/tables"
SUMMARY, DETAIL = [], []
for f in sorted(os.listdir(V8)):
    if f.endswith(".csv"):
        (DETAIL if os.path.getsize(os.path.join(V8, f)) > 60000 else SUMMARY).append(os.path.join(V8, f))
SUMMARY += [os.path.join(TAB2, "Table_S6.csv"), os.path.join(TAB2, "Table_S7.csv")]
for f in ["D_crosschemistry/tables_log1p_20261104/fig6D_delta1_74traits_per_tech.csv",
          "D_crosschemistry/tables_log1p_20261104/fig6D_crosschem_ladder.csv",
          "O1_own_multiome_paired/tables/T4b_cross_channel_rho_variants.csv",
          "O1_own_multiome_paired/tables/T5_split_half_reliability.csv",
          "O1_own_multiome_paired/tables/T2_MS_per_nucleus_rna_vs_atac.csv",
          "O1_own_multiome_paired/tables/L1_by_library.csv",
          "M2b_authoranno_v1/tables/M2b_per_mouse.csv"]:
    SUMMARY.append(os.path.join(R, f))
SUMMARY.append(WORK + "/FY/oligo/scdrs_ms_out/ms_locus_overlap_table.csv")
DETAIL.append(os.path.join(R, "O1_own_multiome_paired/tables/cell_level_table.csv"))


def load(paths):
    vals, src_of = set(), {}
    for p in paths:
        if not os.path.exists(p):
            continue
        try:
            df = pd.read_csv(p)
        except Exception as e:
            print("  [warn]", p, e)
            continue
        for c in df.columns:
            for v in df[c].astype(str):
                for t in toks(v):
                    try:
                        f = float(t)
                    except ValueError:
                        continue
                    vals.add(f)
                    src_of.setdefault(f, set()).add(os.path.basename(p))
    return vals, src_of


VS, SS = load(SUMMARY)
VD, SD = load(DETAIL)
print("汇总表 %d 张 -> %d 个不同数值；明细表 %d 张 -> %d 个不同数值"
      % (len(SUMMARY), len(VS), len(DETAIL), len(VD)))
BIG = sorted(set(src for s in SD.values() for src in s))
print("明细表：", ", ".join(BIG))


def index(vals, decs=range(0, 8)):
    """按小数位建 round(a,dec) 索引，模拟门的宽松判定（不含百分比写法）"""
    return {d: {round(a, d) for a in vals} for d in decs}


IS, ID = index(VS), index(VD)
IALL = {d: IS[d] | ID[d] for d in IS}


def passes(f, IDX, dec):
    return abs(round(f, dec) - f) < 1e-9 and f in IDX[dec] if False else (round(f, dec) in IDX[min(dec, 7)])


rng = np.random.default_rng(7)
print("\n== 编造数字的通过率（每种 500 个）==")
cases = [("rho 3 位小数 U(0,1)", [round(x, 3) for x in rng.random(500)], 3),
         ("P  4 位小数 U(0,0.01)", [round(x, 4) for x in rng.random(500) * 0.01], 4),
         ("rho 2 位小数 U(0,1)", [round(x, 2) for x in rng.random(500)], 2)]
for lab, gen, dec in cases:
    ps = sum(round(v, dec) in IS[dec] for v in gen)
    pd_ = sum(round(v, dec) in IALL[dec] for v in gen)
    print("  %-24s 仅汇总表 %5.1f%%    汇总表+明细表 %5.1f%%" % (lab, 100.0 * ps / len(gen), 100.0 * pd_ / len(gen)))

print("\n== 1e-9 绝对容差：极小数自动放行 ==")
for v in [1e-4, 1e-6, 1e-9, 1e-10, 3.3e-31, 9.9e-77]:
    inS = any(abs(round(a, 6) - v) < 1e-9 for a in VS)
    print("   P = %-9g 在汇总表里能按 1e-9 容差找到出处：%s" % (v, inS))

print("\n== 具体编造值撞在哪张表（宽松语义）==")
for t in ["0.817", "0.5123", "0.9999", "0.123", "0.819"]:
    f = float(t)
    dec = len(t.split(".")[1])
    hs = sorted({s for a in VS if round(a, dec) == f for s in SS.get(a, [])})[:3]
    hd = sorted({s for a in VD if round(a, dec) == f for s in SD.get(a, [])})[:3]
    print("   %-8s 汇总: %-42s 明细: %s" % (t, ", ".join(hs) or "-", ", ".join(hd) or "-"))
