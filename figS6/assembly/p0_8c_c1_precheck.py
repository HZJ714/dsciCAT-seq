import os
# -*- coding: utf-8 -*-
"""p0_8c_c1_precheck.py — Fig6C 定案 C1 的渲染前置检查（NARRATIVE_DESIGN §5 P0-a）。

三步：
 1) 38 风险基因（Table 2，人源符号）→ scDRS 内置同源映射 → RNA/ATAC 矩阵命中数（溯源透明性）
 2) P0-a：L2「OL 亚型特异表达」在 38 行上的命中数（三档定义都给），按 §3 双边规则判定能否占列
 3) 出 Ca/Cb 所需的全部数值表（每亚型 RNA 均值/检出率、ATAC 峰值汇总均值、类水平效应量）

内存：RNA 用分块累加（直接布尔切片会稠密化到 18 GB+，已踩坑）；ATAC 逐列取。
统计口径提醒：RNA 的 X 已是 log2(CP10k+1)，不得再算 CP10k。
"""
import re
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, scipy.stats as st, os
import warnings; warnings.filterwarnings("ignore")
from scdrs.util import load_homolog_mapping

W = os.environ.get('DSCI_WORK_ROOT', '')
U = os.environ.get('DSCI_WORK_ROOT', '')
OUT = f"{W}/figures/routeA_out/P0_8C/tables"; os.makedirs(OUT, exist_ok=True)
SUBS = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]      # 左→右按成熟度轴
ENA, ENB = "MOL3/6", "MOL6"                                        # Fig6A group_q<0.05 的富集亚型
g38 = pd.read_csv(f"{U}/figures/CIRCLE38_20261116/gene38.csv")
def _p(s):                                   # "5 × 10−8" -> 5e-8
    m = re.match(r"\s*(\d+(?:\.\d+)?)\s*[×*]\s*10\s*[−-](\d+)\s*", str(s))
    return float(m.group(1)) * 10 ** (-int(m.group(2))) if m else np.nan
g38["p_num"] = g38.p_value.map(_p)
print("Table 2 行数 = %d | 唯一基因 = %d | p 值域 %.1e–%.1e" % (len(g38), g38.gene.nunique(), g38.p_num.min(), g38.p_num.max()))

# ---- 1) 同源映射 + 矩阵命中
H2M = load_homolog_mapping("hsapiens", "mmusculus")       # dict: 人源符号 -> 鼠源符号列表
rn_all = set(map(str, ad.read_h5ad(f"{W}/FY/oligo/3k_dsci_oligo_subset.h5ad", backed="r").var_names))
def as_list(v):                                    # 坑：scDRS 同源表 value 是字符串，list() 会拆成单字符
    if v is None: return []
    return list(v) if isinstance(v, (list, tuple, set)) else [str(v)]
def pick(g):
    lst = as_list(H2M.get(g))
    hit = [m for m in lst if m in rn_all]
    if hit: return hit[0], "homolog+matrix"
    if lst: return lst[0], "homolog only"
    return g, "passthrough"                              # 无同源条目：按原名试（鼠常用首字母大写）
g38[["mouse", "map_src"]] = g38.gene.apply(lambda g: pd.Series(pick(g)))
print("同源表有条目 %d/38 | 其中命中矩阵 %d | 仅同源无矩阵 %d | 无条目原名试 %d"
      % (sum(len(H2M.get(g, [])) > 0 for g in g38.gene), (g38.map_src == "homolog+matrix").sum(),
         (g38.map_src == "homolog only").sum(), (g38.map_src == "passthrough").sum()))
print("  无同源条目:", list(g38.loc[g38.map_src == "passthrough", "gene"]))

ar = ad.read_h5ad(f"{W}/FY/oligo/3k_dsci_oligo_subset.h5ad", backed="r")
key = ar.obs["cell_type"].astype(str).values
rows = {s: np.flatnonzero(key == s) for s in SUBS}
rn = set(map(str, ar.var_names))
aa = ad.read_h5ad(f"{W}/FY/oligo/fig6_rerun_20261019/dualomics_probe/genescore_oligo_5977.h5ad", backed="r")
an = set(map(str, aa.var_names))
g38["in_rna"] = g38.mouse.isin(rn); g38["in_atac"] = g38.mouse.isin(an)
g38.loc[~g38.in_rna, ["mouse","map_src"]].to_csv("/dev/null")
print("RNA 矩阵命中 %d/38 | ATAC 矩阵命中 %d/38 | 两者都命中 %d/38"
      % (g38.in_rna.sum(), g38.in_atac.sum(), (g38.in_rna & g38.in_atac).sum()))

# ---- 2) 分块累加 RNA 每亚型均值与检出率（只留需要的列，降内存）
need = list(dict.fromkeys(g38.loc[g38.in_rna, "mouse"]))
cols = np.array([list(map(str, ar.var_names)).index(g) for g in need])
Xs = sp.csr_matrix(ar.X[:, cols])
Xg = np.asarray(Xs.todense(), dtype=np.float32)                    # 5962×n，很小
del Xs
Mn, Pc = np.zeros((len(need), len(SUBS))), np.zeros((len(need), len(SUBS)))
for j, s in enumerate(SUBS):
    b = Xg[rows[s]]
    Mn[:, j] = b.mean(0); Pc[:, j] = (b > 0).mean(0) * 100
sub = pd.DataFrame({"gene": g38.loc[g38.in_rna, "gene"].values, "mouse": need,
                    **{f"rna_{s}": Mn[:, j] for j, s in enumerate(SUBS)},
                    **{f"pct_{s}": Pc[:, j] for j, s in enumerate(SUBS)}})
# ATAC 每亚型均值（逐列取，38 列）
al = aa.obs["cell_type"].astype(str).values
av = np.array(list(map(str, aa.var_names)))
ata = {}
for g in sub.mouse.values:
    if g in an:
        x = np.asarray(aa.X[:, av == g].todense()).ravel()
        ata[g] = [np.mean(x[al == s]) for s in SUBS]
print("  ATAC 逐列均值完成：%d 基因" % len(ata))
for j, s in enumerate(SUBS):
    sub[f"atac_{s}"] = sub.mouse.map(lambda g: ata[g][j] if g in ata else np.nan)

# ---- 3) P0-a：L2 三档命中数（在 RNA 命中的 38 子集上算）
mx = sub[[f"rna_{s}" for s in SUBS]].values
pmx = sub[[f"pct_{s}" for s in SUBS]].values
top = mx.argmax(1)
second = np.sort(mx, 1)[:, -2]
lfc_top_vs_rest = mx[np.arange(len(sub)), top] - np.delete(mx, top, axis=1).mean(1)
sub["top_subtype"] = [SUBS[i] for i in top]
L2a = (pmx.max(1) >= 20)                                     # 在 OL 里可检出
L2b = L2a & (lfc_top_vs_rest >= 1.0)                        # 最强亚型比其余均值高 2 倍
L2c = (pmx.max(1) >= 30) & (lfc_top_vs_rest >= 1.585)       # 严格 marker 口径
n = len(sub)
for nm, v in [("L2a 可检出(pct>=20%)", L2a), ("L2b 亚型特异(log2FC>=1)", L2b), ("L2c 严格marker(pct>=30&,lfc>=1.585)", L2c)]:
    print("  %-40s %2d/%d = %.1f%%" % (nm, v.sum(), n, v.mean() * 100))
print("  富集亚型(MOL3/6+MOL6)为 top 的行: %d/%d = %.1f%%"
      % (np.isin(sub.top_subtype, [ENA, ENB]).sum(), n, np.isin(sub.top_subtype, [ENA, ENB]).mean() * 100))
print("  top 亚型分布:", sub.top_subtype.value_counts().to_dict())

# ---- 4) Cb 用的效应量（RNA 与 ATAC 各一套），富集池 vs COP/NFOL
def pool(M_, w):
    num = sum(M_[:, SUBS.index(s)] * w[s] for s in w); den = sum(w.values())
    return num / den
wE = {ENA: len(rows[ENA]), ENB: len(rows[ENB])}
wB = {SUBS[0]: len(rows[SUBS[0]])}
sub["rna_delta_enrich_vs_cop"] = pool(Mn, wE) - pool(Mn, wB)
sub["atac_delta_enrich_vs_cop"] = pool(sub[[f"atac_{s}" for s in SUBS]].values, wE) - pool(sub[[f"atac_{s}" for s in SUBS]].values, wB)
sub["rna_delta_mat"] = sub[f"rna_{ENB}"] - sub[f"rna_{SUBS[0]}"]
sub.to_csv(f"{OUT}/p08c_c1_rna_atac_by_subtype.csv", index=False)
g38.to_csv(f"{OUT}/p08c_c1_gene38_map.csv", index=False)
print("[写出] p08c_c1_rna_atac_by_subtype.csv / p08c_c1_gene38_map.csv")
print(sub[["gene", "mouse", "top_subtype", "rna_delta_enrich_vs_cop"]].sort_values("rna_delta_enrich_vs_cop", ascending=False).to_string(index=False))
