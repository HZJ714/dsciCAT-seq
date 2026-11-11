# -*- coding: utf-8 -*-
"""b1c_riskgene_by_magma.py — 把 ① 面板（Fig6 风险基因 × 亚型）的人→鼠基因匹配
**统一到 scDRS 内置同源表**，与产生分数的口径完全一致。

为什么要有这个脚本（b1 的口径缺陷）：
  b1 用「人类符号大写后与小矩阵符号大写相等」匹配（`up={v.upper():i}`）。
  而 scDRS 打分时走的是包内置同源表：
      scdrs/util.py: load_homolog_mapping("hsapiens","mmusculus")
      scdrs/data/mouse_human_homologs.txt （16,466 对，HUMAN_GENE_SYM -> MOUSE_GENE_SYM）
      compute-score --gs-species hsapiens --h5ad-species mmusculus
  两套口径不一致会被审稿人问倒（"为什么图上的基因集和打分的基因集不同？"）。
  大写同名匹配会把小鼠里恰好同名的基因当成人源直系同源收进来：
      PRM1 (z=4.523) / PRM2 (z=4.695)  小鼠睾丸特异精子蛋白，内置表判无同源 → 假阳性
  本脚本改用内置表映射，并且**用映射后的小鼠名义精确匹配矩阵 var_names**（不再 upper）。

审计结论（写死在这里，免得后人再查）：
  top60 候选 → 内置表可映射 54；这 54 个小鼠名义**精确存在**于 RNA(31,263) 与 ATAC(24,333) 矩阵；
  且 54 个映射结果**全部是同名大小写变体**（CLEC16A->Clec16a），无任何别名映射，
  因此统一口径**不改变任何表达数值**，只把 56 行减到 54 行。
  被剔除 6 个：FAM213B, CD58, FAM69A, FCRL3, PRM2, PRM1
  （前 4 个 b1 也没命中矩阵，所以行数只减 2）。

内置表的已知局限（必须进图注，不能隐瞒）：
  1000 基因 MS 集里 109 个映射不到，其中 76 个小鼠矩阵也无同名（人特异/旧命名），
  33 个矩阵有同名但表判无（WARS 只收现代名 WARS1；FOXP1/TCF7/PVR/ESRRA/CSF2RB/IFITM1/SETD1B
  在表的两列都不存在）→ 属表收录缺口。本脚本**不擅自修补**，因为一修补就得重跑 scDRS，
  已签的交付分数会全部作废。口径以"与打分一致"为准。

输出：
  NEWPANELS_v3/tables/panel6D_riskgene_subtype_RNA.csv
  NEWPANELS_v3/tables/panel6D_riskgene_subtype_ATAC_genescore.csv
  NEWPANELS_v3/tables/homolog_mapping_audit.csv   60 行逐基因映射与去留记录（可交付证据）
  NEWPANELS_v3/tables/panel6D_meta.json

跑法（scdrs 装在 python3.9 用户站点，必须用 3.9 环境）：
  $DSCI_PY b1c_riskgene_by_homolog.py          # python 3.9 + scdrs 1.0.2
"""
import os, sys, json
import numpy as np, pandas as pd
import anndata as ad, scipy.sparse as sp

sys.path.insert(0, os.path.expanduser('"~/.local/lib/python3.9/site-packages'))
from scdrs.util import load_homolog_mapping            # 权威口径：scDRS 包内置同源表

WORK = os.environ.get('DSCI_WORK_ROOT', '')
OUT = WORK + "/figures/routeA_out/NEWPANELS_v3/tables"
os.makedirs(OUT, exist_ok=True)
RNA_H5 = WORK + "/FY/oligo/3k_dsci_oligo_subset.h5ad"
ATAC_H5 = WORK + "/FY/oligo/fig6_rerun_20261019/dualomics_probe/genescore_oligo_5977.h5ad"
RV1 = WORK + "/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs"
TRAIT = "PASS_Multiple_sclerosis"
SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]      # 成熟轴顺序，与 Fig6A 一致
TOPN = 60

H2M = load_homolog_mapping("hsapiens", "mmusculus")
print("scDRS 内置同源表: %d 对（人->鼠）" % len(H2M))


def magma_top(n):
    """rv1.gs 的 MS 行：gene:z 形式，按 z 降序取前 n（人类名义，未转换）。"""
    line = [l.rstrip("\n") for l in open(RV1) if l.startswith(TRAIT + "\t")][0]
    z = {}
    for t in line.split("\t")[1].split(","):
        if ":" in t:
            k, v = t.rsplit(":", 1)
            try:
                z[k.strip()] = float(v)
            except ValueError:
                pass
    return sorted(z.items(), key=lambda x: -x[1])[:n]


def subtype_block(h5, pairs):
    """pairs = [(human, mouse), ...]；用**小鼠名义精确匹配** var_names。
    只切片需要的列，避免 5,977 x 31,263 整矩阵稠密化爆内存。"""
    a = ad.read_h5ad(h5, backed="r")
    ct = a.obs["cell_type"].astype(str)
    pos = {str(v): i for i, v in enumerate(a.var_names)}      # 精确、区分大小写
    idx, used = [], []
    for hu, mo in pairs:
        if mo in pos:
            idx.append(pos[mo]); used.append((hu, mo))
    X = a.X[:, idx]
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    print("  %s: 矩阵 %s，精确命中 %d/%d 基因" % (os.path.basename(h5), a.shape, len(used), len(pairs)))
    det = {s: (ct == s).values for s in SUB}
    rows = []
    for j, (hu, mo) in enumerate(used):
        v = X[:, j]
        rec = dict(gene=mo, gene_human=hu)
        dcols = []
        for s in SUB:
            seg = v[det[s]]
            rec[s] = float(np.mean(seg))
            dcols.append(float(np.mean(seg > 0)) * 100.0)
        rec["det_pct_mean"] = float(np.mean(dcols))            # 六亚型平均检测率，供②筛行用
        rec["MOL6_minus_COP_NFOL"] = rec["MOL6"] - rec["COP/NFOL"]
        rows.append(rec)
    return pd.DataFrame(rows), {"n_cells": int(a.shape[0]), "n_genes": int(a.shape[1]),
                                "genes_hit": len(used),
                                "nuclei_per_subtype": {s: int(det[s].sum()) for s in SUB}}


top = magma_top(TOPN)
audit = []
pairs = []
for hu, z in top:
    mo = H2M.get(hu)
    if mo is None:
        audit.append(dict(gene_human=hu, magma_z=z, gene_mouse="", kept=False,
                          reason="no ortholog in scDRS built-in table"))
    else:
        audit.append(dict(gene_human=hu, magma_z=z, gene_mouse=mo, kept=True,
                          reason="" if mo.lower() == hu.lower() else "non-eponymous ortholog symbol"))
        pairs.append((hu, mo))
A = pd.DataFrame(audit)
print("[1/4] 同源映射审计: top%d -> 可映射 %d，剔除 %d" % (TOPN, A.kept.sum(), (~A.kept).sum()))
print("      剔除名单: %s" % ", ".join(A[~A.kept].gene_human))
print("      非同名映射: %s" % (", ".join(A[A.reason == "non-eponymous ortholog symbol"].gene_human)
                              or "无（全部为同名大小写变体）"))

print("[2/4] RNA 块")
rna, meta_r = subtype_block(RNA_H5, pairs)
print("[3/4] ATAC gene-score 块")
atac, meta_a = subtype_block(ATAC_H5, pairs)

# 行 = 两块都命中的基因，按 MAGMA z 降序
zmap = dict(top)
common = [(h, m) for h, m in pairs if h in set(rna.gene_human) and h in set(atac.gene_human)]
rna = rna[rna.gene_human.isin([h for h, m in common])].copy()
atac = atac[atac.gene_human.isin([h for h, m in common])].copy()
for d in (rna, atac):
    d["magma_z"] = d.gene_human.map(zmap)
    d.sort_values("magma_z", ascending=False, inplace=True)
    d.reset_index(drop=True, inplace=True)
    d.insert(0, "rank_magma", range(1, len(d) + 1))
assert list(rna.gene) == list(atac.gene) and rna.magma_z.is_monotonic_decreasing, "行序/行集合必须一致"

note = ("cross-species matching: human MAGMA symbols converted to mouse symbols by the scDRS "
        "built-in table scdrs/data/mouse_human_homologs.txt (16,466 pairs), i.e. the identical "
        "rule used for scDRS scoring (--gs-species hsapiens --h5ad-species mmusculus)")
for d in (rna, atac):
    d["homolog_source"] = "scDRS built-in mouse_human_homologs.txt"
A.to_csv(OUT + "/homolog_mapping_audit.csv", index=False)
rna.to_csv(OUT + "/panel6D_riskgene_subtype_RNA.csv", index=False)
atac.to_csv(OUT + "/panel6D_riskgene_subtype_ATAC_genescore.csv", index=False)

json.dump(dict(top_n_requested=TOPN, ortholog_mapped=int(A.kept.sum()),
               dropped_no_ortholog=list(A[~A.kept].gene_human), rows_delivered=int(len(rna)),
               homolog_table="scdrs/data/mouse_human_homologs.txt", homolog_pairs_in_table=len(H2M),
               ms_gene_set_total=1000, ms_gene_set_used_by_scdrs=891, ms_gene_set_dropped=109,
               magma_source=RV1.split("/")[-1], rna_matrix=RNA_H5.split("/")[-1],
               atac_matrix=ATAC_H5.split("/")[-1], rna_meta=meta_r, atac_meta=meta_a,
               species_note=note,
               value_change_vs_v2="none (all 54 mappings are case-only variants of the same symbol)",
               replaced="NEWPANELS_v2 tables (symbol-uppercase match; included PRM1/PRM2 false positives)"),
          open(OUT + "/panel6D_meta.json", "w"), indent=1)
print("[4/4] 写出 %d 行 × 两块 | RNA 命中 %d | ATAC 命中 %d" % (len(rna), meta_r["genes_hit"], meta_a["genes_hit"]))
print(rna[["rank_magma", "gene", "gene_human", "magma_z", "det_pct_mean", "COP/NFOL", "MOL6"]].head(8).to_string(index=False))
