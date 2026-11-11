# -*- coding: utf-8 -*-
"""b1_riskgene_by_magma.py — ① 面板的正确重算（替掉 panel6D_*_v1 两处问题）

为什么要重算：
1. v1 的行是"按 |MOL6 − COP/NFOL| 表达差"挑的 48 个基因，**不是 MAGMA 权重集**
   （与 rv1.gs 的 1,000 基因只重合 2/48）→ REASSESS P0-8a 要求"行 = MAGMA 权重降序"，v1 不满足。
2. v1 的 ATAC 块是 `peaks/by_celltype/*.bed` ∩ **人 gencode.v45 的 TSS±100 kb**。
   实测 bed 是**小鼠坐标**（chr1 最大 195,371,543 / chr19 61,330,903 = mm10/mm39），
   拿人的 TSS 去交鼠的 peak = 坐标物种错配，那些计数没有生物学含义 → **整块作废**。
   正解 = 直接用同物种、同 5,977 核的 ATAC gene-score 矩阵（genescore_oligo_5977.h5ad，
   它就是 scDRS ATAC 通道打分的输入，peak→gene 链接已在上游按小鼠注释完成）。

跨物种口径（必须进图注）：行来自**人** MS GWAS 的 MAGMA 基因级 z，数据是**小鼠**脑核，
匹配靠符号大写相等（top60 命中 58/60），同源关系未重新验证。

输出 NEWPANELS_v2/tables/panel6D_riskgene_subtype_{RNA,ATAC}.csv + 作废记录。
"""
import os, sys, json
import numpy as np, pandas as pd
import anndata as ad

WORK = os.environ.get('DSCI_WORK_ROOT', '')
OUT = WORK + "/figures/routeA_out/NEWPANELS_v2/tables"
os.makedirs(OUT, exist_ok=True)
RNA_H5 = WORK + "/FY/oligo/3k_dsci_oligo_subset.h5ad"
ATAC_H5 = WORK + "/FY/oligo/fig6_rerun_20261019/dualomics_probe/genescore_oligo_5977.h5ad"
RV1 = WORK + "/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs"
SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]      # 成熟轴顺序，与 Fig6A 一致
TOPN = 60
np.random.seed(0)


def magma_top(n):
    line = [l.rstrip("\n") for l in open(RV1) if l.startswith("PASS_Multiple_sclerosis\t")][0]
    z = {}
    for t in line.split("\t")[1].split(","):
        if ":" in t:
            k, v = t.rsplit(":", 1)
            try:
                z[k.strip().upper()] = float(v)
            except ValueError:
                pass
    return sorted(z.items(), key=lambda x: -x[1])[:n]


def subtype_mean(h5, genes):
    """按亚型取平均（稀疏矩阵只切片需要的列，避免整矩阵稠密化爆内存）"""
    a = ad.read_h5ad(h5, backed="r")
    ct = a.obs["cell_type"].astype(str)
    up = {str(v).upper(): i for i, v in enumerate(a.var_names)}
    idx, used = [], []
    for g in genes:
        if g in up:
            idx.append(up[g]); used.append(g)
    X = a.X[:, idx]
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    print("  %s: 矩阵 %s，命中 %d/%d 基因" % (os.path.basename(h5), a.shape, len(used), len(genes)))
    rows = []
    for g, j in zip(used, range(len(used))):
        v = X[:, j]
        rec = dict(gene=g)
        for s in SUB:
            rec[s] = float(np.mean(v[(ct == s).values]))
        rec["MOL6_minus_COP_NFOL"] = rec["MOL6"] - rec["COP/NFOL"]
        rows.append(rec)
    return pd.DataFrame(rows), {"n_cells": int(a.shape[0]), "genes_hit": len(used),
                                "nuclei_per_subtype": {s: int((ct == s).sum()) for s in SUB}}


top = magma_top(TOPN)
genes = [g for g, z in top]
zmap = dict(top)

print("[1/3] RNA 块")
rna, meta_r = subtype_mean(RNA_H5, genes)
print("[2/3] ATAC gene-score 块")
atac, meta_a = subtype_mean(ATAC_H5, genes)

# 行 = 两块都命中的基因，按 MAGMA z 降序
common = [g for g in genes if g in set(rna.gene) and g in set(atac.gene)]
rna = rna[rna.gene.isin(common)].copy()
atac = atac[atac.gene.isin(common)].copy()
for d in (rna, atac):
    d["magma_z"] = d.gene.map(zmap)
rna = rna.sort_values("magma_z", ascending=False).reset_index(drop=True)
atac = atac.sort_values("magma_z", ascending=False).reset_index(drop=True)
assert rna.magma_z.is_monotonic_decreasing and list(rna.gene) == list(atac.gene), "行序/行集合必须一致"

rna.insert(len(rna.columns) - 1, "nuclei_per_subtype", "; ".join("%s=%d" % (k, v) for k, v in meta_r["nuclei_per_subtype"].items()))
atac.insert(len(atac.columns) - 1, "nuclei_per_subtype", "; ".join("%s=%d" % (k, v) for k, v in meta_a["nuclei_per_subtype"].items()))
rna.to_csv(OUT + "/panel6D_riskgene_subtype_RNA.csv", index=False)
atac.to_csv(OUT + "/panel6D_riskgene_subtype_ATAC_genescore.csv", index=False)

json.dump(dict(top_n_requested=TOPN, rows_delivered=int(len(rna)), magma_source=RV1.split("/")[-1],
               rna_matrix=RNA_H5.split("/")[-1], atac_matrix=ATAC_H5.split("/")[-1],
               rna_meta=meta_r, atac_meta=meta_a,
               species_note="rows = human MAGMA gene-level z; data = mouse brain nuclei; symbol-uppercase match, orthology not re-verified",
               voided="panel6D_riskgene_subtype_ATAC_peakcount.csv (v1) is void: human GENCODE v45 TSS intersected with mouse-coordinate peak beds"),
          open(OUT + "/panel6D_meta.json", "w"), indent=1)
print("[3/3] 写出 %d 行 × 两块；RNA 命中 %d，ATAC 命中 %d" % (len(rna), meta_r["genes_hit"], meta_a["genes_hit"]))
print(rna[["gene", "magma_z", "COP/NFOL", "MOL6", "MOL6_minus_COP_NFOL"]].head(8).to_string(index=False))
