#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p0_8c_logo_genelevel.py —— Fig6C 基因级检验的 self-inclusion 精确消融（P0-8c 增补）

为什么需要这一步
----------------
scDRS 的 norm_score 本身就是 MS 基因集表达的加权和，所以 Fig6C 的
"逐基因表达 vs norm_score 相关"里，被测基因同时是**被相关对象的成分**：

    cov(x_g, s) = w_g·var(x_g) + Σ_{j≠g} w_j·cov(x_g, x_j)
                  └ 机械项 ┘      └ 真实共变项 ┘

而零分布（1000 个随机 GWAS 集的 ctrl_norm_score）几乎不含该基因，机械项
没有被零分布吸收 ⇒ MC P / FDR 会把一部分机械项当成信号。

本脚本用**项目自带的 c2_logo.py 消融产物**（figures/routeA_out/GENELEVEL_v1/out_logo/，
每个基因一份"去掉该基因、同一 scdrs compute-score 参数重打"的分核分数）直接做
leave-one-gene-out 消融，不依赖任何权重重构：

    r_delivery(g) = Pearson(x_g, 交付 norm_score)                    ← 现行 Fig6C 统计量
    r_FULLbatch(g)= Pearson(x_g, 同批 FULL_<batch> 分数)             ← 消融的配对对照
    r_LOO(g)      = Pearson(x_g, 同批 LOO_g 分数)                    ← 去掉自己之后
    机械项 = r_FULLbatch − r_LOO（同批配对，抵消 ctrl 抽样的批间抖动）

零分布沿用交付件的 1000 个 ctrl_norm_score（LOO 的 .score.gz 不含 ctrl 列，
`--flag-return-ctrl-norm-score True` 只写进 full_score 且未落盘），
即 LOO 统计量与交付统计量共用同一零分布 ⇒ 两套 P 可横向比较；此近似的方向是
**保守**的（ctrl 集与 999 基因集的差异 < 与 1000 基因集的差异）。

口径纪律：细胞池 = all6_pool（与 Fig6A/Fig6C 主口径一致），基因 = 223 个被检基因。
"""
import os, sys, glob, gzip
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
import scipy.stats as st

WORK = os.environ.get('DSCI_WORK_ROOT', '')
OUT = WORK + "/figures/routeA_out/P0_8C"
RNA_H5 = WORK + "/FY/oligo/3k_dsci_oligo_subset.h5ad"
RROOT = WORK + "/FY/oligo/scdrs_traits_draft/full74"
TRAIT = "PASS_Multiple_sclerosis"
LOGODIR = WORK + "/figures/routeA_out/GENELEVEL_v1"
SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
EXPECT = "203/490/1016/688/991/2574"          # Fig6A 细胞口径，变了就停

S = pd.read_csv(sorted(glob.glob(RROOT + "/*/" + TRAIT + ".full_score.gz"))[0],
                sep="\t", index_col=0)
ctrl_cols = [c for c in S.columns if c.startswith("ctrl_norm_score_")]
assert len(ctrl_cols) == 1000, "ctrl 列数 %d != 1000" % len(ctrl_cols)

a = ad.read_h5ad(RNA_H5)
key = (pd.Series(a.obs["cell_type"].astype(str)).str.upper()
       .str.replace("/", "", regex=False).str.replace("-", "", regex=False)).values
mask = a.obs_names.isin(set(S.index))
grp = {s: np.where(mask & (key == s.replace("/", "")))[0] for s in SUB}
print("[1] 细胞口径 %s（须等于 %s）" % ("/".join(str(len(grp[s])) for s in SUB), EXPECT))
assert "/".join(str(len(grp[s])) for s in SUB) == EXPECT, "细胞口径与 Fig6A 不同批！"
rows = np.concatenate([grp[s] for s in SUB])
names = a.obs_names[rows]

K = pd.read_csv(OUT + "/tables/p08c_genelevel_all6_pool.csv")
pos = {str(v): i for i, v in enumerate(a.var_names)}
cidx = np.array([pos[g] for g in K.gene])
X = a.X[np.ix_(rows, cidx)]
X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
X = X.astype(np.float32)
y = S.loc[names, "norm_score"].values.astype(np.float32)
C = S.loc[names, ctrl_cols].values.astype(np.float32)

Ec = X - X.mean(0)
yc = y - y.mean()
Cc = C - C.mean(0)
sdE = np.sqrt((Ec ** 2).sum(0))
sdC = np.sqrt((Cc ** 2).sum(0))
sdy = np.sqrt((yc ** 2).sum())
r_del = (Ec * yc[:, None]).sum(0) / (sdE * sdy + 1e-12)
null = (Ec.T @ Cc) / (sdE[:, None] * sdC[None, :] + 1e-12)      # 223 x 1000
print("[2] 交付口径复现：中位 r=%.4f（表里 %.4f），最大偏差 %.1e"
      % (np.median(r_del), K.r.median(), np.abs(r_del - K.r.values).max()))

L = pd.read_csv(LOGODIR + "/logo_delta1_contributions.csv")
L["gkey"] = L.gene.str.upper()
Lm = L.set_index("gkey")


def load_score(path):
    """只读 norm_score 列（gzip 流式，避免整表进内存）"""
    with gzip.open(path, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        j = head.index("norm_score")
        cells = []
        vals = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            cells.append(f[0]); vals.append(float(f[j]))
    return pd.Series(vals, index=cells)


cache, rec = {}, []
for i, g in enumerate(K.gene):
    G = str(K.gene_human.iloc[i]).upper()
    if G not in Lm.index:
        rec.append((g, np.nan, np.nan, np.nan, "")); continue
    row = Lm.loc[G]
    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]
    tag, gi = str(row["batch"]), int(row["gene_index"])
    flo = glob.glob("%s/%s/LOO_%04d_*.score.gz" % (LOGODIR + "/out_logo", tag, gi))
    ffl = glob.glob("%s/%s/FULL_%s.score.gz" % (LOGODIR + "/out_logo", tag, tag))
    if not flo or not ffl:
        rec.append((g, np.nan, np.nan, np.nan, tag)); continue
    if tag not in cache:
        cache[tag] = load_score(ffl[0])
    sL = load_score(flo[0]).reindex(names)
    sF = cache[tag].reindex(names)
    if sL.isna().any() or sF.isna().any():
        rec.append((g, np.nan, np.nan, np.nan, tag)); continue
    xg = Ec[:, i]
    loo = sL.values.astype(np.float32); loo = loo - loo.mean()
    full = sF.values.astype(np.float32); full = full - full.mean()
    r_loo = float(xg @ loo / (np.sqrt((xg ** 2).sum()) * np.sqrt((loo ** 2).sum()) + 1e-12))
    r_fb = float(xg @ full / (np.sqrt((xg ** 2).sum()) * np.sqrt((full ** 2).sum()) + 1e-12))
    rec.append((g, r_loo, r_fb, float(np.corrcoef(loo, full)[0, 1]), tag))
R = pd.DataFrame(rec, columns=["gene", "r_LOO", "r_FULLbatch", "r_LOOvsFULL", "batch"])
R = R.merge(K[["gene", "gene_human", "r", "mc_p", "mc_z", "fdr_bh", "mean_pct"]], on="gene")
R["delta_self"] = R.r_FULLbatch - R.r_LOO
R["null_med"] = np.median(null, 1)
R["z_LOO"] = (R.r_LOO - null.mean(1)) / (null.std(1, ddof=1) + 1e-12)
mp = np.full(len(R), np.nan)
for i in range(len(R)):
    if np.isfinite(R.r_LOO.values[i]):
        mp[i] = (np.sum(null[i] >= R.r_LOO.values[i]) + 1) / (null.shape[1] + 1)
R["mc_p_LOO"] = mp


def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); n = len(p); q = np.empty(n)
    qs = p[o] * n / np.arange(1, n + 1)
    q[o] = np.minimum.accumulate(qs[::-1])[::-1]
    return np.clip(q, 0, 1)


R["fdr_LOO"] = bh(np.where(np.isfinite(R.mc_p_LOO), R.mc_p_LOO, 1.0))
ok = R.dropna(subset=["r_LOO"])
print("[3] 可消融基因 %d/%d（c2_logo 覆盖 %d 个 MS 基因）" % (len(ok), len(R), L.gene.nunique()))
print("    LOO 分数与同批 FULL 分数（逐细胞）中位 Pearson=%.4f" % ok.r_LOOvsFULL.median())
print("    中位 r: 交付 %+.4f | 同批 FULL %+.4f | LOO %+.4f ⇒ 机械份额 %.0f%%"
      % (R.r.median(), ok.r_FULLbatch.median(), ok.r_LOO.median(),
         100 * (1 - ok.r_LOO.median() / ok.r_FULLbatch.median())))
print("    LOO 后：r>0 %d/%d | MC P<0.05 %d | FDR<0.05 %d | FDR<0.10 %d"
      % (int((ok.r_LOO > 0).sum()), len(ok), int((ok.mc_p_LOO < 0.05).sum()),
         int((ok.fdr_LOO < 0.05).sum()), int((ok.fdr_LOO < 0.10).sum())))
old = list(K.sort_values(["mc_p", "r"], ascending=[True, False]).gene[:10])
new = list(ok.sort_values(["mc_p_LOO", "r_LOO"], ascending=[True, False]).gene[:10])
print("    新旧 r 排名 Spearman=%.3f | top10 交集 %d/10"
      % (st.spearmanr(ok.r_LOO, ok.r)[0], len(set(new) & set(old))))
print("    LOO top10: " + ", ".join(ok.sort_values("r_LOO", ascending=False).gene_human[:10]))
R.to_csv(OUT + "/tables/p08c_genelevel_LOOexact.csv", index=False)
print("[4] 写出 tables/p08c_genelevel_LOOexact.csv")
