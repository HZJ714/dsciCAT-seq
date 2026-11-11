#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p0_8c_s6c_logo_audit.py —— 把 S6C 的三层稳健性口径全部在"消融分数"下重算一遍

要回答的问题（回执 §8.7 之后遗留的）：S6C 里那些 FDR 数字是不是和主图一样被 self-inclusion 污染？
污染有多重？换池/换层能不能救回来？

三个层（都是逐基因 Pearson vs 分数，零分布都用交付件同一批 1000 个 ctrl_norm_score）：
  L1 all6_pool              —— 主图口径（5962 细胞）
  L2 enriched_pool          —— S6C 的"富集亚型池"稳健性（Fig6A group_q<0.05 = MOL3/6 + MOL6）
  L3 all6 + 亚型内去均值     —— S6C 的"去掉成熟度梯度"稳健性（回执里曾报 85/223 FDR<0.05）
每个层各算两版分数：
  DEL = 交付 norm_score（含被检基因自己 ⇒ 机械项在分子里）
  LOO = 把该基因从 1000 基因集删掉、同一 scdrs compute-score 参数重打分（c2_logo 产物）

关键机制（为什么"换层"救不了）：亚型内去均值只缩小分母（去掉分数在亚型间的均值差），
分子里的自协变项 w_g·var(x_g) 一点没少 ⇒ 机械占比不降反升。这里用数字直接验证。
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
EXPECT = "203/490/1016/688/991/2574"

S = pd.read_csv(sorted(glob.glob(RROOT + "/*/" + TRAIT + ".full_score.gz"))[0], sep="\t", index_col=0)
ctrl_cols = [c for c in S.columns if c.startswith("ctrl_norm_score_")]
a = ad.read_h5ad(RNA_H5)
key = (pd.Series(a.obs["cell_type"].astype(str)).str.upper()
       .str.replace("/", "", regex=False).str.replace("-", "", regex=False)).values
mask = a.obs_names.isin(set(S.index))
grp = {s: np.where(mask & (key == s.replace("/", "")))[0] for s in SUB}
assert "/".join(str(len(grp[s])) for s in SUB) == EXPECT, "细胞口径与 Fig6A 不同批！"
ENRICH = ["MOL3/6", "MOL6"]   # Fig6A group_q<0.05（用 SUB 原标签，grp 的键）
K = pd.read_csv(OUT + "/tables/p08c_genelevel_all6_pool.csv")
pos = {str(v): i for i, v in enumerate(a.var_names)}
cidx = np.array([pos[g] for g in K.gene])
Xall = a.X[np.ix_(np.arange(a.shape[0]), cidx)]
Xall = np.asarray(Xall.todense()) if sp.issparse(Xall) else np.asarray(Xall)
Xall = Xall.astype(np.float32)
yall = S.reindex(a.obs_names)["norm_score"].values.astype(np.float32)
Call = S.reindex(a.obs_names)[ctrl_cols].values.astype(np.float32)
names = list(a.obs_names)
print("[1] 细胞 %d | 基因 %d | ctrl %d" % (mask.sum(), len(K), Call.shape[1]))

L = pd.read_csv(LOGODIR + "/logo_delta1_contributions.csv")
L["gkey"] = L.gene.str.upper()
Lm = L.set_index("gkey")


def load_norm(path):
    with gzip.open(path, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        j = head.index("norm_score")
        cells, vals = [], []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            cells.append(f[0]); vals.append(float(f[j]))
    return pd.Series(vals, index=cells)


def center(M, gid):
    """按亚型去均值（gid: 每行属于哪个亚型）"""
    out = M.copy()
    for g in np.unique(gid):
        m = gid == g
        out[m] -= out[m].mean(0)
    return out


def pearson(A, b):
    """A (n x p) 已中心化列；b (n,) 已中心化 -> (p,) 相关"""
    b = b - b.mean()
    num = A.T @ b
    den = np.sqrt((A ** 2).sum(0)) * np.sqrt((b ** 2).sum())
    return num / (den + 1e-12)


def layer(idx, gids, tag):
    """该层的：中心化表达矩阵、中心化交付分数、223x1000 零分布（全部同口径处理）"""
    Xc = Xall[idx]; Cc = Call[idx]; yv = yall[idx].copy()
    if gids is not None:
        Xc = center(Xc, gids); Cc = center(Cc, gids); yv = center(yv[:, None], gids)[:, 0]
    Xc = Xc - Xc.mean(0); Cc = Cc - Cc.mean(0); yv = yv - yv.mean()
    sdX = np.sqrt((Xc ** 2).sum(0)); sdC = np.sqrt((Cc ** 2).sum(0))
    null = (Xc.T @ Cc) / (sdX[:, None] * sdC[None, :] + 1e-12)
    print("   [%s] 细胞 %d 零分布 %s" % (tag, len(idx), null.shape))
    return Xc, Cc, sdX, null, yv


def mcp(null, r):
    """add-one MC P（单侧），与交付件同式"""
    return (np.sum(null >= r[:, None], 1) + 1) / (null.shape[1] + 1)


def bh(p):
    o = np.argsort(p); q = np.empty_like(p); m = len(p)
    q[o] = np.minimum.accumulate((p[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    return q


LAY = {}
idx_all = np.concatenate([grp[s] for s in SUB]); gid_all = np.concatenate([[i] * len(grp[s]) for i, s in enumerate(SUB)])
idx_enr = np.concatenate([grp[s] for s in ENRICH]); gid_enr = np.concatenate([[i] * len(grp[s]) for i, s in enumerate(ENRICH)])
LAY["L1 all6"] = layer(idx_all, None, "L1 all6 不去均值")
LAY["L2 enriched"] = layer(idx_enr, None, "L2 富集池")
LAY["L3 all6+demean"] = layer(idx_all, gid_all, "L3 亚型内去均值")
LAY["L4 enriched+demean"] = layer(idx_enr, gid_enr, "L4 富集池+去均值")

# 每层两套观察统计量：DEL（同一层中心化后的分数向量）与 LOO（逐基因消融分数）
res = {k: {"r_del": np.zeros(len(K)), "r_loo": np.zeros(len(K))} for k in LAY}
for k, (Xc, Cc, sdX, null, yv) in LAY.items():
    res[k]["r_del"] = (Xc.T @ yv) / (sdX * np.sqrt((yv ** 2).sum()) + 1e-12)

print("[2] 读 223 个消融分数 …")
loo = {}
for i, g in enumerate(K.gene):
    G = str(K.gene_human.iloc[i]).upper()
    if G not in Lm.index:
        loo[i] = None; continue
    row = Lm.loc[G]
    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]
    fl = glob.glob("%s/%s/LOO_%04d_*.score.gz" % (LOGODIR + "/out_logo", str(row["batch"]), int(row["gene_index"])))
    loo[i] = load_norm(fl[0]) if fl else None
print("   命中 %d/%d" % (sum(v is not None for v in loo.values()), len(K)))

for k, (Xc, Cc, sdX, null, yv) in LAY.items():
    idxs = idx_all if ("L1" in k or "L3" in k) else idx_enr
    nm = np.array([names[j] for j in idxs])
    rl = np.full(len(K), np.nan)
    for i in range(len(K)):
        s = loo[i]
        if s is None:
            continue
        v = s.reindex(nm).values.astype(np.float32)
        if np.isnan(v).any():
            continue
        if "demean" in k:
            gid = gid_all if "L3" in k else gid_enr
            for g in np.unique(gid):
                m = gid == g
                v[m] -= v[m].mean()
        v = v - v.mean()
        rl[i] = (Xc[:, i] @ v) / (np.linalg.norm(Xc[:, i]) * np.linalg.norm(v) + 1e-12)
    res[k]["r_loo"] = rl

rows = []
for k, (Xc, Cc, sdX, null, yv) in LAY.items():
    for ver in ("r_del", "r_loo"):
        r = res[k][ver]
        ok = np.isfinite(r)
        p = np.full(len(K), np.nan); p[ok] = mcp(null[ok], r[ok])
        q = bh(np.where(np.isnan(p), 1.0, p))
        setmed = np.nanmedian(r)
        nullset = np.median(null, 0)
        z = (setmed - nullset.mean()) / nullset.std()
        P = (np.sum(nullset >= setmed) + 1) / (len(nullset) + 1)
        rows.append(dict(layer=k, score="DEL(含自己)" if ver == "r_del" else "LOO(消融)",
                         median_r=setmed, nominal_lt05=int(np.nansum(p < 0.05)),
                         FDR_lt05=int(np.nansum(q < 0.05)), set_z=z, set_MC_P=P))
T = pd.DataFrame(rows)
print("\n" + "=" * 104)
print("%-20s %-14s %9s %11s %9s %8s %9s" % ("层", "分数", "中位 r", "名义P<.05", "FDR<.05", "集合z", "集合MCP"))
for _, x in T.iterrows():
    print("%-20s %-14s %+9.4f %7d/223 %7d/223 %8.2f %9.4f"
          % (x.layer, x.score, x.median_r, x.nominal_lt05, x.FDR_lt05, x.set_z, x.set_MC_P))
T.to_csv(OUT + "/tables/p08c_s6c_layers_del_vs_loo.csv", index=False)
print("\n[写出] tables/p08c_s6c_layers_del_vs_loo.csv")
# 池间/层间秩相关（LOO 版）——S6C 真正能留的东西
def top20(v):
    o = np.argsort(-np.nan_to_num(v, nan=-9))[:20]
    return set(K.gene.iloc[o])
out = pd.DataFrame({"gene": K.gene, "gene_human": K.gene_human})
for k in LAY:
    out[k + "_DEL"] = res[k]["r_del"]; out[k + "_LOO"] = res[k]["r_loo"]
out.to_csv(OUT + "/tables/p08c_s6c_layers_per_gene.csv", index=False)
print("[写出] tables/p08c_s6c_layers_per_gene.csv")
for tag, pre in (("DEL", "DEL"), ("LOO", "LOO")):
    r12 = st.spearmanr(out["L1 all6_" + pre], out["L2 enriched_" + pre]).statistic
    r13 = st.spearmanr(out["L1 all6_" + pre], out["L3 all6+demean_" + pre]).statistic
    ov = len(top20(out["L1 all6_" + pre].values) & top20(out["L2 enriched_" + pre].values))
    print("%s 口径：L1×L2 rho=%.3f L1×L3 rho=%.3f top20 交集 %d/20" % (tag, r12, r13, ov))
print("DEL×LOO（同层 L1）rho=%.3f ⇒ 排序保住、水平塌掉" % st.spearmanr(out["L1 all6_DEL"], out["L1 all6_LOO"]).statistic)
print("机械偏移（DEL-LOO）各层中位：%s" % ", ".join("%s=%.4f" % (k, np.nanmedian(out[k + "_DEL"] - out[k + "_LOO"])) for k in LAY))
