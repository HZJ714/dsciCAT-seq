import os
# -*- coding: utf-8 -*-
"""p0_8c_dual_setlevel.py —— Fig6C 到底该用哪个基因集：top1000 的基因层"着陆"集合级检验（双组学）

背景（作者指出 v31 偏航）：
  Fig6A/B 用的是 MS GWAS 的 **MAGMA top1000 基因集**（scDRS 打分集，meta.json: ms_set_size=1000），
  而 38/331 是 **FigS6A/B/E** 的 *cis*-window 位点解析，不是同一个问题。
  ⇒ 主图 C 若存在，域必须与 A/B 一致 = top1000。

本脚本做**非循环**的基因层检验（完全不碰 scDRS 分数，故无 self-inclusion）：
  显示量 = 被检验量：
      Δ[g] = (1006·mean[MOL3/6] + 2574·mean[MOL6])/3580 − mean[COP/NFOL]
  对照 = **检测率匹配**的非 MS 背景基因（六亚型平均阳性率十分位分层，同层内随机抽一个对照，10,000 次置换）
  通道：RNA（log2(CP10k+1)，X 已是该尺度）与 ATAC（gene score）分别做，并给联合标准化统计量。

内存：分块累加列和（1024 行/块），ATAC 稠密块 ~200 MB，峰值 < 0.6 GB。
"""
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, os
OUT = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'routeA_out/P0_8C/tables')
SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
ENA, ENB = "MOL3/6", "MOL6"
W = np.array([203, 490, 1016, 688, 1006, 2574], float)
fA = W[4] / (W[4] + W[5])                      # 0.2810
NPERM, SEED = 10000, 20261116

# ---------------------------------------------------------------- 1. 两个通道的 6 亚型均值 + 检测率
def channel(path, tag):
    cache = OUT + "/p08c_dual_%s.npz" % tag.lower()
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        print("%s: 用缓存 %s" % (tag, os.path.basename(cache)), flush=True)
        return z["vn"], z["MEAN"], z["PCT"], z["d"]
    a = ad.read_h5ad(path, backed="r")
    key = pd.Series(a.obs["cell_type"].astype(str)).str.upper().str.replace("/", "", regex=False).values
    vn = np.array(list(map(str, a.var_names)))
    rows = {s: np.flatnonzero(key == s.upper().replace("/", "")) for s in SUB}
    G = vn.size
    sm = np.zeros((6, G)); pos = np.zeros((6, G))
    for i in range(0, key.size, 1024):
        blk = np.arange(i, min(i + 1024, key.size))
        X = a.X[blk]
        X = X.toarray() if sp.issparse(X) else np.asarray(X)
        l = key[blk]
        for j, s in enumerate(SUB):
            m = l == s.upper().replace("/", "")
            if m.any():
                sm[j] += X[m].sum(0); pos[j] += (X[m] > 0).sum(0)
        if i and i % 4096 == 0:
            rss = int(open("/proc/self/status").read().split("VmRSS:")[1].split("\n")[0].split()[0]) / 1024
            print("   %s %d/%d 内存 %.0f MB" % (tag, blk[-1], key.size, rss), flush=True)
    MEAN = sm / W[:, None]
    PCT = pos / W[:, None] * 100
    d = fA * MEAN[4] + (1 - fA) * MEAN[5] - MEAN[0]        # 显示量 = 被检验量
    np.savez_compressed(cache, vn=vn, MEAN=MEAN, PCT=PCT, d=d)
    print("%s: 细胞 %d 基因 %d | Δ 中位 %+.4f" % (tag, key.size, G, np.median(d)), flush=True)
    return vn, MEAN, PCT, d

vn_r, Mr, Pr, dr = channel(os.path.join(os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')), 'oligo/3k_dsci_oligo_subset.h5ad'), "RNA")
vn_a, Ma, Pa, da = channel(os.path.join(os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')), 'oligo/fig6_rerun_20261019/dualomics_probe/genescore_oligo_5977.h5ad'), "ATAC")

# ---------------------------------------------------------------- 2. 基因集
L = pd.read_csv(OUT + "/p08c_genelevel_all6_pool_allgenes.csv")           # 875 映射（含 near_zero）
ms223 = set(L[~L.near_zero].gene.astype(str))
ms875 = set(L.gene.astype(str))
cis38 = set(pd.read_csv(OUT + "/p08c_c1_rna_atac_by_subtype.csv").mouse.astype(str))
anch = ["Mbp", "Plp1", "Mog", "Mag", "Cnp", "Mobp", "Sox10", "Olig2"]

def idx_of(vn, gs):
    s = {str(g) for g in gs}
    return np.flatnonzero(np.isin(vn, list(s)))

# ---------------------------------------------------------------- 3. 检测率匹配的置换检验
def matched_test(vn, PCT, d, gs, name, minpct=10.0, res=None):
    """同层（六亚型平均阳性率十分位）内随机抽一个对照 ⇒ 对照集均值分布 ⇒ 观测值分位。"""
    gi = idx_of(vn, gs)
    if res is not None:
        gi = np.array([i for i in gi if res[i]])
    mp = PCT.mean(0)
    pool = mp >= minpct
    ok = np.array([pool[i] for i in gi]) if gi.size else np.array([], bool)
    gi = gi[ok]
    if gi.size < 5:
        print("%-12s 入池 %d 太少" % (name, gi.size)); return None
    bg = np.flatnonzero(pool & ~np.isin(vn, list(ms875)))
    qs = np.quantile(mp[bg], np.linspace(0, 1, 11)); qs[0], qs[-1] = -1, 1e9
    bdec = np.digitize(mp[bg], qs[1:-1])
    bins = {b: bg[bdec == b] for b in range(10)}
    gdec = np.digitize(mp[gi], qs[1:-1])
    rng = np.random.default_rng(SEED)
    obs_r = d[gi].mean(); obs_m = np.median(d[gi])
    null_r = np.empty(NPERM); null_m = np.empty(NPERM)
    for k in range(NPERM):
        pick = np.concatenate([rng.choice(bins[b], size=int((gdec == b).sum()), replace=True)
                               for b in range(10) if (gdec == b).sum() and len(bins[b])])
        null_r[k] = d[pick].mean(); null_m[k] = np.median(d[pick])
    p_hi = (np.sum(null_r >= obs_r) + 1) / (NPERM + 1); p_lo = (np.sum(null_r <= obs_r) + 1) / (NPERM + 1)
    z = (obs_r - null_r.mean()) / null_r.std()
    print("%-12s n=%4d | 均值Δ%+.4f 匹配零%+.4f±%.4f z=%+.2f P(高)=%.4f P(低)=%.4f | 中位Δ%+.4f"
          % (name, gi.size, obs_r, null_r.mean(), null_r.std(), z, p_hi, p_lo, obs_m), flush=True)
    return dict(set=name, n=gi.size, mean_d=obs_r, null_mean=null_r.mean(), null_sd=null_r.std(),
                z=z, p_greater=p_hi, p_less=p_lo, median_d=obs_m, bg_pool=int(bg.size))

rep = []
for tag, vn, PCT, d in (("RNA", vn_r, Pr, dr), ("ATAC", vn_a, Pa, da)):
    print("\n=== %s 通道（显示量 Δ = %.4f·MOL3/6 + %.4f·MOL6 − COP/NFOL）===" % (tag, fA, 1 - fA), flush=True)
    for nm, gs in (("髓鞘锚(正对照)", anch), ("MS_top1000可测", ms223), ("MS_875映射", ms875), ("cis38", cis38)):
        r = matched_test(vn, PCT, d, gs, tag + ":" + nm)
        if r: r["channel"] = tag; rep.append(r)

pd.DataFrame(rep).to_csv(OUT + "/p08c_dual_setlevel_matched.csv", index=False)

# ---------------------------------------------------------------- 4. 给热图用的 top1000 着陆矩阵
pos_a = {g: i for i, g in enumerate(vn_a)}
ia = np.array([pos_a.get(g, -1) for g in vn_r])
M = pd.DataFrame({"gene": vn_r,
                  **{"rna_" + s: Mr[j] for j, s in enumerate(SUB)},
                  **{"atac_" + s: np.where(ia >= 0, Ma[j][np.clip(ia, 0, None)], np.nan) for j, s in enumerate(SUB)},
                  "rna_pct": Pr.mean(0), "rna_d": dr,
                  "atac_pct": np.where(ia >= 0, Pa.mean(0)[np.clip(ia, 0, None)], np.nan),
                  "atac_d": np.where(ia >= 0, da[np.clip(ia, 0, None)], np.nan)})
M["set223"] = M.gene.isin(ms223); M["set875"] = M.gene.isin(ms875); M["set38"] = M.gene.isin(cis38)
M.to_csv(OUT + "/p08c_dual_landing_genomewide.csv", index=False)
sub = M[M.set223].sort_values("rna_d", ascending=False)
print("\ntop1000 可测 %d 基因：RNA Δ>0 %d, <0 %d；ATAC Δ>0 %d | |Δ| p50/p90 RNA %.2f/%.2f ATAC %.2f/%.2f"
      % (len(sub), (sub.rna_d > 0).sum(), (sub.rna_d < 0).sum(), (sub.atac_d > 0).sum(),
         sub.rna_d.abs().median(), sub.rna_d.abs().quantile(.9), sub.atac_d.abs().median(), sub.atac_d.abs().quantile(.9)))
print("RNA×ATAC Δ Spearman (223) = %s" % sub[["rna_d", "atac_d"]].corr(method="spearman").iloc[0, 1].round(3))
print("\n检测率阈值敏感性（RNA 显示量 Δ 均值 vs 匹配零）：")
for q in (10, 20, 30, 40, 50):
    r = matched_test(vn_r, Pr, dr, ms223, "RNA pct>=%d" % q, minpct=q)
print("[写出] p08c_dual_setlevel_matched.csv, p08c_dual_landing_genomewide.csv", flush=True)
