import os
# -*- coding: utf-8 -*-
"""p0_8c_handpick_audit.py —— "手动挑几个 MS 基因"能不能成立？三种挑法的合法性与实测

作者提议：top1000 阴性，那就手动挑几个 MS 相关基因展示亚型特异性表达，与 A/B 呼应。
问题不在"挑"，而在**挑的依据是否独立于要展示的那个量**。这里把三种挑法全部实测：

  P1 依据 = MAGMA z（GWAS 侧连续轴，与表达数据无关）→ 剂量–反应检验：z 越高，着陆越偏富集亚型？
  P2 依据 = 位置（38 个 *cis*-window 基因，Table 2 预定义）→ **对照换成同一 331 OL 面板里剩下的 293 个基因**
     （这是唯一能回答"是 MS 特异还是 OL 通用"的比较；之前用全基因组做对照，阳性被"分母本身是 OL 基因"污染）
  P3 依据 = 对 A/B 统计量的贡献（S2_exact_gene_contribution）→ 只能描述"信号弥散度"，不能当富集证据
  P4 依据 = 我们自己的着陆量（即真正意义的"手动挑好看的"）→ 循环，直接判死，仅作对照展示其虚假显著性
"""
import numpy as np, pandas as pd, scipy.stats as sp
TAB = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'routeA_out/P0_8C/tables')
GS = os.path.join(os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')), 'oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs')
PLAT = os.path.join(os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')), 'oligo/scdrs_ms_out/genes_of_interest.txt')
SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
NPERM, SEED = 10000, 20261116

M = pd.read_csv(TAB + "/p08c_dual_landing_genomewide.csv")
L = pd.read_csv(TAB + "/p08c_genelevel_all6_pool_allgenes.csv")
ms223 = set(L[~L.near_zero].gene.astype(str))
P331 = set(x.strip().upper() for x in open(PLAT) if x.strip())
C38 = set(str(x).upper() for x in pd.read_csv(TAB + "/p08c_c1_rna_atac_by_subtype.csv").mouse.astype(str))
print("331 面板 %d | 38 命中 %d | 面板内未命中 %d | 223 可测 %d"
      % (len(P331), len(C38), len(P331 - C38), len(ms223)))

# ---- P1: MAGMA z 连续轴（GWAS 侧，与表达无关）
zmap = {}
for line in open(GS):
    p = line.rstrip("\n").split("\t")
    if len(p) == 2 and "Multiple_sclerosis" in p[0]:
        for tok in p[1].split(","):
            if ":" in tok:
                g, z = tok.rsplit(":", 1)
                try:
                    zmap[g.upper()] = float(z)
                except ValueError:
                    pass
print("MS MAGMA z 可用基因 %d" % len(zmap))
h2m = dict(zip(L.gene_human.astype(str).str.upper(), L.gene.astype(str)))
M["magma_z"] = [zmap.get(str(g).upper(), np.nan) for g in M.gene.astype(str).str.upper()]
M["gene_u"] = M.gene.astype(str).str.upper()
M["magma_z"] = [zmap.get(g, np.nan) for g in M.gene_u]
out = []
for tag in ("rna", "atac"):
    d, pct = M[tag + "_d"].values, M[tag + "_pct"].values
    mean = M[[tag + "_" + s for s in SUB]].mean(1).values
    pool = (pct >= 10) & np.isfinite(d) & np.isfinite(mean)
    zz = M.magma_z.values
    # (a) 风险集内部：z 越高，着陆是否越高？
    sel = pool & np.isfinite(zz)
    rho, pv = sp.spearmanr(zz[sel], d[sel])
    print("\n=== %s ===" % tag.upper())
    print("P1 风险集内剂量-反应：n=%d Spearman(z, Δ) = %+.3f  P=%.3f" % (sel.sum(), rho, pv))
    # (b) 全基因组（含非风险基因）：偏相关，控制检测率与平均表达
    import statsmodels.api as sm
    X = np.column_stack([np.log1p(pct[pool]), mean[pool], np.isfinite(zz)[pool].astype(float),
                         np.nan_to_num(zz[pool])])
    ym = d[pool]
    fit = sm.OLS(ym, sm.add_constant(X)).fit()
    print("P1b 全基因组控制(检测率,平均表达,是否入top1000)后 z 的偏效应 β=%.5f  SE=%.5f  P=%.2f"
          % (fit.params[3], fit.bse[3], fit.pvalues[3]))
    out.append(dict(channel=tag, test="P1 dose-response within risk set", n=int(sel.sum()), stat=rho, p=pv))

# ---- P2: 38 vs 同面板 293（唯一能区分 MS 特异 / OL 通用的对照）
res2 = []
for tag in ("rna", "atac"):
    d, pct = M[tag + "_d"].values, M[tag + "_pct"].values
    mean = M[[tag + "_" + s for s in SUB]].mean(1).values
    pool = (pct >= 10) & np.isfinite(d) & np.isfinite(mean)
    gu = M.gene_u.values
    is38 = pool & np.array([g in C38 for g in gu])
    is293 = pool & np.array([g in (P331 - C38) for g in gu])
    if is38.sum() < 5 or is293.sum() < 20:
        print("P2 %s 样本不足 38=%d 293=%d" % (tag, is38.sum(), is293.sum())); continue
    obs = d[is38].mean() - d[is293].mean()
    # 二维匹配：在 38∪293 合池里按 (检测率, 平均表达) 5×5 分层，层内置换标签
    sub_i = np.flatnonzero(is38 | is293)
    qb = np.quantile(pct[sub_i], [.2, .4, .6, .8]); qm = np.quantile(mean[sub_i], [.2, .4, .6, .8])
    cellb = np.digitize(pct[sub_i], qb) * 6 + np.digitize(mean[sub_i], qm)
    lab = np.array([gu[i] in C38 for i in sub_i])
    cells = {c: np.flatnonzero(cellb == c) for c in np.unique(cellb)}
    rng = np.random.default_rng(SEED + 3)
    null = np.empty(NPERM)
    for k in range(NPERM):
        pl = lab.copy()
        for c, idx in cells.items():
            if len(idx) > 1:
                pl[idx] = rng.permutation(pl[idx])
        null[k] = d[sub_i][pl].mean() - d[sub_i][~pl].mean()
    p_hi = (np.sum(null >= obs) + 1) / (NPERM + 1)
    zsc = (obs - null.mean()) / null.std()
    print("P2 %s：38(n=%d) Δ̄%+.4f vs 同面板293(n=%d) Δ̄%+.4f ⇒ 差 %+.4f  层内置换零%+.4f±%.4f  z=%+.2f P=%.4f"
          % (tag.upper(), is38.sum(), d[is38].mean(), is293.sum(), d[is293].mean(), obs,
             null.mean(), null.std(), zsc, p_hi))
    res2.append(dict(channel=tag, n38=int(is38.sum()), n293=int(is293.sum()), d38=d[is38].mean(),
                     d293=d[is293].mean(), diff=obs, null_mean=null.mean(), null_sd=null.std(),
                     z=zsc, p=p_hi))
pd.DataFrame(res2).to_csv(TAB + "/p08c_p2_38vs293.csv", index=False)
print("\n[写出] p08c_p2_38vs293.csv")

# ---- P3: 对 A/B 统计量的贡献集中度（描述性）
C = pd.read_csv(os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'SUPPL_ANALYSES_20261112/tables/S2_exact_gene_contribution.csv'))
C = C.sort_values("share_of_delta_gene_term", key=np.abs, ascending=False)
tot = C.share_of_delta_gene_term.abs().sum()
for k in (1, 5, 10, 20):
    print("P3 前 %2d 个基因占 |贡献| 的 %.1f%%（该表共 %d 行，|贡献|合计 %.2f）"
          % (k, 100 * C.share_of_delta_gene_term.abs().head(k).sum() / tot, len(C), tot))
print("P3 top5:", ", ".join("%s %.1f%%" % (r.gene, 100 * r.share_of_delta_gene_term / tot)
                            for r in C.head(5).itertuples()))

# ---- P4: 反面示范——按着陆量挑 top10，再看它在同框架下"多显著"（说明为何不能用）
for tag in ("rna", "atac"):
    d, pct = M[tag + "_d"].values, M[tag + "_pct"].values
    mean = M[[tag + "_" + s for s in SUB]].mean(1).values
    pool = (pct >= 10) & np.isfinite(d) & np.isfinite(mean)
    idx = np.flatnonzero(pool & M.gene.astype(str).isin(ms223).values)
    top = idx[np.argsort(-np.abs(d[idx]))][:10]
    bg = np.flatnonzero(pool & ~M.gene.astype(str).isin(set(L.gene.astype(str))).values)
    qb = np.quantile(pct[bg], [.2, .4, .6, .8]); qm = np.quantile(mean[bg], [.2, .4, .6, .8])
    cbg = np.digitize(pct[bg], qb) * 6 + np.digitize(mean[bg], qm)
    cells = {c: bg[np.flatnonzero(cbg == c)] for c in np.unique(cbg)}
    csel = np.digitize(pct[top], qb) * 6 + np.digitize(mean[top], qm)
    rng = np.random.default_rng(SEED + 7)
    null = np.empty(NPERM)
    for k in range(NPERM):
        pick = np.array([cells[c][rng.integers(len(cells[c]))] if c in cells else -1 for c in csel])
        null[k] = d[pick[pick >= 0]].mean()
    obs = d[top].mean()
    print("P4 %s：按 |Δ| 挑 top10 → Δ̄=%+.4f vs 匹配零 %+.4f±%.4f  z=%+.1f  P=%.4f（**这是挑出来的，无效**）"
          % (tag.upper(), obs, np.nanmean(null), np.nanstd(null), (obs - np.nanmean(null)) / np.nanstd(null),
             (np.sum(null >= obs) + 1) / (NPERM + 1)))
