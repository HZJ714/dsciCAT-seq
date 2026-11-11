import os
# -*- coding: utf-8 -*-
"""p0_8c_rna_setlevel.py — Fig6C 基因层的"非循环"版本：亚型间表达差的集合级检验。

统计量完全不碰 scDRS 分数 ⇒ 没有自包含（self-inclusion）问题。
  每个基因 d = mean(MOL6) − mean(COP/NFOL)（log2(CP10k+1) 尺度，X 已是该尺度，勿再归一化）
  检验：MS 风险基因集合（223 / 875）与"检测率匹配的非 MS 背景基因"的 d 分布比较（MWU，单侧）
  正对照：8 个髓鞘/OL 锚基因（若它们不显著，说明检验本身没灵敏度）

内存：分块累加列和，峰值 ~0.4 GB（直接 a.X[布尔] 会稠密化到 18 GB+，已踩过坑）。
"""
import numpy as np, pandas as pd, anndata as ad, glob, scipy.stats as st, scipy.sparse as sp
OUT = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'routeA_out/P0_8C/tables')
SUBS = ["COPNFOL", "MOL3", "MOL36", "MOL2", "MFOL", "MOL6"]   # 6 亚型标签（已去 '/'）
a = ad.read_h5ad(os.path.join(os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')), 'oligo/3k_dsci_oligo_subset.h5ad'), backed="r")
S = pd.read_csv(sorted(glob.glob(os.path.join(os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')), 'oligo/scdrs_traits_draft/full74/*/PASS_Multiple_sclerosis.full_score.gz')))[0], sep="\t", index_col=0)
key = pd.Series(a.obs["cell_type"].astype(str)).str.upper().str.replace("/", "", regex=False).values
keep = a.obs_names.isin(set(S.index))
lab = key[keep]
rows = np.flatnonzero(keep)
vn = np.array(list(map(str, a.var_names)))
G = len(vn)
grp = {s: rows[lab == s] for s in SUBS}
print("细胞 %d | 基因 %d | 各组 %s" % (rows.size, G, {s: len(v) for s, v in grp.items()}), flush=True)

sm = np.zeros((len(SUBS), G)); pos = np.zeros((len(SUBS), G))
CH = 1024
for i in range(0, rows.size, CH):
    blk = rows[i:i + CH]
    X = a.X[blk]
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    l = key[blk]
    for j, s in enumerate(SUBS):
        m = l == s
        if m.any():
            sm[j] += X[m].sum(0); pos[j] += (X[m] > 0).sum(0)
    if i % (CH * 8) == 0:
        print("  %d/%d" % (i + CH, rows.size), flush=True)
n = np.array([len(grp[s]) for s in SUBS])[:, None]
MEAN, PCT = sm / n, pos / n * 100
d = MEAN[SUBS.index("MOL6")] - MEAN[SUBS.index("COPNFOL")]
pct2 = PCT[SUBS.index("MOL6")] + PCT[SUBS.index("COPNFOL")]
mean_all = PCT.mean(0)
pd.DataFrame({"gene": vn, "d_mol6_cop": d, "pct_both": pct2, "mean_pct": mean_all,
              **{s: MEAN[SUBS.index(s)] for s in SUBS}}).to_csv(OUT + "/p08c_rna_delta_genomewide.csv", index=False)

L = pd.read_csv(OUT + "/p08c_subtype_stats_long.csv")
hit = set(L[L.subtype == "MOL6"].gene.astype(str))
ms223 = set(pd.read_csv(OUT + "/p08c_genelevel_all6_pool.csv").gene.astype(str))
anch = ["Mbp", "Plp1", "Mog", "Mag", "Cnp", "Mobp", "Sox10", "Olig2"]
inmat = lambda ss: np.array(sorted(set(vn) & {str(g) for g in ss}))
ok = (pct2 >= 20) & np.isfinite(d)
bg = ok & ~np.isin(vn, list(hit))
print("\n背景（非 MS 映射基因、MOL6+COP 合计阳性率>=20%%）= %d | MS223 入池 %d | 锚块入池 %d"
      % (bg.sum(), np.isin(vn, list(ms223))[ok].sum(), np.isin(vn, anch)[ok].sum()), flush=True)
rep = []
def t(name, gs):
    idx = ok & np.isin(vn, list(gs))
    v = d[idx]
    if v.size < 5:
        print("%-16s n=%d 太少" % (name, v.size)); return
    u = st.mannwhitneyu(v, d[bg], alternative="greater")
    dk = (np.mean(v) - np.mean(d[bg])) / np.std(d[bg])
    print("%-16s n=%5d 中位d%+.4f 背景%+.4f 差%+.4f MWU P=%.2g CohenD=%.3f" % (name, v.size, np.median(v), np.median(d[bg]), np.median(v) - np.median(d[bg]), u.pvalue, dk), flush=True)
    rep.append(dict(set=name, n=v.size, med_d=np.median(v), bg_d=np.median(d[bg]), diff=np.median(v) - np.median(d[bg]), mwu_p=u.pvalue, cohend=dk))
t("髓鞘锚块(正对照)", anch); t("MS 223", ms223); t("MS 875 映射", hit)
for q in (10, 20, 30, 40, 50, 60, 70, 80):
    sel = np.array(sorted(set(vn) & ms223 & set(vn[(pct2 >= q)])))
    if sel.size >= 15: t("MS pct>=%d" % q, sel)
pd.DataFrame(rep).to_csv(OUT + "/p08c_rna_setlevel.csv", index=False)
print("\n[写出] p08c_rna_delta_genomewide.csv, p08c_rna_setlevel.csv", flush=True)
