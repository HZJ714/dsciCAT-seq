import os
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""s3c_depth_matched_recheck.py — 「深度匹配重检」（cell 2b）
=================================================================================
背景：S3 已判 FAIL —— 两个亚簇的文库深度不同（sub-0 中位 2,682 UMI / 1,740 基因，
sub-1 中位 4,292 UMI / 2,322 基因，BH p<0.003、Cliff |δ|≈0.4）。146 个"独立差异基因"
很可能只是深度的影子。本脚本把**两组一起下采样到同一目标深度**后重跑同一条 S1 判据，
看差异基因是否还在。

⚠ 这是看到 S3 结果之后加的敏感性检验，**只能削弱结论、不能拯救结论**（判据方向写死：
   若匹配后合格基因 < 3 ⇒ S1 视为不稳健）。阈值与 S1 完全一致，不另立标准。
"""
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

RNG = np.random.default_rng(20261029)
D = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'u13_cop_4k/s3c')
CALL = {'Bmp4', 'Gpr17', 'C1ql1', 'Ptprz1'}

cnt = pd.read_csv(f'{D}/s3c_counts.csv')
genes = cnt['gene'].values
X = cnt.drop(columns='gene').values.astype(np.float64)      # genes × cells
meta = pd.read_csv(f'{D}/s3c_meta.csv')
lab = meta['cluster'].values                                 # 0/1
assert list(meta['cell']) == list(cnt.columns[1:])

TARGET = int(min(np.median(meta.loc[lab == 0, 'nCount_RNA']), np.median(meta.loc[lab == 1, 'nCount_RNA'])))
print('共同目标深度 T = %d UMI（两组中位数取小）; 细胞数 %d/%d' % (TARGET, (lab == 0).sum(), (lab == 1).sum()))


def downsample_to(X, col_totals, T):
    """每列多项式下采样到恰好 T；原本 < T 的列保持原样（无法上调），返回 (矩阵, 未动的列数)。"""
    Y = X.copy()
    untouched = 0
    for j in range(X.shape[1]):
        n = col_totals[j]
        if n <= T:
            untouched += 1
            continue
        p = X[:, j] / n
        Y[:, j] = RNG.multinomial(T, p)
    return Y, untouched


totals = meta['nCount_RNA'].values.astype(float)
Y, keep_low = downsample_to(X, totals, TARGET)
print('原本身低于 T、未做下采样的细胞 = %d / %d' % (keep_low, X.shape[1]))


def bh(p):
    """BH-FDR。⚠ 上一版把「排序后的 p」和「未排序的秩」直接相乘 ⇒ 秩与 p 错配，
       等于随机打乱的假 FDR（Prom1 真 p=4.3e-21 却算出 q=0.099）。现在按标准写法。"""
    p = np.asarray(p, dtype=float)
    n = p.size
    o = np.argsort(p, kind='stable')
    ps = p[o]
    qs = np.minimum.accumulate((ps * n / np.arange(1, n + 1))[::-1])[::-1]
    q = np.empty(n)
    q[o] = np.clip(qs, 0.0, 1.0)
    return q


def de_lognorm(M):
    """Seurat 同口径：CPM×1e4 再 log1p；两簇均值差 + Mann-Whitney + BH。"""
    lib = M.sum(axis=0, keepdims=True)
    L = np.log1p(M / np.where(lib == 0, 1, lib) * 1e4)
    d = L[:, lab == 0].mean(1) - L[:, lab == 1].mean(1)
    p = np.full(M.shape[0], np.nan)
    for i in range(M.shape[0]):
        a, b = L[i, lab == 0], L[i, lab == 1]
        if a.std() == 0 and b.std() == 0:
            p[i] = 1.0
            continue
        p[i] = mannwhitneyu(a, b, alternative='two-sided').pvalue
    q = bh(np.nan_to_num(p, nan=1.0))
    return d, p, q


# BH 自检（秩错配那个 bug 的回归测试）：与手算值比对
_chk = np.array([1e-12, 2e-6, 3e-4, 0.2, 0.9])
assert np.allclose(bh(_chk), np.array([5e-12, 5e-6, 5e-4, 0.25, 0.9]), atol=1e-12), 'BH 自检失败'

for tag, M in (('asis', X), ('matched', Y)):
    d, p, q = de_lognorm(M)
    not_call = ~np.isin(genes, list(CALL))
    qual = not_call & (np.abs(d) >= 0.5) & (q < 0.05)
    print('\n== %s ==' % tag)
    print('  可测基因 %d | 全部基因 FDR<0.05 %d | 非判据 & |Δ|>=0.5 & FDR<0.05 = %d'
          % (len(genes), int((q < 0.05).sum()), int(qual.sum())))
    idx = np.argsort(-np.abs(np.where(qual, d, 0)))[:12]
    print('  top: ' + ', '.join('%s(%+.2f)' % (genes[i], d[i]) for i in idx if qual[i]))
    
    pd.DataFrame({'gene': genes, 'mean_diff': d, 'p': p, 'p_bh': q, 'qualifies': qual}).to_csv(
        f'{D}/s3c_de_{tag}.csv', index=False)
