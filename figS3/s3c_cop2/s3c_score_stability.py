import os
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""s3c_score_stability.py — 给扰动网格算「两簇分裂复现率 + ARI」（判据预先冻结，见
   S3C_COP2_FROZEN_CRITERIA_20261029.md 的 S2：恰好 2 簇且各 ≥10 核 记为复现；
   要求复现率 ≥70% 且 ARI 均值 ≥0.70）。ARI 用 sklearn（R 端没装 mclust）。
输出 s3c_stability_scored.csv / s3c_S2.json
"""
import json, os
import numpy as np, pandas as pd
from sklearn.metrics import adjusted_rand_score as ARI

D = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'u13_cop_4k/s3c')
long = pd.read_csv(f'{D}/s3c_stability_long.csv')
asg = pd.read_csv(f'{D}/s3c_assignments.csv')

# 逐次运行的 ARI（与基线 cluster_id 比较，只用该次运行里存在的细胞）
a = {}
for rid, g in asg.groupby('run'):
    a[int(rid)] = ARI(g['baseline'], g['cluster'])
long['ari'] = long['run'].map(a)
long['split2'] = (long['n_clusters'] == 2) & (long['min_size'] >= 10)
long['kind'] = np.where(long['cond'] == 'full', 'full', 'downsampled')
long.to_csv(f'{D}/s3c_stability_scored.csv', index=False)

def rate(sel):
    g = long[sel]
    return dict(n=len(g), rec=int(g.split2.sum()), pct=100.0 * g.split2.mean(),
                ari=float(g.ari.mean()), ari_med=float(g.ari.median()))

grid = {f'r{r}|{k}': rate((long.res == r) & (long.kind == k))
        for r in sorted(long.res.unique()) for k in ('full', 'downsampled')}
S2 = dict(primary_frozen_grid=rate(pd.Series(True, index=long.index)),   # 冻结原文：整个扰动网格
          secondary_res_matched=rate(long.res == 0.5),                  # 只看与图同分辨率的那一档
          grid=grid)
# **主口径按冻结原文 = 整个网格**。res=0.1/1.0 本来就不该给出 2 簇，这会压低复现率——
# 这一点如实写进图注，不做事后改判据；res 匹配档只作次要参考值一并报出。
S2['pass'] = bool(S2['primary_frozen_grid']['pct'] >= 70 and S2['primary_frozen_grid']['ari'] >= 0.70)
json.dump(S2, open(f'{D}/s3c_S2.json', 'w'), indent=1)

for key in ('primary_frozen_grid', 'secondary_res_matched'):
    v = S2[key]
    print('%-22s %3d/%3d 复现 = %5.1f%% | ARI 均值 %.3f 中位 %.3f'
          % (key, v['rec'], v['n'], v['pct'], v['ari'], v['ari_med']))
print('S2 判定（主口径 = 冻结原文的整个网格）: %s' % ('PASS' if S2['pass'] else 'FAIL'))
for k in sorted(grid, key=lambda s: (float(s.split('|')[0][1:]), s.split('|')[1])):
    v = grid[k]
    print('  %-18s %2d/%2d = %5.1f%%  ARI %.3f' % (k, v['rec'], v['n'], v['pct'], v['ari']))
print('\n簇数分布（res=0.5）:', long[long.res == 0.5].n_clusters.value_counts().sort_index().to_dict())
