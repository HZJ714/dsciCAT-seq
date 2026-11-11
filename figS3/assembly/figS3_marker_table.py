#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""FigS3 的"这些簇是怎么注释的"证据 —— 出**表**不出图 (2026-10-13 取舍, 见 README_Fig3 §15)。

口径与 figS3_markerspec.py 完全一致, 全部现算, 不手写 marker 字典:
  矩阵 = scripts/cache/avgexpr_{EL14,FY13,FY14}.csv (每簇平均表达, 5,597 个三库共有 HVG)
  特异度 = 该簇均值 − 其它簇最大值  (drop-out 稳健的"只在这个簇高"的度量)
  每簇取特异度前 8 个基因 = 该簇的注释依据; n / % 来自 counts_*.csv

输出: ../tables/TableS3_cluster_markers.{csv,md}
用法: python3 figS3_marker_table.py
"""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
CACHE, OUT = f'{B}/scripts/cache', f'{B}/tables'
LIBS = [('HT-scCAT-seq', 'EL14'), ('dsciCAT-seq-50k', 'FY13'), ('dsciCAT-seq-100k', 'FY14')]
TOPN = 8


def markers(tag):
    ae = pd.read_csv(f'{CACHE}/avgexpr_{tag}.csv', index_col=0)
    cn = pd.read_csv(f'{CACHE}/counts_{tag}.csv').set_index('Var1')['Freq']  # Var1=簇名, Freq=n
    v = ae.values
    rows = []
    for i, cl in enumerate(ae.index):
        oth = np.delete(v, i, axis=0).max(0)                  # 其它簇的最大值 (逐基因)
        d = pd.DataFrame({'g': ae.columns, 'm': v[i], 'o': oth})
        d['spec'] = d['m'] - d['o']
        top = d.sort_values('spec', ascending=False).head(TOPN)
        n = int(cn.get(cl, 0))
        rows.append(dict(library=LIBS[[t for _, t in LIBS].index(tag)][0], cluster=cl,
                         n_cells=n, pct=round(100.0 * n / cn.sum(), 2),
                         markers='; '.join(f"{r.g}({r.m:.1f})" for r in top.itertuples()),
                         marker_genes='; '.join(top['g'])))
    return pd.DataFrame(rows)


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    df = pd.concat([markers(t) for _, t in LIBS], ignore_index=True)
    df.to_csv(f'{OUT}/TableS3_cluster_markers.csv', index=False)
    with open(f'{OUT}/TableS3_cluster_markers.md', 'w') as fh:
        fh.write('# Table S3 | 簇注释依据（每簇特异度前 8 个基因）\n\n'
                 '基因 = 三库共有高变基因 (n = 5,597)；特异度 = 该簇平均表达 − 其它簇最大值；'
                 '括号内为该簇平均表达 (log-normalised)。全部由 `figS3_marker_table.py` 现算。\n')
        for lib, g in df.groupby('library', sort=False):
            fh.write(f'\n## {lib}（{len(g)} 簇，{g.n_cells.sum():,} cells）\n\n'
                     '| Cluster | n | % | Top-8 cluster-specific genes (mean) |\n'
                     '|---|---|---|---|\n')
            for r in g.itertuples():
                fh.write(f'| {r.cluster} | {r.n_cells:,} | {r.pct}% | {r.markers} |\n')
    print(f'{OUT}/TableS3_cluster_markers.csv   {df.shape[0]} 行 (87 簇)')
    for lib, g in df.groupby('library', sort=False):
        print(f'  {lib:<18} {len(g):>2} 簇  n={g.n_cells.sum():>6,}  '
              f'例: {g.iloc[0].cluster} -> {g.iloc[0].marker_genes[:60]}')
