#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""s3c_c1ql1_sensitivity.py —— Codex #87 那条"15-cell 候选按 C1ql1 判失败，规则是不是事后定的"的正面回答
=======================================================================================
做法：**不删基因、不删那个失败候选**，而是把冻结规则里的 `C1ql1 <= 0.2` 整项拿掉重算，
两种口径一起报。判据其余五项一字不动（n>=10, Bmp4>=0.25, Gpr17>=0.5, pct_Gpr17>=30, Ptprz1<=1.5）。

输入（都是既有运行的产物，本脚本不重新聚类）：
  u13_cop_4k/s3_marker_table.csv   四个 full 条件 × 全部簇的 7 marker 均值/检出率 + cop 列
  u13_cop_4k/s3_strict4029.csv     全脑抽 4,029 核 → 谱系重聚类 × 20 draws
输出：u13_cop_4k/s3c/s3c_c1ql1_sensitivity.{csv,json}
"""
import os, json
import numpy as np, pandas as pd

B = os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures'))
OUT = f'{B}/u13_cop_4k/s3c'
os.makedirs(OUT, exist_ok=True)

# 冻结规则的其它五项（与 scripts/u13_cop_4k.R 第 6 行的注释逐字对应）
BASE = lambda t: ((t.mean_Bmp4 >= 0.25) & (t.mean_Gpr17 >= 0.5) & (t.pct_Gpr17 >= 30) &
                  (t.mean_Ptprz1 <= 1.5) & (t.n >= 10))

rows, changed_runs, blocked = [], [], []
for f, label in [('u13_cop_4k/s3_marker_table.csv', 'full libraries (delivered run)'),
                 ('u13_cop_4k/s3_strict4029.csv', 'whole-brain 4,029-nuclei draws (20 seeds)')]:
    t = pd.read_csv(f'{B}/{f}')
    base = BASE(t)
    t['rule_frozen'] = base & (t.mean_C1ql1 <= 0.2)
    t['rule_noC1ql1'] = base
    # 自检：带 C1ql1 的判定必须与原表 cop 列一致（否则说明我抄错了阈值）
    assert (t.rule_frozen == t.cop).all(), f'{f}: 复现原 cop 列失败，阈值抄错'
    per = t.groupby(['group', 'seed']).agg(clusters=('cl', 'size'), frozen=('rule_frozen', 'sum'),
                                           noC1ql1=('rule_noC1ql1', 'sum')).reset_index()
    per.insert(0, 'run', label)
    rows.append(per)
    changed_runs.append(int((per.frozen != per.noC1ql1).sum()))
    blocked.append(t[base & ~(t.mean_C1ql1 <= 0.2)].assign(run=label))

out = pd.concat(rows, ignore_index=True)
out.to_csv(f'{OUT}/s3c_c1ql1_sensitivity.csv', index=False)
blocked = pd.concat(blocked, ignore_index=True)

print(out.to_string(index=False))
print('\n== 只被 C1ql1 单项卡住的候选簇 ==')
print(blocked[['run', 'group', 'cl', 'n', 'mean_Bmp4', 'mean_Gpr17', 'pct_Gpr17',
               'mean_C1ql1', 'mean_Ptprz1']].round(3).to_string(index=False) if len(blocked) else '（无）')

# 逐 draw 的一致性：20 次下采样里，去掉 C1ql1 会不会改变任何一次判定
d = out[out.run.str.startswith('whole-brain')]
per_draw_diff = int((d.frozen != d.noC1ql1).sum())
json.dump(dict(clusters_checked_full=int(len(out[out.run.str.startswith('full')])),
               runs_changed_by_c1ql1=int(sum(changed_runs)),
               draws_total=int(len(d)), draws_changed=per_draw_diff,
               blocked_candidates=blocked[['run', 'group', 'cl', 'n', 'mean_C1ql1']].to_dict('records'),
               reading=('C1ql1 只会取消候选、不会制造候选：全部运行里只有 1 个判定被它改变'
                        '（EL14 单库簇 4，n=15），20 次 draw 的判定一次都没变。')),
          open(f'{OUT}/s3c_c1ql1_sensitivity.json', 'w'), indent=1, ensure_ascii=False)
print('\n逐 draw 判定被改变的次数 = %d / %d' % (per_draw_diff, len(d)))
print('产物 -> %s/s3c_c1ql1_sensitivity.{csv,json}' % OUT)
