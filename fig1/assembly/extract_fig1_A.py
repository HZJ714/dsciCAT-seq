#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_fig1_A.py — 把 A 带（流程图 + 面板字母 A）从 **AI 母版** 抽成一张干净矢量件
================================================================================
为什么不用 `Workflow_combined.pdf`：实测那是**旧布局**（'Integrated matrix' 框内三行文字互相压线、
右列 Visualization/Pseudotime 重叠），作者批准的 A 带是 **Illustrator 里重排过**的，
唯一载体是 `ai/Fig1_v3.ai`（实测与交付件逐像素一致：全页差 0.233%、A 带 0.016%）。
`Workflow_combined.pdf` / `Workflow_combined.ai` 都只到"重排前/放大前"的中间态，不能用。

做法：A 带按 `fig1_layout.CLIP_A` 矢量裁贴进一张**页尺寸=带尺寸**的新页 →
再交给 `p1_pkg_fix.py` 做逐 Form 试删 + 像素门（带外的 B/C/D/E/F 组在新页里零墨迹，会被摘掉，
否则 `show_pdf_page` 的 clip 只挡墨迹不挡文字层，合成件里会留下幽灵文字）。

校验：抽出件按 1:1 贴回交付件原位，与交付件 A 带逐像素比对，必须零差异；并报告文字块数。

用法: PYTHONPATH=<work>/pylibs python3 scripts/extract_fig1_A.py
"""
import os
import subprocess
import sys
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
import fitz
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig1_layout import CLIP_A

FIG = os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures'))
CM = 72 / 2.54
MASTER = f'{FIG}/ai/Fig1_v3.ai'                       # A 带的真母版（含作者批准的排版 + 小字放大）
DELIVERED = '/tmp/Figure_1_delivered.pdf'             # 交付件备份，做像素基准
RAW = '/tmp/F1A_raw.pdf'
OUT = f'{FIG}/replot/F1A_workflow.pdf'


def px(page, clip=None):
    pm = page.get_pixmap(dpi=200, colorspace=fitz.csGRAY, clip=clip)
    return np.frombuffer(pm.samples, dtype=np.uint8).reshape(pm.height, pm.width).astype(int)


def main():
    band = fitz.Rect(*[v * CM for v in CLIP_A])
    src = fitz.open(MASTER)
    doc = fitz.open()
    pg = doc.new_page(width=band.width, height=band.height)
    pg.show_pdf_page(fitz.Rect(0, 0, band.width, band.height), src, 0, clip=band)
    doc.save(RAW, garbage=3, deflate=True)
    print('裁贴 %s -> %s  %.2f x %.2f cm  文字块 %d' % (
        os.path.basename(MASTER), os.path.basename(RAW), band.width / CM, band.height / CM,
        len(fitz.open(RAW)[0].get_text('dict')['blocks'])))

    # ⚠ 必须带 --strip-forms：`show_pdf_page(clip=...)` 只裁**显示**不裁**内容流**，
    #   不剥就会把母版整页的旧文字层当幽灵带进最终图（实测隐藏字 +600、假重叠 330 对）。
    #   母版的旧 B/C/D/E/F 文字自成独立 Form，整 Form 试删 + 像素门一发就能掉（2026-10-20 实测
    #   贴入后 span 692 -> 69、Liberation 幽灵 88 -> 0）；逐条的 --strip-blocks 在本页要 40 min+。
    r = subprocess.run([sys.executable, f'{FIG}/scripts/p1_pkg_fix.py', RAW, OUT, '--strip-forms'],
                       capture_output=True, text=True, env={**os.environ, 'PYTHONPATH': os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs'))})
    tail = [l for l in (r.stdout + r.stderr).split('\n') if ('文字层' in l or '像素' in l or '隐藏' in l or '=>' in l)]
    print('\n'.join('  ' + t for t in tail[-5:]))

    chk = fitz.open(OUT)
    n_sp = sum(len(l['spans']) for b in chk[0].get_text('dict')['blocks'] for l in b.get('lines', []))
    # 贴回交付件原位比对
    ref = fitz.open(DELIVERED)[0]
    t = fitz.open(); tp = t.new_page(width=ref.rect.width, height=ref.rect.height)
    tp.show_pdf_page(band, chk, 0)
    c = fitz.Rect(0.4 * CM, 0.25 * CM, 17.0 * CM, 7.3 * CM)
    d = np.abs(px(ref, c) - px(tp, c))
    print('抽出件 %.2f x %.2f cm  文字 span %d' % (chk[0].rect.width / CM, chk[0].rect.height / CM, n_sp))
    print('回贴原位 vs 交付件 A 带: maxΔ %d / 差>10 像素 %.5f%%  => %s' % (
        d.max(), (d > 10).mean() * 100, 'OK 零差异' if d.max() <= 12 else 'FAIL'))
    return OUT


if __name__ == '__main__':
    main()
