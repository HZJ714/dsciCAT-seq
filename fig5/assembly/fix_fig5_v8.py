#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fix_fig5_v8.py — Figure 5 v7 → v8：B 面板加**网络边方向/权重图例**（Codex #91 / G4）
==================================================================================
只加一个图例块，**不重画网络**：B 面板的 44 个节点、44 条边、全部标签与坐标一字不动
（自检 = 除图例盒以外与交付版逐像素 0 可见差）。

为什么图例放左上角：实测交付版 `Figure_5.pdf` 的 B 面板内最大空白矩形 =
5.25 × 2.00 cm @ x[6.60,11.85] y[0.45,2.45]（网络是放射状，左上角天然空）；
图例块实测占 4.94 × 1.25 cm，落在盒内且离最近标签 >0.1 cm。

图例内容全部**实测自成品矢量**（见 `scripts/network_key.py` 模块头）：
边宽 1.00–4.00 pt（原图）→ 印刷 0.36–1.44 pt，线性映射 TF2G importance；
边色 = 源 TF（MOL3 实测 5 种 tab10 色 = 5 个 TF）；无 alpha 编码（配置里的
transparency_edge_by 在成品里没有生效 ⇒ 图例不写"深浅"）；
regulation 列全为 1 ⇒ "全部为预测的激活连接"。

⚠ 加完必须跑 `finish()`：figkit 写 Arial 会把 526 KB 整字库带进来，
  `doc.subset_fonts()` 裁回 30 KB（实测逐像素 0 差、文字层逐字符相同），
  否则本图会从"全门通过"退化成带【未子集】豁免的件。

用法: PYTHONPATH=<work>/pylibs python3 scripts/fix_fig5_v8.py
输出: <work>/figures/Fig5_v8.pdf   （交付块替换由 finalize 步骤另行确认）
"""
import os, sys, subprocess
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fitz
import numpy as np
from figkit import CM, save_doc
from network_key import draw_key

B = os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures'))
SRC = f'{B}/FINAL_SET_20261021/Figure_5.pdf'      # 交付版 = Fig5_v7 + 字母常规补丁
OUT = f'{B}/Fig5_v8.pdf'
KEY_AT = (6.75, 0.60)
KEY_MAXW = 5.00
FREE_BOX = (6.60, 0.45, 11.72, 1.90)              # 400 dpi 实测 0 墨像素的盒（原 5.25x2.0 盒右下角有一条边线穿过）


def main():
    doc = fitz.open(SRC)
    pg = doc[0]
    # 0) 落笔前实测盒子是空的
    x0, y0, x1, y1 = FREE_BOX
    pm = pg.get_pixmap(dpi=200, clip=fitz.Rect(x0 * CM, y0 * CM, x1 * CM, y1 * CM), alpha=False)
    a = np.frombuffer(pm.samples, np.uint8).reshape(pm.height, pm.width, pm.n)[:, :, :3]
    print('空白盒 x[%.2f,%.2f] y[%.2f,%.2f] 墨占比 = %.4f%%' % (
        x0, x1, y0, y1, 100 * (a.min(2) < 235).mean()))
    assert (a.min(2) < 235).mean() < 0.0005, '盒子不空'
    w, h = draw_key(pg, KEY_AT[0], KEY_AT[1], KEY_MAXW, imp=(0.6, 4.2))
    assert KEY_AT[0] + w <= x1 + 0.02 and KEY_AT[1] + h <= y1 + 0.02, '图例越界'
    save_doc(doc, OUT, clean=False)

    # 1) 自检：图例盒以外必须 0 可见差
    old = fitz.open(SRC)[0]
    new = fitz.open(OUT)[0]
    for dpi in (150, 300):
        A = new.get_pixmap(dpi=dpi, alpha=False)
        Bm = old.get_pixmap(dpi=dpi, alpha=False)
        x = np.frombuffer(A.samples, np.uint8).reshape(A.height, A.width, A.n)[:, :, :3].astype(int)
        y = np.frombuffer(Bm.samples, np.uint8).reshape(Bm.height, Bm.width, Bm.n)[:, :, :3].astype(int)
        d = np.abs(x - y).max(2) > 24
        box = np.zeros_like(d)
        box[int(y0 * dpi / 2.54):int(y1 * dpi / 2.54), int(x0 * dpi / 2.54):int(x1 * dpi / 2.54)] = True
        print('  %d dpi: 图例盒内 %d px，盒外 %d px' % (dpi, int((d & box).sum()), int((d & ~box).sum())))
        assert (d & ~box).sum() == 0, '图例以外有改动'
    return OUT


def finish(path=OUT):
    S = os.path.dirname(os.path.abspath(__file__))
    env = dict(os.environ, PYTHONPATH=os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
    subprocess.run([sys.executable, f'{S}/arialize.py', '--inplace', path], env=env,
                   capture_output=True, text=True)
    d = fitz.open(path)
    d.subset_fonts(verbose=False)
    tmp = path + '.sub.tmp'
    d.save(tmp, garbage=4, deflate=True)
    d.close()
    os.replace(tmp, path)
    print('finish(): arialize + subset_fonts 完成')
    return path


if __name__ == '__main__':
    main()
    finish()
