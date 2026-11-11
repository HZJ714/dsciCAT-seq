#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fix_fig5_v9.py — Figure 5 的网络图例改成**单行**，放到 MOL3 上面、与 B 平齐
==========================================================================
作者第二轮意见（针对 Fig5 的 G4 图例）：
  「图例确实要加，但是下面的小字是不是有点多余，感觉可以写图注里，影响观感；
    然后把 grn 的图例一排放 MOL3 的上面，和 B 平齐，是不是会好点。
    figs5 因为空间多，倒是可以所有信息保留。」

⇒ 只动 **Figure 5**：
  ① 四行说明（TF2G importance / edge colour / all edges activating / gene node colour）
     从图面上删掉，全部进图注（`edits6c_fig5_key_legend.py`）；
  ② 剩下的图形符号排成一行 `◇ TF ──▸ ◇ target gene  ─ 1.0×  ═ 4.0×`（实测宽 4.84 cm），
     放在 y 0.03–0.40（与 B 字母同一行）、x 8.35 起 —— 实测该带（x8.1–18, y0–0.55）
     在 10-21 原件里**暗像素占比 0.0000%**，MOL3 标题从 y 0.92 才开始，不会压到任何东西。
  ③ Figure S5 的三行版图例**不动**（空间够，信息全留）。

起点选 **10-21 原始交付件**（`_superseded/Figure_5_20261021_deprecated.pdf`）而不是当前
交付件：那件**从来没有图例**，所以"删掉旧图例"这件事根本不需要做——不用 redaction，
就不会碰到 `apply_redactions` 在带嵌套 Form 的页面上可能连矢量一起吃掉的雷。
它唯一的旧毛病是【未子集】（整份 Arial 内嵌），用 `doc.subset_fonts()` 收掉，
自检 = 子集化前后逐像素 0 差、文字层逐字符不变（S5 上已实测过这条性质）。

用法: PYTHONPATH=<work>/pylibs python3 scripts/fix_fig5_v9.py [--apply]
产物  Fig5_v9.pdf（18.00 × 18.46 cm，与 10-21 交付件同尺寸）
"""
import os, sys, hashlib, shutil
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fitz
import numpy as np
from figkit import CM
from network_key import draw_row

FIG = os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures'))
SRC = f'{FIG}/FINAL_SET_20261021/_superseded/Figure_5_20261021_deprecated.pdf'
DST = f'{FIG}/Fig5_v9.pdf'
X, Y = 8.35, 0.03                      # 行左上角（cm）
ROW = (8.10, 0.00, 13.45, 0.58)        # 自检允许改动区 = 单行图例足迹（含余量）
NEW_TOKENS = ['TF', 'target gene', '1.0×', '4.0×']


def md5s(p):
    return hashlib.md5(open(p, 'rb').read()).hexdigest()


def main(apply=False):
    doc = fitz.open(SRC)
    pg = doc[0]
    ref = fitz.open(SRC)[0]          # ⚠ 基准必须是**没画过图例**的那一份：pg 会被就地画上图例
    assert abs(pg.rect.width / CM - 18.0) < 0.01 and abs(pg.rect.height / CM - 18.46) < 0.02
    base_lines = pg.get_text().splitlines()
    w = draw_row(pg, X, Y)
    assert X + w < 17.8, '单行图例越页'
    doc.subset_fonts()
    doc.save(DST, garbage=4, deflate=True)
    new = fitz.open(DST)[0]

    # ---- 1. 文字层：只多了这 4 个词，且没有 nbsp / 软连字符 ----
    t = new.get_text()
    assert '\xa0' not in t and '\xad' not in t, '文字层出现 \\xa0/\\xad'
    extra = sorted(set(t.splitlines()) - set(base_lines))
    print('文字层新增行: %s' % extra)
    assert set(NEW_TOKENS) <= set(extra), '四个标签没都写上'
    assert set(extra) <= set(NEW_TOKENS), '除图例标签外还多了别的字: %s' % (set(extra) - set(NEW_TOKENS))

    # ---- 2. 像素：改动必须全在单行图例足迹内 ----
    for dpi in (150, 300):
        A = new.get_pixmap(dpi=dpi, alpha=False)
        B = ref.get_pixmap(dpi=dpi, alpha=False)
        a = np.frombuffer(A.samples, np.uint8).reshape(A.height, A.width, A.n)[:, :, :3].astype(int)
        b = np.frombuffer(B.samples, np.uint8).reshape(B.height, B.width, B.n)[:, :, :3].astype(int)
        d = np.abs(a - b).max(2) > 24
        m = np.zeros_like(d)
        x0, y0, x1, y1 = ROW
        m[int(y0 * dpi / 2.54):int(y1 * dpi / 2.54), int(x0 * dpi / 2.54):int(x1 * dpi / 2.54)] = True
        print('  %3d dpi: 总变化 %d px，单行图例足迹以外 %d px' % (dpi, int(d.sum()), int((d & ~m).sum())))
        assert (d & ~m).sum() == 0, '图例以外有改动'
    print('  尺寸 %.2f × %.2f cm  体积 %.2f → %.2f MB' % (
        new.rect.width / CM, new.rect.height / CM,
        os.path.getsize(SRC) / 1e6, os.path.getsize(DST) / 1e6))
    new.get_pixmap(dpi=300, clip=fitz.Rect(7.6 * CM, 0, 18 * CM, 2.6 * CM)).save('/tmp/f5_v9_head.png')
    if apply:
        dst = f'{FIG}/FINAL_SET_20261021/Figure_5.pdf'
        print('  交付件 md5 %s（应为 2302687ee28b）' % md5s(dst)[:12])
        assert md5s(dst).startswith('2302687ee28b'), '要覆盖的不是本轮 G4 版'
        shutil.copy2(DST, dst)
        print('  安装 Figure_5.pdf  md5 %s' % md5s(dst)[:12])
    else:
        print('（试运行）产物', DST)


if __name__ == '__main__':
    main(apply='--apply' in sys.argv)
