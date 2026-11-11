#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""点径普查（作者 2026-10-21 定案：UMAP 点 0.45 pt / feature plot 1.60 pt / dotplot 大小是数据不动）
用法: python3 diaudit.py FILE[:区域] ...   区域 = all | top(y<5.77) | dpanel(x<4.65,y>5.77) | cm
口径与 scripts/make_word_preview_derived_20261020.py 的 collect() 完全一致（同一套判定条件），
只是不改文件、只报数。
"""
import collections, os, sys
import fitz
CM = 72 / 2.54


def census(page, region='all', lo=0.02, hi=4.0):
    """单位与 make_word_preview_derived_20261020.collect() 一致：径长用 **pt**（=PDF 点=1/72 in），
    区域框用 **cm**（与作者定案里写的 0.45 pt / 1.60 pt 同一把尺）。"""
    d = collections.Counter()
    tot = 0
    for dr in page.get_drawings():
        r = dr['rect']
        w_pt = max(r.width, r.height)
        if not (lo < w_pt < hi) or not (0.7 < r.width / max(r.height, 1e-6) < 1.4):
            continue
        if dr.get('fill') is None:
            continue
        cx, cy = (r.x0 + r.x1) / 2 / CM, (r.y0 + r.y1) / 2 / CM
        if region == 'top' and not cy < 5.77:
            continue
        if region == 'dpanel' and not (cx < 4.65 and cy > 5.77):
            continue
        d[round(w_pt, 3)] += 1
        tot += 1
    return d, tot


def report(path, regions=('all',)):
    if not os.path.exists(path):
        print('缺失 %s' % path); return
    doc = fitz.open(path); p = doc[0]
    print('\n== %s  (%.2f x %.2f cm)' % (os.path.basename(path), p.rect.width / CM, p.rect.height / CM))
    for reg in regions:
        d, tot = census(p, reg)
        if not tot:
            print('   [%s] 圆点 0 个' % reg); continue
        top = sorted(d.items(), key=lambda kv: -kv[1])[:6]
        modal = top[0]
        print('   [%s] 圆点 %d 个 | 众数径 %.3f pt (x%d, %.1f%%) | 前几名 %s'
              % (reg, tot, modal[0], modal[1], 100 * modal[1] / tot,
                 ', '.join('%.3f×%d' % (k, v) for k, v in top[:5])))
        bands = {'<0.20 pt(糊)': sum(v for k, v in d.items() if k < 0.20),
                 '0.45±0.03': sum(v for k, v in d.items() if abs(k - 0.45) <= 0.03),
                 '1.6±0.05': sum(v for k, v in d.items() if abs(k - 1.60) <= 0.05)}
        print('        分档: ' + ' | '.join('%s %d (%.1f%%)' % (b, c, 100 * c / tot) for b, c in bands.items()))
    doc.close()


if __name__ == '__main__':
    for spec in sys.argv[1:]:
        if ':' in spec:
            f, reg = spec.split(':', 1)
            report(f, tuple(r.strip() for r in reg.split(',')))
        else:
            report(spec)
