import os
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""inkcheck.py — barcol v3 的第二遍：600 dpi「紧墨盒」重叠候选清单
================================================================================
为什么需要第二遍：barcol v3 的判据是 4-邻域连通墨块「真的贴上」(≤0.24pt)，对 6pt 小字
太严 —— FigS5_noc2 的 A 带实测 `196.00` 压在 `200` 刻度上 2.2×5.1pt，v3 只报 1 对。
本脚本把每个 span 渲染后取**紧墨 bounding box**，列出两两重叠 >tol(默认 0.3pt) 的候选。
注意：紧墨盒重叠 ≠ 笔画相接（字母 x 高度不同也可能盒套盒），所以这是**候选清单**，
最终判定仍回 barcol 的连通墨块或肉眼。
用法: python3 inkcheck.py a.pdf [b.pdf ...] [--tol 0.3] [--clip x0,y0,x1,y1(cm)]
"""
import sys, os
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
import fitz, numpy as np

DPI = 600
PT = 72 / 25.4          # pt/cm
SC = None                # px/pt = DPI/72，见 main


def spans_of(page, clip=None):
    out = []
    for b in page.get_text('dict', clip=clip)['blocks']:
        for l in b.get('lines', []):
            for s in l['spans']:
                t = s['text'].strip()
                if t:
                    out.append((t, fitz.Rect(s['bbox']), round(s['size'], 1), s['font']))
    return out


def tight_boxes(page, sps, sc):
    """逐 span 取紧墨盒：只点亮该 span 的着色路径不可行，改为整页渲染 + 按 span bbox 区域裁切取墨。"""
    pm = page.get_pixmap(dpi=DPI, colorspace=fitz.csGRAY)
    g = np.frombuffer(pm.samples, np.uint8).reshape(pm.height, pm.width)
    ink = g < 200
    boxes = []
    for t, r, z, f in sps:
        x0, y0 = r.x0 * sc, r.y0 * sc
        x1, y1 = r.x1 * sc, r.y1 * sc
        X0, Y0 = max(0, int(x0) - 3), max(0, int(y0) - 3)
        X1, Y1 = min(pm.width, int(x1) + 4), min(pm.height, int(y1) + 4)
        sub = ink[Y0:Y1, X0:X1]
        ys, xs = np.nonzero(sub)
        if len(xs) == 0:
            continue
        # 换算回 pt（页面坐标），与 PDF bbox 同系
        boxes.append((t, fitz.Rect((X0 + xs.min()) / sc, (Y0 + ys.min()) / sc,
                                   (X0 + xs.max() + 1) / sc, (Y0 + ys.max() + 1) / sc), z, f))
    return boxes


def main():
    files = [a for a in sys.argv[1:] if not a.startswith('--')]
    tol = 0.3
    clip = None
    for a in sys.argv[1:]:
        if a.startswith('--tol'):
            tol = float(a.split('=')[1]) if '=' in a else float(sys.argv[sys.argv.index(a) + 1])
        if a.startswith('--clip='):
            v = [float(x) * PT for x in a.split('=')[1].split(',')]
            clip = fitz.Rect(*v)
    for f in files:
        d = fitz.open(f)
        pg = d[0]
        sc = DPI / 72.0                # 页面坐标是 pt，别乘 25.4
        bx = tight_boxes(pg, spans_of(pg, clip), sc)
        hits = []
        for i in range(len(bx)):
            for j in range(i + 1, len(bx)):
                a, b = bx[i], bx[j]
                r = a[1] & b[1]
                if r.is_valid and r.width > tol and r.height > tol:
                    hits.append((r.width, r.height, a[0], b[0], round(a[2], 1), round(b[2], 1)))
        print('%-28s span %3d | 紧墨盒重叠>%0.1fpt: %d 对' % (os.path.basename(f), len(bx), tol, len(hits)))
        for w, h, t1, t2, z1, z2 in sorted(hits, key=lambda x: -x[0] * x[1]):
            print('    %-26r ⨯ %-26r  %.2f × %.2f pt  (%.1f/%.1fpt)' % (t1[:26], t2[:26], w, h, z1, z2))
        # 顺带报所有文字行的实际墨间距
        lines = sorted({round(b[1].y0, 1): b[0] for b in bx}.items())
        prev = None
        for y, t in lines:
            if prev is not None and y - prev < 6.0:
                pass
            prev = y


if __name__ == '__main__':
    main()
