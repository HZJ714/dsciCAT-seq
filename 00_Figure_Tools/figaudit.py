import os
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
figaudit.py — 成品 PDF 的文字体检: 字号窗口 (6-8pt) + 真实墨迹重叠
================================================================================
为什么不用 figkit.overlaps(): 它用 PyMuPDF 的 span bbox 中心推基线, 误差可达 0.05cm,
会把"只差一点"的标签判成重叠; 而 PyMuPDF 又会把同一基线上相邻的两个 text 合并成一个
span (如 "Zfp354c ext" + "Arnt2"), 得到横跨两者的假 bbox -> 假阳性。
本模块直接读 rawdict 的**逐字符 bbox + origin(基线)**, 按空格/大间隙切成词,
用 Liberation Sans/Arial 的 cap/x-height/descender 比例算真实墨迹盒 (cm), 再判碰。
用法: python3 figaudit.py a.pdf [b.pdf ...]
"""
import sys, os, collections
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
import fitz

CM = 72 / 2.54
CAP, XH, DESC = 0.716, 0.519, 0.212          # Arial/Liberation Sans em 比例
UPCH = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789bdfhklt/−+%(),.#')
DNCH = set('gjpqy_')


def words_of(page, clip=None):
    """page -> [{t, x0,y0,x1,y1 (cm, 墨迹盒), z (pt), c, bold}]  按词切分"""
    out = []
    d = page.get_text('rawdict', flags=fitz.TEXT_PRESERVE_WHITESPACE)
    for b in d['blocks']:
        if b['type'] != 0:
            continue
        for l in b['lines']:
            dirv = l.get('dir', (1, 0))
            cur = []
            for s in l['spans']:
                for ch in s['chars']:
                    c = ch['c']
                    if c.strip() == '':
                        if cur:
                            out.append(_mk(cur, dirv)); cur = []
                        continue
                    if cur:
                        gap = ch['bbox'][0] - cur[-1]['bbox'][2]
                        if gap > 0.35 * s['size']:
                            out.append(_mk(cur, dirv)); cur = []
                    cur.append(dict(c=c, bbox=ch['bbox'], size=s['size'],
                                    origin=ch['origin'], color=s['color'],
                                    font=s['font']))
            if cur:
                out.append(_mk(cur, dirv))
    if clip:
        out = [w for w in out if fitz.Rect(w['x0'] * CM, w['y0'] * CM,
                                           w['x1'] * CM, w['y1'] * CM).intersects(clip)]
    return out


def _mk(chars, dirv):
    t = ''.join(c['c'] for c in chars)
    col = ((chars[0]['color'] >> 16 & 255) / 255, (chars[0]['color'] >> 8 & 255) / 255,
           (chars[0]['color'] & 255) / 255)
    bold = 'bold' in chars[0]['font'].lower()
    z = max(c['size'] for c in chars)
    if abs(dirv[1]) > 0.5:                       # 竖排: 直接用旋转后的 char 外接盒
        return dict(t=t, x0=min(c['bbox'][0] for c in chars) / CM,
                    x1=max(c['bbox'][2] for c in chars) / CM,
                    y0=min(c['bbox'][1] for c in chars) / CM,
                    y1=max(c['bbox'][3] for c in chars) / CM,
                    z=z, vert=True, c=col, bold=bold)
    x0 = min(c['bbox'][0] for c in chars) / CM
    x1 = max(c['bbox'][2] for c in chars) / CM
    base = min(c['origin'][1] for c in chars) / CM        # 基线 (y 向下)
    em = z / 72 * 2.54
    up = CAP if any(c in UPCH for c in t) else XH
    dn = DESC if any(c in DNCH for c in t) else 0.0
    return dict(t=t, x0=x0, x1=x1, y0=base - up * em, y1=base + dn * em, z=z,
                vert=False, c=col, bold=bold)


def collide(ws, pad=0.0, minarea=0.0008):
    bad = []
    for i in range(len(ws)):
        for j in range(i + 1, len(ws)):
            a, b = ws[i], ws[j]
            ox = min(a['x1'], b['x1']) - max(a['x0'], b['x0']) - 2 * pad
            oy = min(a['y1'], b['y1']) - max(a['y0'], b['y0']) - 2 * pad
            if ox > 0 and oy > 0 and ox * oy > minarea:
                bad.append((a, b, ox * oy))
    return bad


def check(path, lo=6.0, hi=8.0, verbose=True):
    doc = fitz.open(path)
    sizes, allbad, nword = collections.Counter(), [], 0
    for pi, pg in enumerate(doc):
        ws = words_of(pg)
        nword += len(ws)
        for w in ws:
            sizes[round(w['z'], 1)] += len(w['t'])
        for a, b, ar in collide(ws):
            allbad.append((pi, a['t'], b['t'], round(ar, 4)))
    out_of_range = {z: n for z, n in sizes.items() if z < lo - 1e-6 or z > hi + 1e-6}
    if verbose:
        print(f'{os.path.basename(path):38s} 词 {nword:4d} | 字号 {dict(sorted(sizes.items()))}'
              f' | 越界 {out_of_range or "无"} | 重叠 {len(allbad)}')
        for pi, a, b, ar in allbad[:15]:
            print(f'    p{pi} {a!r} <> {b!r}  {ar} cm2')
    return sizes, out_of_range, allbad


if __name__ == '__main__':
    for f in sys.argv[1:]:
        check(f)
