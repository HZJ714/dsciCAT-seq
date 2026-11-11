#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pdfix.py — 成品 PDF 内的「墨迹空间」标签消解 + 精确重写
================================================================================
为什么需要它: matplotlib 的 Agg 渲染器算出的文本行盒, 与 PDF 后端真正落笔的位置差
~1.3-1.7pt (且随 va/字号变化), 所以在 matplotlib 里"消解到 0 重叠"的标签, 存成 PDF
后仍会残留 2mm² 级别的擦碰。**唯一可靠的坐标系是成品 PDF 本身** —— 这里读它的逐字符
墨迹盒, 在 cm 空间推开, 再按原基线 origin 精确写回。

流程: words() 读词 (精确墨迹盒 + 基线 origin + 方向/色/字号)
      -> resolve() 分组消解 (可动/固定分开, 软边界, 可加障碍盒)
      -> rewrite() 整页只删文字 (figkit.strip_text_stream, 矢量与图像全保留)
                 + 按新 origin 精确重写 + 需要时补画导线 (锚点 -> 文字盒边缘)。
单位: 一律 cm, y 向下, 原点 = 页面左上角。
"""
import os, sys, math
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fitz
from figkit import CM, FONT, strip_text_stream, save_doc
from figaudit import words_of, collide, CAP, XH, DESC, UPCH, DNCH


def rgb(c):
    return ((c >> 16 & 255) / 255, (c >> 8 & 255) / 255, (c & 255) / 255)


def words(page):
    """page -> 词表 [{t,x0,y0,x1,y1,ox,oy,z,c,bold,dir,vert,lid}]
    lid = 行号 (同一个 matplotlib Text 的多段/多行共用一个 lid, 消解时整体移动)"""
    out = []
    d = page.get_text('rawdict', flags=fitz.TEXT_PRESERVE_WHITESPACE)
    lid = -1
    for b in d['blocks']:
        if b['type'] != 0:
            continue
        for l in b['lines']:
            lid += 1
            dirv = tuple(round(v, 3) for v in l.get('dir', (1, 0)))
            vert = abs(dirv[1]) > 0.5
            cur = []
            def flush():
                if not cur:
                    return
                t = ''.join(g['c'] for g in cur)
                xs = [g['bbox'][0] for g in cur] + [g['bbox'][2] for g in cur]
                ys = [g['bbox'][1] for g in cur] + [g['bbox'][3] for g in cur]
                z = max(g['size'] for g in cur)
                em = z / 72 * 2.54
                up = (CAP if any(c in UPCH for c in t) else XH) if not vert else 0
                dn = (DESC if any(c in DNCH for c in t) else 0.0) if not vert else 0
                base = min(g['origin'][1] for g in cur) / CM
                w = dict(t=t, z=z, c=rgb(cur[0]['color']), lid=lid,
                         bold='bold' in cur[0]['font'].lower(), dir=dirv, vert=vert,
                         ox=min(g['origin'][0] for g in cur) / CM, oy=base)
                if vert:                            # 竖排: 直接用旋转后的外接盒
                    w.update(x0=min(xs) / CM, x1=max(xs) / CM,
                             y0=min(ys) / CM, y1=max(ys) / CM)
                else:
                    w.update(x0=min(xs) / CM, x1=max(xs) / CM,
                             y0=base - up * em, y1=base + dn * em)
                out.append(w)
                cur.clear()
            for s in l['spans']:
                for ch in s['chars']:
                    if ch['c'].strip() == '':
                        flush(); continue
                    if cur and ch['bbox'][0] - cur[-1]['bbox'][2] > 0.35 * s['size']:
                        flush()
                    cur.append(dict(c=ch['c'], bbox=ch['bbox'], origin=ch['origin'],
                                    size=s['size'], color=s['color'], font=s['font']))
            flush()
    return out


def resolve(ws, movable, bounds=None, obstacles=(), gap=0.02, iters=1200, push=0.6):
    """原地给 ws[i] 写 dx/dy。movable=可动下标; obstacles=[(x0,y0,x1,y1)] 额外固定盒。"""
    movable = set(movable)
    for w in ws:
        w.setdefault('dx', 0.0); w.setdefault('dy', 0.0)
    obs = [tuple(o) for o in obstacles]

    def bx(w):
        return (w['x0'] + w['dx'] - gap, w['y0'] + w['dy'] - gap,
                w['x1'] + w['dx'] + gap, w['y1'] + w['dy'] + gap)

    def clamp(w):
        if not bounds:
            return
        x0, y0, x1, y1 = bounds
        bw, bh = w['x1'] - w['x0'], w['y1'] - w['y0']
        if bw > x1 - x0 or bh > y1 - y0:
            return
        w['dx'] = max(w['dx'], x0 - w['x0']); w['dx'] = min(w['dx'], x1 - w['x1'])
        w['dy'] = max(w['dy'], y0 - w['y0']); w['dy'] = min(w['dy'], y1 - w['y1'])

    left = -1
    for _ in range(iters):
        for i in movable:
            clamp(ws[i])
        boxes = [bx(w) for w in ws] + obs
        bad = []
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                a, b = boxes[i], boxes[j]
                if not (i in movable or j in movable):
                    continue                       # 固定-固定: 不是本算法能解决的
                ox = min(a[2], b[2]) - max(a[0], b[0])
                oy = min(a[3], b[3]) - max(a[1], b[1])
                if ox > 0 and oy > 0:
                    bad.append((i, j, ox, oy))
        if not bad:
            left = 0
            break
        for i, j, ox, oy in bad:
            ai, bj = i in movable, j in movable
            d = max(min(ox, oy) * push, 0.006)
            A, B = boxes[i], boxes[j]
            horiz = ox <= oy
            sg = 1 if ((A[0] + A[2]) >= (B[0] + B[2]) if horiz
                       else (A[1] + A[3]) >= (B[1] + B[3])) else -1
            key = 'dx' if horiz else 'dy'
            if ai and bj:
                ws[i][key] += sg * d / 2
                ws[j][key] -= sg * d / 2
            elif ai:
                ws[i][key] += sg * d
            else:
                ws[j][key] -= sg * d
        return left


def rewrite(src, dst, ws, page_index=0, leaders=None, lcolor=(0.45, 0.45, 0.45),
            lw=0.25, lthresh=0.035, lmove_only=True):
    """按 ws[i]['dx'/'dy'] 重写文字层 (整页只删文字, 图形保留)。
    leaders: [(锚点x, 锚点y, 词下标)] -> 从词盒边缘画到锚点的细线。"""
    doc = fitz.open(src)
    page = doc[page_index]
    strip_text_stream(doc, page)
    hor = [w for w in ws if not w['vert']]
    ver = [w for w in ws if w['vert']]
    groups = {}
    for w in hor:
        groups.setdefault((w['c'], w['bold']), []).append(w)
    for (col, bold), gs in groups.items():
        f = FONT['b' if bold else 'r']
        tw = fitz.TextWriter(page.rect, color=col)
        for w in gs:
            tw.append(fitz.Point((w['ox'] + w.get('dx', 0.0)) * CM,
                                 (w['oy'] + w.get('dy', 0.0)) * CM), w['t'],
                      font=f, fontsize=w['z'])
        tw.write_text(page, color=col)
    for w in ver:                                   # 竖排: 横排写入后绕 origin 旋转
        f = FONT['b' if w['bold'] else 'r']
        ang = round(math.degrees(math.atan2(w['dir'][1], w['dir'][0])))
        tw = fitz.TextWriter(page.rect, color=w['c'])
        p0 = fitz.Point((w['ox'] + w.get('dx', 0.0)) * CM,
                        (w['oy'] + w.get('dy', 0.0)) * CM)
        tw.append(p0, w['t'], font=f, fontsize=w['z'])
        tw.write_text(page, color=w['c'], morph=(p0, fitz.Matrix(ang)))
    if leaders:
        sh = page.new_shape()
        for axp, ayp, i in leaders:
            w = ws[i]
            if lmove_only and abs(w.get('dx', 0)) < lthresh and abs(w.get('dy', 0)) < lthresh:
                continue
            x0, y0 = w['x0'] + w['dx'], w['y0'] + w['dy']
            x1, y1 = w['x1'] + w['dx'], w['y1'] + w['dy']
            cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
            ddx, ddy = axp - cx, ayp - cy
            if (ddx * ddx + ddy * ddy) ** 0.5 < lthresh:
                continue
            ts = []
            for ex in (x1, x0):
                if ddx:
                    t = (ex - cx) / ddx
                    if t > 0 and y0 - 0.02 <= cy + ddy * t <= y1 + 0.02:
                        ts.append(t)
            for ey in (y0, y1):
                if ddy:
                    t = (ey - cy) / ddy
                    if t > 0 and x0 - 0.02 <= cx + ddx * t <= x1 + 0.02:
                        ts.append(t)
            t = min(ts) if ts else 0.0
            sh.draw_line(fitz.Point((cx + ddx * t) * CM, (cy + ddy * t) * CM),
                         fitz.Point(axp * CM, ayp * CM))
        sh.finish(color=lcolor, width=lw)
        sh.commit()
    save_doc(doc, dst)
    return dst


def check(path, minarea=0.0008):
    """自检: (词数, 残余重叠对数) —— 用 figaudit 的独立墨迹口径"""
    doc = fitz.open(path)
    n, bad = 0, []
    for pg in doc:
        ws = words_of(pg)
        n += len(ws)
        bad += collide(ws, minarea=minarea)
    doc.close()
    return n, bad


# ─────────────────────────── 按"行单元"消解 ───────────────────────────
def units_of(ws, lids=None):
    """ws -> [{'idx':[词下标], 'box':(x0,y0,x1,y1)}]  (同一 lid 合成一个整体单元)"""
    byl = {}
    for i, w in enumerate(ws):
        byl.setdefault(w.get('lid', i), []).append(i)
    us = []
    for lid, idx in sorted(byl.items()):
        us.append(dict(idx=idx, lid=lid,
                       box=(min(w['x0'] for w in [ws[i] for i in idx]),
                            min(w['y0'] for w in [ws[i] for i in idx]),
                            max(w['x1'] for w in [ws[i] for i in idx]),
                            max(w['y1'] for w in [ws[i] for i in idx])),
                       dx=0.0, dy=0.0))
    return us


def resolve_units(us, movable, bounds=None, obstacles=(), gap=0.02, iters=1500,
                  push=0.6):
    """us: units_of() 结果 (原地改 dx/dy); movable: 可动单元下标集合。
    obstacles: [(x0,y0,x1,y1)] 或 [(x0,y0,x1,y1,owner)] —— owner 是"拥有"这个障碍的
    单元下标 (例如标签自己的节点圆点), 该单元可以压在它上面, 别的单元不行。"""
    movable = set(movable)
    obs, own = [], {}
    for o in obstacles:
        o = tuple(o)
        obs.append(o[:4])
        own[len(obs) - 1] = o[4] if len(o) > 4 else None

    def bx(u):
        return (u['box'][0] + u['dx'] - gap, u['box'][1] + u['dy'] - gap,
                u['box'][2] + u['dx'] + gap, u['box'][3] + u['dy'] + gap)

    def clamp(u):
        if not bounds:
            return
        x0, y0, x1, y1 = bounds
        w, h = u['box'][2] - u['box'][0], u['box'][3] - u['box'][1]
        if w > x1 - x0 or h > y1 - y0:
            return
        u['dx'] = max(u['dx'], x0 - u['box'][0]); u['dx'] = min(u['dx'], x1 - u['box'][2])
        u['dy'] = max(u['dy'], y0 - u['box'][1]); u['dy'] = min(u['dy'], y1 - u['box'][3])

    left = -1
    for _ in range(iters):
        for i in movable:
            clamp(us[i])
        boxes = [bx(u) for u in us] + obs
        nunit = len(us)
        bad = []
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                a, b = boxes[i], boxes[j]
                if not (i in movable or j in movable):
                    continue                       # 固定-固定: 不是本算法能解决的
                if i >= nunit and own.get(i) == j:      # j 自己的节点
                    continue
                if j >= nunit and own.get(j) == i:      # i 自己的节点
                    continue
                ox = min(a[2], b[2]) - max(a[0], b[0])
                oy = min(a[3], b[3]) - max(a[1], b[1])
                if ox > 0 and oy > 0:
                    bad.append((i, j, ox, oy))
        if not bad:
            left = 0
            break
        for i, j, ox, oy in bad:
            ai, bj = i in movable, j in movable
            d = max(min(ox, oy) * push, 0.006)
            A, B = boxes[i], boxes[j]
            horiz = ox <= oy
            sg = 1 if ((A[0] + A[2]) >= (B[0] + B[2]) if horiz
                       else (A[1] + A[3]) >= (B[1] + B[3])) else -1
            key = 'dx' if horiz else 'dy'
            if ai and bj:
                us[i][key] += sg * d / 2; us[j][key] -= sg * d / 2
            elif ai:
                us[i][key] += sg * d
            else:
                us[j][key] -= sg * d
    return left


def apply_units(ws, us):
    for u in us:
        for i in u['idx']:
            ws[i]['dx'], ws[i]['dy'] = u['dx'], u['dy']
