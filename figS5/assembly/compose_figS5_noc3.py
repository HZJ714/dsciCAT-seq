#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compose_figS5_noc3.py — FigS5 第三版：**面板顺序改成与主图 Fig 5 逐位对齐**
==============================================================================
为什么出第三版 (2026-10-16, 详证 `S5_PANEL_ORDER_REVIEW_20261016.md`):
  noc2 的 A/B/C/D 与主图 5 的 A/B/C/D **一格都没对上** —— 读者按字母找支撑会错位:
    主图 A 状态特异 eRegulon  <- 实际支撑是 noc2 的 B (eRSS top-10 排序)
    主图 B MOL3 GRN           <- 实际支撑是 noc2 的 D (另三态网络)
    主图 C marker 闭合         <- 实际支撑是 noc2 的 C (Nrxn3 位点轨道, 5C 的 example)
    主图 D 跨技术可复现        <- noc2 里没有对应块 (旧 c/d 已删)
  本版排法 = **B -> A, D -> B, C -> C, A -> D**:
    新 A  eRSS top-10 ×6 状态      <-> 5A
    新 B  三张状态 GRN 网络         <-> 5B
    新 C  Nrxn3 位点轨道 + 表达      <-> 5C (正文点名的 example)
    新 D  降采样稳定性 (QC)          <-> 5D —— 5D 问"换个平台还一样吗",
        本块问"少一半细胞还一样吗", 同一问句的两个轴, 不是垃圾桶。

相对 noc2 的三处实质改动:
  1) 顺序重排 (同上), 页面尺寸**不变** 18.00 × 46.75 cm;
  2) 稳定性块不再从 `FigS5.pdf` 整页裁贴, 改由 `figs5_stability.py` 按槽位真实尺寸
     重画后原尺寸贴入 -> 一次消掉 ①**5,057 个隐藏字** (裁贴带进来的源页文字算子,
     含两份面板字母) 与 ②三类真压字 (196.00⨯刻度200 = 2.2×5.1pt、294.00⨯子图标题
     = 10.9×2.6pt、15 个上标 3 ⨯ 10 的 0 = 1.9–3.1pt —— 根因是源图缩 0.477 而文字
     被抬回 6pt);
  3) 网络块第二行的 MOL6 由左对齐改**居中** (换序后它落在页面中段, 左上+右上+左下
     的 Z 形空角会显得像缺一块)。
字体/字号/槽位宽/网络消解逻辑与 noc2 完全一致 (可见字 6–8pt Arial, 面板字母 8pt Bold)。
输出: <work>/figures/FigS5_noc3.pdf (中间产物 replot/S5d_*_bare.pdf 与 noc2 同名共用)
用法: 先跑 figs5_stability.py, 再 python3 compose_figS5_noc3.py
"""
import os, sys, collections
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
for _d in ('00_Figure_Tools', os.path.join('fig5', 'assembly')):
    _p = os.path.join(_ROOT, _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)
import fitz
from figkit import (CM, FONT, cover, letter, save_doc, write_spans, strip_text_stream,
                    CAP, XH, DESC, UPCH, DNCH)
import pdfix

B = os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures'))
REP = f'{B}/replot'
PSC = os.environ.get('DSCI_SCENIC_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pi/plots-scenic+'))
OUT = f'{B}/FigS5_noc3.pdf'
W = 18.0                                        # 版面宽 (cm)

# ---- 纵向流水 (与 noc2 同一批块高, 只是换了先后) ----
A_Y0 = 0.34                                     # A: eRSS top-10 网格 (17.60 × 11.20)
B_T1 = 11.98                                    # B: 第一行状态名标题墨迹中心
B_N1 = 12.20                                    # B: 第一行图顶
B_T2 = 20.90                                    # B: 第二行标题
B_N2 = 21.12
C_Y0 = 30.02                                    # C: Nrxn3 track (18.00 × 11.20)
D_Y0 = 41.82                                    # D: 降采样稳定性 (18.00 × 4.90)
LETTER_SZ = 8.0
NET_W = 8.70                                    # 网络槽位宽
NET_CROP = (2.95, 3.55, 27.10, 26.60)           # 原图内容区 (裁掉四周留白)
NET_H = (NET_CROP[3] - NET_CROP[1]) * NET_W / (NET_CROP[2] - NET_CROP[0])
NETS = [('grn_mol2_network.pdf', 'MOL2', 0.10, B_T1, B_N1),
        ('grn_mol3_6_network.pdf', 'MOL3_6', 9.20, B_T1, B_N1),
        ('grn_mol6_network.pdf', 'MOL6', (W - NET_W) / 2, B_T2, B_N2)]   # 第二行居中
Z_GENE, Z_TF = 6.0, 8.0                         # 基因标签 / TF(hub) 标签字号
LEAD_TH = 0.06                                  # 位移超过 0.06cm 才画导线


def tw_cm(t, z, bold=False):
    return fitz.Font.text_length(FONT['b' if bold else 'r'], t, fontsize=z) / CM


def ink_h_of(t, z):
    em = z / 72 * 2.54
    return ((CAP if any(c in UPCH for c in t) else XH)
            + (DESC if any(c in DNCH for c in t) else 0.0)) * em


def box_edge_t(cx, cy, x0, y0, x1, y1, tx, ty):
    """从盒中心 (cx,cy) 朝目标 (tx,ty) 出发, 返回盒边上的参数 t (0=中心, 1=目标)"""
    dx, dy = tx - cx, ty - cy
    ts = []
    for ex in (x1, x0):
        if dx:
            t = (ex - cx) / dx
            if t > 0 and min(y0, y1) - 0.02 <= cy + dy * t <= max(y0, y1) + 0.02:
                ts.append(t)
    for ey in (y1, y0):
        if dy:
            t = (ey - cy) / dy
            if t > 0 and min(x0, x1) - 0.02 <= cx + dx * t <= max(x0, x1) + 0.02:
                ts.append(t)
    return min(ts) if ts else 0.0


def network_layer(pg, src, x0, y0, w_slot, ttl, ttl_y):
    """贴原网络的**图形** (去文字) + 6/8pt 重写文字 (墨迹消解) + 状态名标题"""
    s = w_slot / (NET_CROP[2] - NET_CROP[0])
    h_slot = (NET_CROP[3] - NET_CROP[1]) * s
    doc = fitz.open(src)
    ps = doc[0]
    ws = pdfix.words(ps)
    items, anchors = [], []
    for w in ws:
        if not (NET_CROP[0] <= w['x0'] and w['x1'] <= NET_CROP[2]
                and NET_CROP[1] <= w['y0'] and w['y1'] <= NET_CROP[3]):
            continue
        tf = w['z'] > 15                        # 原图: TF 名 22pt, 基因名 11pt
        z = Z_TF if tf else Z_GENE
        cx = x0 + ((w['x0'] + w['x1']) / 2 - NET_CROP[0]) * s
        cy = y0 + ((w['y0'] + w['y1']) / 2 - NET_CROP[1]) * s
        items.append(dict(t=w['t'], x=cx, y=cy, z=z, w=tw_cm(w['t'], z, tf),
                          h=ink_h_of(w['t'], z), c=w['c'], bold=tf))
        anchors.append((cx, cy))
    us = [dict(idx=[i], box=(it['x'] - it['w'] / 2, it['y'] - it['h'] / 2,
                            it['x'] + it['w'] / 2, it['y'] + it['h'] / 2),
               dx=0.0, dy=0.0) for i, it in enumerate(items)]
    left = pdfix.resolve_units(us, list(range(len(us))),
                               bounds=(x0 - 0.06, y0 - 0.06, x0 + w_slot + 0.06,
                                       y0 + h_slot + 0.06), gap=0.012, iters=6000)
    bare = f'{REP}/S5d_{os.path.basename(src).replace(".pdf", "")}_bare.pdf'
    strip_text_stream(doc, ps)
    save_doc(doc, bare)
    pg.show_pdf_page(fitz.Rect(x0 * CM, y0 * CM, (x0 + w_slot) * CM, (y0 + h_slot) * CM),
                     fitz.open(bare), 0, clip=fitz.Rect(*[v * CM for v in NET_CROP]))
    sh = pg.new_shape()
    nlead = 0
    for i, u in enumerate(us):
        items[i]['x'] += u['dx']; items[i]['y'] += u['dy']
        if (u['dx'] ** 2 + u['dy'] ** 2) ** 0.5 < LEAD_TH:
            continue
        nlead += 1
        ax, ay = anchors[i]
        it = items[i]
        t = box_edge_t(it['x'], it['y'], it['x'] - it['w'] / 2, it['y'] - it['h'] / 2,
                       it['x'] + it['w'] / 2, it['y'] + it['h'] / 2, ax, ay)
        sh.draw_line(fitz.Point((it['x'] + (ax - it['x']) * t) * CM,
                                (it['y'] + (ay - it['y']) * t) * CM),
                     fitz.Point(ax * CM, ay * CM))
    sh.finish(color=(0.60, 0.60, 0.60), width=0.2)
    sh.commit()
    write_spans(pg, items)
    write_spans(pg, [dict(t=ttl, x=x0 + w_slot / 2, y=ttl_y, z=7.0,
                          w=tw_cm(ttl, 7.0, True), h=ink_h_of(ttl, 7.0), bold=True)])
    print(f'  {os.path.basename(src):26s} -> {ttl:7s} 标签 {len(items):3d} '
          f'(TF {sum(1 for v in items if v["bold"])}) 移动 {nlead:3d} 残余重叠 {left}')
    return h_slot


def paste(pg, src, y0, x0=None, label=None, letter_dy=-0.32):
    """原尺寸贴入 (源件本身就是槽位尺寸, 不缩放)"""
    d = fitz.open(src)
    bw, bh = d[0].rect.width / CM, d[0].rect.height / CM
    x = (W - bw) / 2 if x0 is None else x0
    pg.show_pdf_page(fitz.Rect(x * CM, y0 * CM, (x + bw) * CM, (y0 + bh) * CM), d, 0)
    if label:
        ly = 0.02 if y0 < 1 else y0 + letter_dy
        letter(pg, 0.02, ly, label, LETTER_SZ)
    print(f'  {os.path.basename(src):30s} -> {label}  y=[{y0:.2f},{y0 + bh:.2f}]cm  {bw:.2f}×{bh:.2f}')
    return y0 + bh


def main():
    out = fitz.open()
    pg = out.new_page(width=W * CM, height=46.75 * CM)
    paste(pg, f'{REP}/S5b_rss_grid_wide.pdf', A_Y0, label='A')          # 新 A: eRSS
    for f, ttl, x0, ty, ny in NETS:                                     # 新 B: 三态网络
        network_layer(pg, f'{PSC}/{f}', x0, ny, NET_W, ttl, ty)
    letter(pg, 0.02, B_T1 - 0.36, 'B', LETTER_SZ)
    paste(pg, f'{REP}/S5c_nrxn3_track_box.pdf', C_Y0, x0=0.0, label='C')  # 新 C: Nrxn3
    paste(pg, f'{REP}/S5_stability_wide.pdf', D_Y0, x0=0.0, label='D')    # 新 D: 降采样
    save_doc(out, OUT, clean=False)
    p0 = fitz.open(OUT)[0]
    sizes = collections.Counter()
    for b in p0.get_text('dict')['blocks']:
        for l in b.get('lines', []):
            for s in l['spans']:
                if s['text'].strip():
                    sizes[round(s['size'], 1)] += len(s['text'].strip())
    print(f'{os.path.basename(OUT)}  {W:.2f} × {p0.rect.height / CM:.2f} cm | 可见字号 {dict(sorted(sizes.items()))}')


if __name__ == '__main__':
    main()
