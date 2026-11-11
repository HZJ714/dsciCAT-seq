#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
network_key.py — SCENIC+ GRN 网络图的**边方向/权重图例块**（Codex #91 / G4）
============================================================================
Codex 原话（`CODEX_COMMENTS_EXTRACT_20261028.txt` #91）：
  「网络标签/连线交叠，**边无方向/权重图例**；中心 Mbnl2 多次出现为 TF 样节点，需优先审计注释。」

本模块只画一个图例块，**不动任何网络几何**（节点坐标、边、标签全部保持交付版原样）——
重画网络的收益是零、风险是版面全变，而审稿人要的是"看懂这张图的钥匙"。

图例内容 = 从成品 PDF 矢量里**实测出来的**编码，不是从代码里读出来的意图：
  · 边宽：原图 1.00–4.00 pt（`width_edge_by = TF2G_importance, min 1 / max 4`），
    落到 8.70 cm 槽位（×0.3602）= **0.36–1.44 pt**；实测四条网络的宽度区间正是 1.00–4.00。
  · 边色：每个源 TF 一个 tab10 色（实测每图 5/5/3/2 种颜色，与 TF 数一致）。
  · 透明度：配置里写了 `transparency_edge_by`，**实测成品里没有任何 ExtGState/alpha**
    ⇒ 图例只说"颜色 = 源 TF、粗细 = importance"，**不写"深浅/透明度也编码"**（那是代码意图，
      不是图面事实；本轮全部改动的纪律就是"图上没有的话不写"）。
  · 节点：大菱形 = TF（粗体名），小菱形 = 靶基因；实测 0.611/0.417 cm（原图）→ 0.22/0.15 cm（槽位）。
    ⚠ `plot_networkx()` 把所有节点都画成 'D'（菱形），`shape_node_by` 的 circle/ellipse **没有生效**
    ⇒ 图例用菱形，不用"圆/椭圆"。
  · 方向：`create_nx_graph` 用 `nx.from_pandas_edgelist` ⇒ **图是无向的、没有箭头**；
    但每条边的两端类型唯一（TF ↔ gene），方向由节点类型决定 ⇒ 图例用"大菱形 → 小菱形"表达。
  · 符号：四张网络的 `regulation` 列实测全为 1 ⇒ "全部为预测的激活连接"。

用法:
  from network_key import draw_key
  draw_key(page, x=6.7, y=0.6, maxw=5.0, imp=(0.6, 4.2))
自检: python3 network_key.py            → /tmp/key_preview.png + 重叠自检
"""
import os, sys
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fitz
from figkit import CM, FONT, write_spans, ink_h

TAB10 = [(0.122, 0.467, 0.706), (1.0, 0.498, 0.055), (0.173, 0.627, 0.173),
         (0.839, 0.153, 0.157), (0.58, 0.404, 0.741)]
GENE_MID = (0.54, 0.54, 1.0)          # bwr 中段偏蓝，实测基因节点常见色
BLK = (0.0, 0.0, 0.0)
TF_D, GENE_D = 0.22, 0.15             # 菱形对角宽 cm（= 原图 0.611/0.417 × 槽位 0.3602）
W_LO, W_HI = 0.36, 1.44               # 印刷后的边宽 pt（原图 1.00–4.00 × 0.3602）


def tw(t, z, bold=False):
    return fitz.Font.text_length(FONT['b' if bold else 'r'], t, fontsize=z) / CM


def _wrap(text, maxw, z):
    out, cur = [], ''
    for word in text.split(' '):
        t = (cur + ' ' + word).strip()
        if tw(t, z) <= maxw or not cur:
            cur = t
        else:
            out.append(cur); cur = word
    if cur:
        out.append(cur)
    return out


def draw_key(page, x, y, maxw, imp=(0.5, 5.0), z=6.0, pitch=0.255, gap_glyph=0.10,
             color=BLK, verbose=True):
    """(x,y) = 图例块左上角 (cm)；块宽不超过 maxw。返回 (实测宽, 高) cm。"""
    items, cmds = [], []
    cy = y + ink_h('H', z) / 2
    # ---- 行 1：方向 = 大菱形（TF）→ 小菱形（靶基因）----
    t1a, t1b = 'TF', 'target gene'
    cx = x
    cmds.append(('dia', cx + TF_D / 2, cy, TF_D, TAB10[0])); cx += TF_D + gap_glyph
    w = tw(t1a, z); items.append(dict(t=t1a, x=cx + w / 2, y=cy, z=z, c=color)); cx += w + 0.14
    cmds.append(('arrow', cx, cx + 0.52, cy)); cx += 0.52 + 0.14
    cmds.append(('dia', cx + GENE_D / 2, cy, GENE_D, GENE_MID)); cx += GENE_D + gap_glyph
    w = tw(t1b, z); items.append(dict(t=t1b, x=cx + w / 2, y=cy, z=z, c=color))
    used = cx + w
    row_w = used - x
    cy += pitch
    # ---- 行 2：权重 = 三条样例线（印刷真实宽度）----
    for i, lw in enumerate([W_LO, (W_LO + W_HI) / 2, W_HI]):
        cmds.append(('line', x + i * 0.48, x + i * 0.48 + 0.42, cy, lw))
    tx = x + 3 * 0.48 + 0.04
    ln2 = _wrap('TF2G importance, low to high', maxw - (tx - x), z)
    items.append(dict(t=ln2[0], x=tx + tw(ln2[0], z) / 2, y=cy, z=z, c=color))
    row_w = max(row_w, tx + tw(ln2[0], z) - x)
    for k, ln in enumerate(ln2[1:], start=1):
        cy += pitch
        items.append(dict(t=ln, x=x + tw(ln, z) / 2, y=cy, z=z, c=color))
        row_w = max(row_w, tw(ln, z))
    cy += pitch
    # ---- 行 3–5：其余编码 ----
    for t in ['edge colour = the source TF (eRegulon)',
              'all edges = predicted activating links',
              'gene node colour = state log2FC (blue low, red high)']:
        for ln in _wrap(t, maxw, z):
            w = tw(ln, z)
            items.append(dict(t=ln, x=x + w / 2, y=cy, z=z, c=color))
            row_w = max(row_w, w)
            cy += pitch
    # ---- 落笔 ----
    sh = page.new_shape()
    for c in cmds:
        if c[0] == 'line':
            _, x0, x1, yy, lw = c
            sh.draw_line(fitz.Point(x0 * CM, yy * CM), fitz.Point(x1 * CM, yy * CM))
            sh.finish(color=TAB10[0], width=lw)
        elif c[0] == 'arrow':
            _, x0, x1, yy = c
            sh.draw_line(fitz.Point(x0 * CM, yy * CM), fitz.Point((x1 - 0.15) * CM, yy * CM))
            sh.finish(color=color, width=0.5)
            sh.draw_polyline([fitz.Point((x1 - 0.16) * CM, (yy - 0.06) * CM),
                              fitz.Point(x1 * CM, yy * CM),
                              fitz.Point((x1 - 0.16) * CM, (yy + 0.06) * CM),
                              fitz.Point((x1 - 0.16) * CM, (yy - 0.06) * CM)])
            sh.finish(color=None, fill=color, width=0)
        else:
            _, ccx, cyy, d, fill = c
            h = d / 2
            sh.draw_polyline([fitz.Point(ccx * CM, (cyy - h) * CM), fitz.Point((ccx + h) * CM, cyy * CM),
                              fitz.Point(ccx * CM, (cyy + h) * CM), fitz.Point((ccx - h) * CM, cyy * CM),
                              fitz.Point(ccx * CM, (cyy - h) * CM)])
            sh.finish(color=None, fill=fill, width=0)
    sh.commit()
    write_spans(page, items)
    h = cy - y - pitch + ink_h('H', z)
    if verbose:
        print('  network_key @ (%.2f,%.2f) maxw=%.2f → 实测 %.2f × %.2f cm，%d 条文字'
              % (x, y, maxw, row_w, h, len(items)))
    return row_w, h


def draw_row(page, x, y, z=6.0, color=BLK, lo=W_LO, hi=W_HI, verbose=True):
    """**单行**版（Figure 5 用）。作者："下面的小字有点多余，写图注里，影响观感；
    把图例一排放 MOL3 的上面、和 B 平齐。" ⇒ 只留图形符号 + 四个短标签，四行说明进图注。
    (x,y) = 行左上角（与 draw_key 同义：y 是墨迹上沿），返回实测宽 cm。
    排布： ◇ TF ──▸ ◇ target gene   ─ 1.0×  ═ 4.0×      （线宽 = 印刷真实宽度）"""
    cy = y + ink_h('H', z) / 2
    items, cmds, cx = [], [], x
    for t, d, fill, gap in (('TF', TF_D, TAB10[0], 0.10), ('target gene', GENE_D, GENE_MID, 0.10)):
        cmds.append(('dia', cx + d / 2, cy, d, fill)); cx += d + gap
        w = tw(t, z); items.append(dict(t=t, x=cx + w / 2, y=cy, z=z, c=color)); cx += w
        if t == 'TF':
            cx += 0.14
            cmds.append(('arrow', cx, cx + 0.52, cy)); cx += 0.52 + 0.14
    cx += 0.26
    for lw, lab in ((lo, '1.0\u00d7'), (hi, '4.0\u00d7')):
        cmds.append(('line', cx, cx + 0.42, cy, lw)); cx += 0.42 + 0.05
        w = tw(lab, z); items.append(dict(t=lab, x=cx + w / 2, y=cy, z=z, c=color)); cx += w + 0.10
    sh = page.new_shape()
    for c in cmds:
        if c[0] == 'line':
            _, x0, x1, yy, lw = c
            sh.draw_line(fitz.Point(x0 * CM, yy * CM), fitz.Point(x1 * CM, yy * CM))
            sh.finish(color=TAB10[0], width=lw)
        elif c[0] == 'arrow':
            _, x0, x1, yy = c
            sh.draw_line(fitz.Point(x0 * CM, yy * CM), fitz.Point((x1 - 0.15) * CM, yy * CM))
            sh.finish(color=color, width=0.5)
            sh.draw_polyline([fitz.Point((x1 - 0.16) * CM, (yy - 0.06) * CM),
                              fitz.Point(x1 * CM, yy * CM),
                              fitz.Point((x1 - 0.16) * CM, (yy + 0.06) * CM),
                              fitz.Point((x1 - 0.16) * CM, (yy - 0.06) * CM)])
            sh.finish(color=None, fill=color, width=0)
        else:
            _, ccx, cyy, d, fill = c
            h = d / 2
            sh.draw_polyline([fitz.Point(ccx * CM, (cyy - h) * CM), fitz.Point((ccx + h) * CM, cyy * CM),
                              fitz.Point(ccx * CM, (cyy + h) * CM), fitz.Point((ccx - h) * CM, cyy * CM),
                              fitz.Point(ccx * CM, (cyy - h) * CM)])
            sh.finish(color=None, fill=fill, width=0)
    sh.commit()
    write_spans(page, items)
    if verbose:
        print('  network_key(单行) @ (%.2f,%.2f) → 宽 %.2f cm，%d 条文字' % (x, y, cx - x - 0.10, len(items)))
    return cx - x - 0.10


# ---------------------------------------------------------------- 自检
if __name__ == '__main__':
    doc = fitz.open()
    pg = doc.new_page(width=18 * CM, height=9 * CM)
    draw_key(pg, 0.3, 0.3, 5.0, imp=(0.6, 4.2))
    draw_key(pg, 6.5, 0.3, 4.4, imp=(0.5, 5.0))
    draw_key(pg, 12.5, 0.3, 3.2, imp=(0.6, 4.2))
    doc.save('/tmp/key_preview.pdf')
    pg.get_pixmap(dpi=350, clip=fitz.Rect(0, 0, 18 * CM, 8 * CM)).save('/tmp/key_preview.png')
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import pdfix
    n, bad = pdfix.check('/tmp/key_preview.pdf')
    print('自检: %d 词, 重叠对 %d' % (n, len(bad)))
    for a, b, ar in bad[:6]:
        print('   %r <> %r  %.2f mm2' % (a['t'], b['t'], ar * 100))
