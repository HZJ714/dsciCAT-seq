#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compose_fig1_v4.py — 组 Fig1_v4.pdf（版式几何全部沿用 scripts/fig1_layout.py）
=============================================================================
与 v3 的唯一结构差别：**面板 C 不再从 `Figure 1.pdf` 矢量裁贴**，而是放
`replot/C_v4_tracks.pdf`（fig1c_v4.py 从原始数据整块重画，11.85 x 6.53 cm = C 槽位 1:1）。
  → 于是 C 字母必须自己画（旧图里那个 C 在裁贴内部，现在裁贴没了）；
  → LTR_C_COVER 的白覆盖也不再需要（不压到任何墨迹）。
其余 B/D/E/F 换成 replot_v4.py 的 *_v4_* 源图（log10 轴标签、矩阵型方法 * 号、
E 右半右脊柱刻度、sci-CAR 鼠侧 ATAC 灰纱斜纹 + ‡）。

用法: <conda-env>/bin/python3 compose_fig1_v4.py
输出: <work>/figures/Fig1_v4.pdf + png/Fig1_v4.png
"""
import os
import sys
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
import fitz
sys.path.insert(0, os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'scripts'))
import figkit as K
import fig1_layout as FL

FIG = os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures'))
REP = os.path.join(FIG, 'replot')
PNG = os.path.join(FIG, 'png')
CM = 72 / 2.54
PW = FL.PW

SRC = {'A': f'{REP}/F1A_workflow.pdf',
       'B': f'{REP}/B_v4_HG14_RNA_ATAC_contamination.pdf',
       'C': f'{REP}/C_v4_tracks.pdf',                                # ← v4：自绘，不再裁贴
       'D': f'{REP}/D_v4_RNA_umis_genes.pdf',
       'E': f'{REP}/E_v4_ATAC_fragments.pdf',
       'F': f'{REP}/F_v4_RNA_similarity_heatmap.pdf'}
SLOTS = dict(FL.SLOTS)
SLOTS['C'] = FL.C_SLOT
LETTERS = FL.LETTERS
LEG = [('dsciCAT-seq-50k', '#AF2337'), ('dsciCAT-seq-100k', '#AF2337'),
       ('HT-scCAT-seq', '#2F3C28'), ('ISSAAC-seq(Droplet)', '#3A7CA5'),
       ('Paired-seq', '#EF476F'), ('sciCAR-seq', '#5F4BB6'),
       ('SUM-seq', '#FF6B35'), ('10xMultiome', '#ECC342')]


def put(page, src, x, y, w, h):
    sd = fitz.open(src); r = sd[0].rect
    sc = min(w * CM / r.width, h * CM / r.height)
    nw, nh = r.width * sc, r.height * sc
    ox, oy = x * CM + (w * CM - nw) / 2, y * CM + (h * CM - nh) / 2
    page.show_pdf_page(fitz.Rect(ox, oy, ox + nw, oy + nh), sd, 0)
    sd.close()
    return sc


def draw_legend(page, y, fs=6.5):
    widths = [K.text_width(t, fs) for t, _ in LEG]
    total = sum(widths) + len(LEG) * 0.45
    gap = (16.70 - total) / (len(LEG) - 1)
    if gap < 0.12:
        return draw_legend(page, y, fs - 0.5)
    x = 0.15
    for (t, c), w in zip(LEG, widths):
        rgb = tuple(int(c[i:i + 2], 16) / 255 for i in (1, 3, 5))
        ym = y - 0.058
        page.draw_line(fitz.Point(x * CM, ym * CM), fitz.Point((x + 0.35) * CM, ym * CM),
                       color=rgb, width=0.9)
        K.put_text(page, fitz.Point((x + 0.45) * CM, y * CM), t, size=fs)
        x += 0.45 + w + gap
    return fs


def letter(page, L, x, y, size=10.0):
    K.put_text(page, fitz.Point(x * CM, (y + size * 0.82 / 72 * 2.54) * CM), L,
               size=size, bold=FL.LTR_BOLD, fname='F1v4Ltr' + L)


def main():
    for f in sorted(SRC.values()):
        if not os.path.exists(f):
            raise SystemExit('缺源图: %s（先跑 replot_v4.py / fig1c_v4.py）' % f)
    doc = fitz.open()
    page = doc.new_page(width=PW * CM, height=FL.PH * CM)
    sc = {}
    sc['A'] = put(page, SRC['A'], 0.0, FL.A_TOP, PW, FL.A_H)
    fs_leg = draw_legend(page, FL.LEG_Y)
    for L in ('B', 'C', 'D', 'E', 'F'):
        x, y, w, h = SLOTS[L]
        sc[L] = put(page, SRC[L], x, y, w, h)
    for L, (x, y) in LETTERS.items():
        if L == 'A':
            continue                                   # A 字母在裁贴带内
        letter(page, L, x, y)
    print('   fit scale:', {k: round(v, 3) for k, v in sc.items()}, 'legend fs=%.1f' % fs_leg)
    out = os.path.join(FIG, 'Fig1_v4.pdf')
    doc.save(out, deflate=True, garbage=3)
    d = fitz.open(out)
    os.makedirs(PNG, exist_ok=True)
    d[0].get_pixmap(dpi=300).save(os.path.join(PNG, 'Fig1_v4.png'))
    print('OK %-14s %.2f x %.2f cm  %.0f KB' % ('Fig1_v4', d[0].rect.width / CM,
          d[0].rect.height / CM, os.path.getsize(out) / 1024))
    d.close()


if __name__ == '__main__':
    main()
