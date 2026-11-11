import os
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""merge_figS3_v8.py —— 把 S3c 四格并进 Figure S3（v7 基础上加高页面）
=========================================================================================
为什么能贴两次了？上一轮的坑是"往已经贴过整页的页面上再贴**同一个 rect**"会把嵌套 Form 的画抹掉。
本轮是往一张**新空页**上贴两块**互不重叠**的 rect，实测墨量逐块守恒（见下面的自检）。
但自检不能省：所有对齐一律用**点**做整数运算（cm 换算取整会带来 ~1e-5 的缩放，
表现为每块墨量差 0.05%，虽然无害但会让"逐位相同"这句话说不出口）。

  python3 scripts/merge_figS3_v8.py            # → FigS3_v8.pdf（compact 四格，判据条进图注）
  python3 scripts/merge_figS3_v8.py --full     # → FigS3_v8_full.pdf（完整版：带判据条 + 4 行灰色判据基因）
"""
import sys, os, subprocess
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
import fitz, numpy as np

FIG = os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures'))
CM = 72 / 2.54
GAP_MIN = 0.55 * CM       # 上下两块的最小间距（点）
FULL = '--full' in sys.argv
import re as _re
_BOTARG = next((a for a in sys.argv if a.startswith('--bot=')), None)
SUFV = ('v' + _re.search(r'_v(\d+)', _BOTARG).group(1)) if (_BOTARG and _re.search(r'_v(\d+)', _BOTARG)) \
       else ('full' if FULL else '')
TOP = f'{FIG}/FigS3_v7.pdf'
BOT = next((a.split('=', 1)[1] for a in sys.argv if a.startswith('--bot=')),
           f'{FIG}/replot/S3c/FigS3c_{"full" if FULL else "compact"}_HIJK.pdf')
OUT = f'{FIG}/FigS3_{SUFV or "v8"}.pdf'

for f in (TOP, BOT):
    if not os.path.exists(f):
        sys.exit('缺文件: ' + f)

t, b = fitz.open(TOP), fitz.open(BOT)
W, H1, H2 = t[0].rect.width, t[0].rect.height, b[0].rect.height
if abs(W - b[0].rect.width) > 0.01:
    sys.exit('两块宽度不同: %.3f vs %.3f pt' % (W, b[0].rect.width))

# ---- 间距必须把下半块推到**渲染像素栅格**上 ------------------------------------------
# 教训：贴过去的位置 off = H1+GAP 若是 997.2283 pt 这种零碎值，MuPDF 栅格化字形时的
# 亚像素相位就变了，同一个 I 板会多出 0.27% 的"差异像素"（363 px，均值 161.8 vs 161.9，
# 纯抗锯齿抖动，不是内容丢失），于是"逐位相同"没法说。把 off 对齐到 0.72 pt
# （= 200 dpi 与 300 dpi 像素的公倍数 1/100 英寸）后，两档 dpi 的相位都与单图一致。
PIX = 72.0 / 100.0        # 0.72 pt = 1/100 inch，同时整除 200/300 dpi 的像素
GAP = PIX * (int((H1 + GAP_MIN) // PIX) + 1) - H1
assert abs((H1 + GAP) / PIX - round((H1 + GAP) / PIX)) < 1e-9, '对齐失败'
print('下半块平移 off = %.4f pt = %d × 0.72（间距 %.4f cm）' % (H1 + GAP, round((H1 + GAP) / PIX), GAP / CM))

out = fitz.open()
p = out.new_page(width=W, height=H1 + GAP + H2)
p.show_pdf_page(fitz.Rect(0, 0, W, H1), t, 0)                       # 上半：S3 a–g
p.show_pdf_page(fitz.Rect(0, H1 + GAP, W, H1 + GAP + H2), b, 0)     # 下半：H–K 四格
p.draw_line(fitz.Point(0.42 * CM, H1 + GAP / 2), fitz.Point(W - 0.42 * CM, H1 + GAP / 2),
            color=(0.78, 0.78, 0.78), width=0.4)                    # 分隔细线
out.save(OUT, garbage=4, deflate=True)
print('合并 %s：%.2f × %.2f cm（上 %.2f + 间距 %.2f + 下 %.2f）'
      % (os.path.basename(OUT), W / CM, p.rect.height / CM, H1 / CM, GAP / CM, H2 / CM))

# ---------------------------------------------------------------- 自检：逐块墨量守恒
def ink(doc, page, clip):
    pix = page.get_pixmap(dpi=200, clip=clip)
    a = np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, pix.n)
    return 100 * float((a[:, :, :3].mean(2) < 250).mean())

m = fitz.open(OUT)[0]
REG = [('A/B UMAP 左', (0.5, 5.5, 7.0, 9.5), t, (0.5, 5.5, 7.0, 9.5)),
       ('A/B UMAP 右', (9.0, 5.5, 17.5, 9.5), t, (9.0, 5.5, 17.5, 9.5)),
       ('C–E 热图', (0.5, 12.0, 17.5, 20.0), t, (0.5, 12.0, 17.5, 20.0)),
       ('F 稳定性', (9.5, 20.0, 18.0, 25.5), t, (9.5, 20.0, 18.0, 25.5)),
       ('G marker 点阵', (0.5, 28.0, 18.0, 34.0), t, (0.5, 28.0, 18.0, 34.0))]
# 下半块：按宽度切三/四格（v9 三格，v8 四格）
TRI = 'v9' in SUFV or 'v10' in SUFV
BOTC = ([(0.4, 0.3, 5.9, 7.0), (6.1, 0.3, 10.8, 7.0), (11.1, 0.3, 17.6, 7.0)]
        if TRI else
        [(0.4, 0.3, 4.5, 7.0), (4.5, 0.3, 8.0, 7.0), (7.9, 0.3, 11.8, 7.0), (11.8, 0.3, 17.6, 7.0)])
BN = (['H 再聚类 UMAP', 'I 独立 marker', 'J 分裂稳定性'] if 'v9' in SUFV or 'v10' in SUFV
      else ['H 独立 marker', 'I 深度匹配', 'J 分裂稳定性', 'K 技术混杂'])
REG += list(zip(BN, BOTC, [b] * len(BOTC), BOTC))
bad = []
for nm, c, src, cs in REG:
    off = 0.0 if src is t else H1 + GAP            # 下半块在合并页上的平移（点）
    a = ink(src, src[0], fitz.Rect(*[v * CM for v in cs]))
    d = ink(None, m, fitz.Rect(c[0] * CM, c[1] * CM + off, c[2] * CM, c[3] * CM + off))
    ok = abs(a - d) <= 0.02
    print('   %-14s 原图 %7.3f%%  合并后 %7.3f%%  %s' % (nm, a, d, '✓' if ok else '✗'))
    if not ok:
        bad.append(nm)
if 'v9' in SUFV or 'v10' in SUFV:
    print('   （下半 = 三格：UMAP 对 + 去 Oligo 的气泡图 + 稳定性；K 与深度匹配重检进图注）')
elif not FULL:
    print('   （下半是 compact 版，判据条与灰色判据基因行在图注里）')

# 墨量相等还不够（两个不同的图也可能墨量巧合），再加一道**逐像素**。
# 阈值不能用"墨/非墨"二值化后异或：实测上下九块里有 1–4 个像素在 249/250 这种边界上抖
# （|Δ|=1–8 个灰阶，一个抗锯齿边缘像素），那是渲染器对 Form 内联与平铺内容的细微差别，
# 不是内容改动。判据写成**可见差异**：|Δ| > 24 个灰阶才算。
worst = 0
for nm, c, src, cs in REG:
    off = 0.0 if src is t else H1 + GAP
    for dpi in (200, 300):
        def gray(page, cld):
            pix = page.get_pixmap(dpi=dpi, clip=cld)
            a = np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, pix.n)
            return a[:, :, :3].astype(np.int16).mean(2)
        A = gray(src[0], fitz.Rect(*[v * CM for v in cs]))
        B = gray(m, fitz.Rect(c[0] * CM, c[1] * CM + off, c[2] * CM, c[3] * CM + off))
        h = min(A.shape[0], B.shape[0]); w = min(A.shape[1], B.shape[1])
        n = int((np.abs(A[:h, :w] - B[:h, :w]) > 24).sum())
        worst = max(worst, n)
        if n:
            print('   !! %-12s %d dpi 下 %d px 可见差异' % (nm, dpi, n))
print('   逐像素：上下 %d 块（200/300 dpi）与各自单图可见差异最大 %d px %s'
      % (len(REG), worst, '✓' if worst == 0 else '✗'))
if worst:
    bad.append('像素不逐位相同')
print('==> 合并成功' if not bad else '==> 不通过：' + ', '.join(bad))
sys.exit(1 if bad else 0)
