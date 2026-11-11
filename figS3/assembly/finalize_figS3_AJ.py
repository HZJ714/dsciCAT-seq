#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""finalize_figS3_AJ.py —— 把 FigS3_v10.pdf（A–J）落进交付块，并同步 S3b / PNG / md5
=========================================================================================
为什么需要这个脚本：交付块 `FINAL_SET_20261021/` 里 S3 有**三种形态**并存 ——
  Figure_S3.pdf   18.00 × 34.56 cm  一页 A–G（**投稿主件**）
  Figure_S3a.pdf  18.00 × 20.18 cm  A–D        ┐ 编辑部若强制"一图一页"时的
  Figure_S3b.pdf  18.00 × 14.84 cm  E–G        ┘ 备用两件（切缝在 D/E 之间，重叠 0.46 cm）
只换主件会让备用件**悄悄丢掉 H–J**，所以三件必须一起动。

做的事（顺序即检查顺序，任何一步不过就退出码 1、不落盘）：
  1. 备份被替换的两件到 `_superseded/`（带面板范围与高度，块内既有命名风格）。
  2. 主件 ← FigS3_v10.pdf。
  3. S3b ← 新主件的 E–J 区（**从主件裁**，不是从旧 S3b 续拼，见下面"旧 S3b 的坑"）。
  4. 逐像素自检四件事（200/300 dpi，可见差 = |Δ| > 24 灰阶）：
       a 新主件上半（原 34.56 cm 区） == FigS3_v7.pdf          —— 保证只加了 H–J，A–G 没被动
       b 新主件 E–J 区               == 新 S3b                 —— 两件形态一致
       c 新 S3b 的 H–J 段            == 单图 FigS3c_v10_HIJK.pdf
       d S3a 的 md5 不变                                        —— A–D 一个像素都没碰
  5. 300 dpi PNG 重出两件（块内约定：PIL 写 dpi=(300,300)，读回 299.9994）。
  6. 重写 MD5_checksums.txt（14 件，格式与顺序照旧）。

**旧 S3b 的坑（本轮查出来，顺手修掉）**：交付块里 `Figure_S3b.pdf` 的 E–G 与主件的 E–G
**并不逐像素相同** —— 对齐到正确切缝后仍有 5.55% 可见差，差异全在文字（热图方格对得上、
字有 2–3 px 垂直错位）。来源是 2026-10-21 那次"字母→常规 + 标签对齐"只重出了主件，
S3b 是从旧块照抄的。所以新 S3b 一律**从主件裁**，两件从此保证同源。
"""
import os, sys, shutil, hashlib, zipfile
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
import fitz
import numpy as np
from PIL import Image

FIG = os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures'))
BLK = f'{FIG}/FINAL_SET_20261021'
SUP = f'{BLK}/_superseded'
SRC = f'{FIG}/FigS3_v10.pdf'
S3C = f'{FIG}/replot/S3c/FigS3c_v10_HIJK.pdf'
V7 = f'{FIG}/FigS3_v7.pdf'
DPI = 300
PT = 72 / 2.54
# 旧主件里 E 的顶沿 = 旧主件高 - 旧 S3b 高 = 559.1055 pt。但**裁切平移必须落在渲染像素格上**：
# 559.1055 不是 0.72 pt（200 与 300 dpi 公倍）的整数倍，直接把 E–J 裁出来会让 MuPDF 的字形
# 亚像素相位改变，自检 b/c 立刻报 3761 / 22866 px "差异"（内容没错，全是抗锯齿抖动）。
# 取最近的格点 k=777 → 559.44 pt，切缝只下移 0.012 cm（S3a 重叠 0.46→0.45 cm），换来逐位相同。
CUT_PT = 777 * 0.72            # = 559.44 pt = 19.7360 cm
bad = []


def md5(p):
    return hashlib.md5(open(p, 'rb').read()).hexdigest()


def gray(page, dpi):
    pm = page.get_pixmap(dpi=dpi, alpha=False)
    return np.frombuffer(pm.samples, np.uint8).reshape(pm.height, pm.width, pm.n)[:, :, :3].mean(2)


def visdiff(pa, pb, dpi=200, off_a=None, off_b=None, size=None):
    """两块区域逐像素比；off_* = (x0, y0) pt 起点，size = (w, h) pt。返回可见差像素数。"""
    A = gray(pa, dpi) if off_a is None else gray_clip(pa, off_a, size, dpi)
    B = gray(pb, dpi) if off_b is None else gray_clip(pb, off_b, size, dpi)
    h = min(A.shape[0], B.shape[0]); w = min(A.shape[1], B.shape[1])
    return int((np.abs(A[:h, :w].astype(int) - B[:h, :w].astype(int)) > 24).sum())


def gray_clip(page, off, size, dpi):
    r = fitz.Rect(off[0], off[1], off[0] + size[0], off[1] + size[1])
    pm = page.get_pixmap(dpi=dpi, clip=r, alpha=False)
    return np.frombuffer(pm.samples, np.uint8).reshape(pm.height, pm.width, pm.n)[:, :, :3].mean(2)


def png_from(pdf, name):
    pm = fitz.open(pdf)[0].get_pixmap(dpi=DPI, alpha=False)
    im = Image.frombytes('RGB', (pm.width, pm.height), pm.samples)
    p = os.path.join(BLK, 'png', name + '.png')
    im.save(p, dpi=(DPI, DPI))
    return p, im.size


# ---------------------------------------------------------------- 0. 前置断言
m_new = fitz.open(SRC)[0]
assert abs(m_new.rect.height / PT - 42.8782) < 0.002, m_new.rect.height
old_s3 = f'{BLK}/Figure_S3.pdf'
s3a_md5_before = md5(f'{BLK}/Figure_S3a.pdf')

# ---------------------------------------------------------------- 1. 备份
os.makedirs(SUP, exist_ok=True)
for f, tag in [(old_s3, 'Figure_S3_A-G_34.56cm_deprecated.pdf'),
               (f'{BLK}/Figure_S3b.pdf', 'Figure_S3b_E-G_14.84cm_deprecated.pdf')]:
    d = os.path.join(SUP, tag)
    if not os.path.exists(d):
        shutil.copy2(f, d)
        print('备份 %s -> _superseded/%s (%.2f MB, md5 %s)'
              % (os.path.basename(f), tag, os.path.getsize(d) / 1e6, md5(d)[:12]))
    else:
        print('已存在备份 _superseded/%s（不动）' % tag)

# 重跑安全：几何自检一律对着 _superseded 里的**旧件**，别管块里现在放的是第几版
BAK_S3 = os.path.join(SUP, 'Figure_S3_A-G_34.56cm_deprecated.pdf')
BAK_B = os.path.join(SUP, 'Figure_S3b_E-G_14.84cm_deprecated.pdf')
old_h = fitz.open(BAK_S3)[0].rect.height
assert abs(old_h - 979.6535) < 0.01, old_h
assert abs(fitz.open(BAK_B)[0].rect.height - 420.5480) < 0.01
assert abs(old_h - fitz.open(BAK_B)[0].rect.height - 559.1055) < 0.01   # 旧切缝

# ---------------------------------------------------------------- 2/3. 落盘
shutil.copy2(SRC, old_s3)
print('主件 Figure_S3.pdf <- FigS3_v10.pdf  %.2f × %.2f cm  %.2f MB'
      % (m_new.rect.width / PT, m_new.rect.height / PT, os.path.getsize(old_s3) / 1e6))

H_B = m_new.rect.height - CUT_PT                      # 新 S3b 高（pt）
d = fitz.open()
pg = d.new_page(width=m_new.rect.width, height=H_B)
pg.show_pdf_page(fitz.Rect(0, 0, m_new.rect.width, H_B), fitz.open(old_s3), 0,
                 clip=fitz.Rect(0, CUT_PT, m_new.rect.width, m_new.rect.height))
d.save(f'{BLK}/Figure_S3b.pdf', garbage=4, deflate=True)
d.close()
print('备用件 Figure_S3b.pdf = 主件 E–J 裁切  %.2f × %.2f cm' % (m_new.rect.width / PT, H_B / PT))

# ---------------------------------------------------------------- 4. 自检
mm = fitz.open(old_s3)[0]
v7 = fitz.open(V7)[0]
s3b = fitz.open(f'{BLK}/Figure_S3b.pdf')[0]
s3c = fitz.open(S3C)[0]
n = visdiff(mm, v7, 200, (0, 0), (0, 0), (510.2362, 979.6535))
print('  a 主件上半(原 34.56cm 区) vs FigS3_v7.pdf   可见差 %d px %s' % (n, '✓' if n == 0 else '✗'))
bad += [] if n == 0 else ['主件上半与 v7 不一致']
n = visdiff(mm, s3b, 200, (0, CUT_PT), (0, 0), (510.2362, H_B))
print('  b 主件 E–J 区 vs 新 S3b                    可见差 %d px %s' % (n, '✓' if n == 0 else '✗'))
bad += [] if n == 0 else ['S3b 与主件 E–J 区不一致']
n = visdiff(s3b, s3c, 200, (0, H_B - s3c.rect.height), (0, 0), (510.2362, s3c.rect.height))
print('  c 新 S3b 的 H–J 段 vs 单图 FigS3c_v10      可见差 %d px %s' % (n, '✓' if n == 0 else '✗'))
bad += [] if n == 0 else ['S3b 的 H–J 段与单图不一致']
same = md5(f'{BLK}/Figure_S3a.pdf') == s3a_md5_before
print('  d Figure_S3a.pdf md5 不变                   %s %s' % (s3a_md5_before[:12], '✓' if same else '✗'))
bad += [] if same else ['S3a 被改动了！']
# 旧 S3b 与旧主件的偏差（记录用，不参与判定）
oldb = fitz.open(BAK_B)[0]
n = visdiff(fitz.open(BAK_S3)[0], oldb, 200,
            (0, 979.6535 - oldb.rect.height), (0, 0), (510.2362, oldb.rect.height))
print('  （存档对照：旧 S3b vs 旧主件同位 = %d px 可见差 —— 这就是本轮顺手修掉的旧偏差）' % n)

# ---------------------------------------------------------------- 5. PNG
for pdf, stem in [(old_s3, 'Figure_S3'), (f'{BLK}/Figure_S3b.pdf', 'Figure_S3b')]:
    p, sz = png_from(pdf, stem)
    print('PNG %-14s %d × %d px = %.2f × %.2f cm' % (stem, sz[0], sz[1], sz[0] / 300 * 2.54, sz[1] / 300 * 2.54))

# ---------------------------------------------------------------- 6. md5 清单
# 行序照旧清单（不然 diff 会满屏重排）；新出现的件追加在末尾
old_list = [l.split()[-1][:-4] for l in open(f'{BLK}/MD5_checksums.txt') if l.strip()]
stems = [f[:-4] for f in sorted(os.listdir(BLK)) if f.endswith('.pdf')]
order = [s for s in old_list if s in stems] + [s for s in stems if s not in old_list]
assert len(order) == len(stems), (order, stems)   # 15 件：14 图 + Graphical_Abstract
with open(f'{BLK}/MD5_checksums.txt', 'w') as f:
    for s in order:
        f.write('%s  %s.pdf\n' % (md5(os.path.join(BLK, s + '.pdf')), s))
print('MD5_checksums.txt 重写（%d 件，行序照旧）' % len(order))

print('==> %s' % ('全部通过' if not bad else '未通过: ' + '; '.join(bad)))
sys.exit(1 if bad else 0)
