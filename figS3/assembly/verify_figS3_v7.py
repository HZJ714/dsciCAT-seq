import os
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify_figS3_v7.py —— FigS3_v7.pdf 的"只改了四行字"验收
=========================================================================================
四道门（任何一道不过就退出非零）：
  1. 逐像素：300 dpi 下所有差异像素必须落在四个标签框（含 0.10 cm 余量）里；
  2. 墨量：分区块（A/B UMAP、C–E 热图、F、G）与母版逐块一致（点径/热图一个点都不许变）；
  3. 文字层：四句新文案在、四句旧文案不在；
  4. 字体硬门：交给 final_audit.py（Arial / 嵌入 / 子集 / ≥5.5 pt / Type3 / 隐藏字）。
"""
import sys, os
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
import fitz, numpy as np

SRC = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'FINAL_SET_20261021/Figure_S3.pdf')
NEW = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'FigS3_v7.pdf')
CM = 72 / 2.54
DPI = 300
# 允许改动区 = 旧标签框 ∪ 新句子实际落墨框（新句普遍比旧句长，只按旧框判会把自家新字当成"越界"，
# 上一版就是这么误判成"redact 把隐藏文字翻出来了"）
BOXES = [(9.70, 20.26, 17.30, 20.78), (9.98, 21.32, 12.42, 21.70),
         (15.10, 21.32, 17.28, 21.70), (13.82, 31.07, 15.16, 31.42)]
PAD = 0.05
REGIONS = {'A/B UMAP 左': (0.5, 5.5, 7.0, 9.5), 'A/B UMAP 右': (9.0, 5.5, 17.5, 9.5),
           'C-E 热图': (0.5, 12.0, 17.5, 20.0), 'F 稳定性': (9.5, 20.0, 18.0, 25.5),
           'G marker 点阵': (0.5, 28.0, 18.0, 34.0)}
OLD = ['COP/NFOL resolvability tracks depth, not chemistry', 'never separable (0/20)',
       'always separable (1/1)', 'fails C1ql1']
NEW_TXT = ['COP/NFOL resolvability improves with the number of lineage nuclei analysed',
           'not separable (0/20 draws)', 'separable (1/1 libraries)', 'C1ql1 not met']
fail = []


def arr(f):
    d = fitz.open(f); p = d[0]
    pix = p.get_pixmap(dpi=DPI)
    a = np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, pix.n)
    return a[:, :, :3].astype(np.int16)


def ink(d, clip, mask=True):
    """区块墨量。**必须把四个标签框挖掉再比**：新标题比旧标题长 2.5 cm，
    整块墨量当然会变（53.250%→53.698%），拿它当"图没被动过"的证据是错的；
    挖掉标签区之后剩下的全是热图/柱子/点，那才是一个都不许变的。"""
    p = d[0]
    pix = p.get_pixmap(dpi=200, clip=fitz.Rect(*[c * CM for c in clip]))
    a = np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, pix.n)[:, :, :3]
    m = a.mean(2) < 250
    if mask:
        sc = 200 / 2.54
        for x0, y0, x1, y1 in BOXES:
            m[max(int((y0 - clip[1]) * sc), 0):int((y1 - clip[1]) * sc),
              max(int((x0 - clip[0]) * sc), 0):int((x1 - clip[0]) * sc)] = False
    return 100 * float(m.mean())


A, B = arr(SRC), arr(NEW)
if A.shape != B.shape:
    print('!! 像素尺寸不同', A.shape, B.shape); sys.exit(2)
D = np.abs(A - B).sum(2) > 24
SC = DPI / 2.54
ys, xs = np.nonzero(D)
inside = np.zeros(len(ys), bool)
for x0, y0, x1, y1 in BOXES:
    inside |= (xs / SC >= x0 - PAD) & (xs / SC <= x1 + PAD) & (ys / SC >= y0 - PAD) & (ys / SC <= y1 + PAD)
print('1) 差异像素 %d 个（%.4f%% 全页），其中落在四个标签框内 %d / 框外 %d'
      % (D.sum(), 100 * D.mean(), inside.sum(), (~inside).sum()))
if (~inside).sum() > 0:
    u = np.unique(np.c_[ys[~inside] // 60, xs[~inside] // 60], axis=0)[:6]
    print('   !! 框外改动位置(cm, 抽样):', [(round(c[1] * 2.54 / 20, 2), round(c[0] * 2.54 / 20, 2)) for c in u])
    fail.append('像素差异越界')

ds, dn = fitz.open(SRC), fitz.open(NEW)
print('2) 分区块墨量（点径/热图不许变）')
for nm, clip in REGIONS.items():
    a, b = ink(ds, clip), ink(dn, clip)
    ok = abs(a - b) <= 0.02
    print('   %-12s 母版 %7.3f%%  新版 %7.3f%%  （已挖掉四个标签区）%s' % (nm, a, b, '✓' if ok else '✗'))
    if not ok:
        fail.append('墨量变了: ' + nm)

t = dn[0].get_text()
print('3) 文字层')
for o, n in zip(OLD, NEW_TXT):
    a, b = o in t, n in t
    print('   新句 %-8s / 旧句残留 %-8s  %s' % (b, a, '✓' if (b and not a) else '✗'))
    if not (b and not a):
        fail.append('文字层: ' + o)

print('4) 大小 %.2f MB（母版 %.2f MB）' % (os.path.getsize(NEW) / 1e6, os.path.getsize(SRC) / 1e6))
print('\n==>', '全部通过 ✓' if not fail else '不通过: ' + '; '.join(fail))
sys.exit(1 if fail else 0)
