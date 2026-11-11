#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""patch_figS3_v7.py —— 在**已放大的交付件** Figure_S3.pdf 上做四行文字外科替换（Codex #87 ①②④）
=========================================================================================
为什么以 `FINAL_SET_20261021/Figure_S3.pdf` 为基准（作者 2026-10-29 定）：
  · 它的 UMAP 点已经是 0.45 pt（点径回退的教训见 FIGURES_CODEX_ROUND2_20261029.md §11/§12）；
  · 作者认可它的观感：图面没有多余文字，v6 那种"加通栏定义带"的版式作废，长解释进图注。
  ⇒ 本轮**只换四行字**，其余一个像素都不许动。

替换表（新句宽度用 figkit.text_width 实测，全部落在原标签右侧的空里）：
  1 F 标题  COP/NFOL resolvability tracks depth, not chemistry
          → COP/NFOL resolvability improves with the number of lineage nuclei analysed   (6.2 pt)
  2         never separable (0/20)  → not separable (0/20 draws)                          (5.5 pt)
  3         always separable (1/1)  → separable (1/1 libraries)                           (5.5 pt)
  4 G       fails C1ql1             → C1ql1 not met                                      (5.5 pt)
     ④ 的措词：批注说"15-cell 候选按 C1ql1 判失败"像事后排除。这里**不删候选、不删基因**，
        只把 fails 换成中性的 not met，图注写明规则是预先冻结的，
        另附敏感性表 scripts/s3c_c1ql1_sensitivity.py（拿掉 C1ql1 后 24 个运行只有 1 个判定改变）。

⚠ 三条本轮实测出来的机制坑（都写进代码，别再试错）：
  (a) `page.show_pdf_page(..., overlay=True)` **不能**用在"已经贴过整页的页"上：
      本页有 37 条内容流 + 嵌套 Form，第二次贴会把已贴的画**整个抹掉**
      （实测 A/B UMAP 墨量 7.243% → 0.000%）。⇒ 新字不能用 matplotlib 叠加层。
  (b) `apply_redactions` 本身是安全的，但必须 `images=NONE, graphics=NONE, text=REMOVE`，
      默认 graphics=1 会连矢量一起削。
  (c) 直接 `insert_text(fontfile=ARIAL.TTF)` 会把 1.01 MB 整份字库嵌进去（"未子集化"硬门 FAIL）。
      ⇒ 先用 fontTools 把用到的 35 个字形子集化（39 KB，name1 仍是 Arial），再嵌这份子集。

用法: PYTHONPATH=<work>/pylibs python3 scripts/patch_figS3_v7.py
验收: 逐像素差异只允许落在四个标签框内；文字层新句在、旧句不在；final_audit 硬门全过；
      UMAP 墨量与母版逐区一致（点径不许动）。
"""
import os, sys, subprocess
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fitz
from figkit import CM, text_width

SRC = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'FINAL_SET_20261021/Figure_S3.pdf')
DST = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'FigS3_v7.pdf')
SUB = '/tmp/ArialLabel_s3v7.ttf'
INK = (0.13, 0.13, 0.13)

# (旧句, 新句, 新句基线 x/y cm, 字号, redact 紧框)
EDIT = [
    ('COP/NFOL resolvability tracks depth, not chemistry',
     'COP/NFOL resolvability improves with the number of lineage nuclei analysed',
     9.77, 20.65, 6.2, (9.77, 20.30, 14.80, 20.72)),
    ('never separable (0/20)', 'not separable (0/20 draws)',
     10.05, 21.55, 5.5, (10.05, 21.38, 12.02, 21.62)),
    ('always separable (1/1)', 'separable (1/1 libraries)',
     15.18, 21.55, 5.5, (15.18, 21.38, 17.15, 21.62)),
    ('fails C1ql1', 'C1ql1 not met',
     13.89, 31.30, 5.5, (13.89, 31.13, 14.82, 31.36)),
]

src = fitz.open(SRC)
pg = src[0]
W_pt, H_pt = pg.rect.width, pg.rect.height
W_cm, H_cm = W_pt / CM, H_pt / CM

# ---------------------------------------------------------------- 1. 定位 + 宽度体检
for old, new, x, y, fs, box in EDIT:
    hits = [s for b in pg.get_text('dict')['blocks'] for l in b.get('lines', [])
            for s in l['spans'] if s['text'] == old]
    if len(hits) != 1:
        raise SystemExit('定位失败（%d 处）: %s' % (len(hits), old))
    s = hits[0]
    r = fitz.Rect(s['bbox']) / CM
    if max(abs(r.x0 - box[0]), abs(r.y1 - box[3])) > 0.06:
        raise SystemExit('记录的框漂了: %s -> %.2f,%.2f,%.2f,%.2f' % (old, *box))
    w_new = text_width(new, fs)
    if x + w_new > W_cm - 0.25:
        raise SystemExit('新句放不下（%.2f cm 越页）: %s' % (x + w_new, new))
    print('%-62s → %s' % (old[:60], new[:64]))
    print('    旧宽 %.2f  新宽 %.2f  右沿 %.2f / 页宽 %.2f cm' % (text_width(old, s['size']), w_new, x + w_new, W_cm))

# ---------------------------------------------------------------- 2. 新页 ← 贴母版（空页才安全，见头注 a）
out = fitz.open()
q = out.new_page(width=W_pt, height=H_pt)
q.show_pdf_page(q.rect, src, 0)

# ---------------------------------------------------------------- 3. 只删文字，保留全部矢量/图像（头注 b）
for old, new, x, y, fs, box in EDIT:
    q.add_redact_annot(fitz.Rect(*[v * CM for v in box]))
q.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                   graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                   text=fitz.PDF_REDACT_TEXT_REMOVE)
t = q.get_text()
if any(old in t for old, *_ in EDIT):
    raise SystemExit('旧字没删干净')

# ---------------------------------------------------------------- 4. 子集 Arial 后直接插字（头注 c）
chars = ''.join(sorted(set(''.join(e[1] for e in EDIT))))
from fontTools.ttLib import TTFont
from fontTools.subset import Subsetter, Options
f = TTFont(os.path.join(os.environ.get('DSCI_FONTS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'fonts')), 'ARIAL.TTF'))
o = Options(); o.name_IDs = ['*']; o.notdef_outline = True; o.drop_tables = ['DSIG']
s = Subsetter(o); s.populate(text=chars); s.subset(f); f.save(SUB)
print('子集字体 %.1f KB（原 %.1f KB）/ 用到 %d 个字形'
      % (os.path.getsize(SUB) / 1024, os.path.getsize(os.path.join(os.environ.get('DSCI_FONTS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'fonts')), 'ARIAL.TTF')) / 1024, len(chars)))
for new, x, y, fs in [(e[1], e[2], e[3], e[4]) for e in EDIT]:
    q.insert_text(fitz.Point(x * CM, y * CM), new, fontsize=fs, fontfile=SUB,
                  fontname='FigLbl', color=INK)

out.save(DST, garbage=4, deflate=True)
out.close(); src.close()
print('OK ->', DST, '%.2f MB' % (os.path.getsize(DST) / 1e6))
