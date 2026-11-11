#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
patch_figS5_c_units.py — Figure S5 的 C 面板：撤掉错误的 "RPM-normalized" + 对齐两个 MWU 口径
============================================================================================
Codex #92 前半（`CODEX_COMMENTS_EXTRACT_20261028.txt` 第 279 行）：
  「C 轨道宣称 RPM-normalized 但 y 轴达 20 万/50 万且各行范围不同，需核实是否实际用 raw counts
    或累计权重；用统一合理量纲、mm10 坐标和对数/固定尺度……图内 Nrxn3 p<10⁻²² 与正文 6.2×10⁻⁸⁰
    需对齐比较。」

**两条都实测确认，属"图面写了错话"，不是读图口径问题：**

1) 量纲。上游 `pi/scripts-scenic+/grn_mef2c_nrxn3_track.py`：
       rpkm = 每底 counts × 1e6 / (该态 chr12 总插入)        # 每底、按各态文库量归一
       y    = 200 底求和 × (1e6 / 200)                        # 画图时又乘了一次 1e6/BIN
   ⇒ y ∝ 每底平均插入数 × 1e12/total_chr —— **不是 RPM**（RPM 量级是个位数到几百）。
   实测缓存 `scripts/cache/nrxn3_tracks.npz`：六态 y 最大 477,918 / 622,101 / 300,806 /
   308,925 / 312,953 / 293,123；交付版刻度 500000 / 500000 / 200000 / 200000 / 200000 / 200000
   ⇒ 每行既按各自文库归一、又各自 autoscale。标题句改成
   "ATAC insertions per 200-bp bin, library-size normalized per row"（新串仅宽 +0.26 cm，
   右端 11.27 cm，离右侧表达图 14.45 cm 还差 3.2 cm）。
2) p 值。用 `scripts/cache/nrxn3_expr.npz` 重算单侧 Mann-Whitney U：
       MOL3 vs 其余五态**合并** (n=688 vs 5,289)  p = 6.15e-80   ← 正文那个数
       逐对：vs COP_NFOL 1.1e-22 / MFOL 1.5e-42 / MOL2 4.82e-55 / MOL3_6 1.67e-51 / MOL6 1.12e-70
   ⇒ 图上 "MWU p<1e-22" 是**最弱逐对**（且是最小的那组 n=203），正文是**合并**口径：
   两个都没算错，但图上没说是哪个比较。改法：图内两处 "p<1e-22" → "p<1e-79"
   （合并口径的保守写法，6.15e-80 < 1e-79 ✓，**字符数完全相同** ⇒ 绿色标注框不必重画）；
   精确值与两个口径同时写进图注。

⚠ 为什么是"成品 PDF 上原位换字"而不是改脚本重渲：实测**盘上的 `figs5_panel_track_box.py`
   已经不能复现交付版 C 面板**（一字未改地重渲，与交付版差 45,311 px = 3.62%，最佳对齐还要
   整体错 1 px）——10-12 之后脚本又动过。重渲会顺手改掉 3.6% 的画面，而本次要改的只有 3 行字。
⚠ 为什么不用 `pdfix.rewrite`：它（以及任何走 PyMuPDF TextWriter 的写法）会把词间写成 U+00A0、
   连字符写成 U+00AD（交付件这两类字符计数都是 0）。改用 Figure S3 v7 那套已验证配方：
   **贴页 → 逐框 redact(images=NONE, graphics=NONE, text=REMOVE) → 用 fontTools 现裁的 Arial
   子集 insert_text**。自检 = 像素差只落在三行字的足迹并集内 + 文字层无 \xa0/\xad + 门全过。

用法: PYTHONPATH=<work>/pylibs python3 scripts/patch_figS5_c_units.py [--apply]
输入  FigS5_noc4.pdf     输出  FigS5_noc4_cpatched.pdf（--apply 时原地替换，自动备份 /tmp）
"""
import os, sys, shutil
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fitz
import numpy as np
from figkit import CM

FIG = os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures'))
SRC = f'{FIG}/FigS5_noc4.pdf'
DST = f'{FIG}/FigS5_noc4_cpatched.pdf'
SUB = f'{FIG}/scripts/cache/Arial-C-units-patch.ttf'
TOL = 0.06                      # 记录框漂移容差 (cm)

L1_OLD = ('Nrxn3 locus (chr12:88,550,000-89,600,000)  |  ATAC insertions per cell state '
          '(RPM-normalized, 200bp bins)')
L1_NEW = ('Nrxn3 locus (chr12:88,550,000-89,600,000)  |  ATAC insertions per 200-bp bin, '
          'library-size normalized per row')
L3_OLD = ('Right: Nrxn3 expression (boxplot + jitter, linearized 2ˣ−1) per state — MOL3-specific '
          '(median ≈7.0 vs ≈3.1–3.8, ≈2×; MWU p<1e-22)')
L3_NEW = L3_OLD.replace('p<1e-22', 'p<1e-79')
BX_OLD = '≈2×; MWU p<1e-22'
BX_NEW = '≈2×; MWU p<1e-79'
EDIT = [(L1_OLD, L1_NEW), (L3_OLD, L3_NEW), (BX_OLD, BX_NEW)]


def spans(page):
    out = []
    for b in page.get_text('dict')['blocks']:
        for l in b.get('lines', []):
            for s in l['spans']:
                if s['text'].strip():
                    out.append(s)
    return out


def main(apply=False):
    src = fitz.open(SRC)
    pg = src[0]
    W_pt, H_pt = pg.rect.width, pg.rect.height
    S = spans(pg)
    plan = []
    for old, new in EDIT:
        hit = [s for s in S if s['text'] == old]
        if len(hit) != 1:
            raise SystemExit('定位失败（%d 处）: %r' % (len(hit), old[:52]))
        s = hit[0]
        r = fitz.Rect(s['bbox'])
        col = tuple(((s['color'] >> 16) & 255) / 255, ) if False else (
            ((s['color'] >> 16) & 255) / 255, ((s['color'] >> 8) & 255) / 255, (s['color'] & 255) / 255)
        f = fitz.Font(fontfile=os.path.join(os.environ.get('DSCI_FONTS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'fonts')), 'ARIAL.TTF'))
        w_old = f.text_length(old, fontsize=s['size']) / CM
        w_new = f.text_length(new, fontsize=s['size']) / CM
        assert abs(w_old - sum(x for x in [0]) - w_old) < 1          # 占位，防 lint
        ox, oy = s['origin'][0] / CM, s['origin'][1] / CM
        print('  %r' % old[:46])
        print('    框 x %.2f–%.2f  y %.2f–%.2f  %.1f pt  色 %s  旧宽 %.2f 新宽 %.2f cm  右沿 %.2f'
              % (r.x0 / CM, r.x1 / CM, r.y0 / CM, r.y1 / CM, s['size'],
                 tuple(round(c, 3) for c in col), w_old, w_new, ox + w_new))
        assert ox + w_new < (W_pt / CM) - 0.2, '新句越页'
        # ⚠ redact 框只留 0.004 cm 余量：给到 0.02 就会把**上一行**一起吃掉
        #   （第一次跑就是这样丢了绿框里的 "vs ≈3.1–3.8"，MuPDF 的 redaction 判的是
        #    bbox 相交，不是包含）。少掉的行由下面的"逐行集合必须只少这 3 行"拦住。
        plan.append(dict(old=old, new=new, rect=(r.x0 / CM - 0.004, r.y0 / CM - 0.004,
                                                 r.x1 / CM + 0.004, r.y1 / CM + 0.004),
                         x=ox, y=oy, fs=s['size'], col=col, bold='bold' in s['font'].lower()))
    # ---- 贴页 → 只删文字 → 插新字 ----
    out = fitz.open()
    q = out.new_page(width=W_pt, height=H_pt)
    q.show_pdf_page(q.rect, src, 0)
    for p in plan:
        q.add_redact_annot(fitz.Rect(*[v * CM for v in p['rect']]))
    q.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                       graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                       text=fitz.PDF_REDACT_TEXT_REMOVE)
    t = q.get_text()
    for p in plan:
        assert p['old'] not in t, '旧字没删干净: %r' % p['old'][:40]
    lost = sorted(set(pg.get_text().splitlines()) - set(t.splitlines()))
    want = sorted(p['old'] for p in plan)
    assert lost == want, 'redaction 多吃/少吃了几行:\n  实际少: %s\n  应为: %s' % (lost, want)
    print('  redaction 后逐行集合 = 原集合 - 恰好这 3 行 ✓')
    from fontTools.ttLib import TTFont
    from fontTools.subset import Subsetter, Options
    chars = ''.join(sorted(set(''.join(p['new'] for p in plan))))
    ft = TTFont(os.path.join(os.environ.get('DSCI_FONTS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'fonts')), 'ARIAL.TTF'))
    o = Options(); o.name_IDs = ['*']; o.notdef_outline = True; o.drop_tables = ['DSIG']
    ss = Subsetter(o); ss.populate(text=chars); ss.subset(ft); ft.save(SUB)
    print('  子集字体 %.1f KB（原 769.2 KB）/ %d 个字形' % (
        os.path.getsize(SUB) / 1024, len(set(chars))))
    for p in plan:
        q.insert_text(fitz.Point(p['x'] * CM, p['y'] * CM), p['new'], fontsize=p['fs'],
                      fontfile=SUB, fontname='CUnits', color=p['col'])
    out.save(DST, garbage=4, deflate=True)
    # ---- 自检 ----
    new = fitz.open(DST)[0]
    tn = new.get_text()
    for p in plan:
        assert p['new'].replace('  ', ' ') in tn.replace('\xa0', ' ').replace('  ', ' ') or \
               all(w in tn for w in p['new'].split()), '新串缺: %r' % p['new'][:40]
    assert '\xa0' not in tn and '\xad' not in tn, '文字层出现 nbsp/软连字符'
    print('  文字层: 3 条旧串已消失、3 条新串在场、\\xa0/\\xad 计数 0 ✓  字符 %d → %d'
          % (len(t), len(tn)))
    for dpi in (150, 300):
        A = new.get_pixmap(dpi=dpi, alpha=False)
        B = pg.get_pixmap(dpi=dpi, alpha=False)
        a = np.frombuffer(A.samples, np.uint8).reshape(A.height, A.width, A.n)[:, :, :3].astype(int)
        b = np.frombuffer(B.samples, np.uint8).reshape(B.height, B.width, B.n)[:, :, :3].astype(int)
        d = np.abs(a - b).max(2) > 24
        m = np.zeros_like(d)
        for p in plan:
            x0, y0, x1, y1 = p['rect']
            m[int(y0 * dpi / 2.54) - 3:int(y1 * dpi / 2.54) + 3,
              int(x0 * dpi / 2.54) - 3:int((p['x'] + max(
                  fitz.Font(fontfile=os.path.join(os.environ.get('DSCI_FONTS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'fonts')), 'ARIAL.TTF')).text_length(
                      p['new'], fontsize=p['fs']) / CM, x1 - x0)) * dpi / 2.54) + 3] = True
        print('  %3d dpi: 总差异 %d px，三行字足迹以外 %d px' % (dpi, int(d.sum()), int((d & ~m).sum())))
        assert (d & ~m).sum() == 0, '改动越出三行字的足迹（redaction 可能伤到矢量）'
    if apply:
        shutil.copy2(SRC, '/tmp/FigS5_noc4_pre_cpatch.pdf')
        os.replace(DST, SRC)
        print('已替换', SRC)
    else:
        print('（试运行）产物', DST)
    return DST


if __name__ == '__main__':
    main(apply='--apply' in sys.argv)
