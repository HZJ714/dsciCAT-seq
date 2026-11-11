#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
finalize_fig1_v4.py — 把 Fig1_v4.pdf 装成交付块的 Figure_1.pdf（作者 2026-10-29 拍板）
=======================================================================================
背景：Figure 1 是块内唯一还挂着 10-21 旧版的一张。Codex G0/G1（= #78 的 C/D/E 三条）
已在 10-28 的 `Fig1_v4.pdf` 里落地：C 面板整块重画（坐标轴给基因组版本 hg38/mm10、
归一化口径写进图内、基因名 `MIR17HG`/`Hes7`、三块 tile）、D/E 轴标签 `log10`、
sci-CAR 鼠列灰纱+斜纹+`‡`、非等深度方法标 `*`。

作者本轮决定：**A 面板（AI 母版流程图）不改**——两轮条码的对应关系图注 A 与正文都写了；
**"图太密"也不改**——污染率即使换成双组学一张，密度也一样。所以 v4 即为定稿。

⚠ 装完之后 Figure_1.pdf **仍不过 house 硬门**，这是知情的豁免，不是漏网：
  · 最小字号 1.84 pt（A 示意图里的 `Biotin`/`S-S`，AI 母版，抬字会撞版）
  · <5.5 pt 的 span 116 个（其中 96 个在 B–F，是 v4 新加的小注记）
  · Helvetica(Type1) 13 span、Arial/Liberation 未子集、Type3 DejaVuSans、47 个隐藏字
  这些与 10-21 交付版同源（旧版 76 个小字、198 个隐藏字），v4 只多不减的是小注记。
  豁免理由与逐条清单写进 GENERATED_TABLE / CHANGELOG，别让它变成"没人看见的 FAIL"。

同步四处：块内 PNG（300 dpi，像素尺寸必须仍是 2008×2975）、`figures.docx` 内嵌图
（rId20.png）、`MD5_checksums.txt`（15 项）、以及图注（见 `edits8_fig1_legend_sync.py`）。
"""
import os, sys, hashlib, zipfile, shutil
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
import fitz

CM = 72 / 2.54
FIG = os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures'))
B = f'{FIG}/FINAL_SET_20261021'
SUP = f'{B}/_superseded'
SRC = f'{FIG}/Fig1_v4.pdf'
DST = f'{B}/Figure_1.pdf'
OLD_MD5 = '8e3d31f59f71'                       # 10-21 交付版
DOCX_IMG = 'word/media/rId20.png'              # figures.docx 里的 Figure 1 页图
WANT_PX = (2008, 2975)                         # 17.00 × 25.19 cm @300 dpi


def md5(p):
    return hashlib.md5(open(p, 'rb').read()).hexdigest()


def main(apply=False):
    cur = md5(DST)[:12]
    assert cur == OLD_MD5, '交付块 Figure_1.pdf 当前 md5 %s，不是登记的 10-21 件 (%s)' % (cur, OLD_MD5)
    print('  交付块现件 md5 %s ✓ 是 10-21 交付件' % cur)
    print('  候选   Fig1_v4.pdf md5 %s  %.2f × %.2f cm' % (
        md5(SRC)[:12], *(lambda r: (r.width / CM, r.height / CM))(fitz.open(SRC)[0].rect)))

    png = f'{B}/png/Figure_1.png'
    if apply:
        os.makedirs(SUP, exist_ok=True)
        bak = f'{SUP}/Figure_1_20261021_deprecated.pdf'
        if not os.path.exists(bak):
            shutil.copy2(DST, bak)
            print('  备份 → %s' % bak)
        shutil.copy2(SRC, DST)
        fitz.open(DST)[0].get_pixmap(dpi=300, alpha=False).save(png)
        patch_docx()
        md5list()
    else:
        print('  (dry-run，未写盘)')

    from PIL import Image
    im = Image.open(png)
    print('  PNG %dx%d px（要求 %dx%d）' % (im.width, im.height, *WANT_PX))
    assert (im.width, im.height) == WANT_PX, '像素尺寸变了 ⇒ docx extent 要同步改，先停下'


def patch_docx():
    from PIL import Image
    import io
    doc = f'{B}/figures.docx'
    zin = zipfile.ZipFile(doc)
    old = zin.read(DOCX_IMG)
    assert Image.open(io.BytesIO(old)).size == WANT_PX, 'docx 内 Figure 1 图尺寸异常'
    new = open(f'{B}/png/Figure_1.png', 'rb').read()
    assert Image.open(io.BytesIO(new)).size == WANT_PX
    tmp = doc + '.tmp'
    zout = zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED)
    for it in zin.infolist():
        zout.writestr(it, new if it.filename == DOCX_IMG else zin.read(it.filename))
    zout.close(); zin.close(); os.replace(tmp, doc)
    z = zipfile.ZipFile(doc)
    n = len([x for x in z.namelist() if x.startswith('word/media/')])
    assert md5_bytes(z.read(DOCX_IMG)) == md5(f'{B}/png/Figure_1.png')
    print('  figures.docx: %s 已替换（%d → %d B），media=%d 张不变' % (
        DOCX_IMG, len(old), len(new), n))
    assert n == 17, 'media 张数变了！'


def md5_bytes(b):
    return hashlib.md5(b).hexdigest()


def md5list():
    import glob
    names = [os.path.basename(x) for x in sorted(glob.glob(f'{B}/Figure_*.pdf'))
             + [f'{B}/Graphical_Abstract.pdf'] if os.path.exists(x)]
    lines = ['%s  %s' % (md5(f'{B}/{n}'), n) for n in names]
    open(f'{B}/MD5_checksums.txt', 'w').write('\n'.join(lines) + '\n')
    print('  MD5_checksums.txt 重算 %d 行' % len(lines))


if __name__ == '__main__':
    main(apply='--apply' in sys.argv)
