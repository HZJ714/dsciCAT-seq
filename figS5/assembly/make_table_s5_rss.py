#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_table_s5_rss.py — 给手稿加 **Table S5：六个状态各自的 top-10 eRSS**（稿件 v16 起的编号；旧脚本名里的 s6 即此表）（Codex #90 / G3 的后半）
================================================================================================
Codex 原话：「各面板 y 轴尺度不统一易误判特异性大小，建议统一坐标**并在同一表给 top
regulon / rank / RSS**」——前半（统一坐标）在 `scripts/figs5_rss_unified.py`，本脚本做后半。

为什么必须有一张表：#45/#90 判错的根因是**右列三格没有 y 刻度**，读者把左列刻度串过去读，
于是"MOL6 只有 0.33 / MFOL 接近 0.6"（正好是 MOL3_6 与 COP_NFOL 的真极值）。统一坐标解决了
"能不能读对"，这张表解决"要不要靠读"：数值直接给出来，且与图、CSV 三方同源。

数值来源 = 交付图用的同一份 `rss_scores_filtered_96.csv`（与笔记本重算的
`S5a_rss_scores_96.csv` 最大绝对差 0.000e+00，本脚本内断言）。

用法: PYTHONPATH=<work>/pylibs python3 scripts/make_table_s5_rss.py [--apply]
      不加 --apply 只打印将要插入的内容与位置（干跑）
输出: Manuscript_CTM_codex_addressed_v7_20261029.docx（--apply 时）
"""
import os, re, sys, zipfile, shutil, collections
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
import lxml.etree as ET
import pandas as pd

B = os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures'))
SRC = f'{B}/Manuscript_CTM_codex_addressed_v6_20261029.docx'
DST = f'{B}/Manuscript_CTM_codex_addressed_v7_20261029.docx'
W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
W14 = 'http://schemas.microsoft.com/office/word/2010/wordml'
ARCH = os.path.join(os.environ.get('DSCI_SCENIC_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pi/plots-scenic+')), 'rss_scores_filtered_96.csv')
NB = os.path.join(os.environ.get('DSCI_REPO_ROOT', '.'), 'figS5/figS5_out/S5a_rss_scores_96.csv')
ORDER = ['COP_NFOL', 'MFOL', 'MOL2', 'MOL3', 'MOL3_6', 'MOL6']
NCOL = ['State', 'Rank', 'eRegulon (motif label)', 'n regions', 'eRSS']
GRID = [1500, 800, 3300, 1300, 2146]           # 合计 9046 twips = 与 Table S1 同宽


def label(n):
    """'Creb5_extended_+_(172r)' -> ('Creb5 ext', '172', '+')"""
    m = re.match(r'^(.*?)_(extended|direct)?_\+_\((\d+)r\)$', n)
    if not m:
        m = re.match(r'^(.*?)_\+_\((\d+)r\)$', n)
        return m.group(1), m.group(2)
    return m.group(1) + (' ext' if m.group(2) == 'extended' else ''), m.group(3)


def build_rows():
    a = pd.read_csv(ARCH, index_col=0)
    b = pd.read_csv(NB, index_col=0)
    assert list(a.columns) == list(b.columns) and list(a.index) == list(b.index)
    d = float((a.values - b.values).__abs__().max())
    assert d == 0.0, '两份 eRSS 表不一致: %g' % d
    rows = []
    for st in ORDER:
        s = a.loc[st].sort_values(ascending=False)
        for r, (n, v) in enumerate(s.head(10).items(), start=1):
            nm, nreg = label(n)
            rows.append([st, str(r), nm, nreg, '%.3f' % v])
    return a, rows


def mk_p(text, bold=False, sz=20, style=None):
    p = ET.Element(W + 'p')
    if style:
        ET.SubElement(ET.SubElement(p, W + 'pPr'), W + 'pStyle').set(W + 'val', style)
    r = ET.SubElement(p, W + 'r')
    if bold:
        ET.SubElement(r, W + 'rPr')
        ET.SubElement(r.find(W + 'rPr'), W + 'b')
        ET.SubElement(r.find(W + 'rPr'), W + 'sz').set(W + 'val', str(sz))
    t = ET.SubElement(r, W + 't')
    t.text = text
    if text != text.strip():
        t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    return p


def mk_cell(text, bold=False):
    tc = ET.Element(W + 'tc')
    tcpr = ET.SubElement(tc, W + 'tcPr')
    ET.SubElement(tcpr, W + 'vAlign').set(W + 'val', 'center')
    p = ET.SubElement(tc, W + 'p')
    ppr = ET.SubElement(p, W + 'pPr')
    sp = ET.SubElement(ppr, W + 'spacing')
    sp.set(W + 'after', '20'); sp.set(W + 'line', '240'); sp.set(W + 'lineRule', 'auto')
    r = ET.SubElement(p, W + 'r')
    rpr = ET.SubElement(r, W + 'rPr')
    if bold:
        ET.SubElement(rpr, W + 'b')
    ET.SubElement(rpr, W + 'sz').set(W + 'val', '20')
    t = ET.SubElement(r, W + 't')
    t.text = text
    return tc


def mk_row(cells, header=False):
    tr = ET.Element(W + 'tr')
    for i, c in enumerate(cells):
        tc = mk_cell(c, bold=header)
        tc.find(W + 'tcPr').insert(0, ET.Element(W + 'tcW'))
        tc.find(W + 'tcPr/' + W + 'tcW').set(W + 'w', str(GRID[i])); tc.find(
            W + 'tcPr/' + W + 'tcW').set(W + 'type', 'dxa')
        tr.append(tc)
    return tr


def build_table(rows):
    tbl = ET.Element(W + 'tbl')
    tblpr = ET.SubElement(tbl, W + 'tblPr')
    ET.SubElement(tblpr, W + 'tblStyle').set(W + 'val', '6')
    ET.SubElement(tblpr, W + 'tblW').set(W + 'w', '0'); ET.SubElement(tblpr, W + 'tblW').set(W + 'type', 'auto')
    borders = ET.SubElement(tblpr, W + 'tblBorders')
    for side in ('top', 'left', 'bottom', 'right', 'insideH', 'insideV'):
        e = ET.SubElement(borders, W + side)
        e.set(W + 'val', 'single'); e.set(W + 'color', '000000'); e.set(W + 'sz', '4'); e.set(W + 'space', '0')
    grid = ET.SubElement(tbl, W + 'tblGrid')
    for g in GRID:
        ET.SubElement(grid, W + 'gridCol').set(W + 'w', str(g))
    tbl.append(mk_row(NCOL, header=True))
    for r in rows:
        tbl.append(mk_row(r))
    return tbl


def main(apply=False):
    rss, rows = build_rows()
    counts = {st: int((rss.loc[st] > 0.5).sum()) for st in ORDER}
    title = ('Table S5. Top-10 eRegulon specificity scores (eRSS) per lineage state, '
             'as ranked in Figure S5A.')
    note = ('eRSS values from the filtered 96-eRegulon specificity matrix, sorted descending within '
            'each state; identical to the ranking plotted in Figure S5A and to '
            '`S5a_rss_scores_96.csv`. "ext" marks extended eRegulons; n regions is the size of the '
            'cistrome. Numbers of eRegulons with eRSS > 0.5: COP_NFOL %d, MFOL %d, MOL2 %d, MOL3 %d, '
            'MOL3_6 %d, MOL6 %d (Creb5 0.513, Tbx2 0.512, Nfia 0.506, Mbnl2 0.505, Nfix 0.502). '
            'Labels are the SCENIC+ motif-to-TF annotations: Mbnl2 is annotated as a motif-driving '
            'factor although MBNL2 encodes an RNA-binding splicing factor rather than a '
            'sequence-specific DNA-binding transcription factor (Methods).') % tuple(counts[s] for s in ORDER)
    z = zipfile.ZipFile(SRC)
    xml = ET.fromstring(z.read('word/document.xml'))
    body = xml.find(W + 'body')
    kids = list(body)
    tbls = [i for i, e in enumerate(kids) if e.tag == W + 'tbl']
    last = tbls[-1]                                     # Table S5 的 tbl
    j = last + 1
    while j < len(kids) and kids[j].tag == W + 'p' and not ''.join(
            t.text or '' for t in kids[j].iter(W + 't')).strip():
        j += 1                                           # 跳过表后空段
    anchor = kids[j]
    print('插入位置: 第 %d 个 body 子元素之前（当前为 %r）' % (
        j, ''.join(t.text or '' for t in anchor.iter(W + 't'))[:48]))
    print('表体 %d 行（含表头）× %d 列；列宽 %s twips' % (len(rows) + 1, len(NCOL), GRID))
    for r in rows[:3] + rows[-2:]:
        print('   ', r)
    if not apply:
        print('\n[干跑] 加 --apply 才写文件。')
        print('脚注: ' + note[:200] + ' …')
        return
    shutil.copy(SRC, DST + '.tmp')
    pos = list(body).index(anchor)
    body.insert(pos, mk_p(note, sz=18))
    body.insert(pos, build_table(rows))
    body.insert(pos, mk_p(title, bold=True))
    for name, data in [('word/document.xml', ET.tostring(xml, xml_declaration=True,
                                                         encoding='UTF-8', standalone=True))]:
        pass
    zin = zipfile.ZipFile(SRC)
    zout = zipfile.ZipFile(DST, 'w', zipfile.ZIP_DEFLATED)
    for it in zin.infolist():
        buf = zin.read(it.filename)
        if it.filename == 'word/document.xml':
            buf = ET.tostring(xml, xml_declaration=True, encoding='UTF-8', standalone=True)
        zout.writestr(it, buf)
    zout.close(); zin.close()
    os.remove(DST + '.tmp')
    print('\n已写出 %s' % DST)
    return DST


if __name__ == '__main__':
    main(apply='--apply' in sys.argv)
