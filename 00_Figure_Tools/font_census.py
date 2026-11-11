import os
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""font_census.py —— 内嵌字体程序的**权威**清点(沿对象链走, 不依赖库内部解析)

为什么需要它 (2026-10-15 第 3 批实测):
    fitz 的 page.get_fonts() / doc.extract_font() 只看到 **Type0 层**。
    `doc.subset_fonts()` 子集化后, 真实的内嵌程序挂在
        Type0 /DescendantFonts[ N 0 R]  ->  CIDFontType2 /FontDescriptor M
        ->  FontDescriptor /FontName + /FontFile2 K
    链尾。第 3 批就是这种情况: **Type0 的 /BaseFont 仍是 `Arial-BoldMT`(无子集前缀),
    但描述符 /FontName 已变成 `DEWZOW+Arial-BoldMT`**, 程序只有 ~10 KB。
    只看 BaseFont 会误判"未子集化"; 只看 extract_font 会读到 0 KB 而**假通过**。
    另外 subset_fonts() 的解压字节可能 > 压缩字节 10 倍以上, 判阈值必须先确认单位:
    figS6 里 23.6 KB 压缩 / 150 KB 解压, 用解压值配 100 KB 阈会把**正常子集判成违规**。

判据 (最终版):
    以**压缩流字节**为准, 单份 > 100 KB(压缩) = 未子集化。
    解压值只作信息打印 (整库 ARIALBD.TTF 压缩后 ~957 KB 才是"未子集"的样子)。
    另报"前缀不一致": FontName 有 XXXXXX+ 而 BaseFont 没有 -> PDF/A 严格校验可能挑刺。

用法:
    python3 font_census.py [PDF ...]        # 不给参数 = 交付四张 + 本地定稿链
"""
import sys, os, re, io, glob, hashlib
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
import fitz
try:
    from fontTools.ttLib import TTFont
except Exception:
    TTFont = None

SUBSET_MAX_KB_RAW = 100          # 压缩流 >100 KB = 未子集化


def _obj(d, x):
    try:
        return d.xref_object(x, compressed=True) or ''
    except Exception:
        return ''


def walk(d, x0):
    """Type0 xref -> (BaseFont, FontName, 压缩KB, 解压KB, 字形数, FontFile 类型)"""
    s = _obj(d, x0)
    bf = (re.search(r'/BaseFont\s*/([A-Za-z0-9+\-_.]+)', s) or [None, '?'])[1] if '/BaseFont' in s else '?'
    m = re.search(r'/DescendantFonts\s*\[?\s*(\d+)\s+0\s+R', s)
    cid = int(m.group(1)) if m else None
    fn, ff_raw, ff_dec, ff_kind, ff_md5 = None, None, None, None, None
    cands = [x0] + ([cid] if cid else [])
    for c in cands:
        t = _obj(d, c)
        fm = re.search(r'/FontDescriptor\s+(\d+)\s+0\s+R', t)
        if fm:
            u = _obj(d, int(fm.group(1)))
            fn = (re.search(r'/FontName\s*/([A-Za-z0-9+\-_.]+)', u) or [None, '?'])[1]
            ff = re.search(r'/(FontFile2|FontFile3)\s+(\d+)\s+0\s+R', u)
            if ff:
                ff_kind = ff.group(1)
                i = int(ff.group(2))
                raw = d.xref_stream_raw(i) or b''
                ff_raw = len(raw)
                ff_dec = len(d.xref_stream(i) or b'')
                ff_md5 = hashlib.md5(raw).hexdigest()[:10]
    return bf, fn, (ff_raw or 0) / 1024.0, (ff_dec or 0) / 1024.0, ff_kind, ff_md5


def glyphs(d, x0):
    s = _obj(d, x0)
    m = re.search(r'/DescendantFonts\s*\[?\s*(\d+)\s+0\s+R', s)
    cands = [x0] + ([int(m.group(1))] if m else [])
    for c in cands:
        t = _obj(d, c)
        fm = re.search(r'/FontDescriptor\s+(\d+)\s+0\s+R', t)
        if not fm:
            continue
        u = _obj(d, int(fm.group(1)))
        ff = re.search(r'/(FontFile2|FontFile3)\s+(\d+)\s+0\s+R', u)
        if not ff:
            continue
        if ff.group(1) != 'FontFile2' or TTFont is None:
            return None
        try:
            return TTFont(io.BytesIO(d.xref_stream(int(ff.group(2))))).getGlyphOrder().__len__()
        except Exception:
            return None
    return None


def census(path):
    d = fitz.open(path)
    rows = []
    seen_x = set()
    for pno in range(len(d)):
        for x0 in range(1, d.xref_length()):
            if x0 in seen_x:
                continue
            s = _obj(d, x0)
            if '/Subtype/Type0' not in s and '/Subtype/CIDFontType2' not in s:
                continue
            if '/BaseFont' not in s:
                continue
            if '/Subtype/CIDFontType2' in s:
                continue          # 只从 Type0 入口走链, 避免重复计
            seen_x.add(x0)
            bf, fn, raw_kb, dec_kb, kind, md5s = walk(d, x0)
            rows.append(dict(base=bf, name=fn, raw=raw_kb, dec=dec_kb, kind=kind,
                             glyphs=glyphs(d, x0), md5=md5s,
                             mismatch=bool(fn and '+' in fn and '+' not in (bf or ''))))
    return rows


def main():
    args = [a for a in sys.argv[1:] if a.endswith('.pdf')]
    if not args:
        O = os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')
        args = [os.path.join(O, f) for f in ['fig4_composed.pdf', 'figS4_composed.pdf',
                                             'fig6_composed.pdf', 'figS6_composed.pdf']]
    bad = []
    for p in args:
        rows = census(p)
        print('=== %s  (%d 份内嵌程序)' % (os.path.basename(p), len(rows)))
        tot = 0
        byname = {}
        for r in sorted(rows, key=lambda r: -r['raw']):
            flag = '   <<< 未子集化' if r['raw'] > SUBSET_MAX_KB_RAW else ''
            mm = '   ⚠前缀不一致(BaseFont 无 XXXXXX+)' if r['mismatch'] else ''
            byname.setdefault(r['base'], set()).add(r.get('md5'))
            print('    %-28s %-30s 压缩 %6.1f KB  解压 %7.1f KB  字形 %s  流md5 %s%s%s'
                  % (r['base'], r['name'] or '?', r['raw'], r['dec'],
                     r['glyphs'] if r['glyphs'] is not None else '?', r.get('md5'), flag, mm))
            tot += r['raw']
            if r['raw'] > SUBSET_MAX_KB_RAW:
                bad.append((os.path.basename(p), r['base'], r['raw']))
            if r['mismatch']:
                bad.append((os.path.basename(p), r['base'], 'MISMATCH'))
        print('    内嵌程序合计(压缩) %.0f KB\n' % tot)
        for nm, hs in sorted(byname.items()):
            if len(hs) > 1:
                bad.append((os.path.basename(p), nm, '同名不同体 %d 份不同流' % len(hs)))
                print('    ⚠ 同名不同体: %s 由 %d 份**不同**程序共用 -> PDF/A §6.3.2 同名一致性会挑刺'
                      % (nm, len(hs)))
    if bad:
        print('FAIL/注意: ' + '; '.join('%s %s %s' % b for b in bad))
        return 1
    print('全部内嵌程序 ≤ %d KB(压缩) 且 BaseFont/FontName 前缀一致 ✓' % SUBSET_MAX_KB_RAW)
    return 0


if __name__ == '__main__':
    sys.exit(main())
