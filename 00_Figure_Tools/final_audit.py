import os
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
final_audit.py — 定稿候选 PDF 的一张体检表 (字体 / 可见最小字号 / 隐藏文字 / 面板字母)
================================================================================
用法:
    python3 final_audit.py                      # 体检 figure_plan 认定的 8 张定稿
    python3 final_audit.py Fig1_v3.pdf Fig2.pdf # 体检指定文件

判读口径 (与 CTM/期刊 preflight 对齐):
  * 字体      : 目标 = 只出现真 Arial (ArialMT / Arial-BoldMT / Arial-ItalicMT /
                Arial-BoldItalicMT)。Liberation* = 真 Arial 到位(10-13 12:27)之前的克隆;
                DejaVu* = matplotlib 默认; Helvetica = base-14 未嵌。**未嵌入 = 直接卡 preflight**。
  * 最小可见字号: 用"该 span 框内实际墨量"判定可见性 (render mode 3 的隐藏层不算),
                目标 >=5.5 pt (投稿硬底线), 理想 6-8 pt。
  * 隐藏字符  : 页内取出的字符里墨量 <1% 的字符数。大量隐藏字符 = 贴外部源 PDF 带进来的
                旧文字层 (FigS5 系列的 A-F 旧字母、Fig3 的重复 UMAP 标签), 期刊转 AI 时会露出来。
  * 面板字母  : 可见的单个大写字母 (>=6.4pt, 粗体) 按位置排序 = 该图的字母序列。
  * **子集化**  : 第三道字体门 (2026-10-15 加, 前两道查不出这类)。"已嵌入"和"是 Arial"都通过
                仍可能整套字库入 PDF: **单份内嵌流 > 100 KB = 未子集化**。
                实例: 面板字母用 8 pt Arial-Bold, 却内嵌 957 KB / 4,503 字形的整字库 (fig4/S4/6/S6
                第 2 批实测, 我方 FigS1/Fig5_v2 同病), fig6 本身才 883 KB。
                修法 = compose 之后、保存之前 `doc.subset_fonts()` (实测省 3.1 MB / 四张)。
                许可面: Arial 通常只允许渲染时子集嵌入。
                ⚠ **判据只用字节, 不用 numGlyphs**: matplotlib 子集会重写 glyf 但把 maxp.numGlyphs
                  留在原值 (实测 Fig2 `TZMXNH+ArialMT` = 77 KB / 4,503 glyphs —— 是正常子集)。
                  用 numGlyphs 会误伤; 只作参考显示。
  * Type3    : 有 Type3 字形时 get_text() 抽不到其文字 -> 文本层 grep 不可信, 需重渲源面板
                (ps.fonttype=42), 不要在成品上补。
  * 字母 house rule: 全稿面板字母 = Arial **常规** (2026-10-21 作者改口径「主图字体尽可能统一」;
                旧口径 8.0 pt Arial-BoldMT 已废)。Figure_1/2 例外: 10 pt。
                识别不再要求粗体 (见 _is_letter_font), 否则改完常规就报「字母 no」。
"""
import sys, os, io, collections
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
import fitz
import numpy as np
try:
    from fontTools.ttLib import TTFont
except Exception:                                        # fontTools 缺失时子集化门降级为未知
    TTFont = None

SUBSET_MAX_KB = 100            # 内嵌流**压缩字节** >100 KB -> 未子集化 (解压值仅信息; 见 font_census.py 说明)
LETTER_HOUSE_PT, LETTER_HOUSE_FONT = 8.0, 'arialmt'   # 全稿面板字母 house rule (2026-10-21 起为常规字重)
# 2026-10-21 作者改口径: **主图字体尽可能统一** -> Figure 2 的 A–E 改成常规字重 ArialMT
# (原位只动字母, 见 scripts/patch_fig2_letters.py)。所以字母识别不再要求粗体;
# 但只认 ArialMT / Arial-BoldMT 这两种"整页只写一个字母"的字重, 避免把数据点标签
# (如 Fig1 工作流程图里的腺嘌呤 A / 鸟嘌呤 G, 6 pt 斜体) 误判成面板字母。
LETTER_FONTS = ('arialmt', 'arial-boldmt', 'arialregular', 'arialbold')


def _is_letter_font(fn):
    return any(k in fn.lower() for k in LETTER_FONTS)

FIG = os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures'))
CM = 72 / 2.54
DPI = 150
FINALS = ['Fig1_v3.pdf', 'Fig2.pdf', 'Fig3_v10.pdf', 'Fig5_v2.pdf',
          'FigS1.pdf', 'FigS1_12lib.pdf', 'FigS2.pdf', 'FigS3_v3.pdf',
          'FigS5_noc2.pdf']

ARIAL = ('arialmt', 'arial-boldmt', 'arial-italicmt', 'arial-bolditalicmt',
         'arialregular', 'arialbold', 'arialitalic', 'arialbolditalic')


def clean(n):
    return n.split('+')[-1].replace(' ', '').replace('_', '')


def audit(path, verbose=True):
    doc = fitz.open(path)
    page = doc[0]
    pm = page.get_pixmap(dpi=DPI)
    arr = np.frombuffer(pm.samples, dtype=np.uint8).reshape(pm.height, pm.width, pm.n)
    gray = arr[:, :, :3].mean(axis=2)
    scale = DPI / 72.0

    emb = {}
    for f in page.get_fonts(full=True):
        xref, ext, ftype, base = f[0], f[1], f[2], f[3]
        emb[clean(base)] = (xref > 0 and ext not in ('n/a', ''), ftype)

    # ---- 子集化体检: 必须沿对象链走 (Type0 -> DescendantFonts -> FontDescriptor -> FontFile2)
    # ⚠ 不能用 fitz 的 extract_font(): 第 3 批实测它对 subset_fonts() 生成的字体返回空/报错,
    #    于是 >100KB 判据永远不触发 = **假通过**。同理只看 Type0 的 /BaseFont 会漏掉前缀改名。
    import font_census as _fc
    rows = _fc.census(path)
    type3 = []
    subset_census = {}
    prefix_mismatch = []
    name_streams = {}          # 硬门 B: 同 BaseFont 名 -> 必须同一份程序
    for r in rows:
        subset_census.setdefault(r['base'], (int(r['raw']), r['glyphs'] if r['glyphs'] is not None else -1,
                                            int(r['dec'])))
        if r.get('md5'):
            name_streams.setdefault(r['base'], set()).add(r['md5'])
        if r['mismatch']:
            prefix_mismatch.append('%s(/FontName=%s)' % (r['base'], r['name']))
    dup_names = sorted('%s(%d 份不同流)' % (k, len(v)) for k, v in name_streams.items() if len(v) > 1)
    for f in page.get_fonts(full=True):
        if f[2] == 'Type3':
            type3.append(clean(f[3]))
    not_subset = sorted(k for k, v in subset_census.items() if v[0] > SUBSET_MAX_KB)

    vis_chars = collections.Counter()      # font -> visible char count
    vis_min = {}                           # font -> min visible size
    hidden_chars = 0
    hidden_samples = []
    letters = []
    sizes = collections.Counter()
    for b in page.get_text('dict')['blocks']:
        if b['type'] != 0:
            continue
        for l in b['lines']:
            for s in l['spans']:
                t = s['text']
                if not t.strip():
                    continue
                fn = clean(s.get('font', '?'))
                x0, y0, x1, y1 = [v * scale for v in s['bbox']]
                sub = gray[max(0, int(y0)):int(y1) + 1, max(0, int(x0)):int(x1) + 1]
                ink = float((sub < 200).mean()) if sub.size else 0.0
                n = len(t.strip())
                if ink < 0.01:                       # 隐藏层
                    hidden_chars += n
                    if len(hidden_samples) < 6:
                        hidden_samples.append((t.strip()[:28], round(s['size'], 1), fn))
                else:
                    vis_chars[fn] += n
                    vis_min[fn] = min(vis_min.get(fn, 99), round(s['size'], 2))
                    sizes[round(s['size'], 2)] += n
                c = t.strip()
                if len(c) == 1 and c.isalpha() and c.isupper() and s['size'] >= 6.4 \
                        and ink >= 0.01 and _is_letter_font(fn):
                    letters.append((c, s['bbox'][0] / CM, s['bbox'][1] / CM, round(s['size'], 1)))

    w, h = page.rect.width / CM, page.rect.height / CM
    letter_sizes = sorted({sz for _, _, _, sz in letters})
    bad_font = sorted(f for f in vis_chars if f.lower() not in ARIAL)
    not_emb = sorted(f for f, (e, _) in emb.items() if not e)
    minvis = min(vis_min.values()) if vis_min else None
    seq = [c for c, _, _, _ in sorted(letters, key=lambda r: (round(r[2], 1), r[1]))]
    if verbose:
        print('=== %s   %.2f x %.2f cm, %d 页' % (os.path.basename(path), w, h, doc.page_count))
        print('    可见最小字号 : %s pt' % minvis)
        print('    可见字号谱   : %s' % '  '.join('%g:%d' % kv for kv in sorted(sizes.items())))
        print('    可见文字字体 : %s' % (', '.join('%s(%d字,最小%gpt)' % (f, vis_chars[f], vis_min[f])
                                                  for f in sorted(vis_chars))))
        print('    非 Arial 痕迹: %s   %s' % (', '.join(bad_font) if bad_font else '无',
                                             '<<< 需 Arial 化' if bad_font else '✓'))
        print('    未嵌入字体   : %s   %s' % (', '.join(not_emb) if not_emb else '无',
                                             '<<< preflight 会卡' if not_emb else '✓'))
        print('    隐藏字符     : %d %s' % (hidden_chars, ('<<< 旧文字层残留 ' + str(hidden_samples))
                                            if hidden_chars > 200 else ''))
        print('    面板字母(可见): %s' % ''.join(dict.fromkeys(seq)))
        if letter_sizes:
            ok = (letter_sizes == [LETTER_HOUSE_PT])
            print('    字母字号     : %s pt  %s  (house rule = %g pt %s)'
                  % (letter_sizes, '✓' if ok else '<<< 与 house rule 不一致',
                     LETTER_HOUSE_PT, LETTER_HOUSE_FONT.upper()))
        print('    内嵌字体     : %s' % ('  '.join('%s=%.0fKB/%dglyphs' % (k, v[0], v[1])
                                                  for k, v in sorted(subset_census.items())) or '无内嵌流'))
        print('    字体程序     : %s' % ', '.join('%s=%.0fKB/%d字' % (k, v[0], v[2]) for k, v in sorted(subset_census.items())))
        print('    未子集化     : %s   %s'
              % (', '.join('%s(%.0fKB)' % (k, subset_census[k][0]) for k in not_subset) if not_subset else '无',
                 '<<< 整字库嵌入 >%dKB; 修: compose 后、保存前 doc.subset_fonts()' % SUBSET_MAX_KB
                 if not_subset else '✓'))
        if prefix_mismatch:
            print('    ⚠前缀不一致 : %s\n                  WARN(不拦交): Type0 /BaseFont 无 XXXXXX+ 标签而描述符 /FontName 有 —— '
                  'subset_fonts() 的已知行为; PDF/A-1 §6.3.2 严格校验可能挑刺, 一行改 BaseFont 即可对齐'
                  % ', '.join(sorted(set(prefix_mismatch))))
        if dup_names:
            print('    ⚠同名不同体 : %s\n                  WARN(不拦交): 同一 BaseFont 名被多份**不同**内嵌程序共用'
                  '(fitz 标签按 文档+源字体 生成, 各容器自裁子集所致); ISO19005 §6.3.2 同名一致性会挑刺, '
                  '某些 RIP 按名去重会拿错字形。修: 改名时按流 md5 生成互不相同标签, 或去重共享对象'
                  % ', '.join(dup_names))
        print('    Type3        : %s   %s' % (', '.join(type3) if type3 else '无',
                                             '<<< 文本层抽不到, grep 不可信; 重渲源面板 fonttype=42'
                                             if type3 else '✓'))
        doc.close()
    return dict(name=os.path.basename(path), w=w, h=h, minvis=minvis, bad=bad_font,
                notemb=not_emb, hidden=hidden_chars, letters=''.join(dict.fromkeys(seq)),
                notsubset=not_subset, type3=type3, letter_sizes=letter_sizes,
                mismatch=prefix_mismatch, dup_names=dup_names, programs=subset_census)


if __name__ == '__main__':
    argv = [a for a in sys.argv[1:] if not a.startswith('-')] or FINALS
    rows = []
    for a in argv:
        p = a if os.path.isabs(a) else os.path.join(FIG, a)
        rows.append(audit(p))
    print('\n---------------- 汇总 ----------------')
    print('%-18s %-15s %6s %-8s %-8s %6s %6s %s' % ('文件', '尺寸(cm)', '最小pt', '字母',
                                                     '非Arial', '未嵌入', '隐藏字', '未子集化/Type3'))
    for r in rows:
        print('%-18s %-15s %6s %-8s %-8s %6s %6s %s' % (
            r['name'], '%.2f x %.2f' % (r['w'], r['h']), r['minvis'], r['letters'],
            'YES' if r['bad'] else 'no', 'YES' if r['notemb'] else 'no', r['hidden'],
            ('未子集:' + ','.join(r['notsubset']) if r['notsubset'] else '')
            + ('  Type3:' + ','.join(r['type3']) if r['type3'] else '')
            + ('  ⚠P×%d' % len(r['mismatch']) if r.get('mismatch') else '')
            + ('✓' if not r['notsubset'] and not r['type3'] and not r.get('mismatch') else '')))
    hard_fail = [r['name'] for r in rows
                 if r['bad'] or r['notemb'] or r['notsubset'] or r['type3']
                 or (r['minvis'] is not None and r['minvis'] < 5.5)]
    print('\n硬门 (非Arial/未嵌入/未子集化/Type3/最小pt<5.5) : %s'
          % ('FAIL -> ' + ', '.join(hard_fail) if hard_fail else '全部通过 ✓'))
