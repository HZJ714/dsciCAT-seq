#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
figkit.py — 投稿图拼版共用工具 (Fig3 / FigS3 用, 亦可复用于其它图)
================================================================================
职责:
  1) 字体: 加载 Liberation Sans (= Arial 度量) 供 fitz TextWriter / matplotlib 使用
  2) 文字层读取: spans -> 去 (geom_text_repel 白色光晕造成的 17 份重复) -> 标签对象
  3) 墨迹框: 按实际字形上下伸出量算 ink bbox (比字体 bbox 严格得多, 用来判真重叠)
  4) resolve_labels(): 确定性文字-文字碰撞消解 (迭代推开, 可限定面板边界)
  5) retext(): 整页 redact (只删文字, 保留矢量/图像) + TextWriter 全量重写
  6) save_doc(): 临时文件 + os.replace (避免 "must be incremental" 报错)
单位: 内部一律用 **cm** 坐标 (x 向右, y 向下, 原点=页面左上角), 输出时乘 CM 转 pt。
"""
import os, sys, collections
# 持久盘依赖: 本环境每轮重置 (/opt/conda 里的 pymupdf、~/.fonts 里的字体都会消失),
# 所以 pymupdf 装在 <work>/pylibs、字体放在 <work>/fonts, 这里自动引导。
for _p in (os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')),):
    if os.path.isdir(_p) and _p not in sys.path:
        sys.path.insert(0, _p)
import fitz
import numpy as np

CM = 72 / 2.54
FDIR = next(d for d in (os.environ.get('DSCI_FONTS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'fonts')), os.path.expanduser('~/.fonts'))
            if os.path.exists(os.path.join(d, 'LiberationSans-Regular.ttf')))
# ── 字面来源 (2026-10-13 起优先真 Arial) ──────────────────────────────────────────
# 用户提供正版 ARIAL*.TTF 放持久盘 <work>/fonts; 没有时退回 Liberation Sans
# (Red Hat 的 Arial 替身)。实测两者 hmtx 推进宽度逐字符相同 (比值 1.0000), 所以
# 任一方缺席都只改字形不改排版, 历史校验数据仍然成立。
#   Arial:      cap 1467/2048=0.7163  xH 1062/2048=0.5186  hhea 1854/-434
#   Liberation: cap 1409/2048=0.6880  xH 1082/2048=0.5283
ARIAL_FILES = {'r': 'ARIAL.TTF', 'b': 'ARIALBD.TTF', 'i': 'ARIALI.TTF', 'bi': 'ARIALBI.TTF'}
LIB_FILES = {'r': 'LiberationSans-Regular.ttf', 'b': 'LiberationSans-Bold.ttf',
             'i': 'LiberationSans-Italic.ttf', 'bi': 'LiberationSans-BoldItalic.ttf'}
HAS_ARIAL = all(os.path.exists(os.path.join(FDIR, v)) for v in ARIAL_FILES.values())
FFILES = {k: os.path.join(FDIR, (ARIAL_FILES if HAS_ARIAL else LIB_FILES)[k])
          for k in ('r', 'b', 'i', 'bi')}
FAM = 'Arial' if HAS_ARIAL else 'Liberation Sans'
F_REG, F_BOLD, F_ITAL, F_BOLDITAL = (FFILES[k] for k in ('r', 'b', 'i', 'bi'))
FONT = {k: fitz.Font(fontfile=v) for k, v in FFILES.items() if os.path.exists(v)}
# 资源名必须与任何已嵌入的 BaseFont 名不同, 否则 PyMuPDF 会复用旧子集导致缺字
FN = {'r': 'FigAriR', 'b': 'FigAriB', 'i': 'FigAriI', 'bi': 'FigAriBI'}
S_MIN, S_MAX = 6.0, 8.0                          # 投稿字号窗口
# Arial / Liberation Sans em 比例
CAP, XH, ASC, DESC = 0.716, 0.519, 0.905, 0.212
UPCH = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789bdfhklt/−+%(),.#')
DNCH = set('gjpqy_')


def use_matplotlib_arial():
    """让 matplotlib 用真 Arial (缺则 Liberation Sans), 并以 TrueType 嵌入

    要点:
      * 四个面 (Regular/Bold/Italic/BoldItalic) 都要 addfont, 否则 bold/italic 会
        回退到 DejaVu, 同一张图里混两种字形;
      * 显式点名 [FAM, 'Liberation Sans', 'DejaVu Sans'] —— 不能只写 'Arial' 就指望
        系统找到, matplotlib 找不到时会静默用内置 DejaVu Sans (比 Arial 宽 7~10%,
        FigS1 当年就是这么变成 DejaVu 的);
      * pdf.fonttype=42 (TrueType 子集, AI 里是活字), 不是默认的 3 (Type3 非嵌入,
        Illustrator 只能当图形, 期刊 preflight 也会挑)。
    """
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import font_manager, rcParams
    # 两套都注册: 首选族由 HAS_ARIAL 决定, 但历史脚本里写死的 'Liberation Sans'
    # 也必须还能解析到真的 Liberation (否则 matplotlib 会悄悄给 DejaVu)
    for k in ('r', 'b', 'i', 'bi'):
        for base in (ARIAL_FILES, LIB_FILES):
            fp = os.path.join(FDIR, base[k])
            if os.path.exists(fp):
                font_manager.fontManager.addfont(fp)
    # 注册成功与否以 fontManager 里的 family 名为准 (matplotlib 的 FontEntry.name)
    fam_ok = any(f.name == FAM for f in font_manager.fontManager.ttflist)
    fam_use = FAM if fam_ok else 'Liberation Sans'
    rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': [fam_use, 'Liberation Sans', 'DejaVu Sans'],
                     'pdf.fonttype': 42, 'ps.fonttype': 42})
    return fam_use


# 旧函数名保留 (Fig3/FigS3/FigS5 等脚本仍在调用), 行为已升级为 Arial 优先
use_matplotlib_liberation = use_matplotlib_arial


def rgb(v):
    """PyMuPDF span color (int sRGB) -> (r,g,b)"""
    return ((v >> 16 & 255) / 255, (v >> 8 & 255) / 255, (v & 255) / 255)


# ────────────────────────────── 文字层读取 ──────────────────────────────
# ─────────────────── 细胞类型显示名改写 (2026-10-19) ───────────────────
# 'Ependyma' 是数据里的 cluster 原名 (FY/3k_dsci_ht_10x_rna_celltype.csv 等), 指的却是
# 室管膜这层**组织**; 而它出现在 "Celltypes" 轴/图例里, 与其余 32 个细胞类型并列
# => 显示名统一改成 'Ependymal cells'。只改**落到图上的文字**: 数据键、配色索引、
# 编号映射 (num_map) 一律不动, 所以它在 33 类里仍是 8 号。
# 全交付集只有 Figure_2/3/S2/S3 含该词 (已扫描确认), 故放全局不影响其它图。
# \b 边界: 'Ependymal' 内部含 'Ependyma', 但后面接字母 => 不会被二次改写。
TYPE_RENAME = __import__('re').compile(r'\bEpendyma\b')
TYPE_WITH = 'Ependymal cells'


def rename_type(t):
    """标签文字 -> 投稿用显示名 (幂等)"""
    return TYPE_RENAME.sub(TYPE_WITH, t)


def label_w_cm(t, t0, z, bold, orig_cm):
    """改名后的标签必须用**写图时同一套度量**重算宽度: labels() 的 w 会被
    ink()/overlaps()/resolve() 拿去做碰撞消解, 沿用旧宽会漏判重叠。未改名则原样返回。"""
    if t == t0:
        return orig_cm
    return fitz.Font.text_length(FONT['b' if bold else 'r'], t, fontsize=z) / CM


def labels(page, sizes=None, clip=None, tol=0.15, keep_dir=False):
    """page -> 去重后的标签列表 [{t,x,y,w,h,z,c,bold,dir}] (cm, 中心坐标)

    sizes: 只取这些字号 (None=全部); clip: fitz.Rect(pt) 限定区域
    geom_text_repel 的白晕会把同一标签画成 ~17 份微偏移 span, 这里按
    (文本, 位置) 聚类去重, 保留 n 作为 halo 份数信息。
    """
    out = []
    for b in page.get_text('dict')['blocks']:
        if b['type'] != 0:
            continue
        for l in b.get('lines', []):
            for s in l['spans']:
                t0 = s['text'].strip()
                if not t0:
                    continue
                t = rename_type(t0)                      # 显示名改写 (Ependyma -> Ependymal cells)
                bold0 = ('bold' in s['font'].lower() or 'black' in s['font'].lower())
                z = round(s['size'], 2)
                if sizes is not None and z not in sizes:
                    continue
                r = fitz.Rect(s['bbox'])
                if clip is not None and not clip.intersects(r):
                    continue
                cx, cy = (r.x0 + r.x1) / 2 / CM, (r.y0 + r.y1) / 2 / CM
                hit = None
                for g in out:
                    if g['t'] == t and abs(g['x'] - cx) < tol and abs(g['y'] - cy) < tol:
                        hit = g; break
                if hit:
                    hit['n'] += 1
                    col = rgb(s['color'])
                    if max(hit['c']) - min(hit['c']) < 0.02 and max(col) - min(col) >= 0.02:
                        hit['c'] = col          # 原先存的是白晕副本 -> 换成正文色
                    continue
                g = dict(t=t, x=cx, y=cy,
                         w=label_w_cm(t, t0, z, bold0, r.width / CM), h=r.height / CM, z=z,
                         c=rgb(s['color']), bold=bold0, n=1)
                if keep_dir:
                    g['dir'] = tuple(round(v, 2) for v in l['dir'])
                out.append(g)
    return out


def ink_h(t, z):
    """真实墨迹高 (cm)。ink()/write_spans() 里的基线由 g['h'] 推, 所以 h 必须填这个值:
       填成 w*0.75 会把文字整体下移 (h - h_true)/2, 且**每个标签下移量不同** (w 不同),
       于是 resolve 算出的安全间距在成品里失效 -> 出现斜角擦碰。"""
    em = z / 72 * 2.54
    top = CAP if any(c in UPCH for c in t) else XH
    bot = DESC if any(c in DNCH for c in t) else 0.0
    return (top + bot) * em


def ink(g, scale=1.0, pad=0.0):
    """标签的真实墨迹框 (cm Rect, pt 单位输出) — 竖排(dir 纵向)时宽高互换"""
    z = g['z'] * scale
    em = z / 72 * 2.54                      # 1em 高 (cm)
    vert = abs(g.get('dir', (1, 0))[1]) > 0.5 if g.get('dir') else False
    up = any(c in UPCH for c in g['t'])
    dn = any(c in DNCH for c in g['t'])
    top = (CAP if up else XH) * em
    bot = (-DESC if dn else 0.0) * em
    if vert:
        # 竖排 (write_spans 里是先横排再绕中心旋转 90°): 基线到 cap/x-height 的墨迹
        # 宽度 = (top-bot)*em 且关于中心对称; 长度方向 = 文字长度 = w/h 里较长的那个
        ln = max(g['w'], g['h']) * scale
        hw = (top - bot) * em / 2
        return fitz.Rect((g['x'] - hw - pad) * CM, (g['y'] - ln / 2 - pad) * CM,
                         (g['x'] + hw + pad) * CM, (g['y'] + ln / 2 + pad) * CM)
    # 水平: 字体 bbox 中心 -> 基线位置 ≈ y + h/2 - DESC*em
    base = g['y'] + g['h'] * scale / 2 - DESC * em
    return fitz.Rect((g['x'] - g['w'] * scale / 2 - pad) * CM, (base - top - pad) * CM,
                     (g['x'] + g['w'] * scale / 2 + pad) * CM, (base - bot + pad) * CM)


def overlaps(gs, scale=1.0, pad=0.0, minarea=0.004):
    """真墨迹重叠对 [(i,j,area_cm2)]"""
    box = [ink(g, scale, pad) for g in gs]
    bad = []
    for i in range(len(gs)):
        for j in range(i + 1, len(gs)):
            ir = box[i] & box[j]
            if ir.is_valid and ir.get_area() > minarea * CM ** 2:
                bad.append((i, j, ir.get_area() / CM ** 2))
    return bad


# ────────────────────────────── 碰撞消解 ──────────────────────────────
def resolve(gs, scale=1.0, gap=0.02, bounds=None, iters=600, push=0.6, mindisp=0.004,
            fixed=(), padbox=None):
    """确定性推开: 每轮取所有真重叠对, 沿最小分离轴移动 (不动 fixed 下标)。
    gs: [{x,y,w,h,z,dir?}] 会被原地修改; bounds=(x0,y0,x1,y1) cm 软边界。
    返回残余重叠数。"""
    fixed = set(fixed)

    def box(g):
        return ink(g, scale, gap)

    def clamp_():
        """把标签夹回 bounds 内。必须在**每轮开头与收尾**都跑:
           早先只在“有重叠才移动”的分支里夹紧, 初始无重叠时直接 break -> 标签停在界外。"""
        if not bounds:
            return
        x0, y0, x1, y1 = bounds
        for g in gs:
            bw, bh = g['w'] * scale / 2, g['h'] * scale / 2
            g['x'] = min(max(g['x'], x0 + bw), max(x0 + bw, x1 - bw))
            g['y'] = min(max(g['y'], y0 + bh), max(y0 + bh, y1 - bh))

    for it in range(iters):
        clamp_()
        bad = []
        bx = [box(g) for g in gs]
        for i in range(len(gs)):
            for j in range(i + 1, len(gs)):
                ir = bx[i] & bx[j]
                if ir.is_valid and ir.get_area() > 1e-6:
                    bad.append((i, j, ir.width / CM, ir.height / CM))
        if not bad:
            break            # 不能直接 return: 下面的统一夹紧必须跑, 否则标签会停在界外
        for i, j, ox, oy in bad:
            if ox <= 0 and oy <= 0:
                continue
            d = max(min(ox, oy) * push, 0.006)
            a, b = gs[i], gs[j]
            if i in fixed and j in fixed:
                continue
            if i in fixed:
                dx = d if a['x'] > b['x'] else -d
                dy = d if a['y'] > b['y'] else -d
                b['x'] += (dx if ox <= oy else 0.0); b['y'] += (dy if ox > oy else 0.0)
            elif j in fixed:
                dx = d if b['x'] > a['x'] else -d
                dy = d if b['y'] > a['y'] else -d
                a['x'] += (dx if ox <= oy else 0.0); a['y'] += (dy if ox > oy else 0.0)
            else:
                sgnx = 1 if a['x'] >= b['x'] else -1
                sgny = 1 if a['y'] >= b['y'] else -1
                if ox <= oy:
                    a['x'] += sgnx * d / 2; b['x'] -= sgnx * d / 2
                else:
                    a['y'] += sgny * d / 2; b['y'] -= sgny * d / 2
    clamp_()
    bx = [box(g) for g in gs]
    n = 0
    for i in range(len(gs)):
        for j in range(i + 1, len(gs)):
            ir = bx[i] & bx[j]
            if ir.is_valid and ir.get_area() > 1e-6:
                n += 1
    return n


# ────────────────────────────── 重写文字层 ──────────────────────────────
def clamp(z, lo=S_MIN, hi=S_MAX):
    return min(max(z, lo), hi)


def strip_page(page, keep_images=True, keep_graphics=True):
    """整页只删文字 (redact 覆盖全页; 矢量与图像按需保留)"""
    page.add_redact_annot(page.rect)
    page.apply_redactions(
        images=fitz.PDF_REDACT_IMAGE_NONE if keep_images else fitz.PDF_REDACT_IMAGE_REMOVE,
        graphics=fitz.PDF_REDACT_LINE_ART_NONE if keep_graphics
        else fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED)
    page.clean_contents()


def write(page, items, scale=1.0, color=(0, 0, 0), anchor='c'):
    """items: [{t,x,y,z,bold}] (cm) -> TextWriter 写入 (x,y=中心)"""
    tw = fitz.TextWriter(page.rect)
    for g in items:
        f = FONT['b' if g.get('bold') else 'r']
        z = g['z'] * scale
        em = z / 72 * 2.54
        w = fitz.Font.text_length(f, g['t'], fontsize=z) / CM
        up = any(c in UPCH for c in g['t'])
        base_y = g['y'] + (CAP if up else XH) * em / 2          # 中心 -> 基线
        x0 = g['x'] - w / 2 if anchor == 'c' else (g['x'] if anchor == 'l' else g['x'] - w)
        col = g.get('c', color)
        tw.append((x0 * CM, base_y * CM), g['t'], font=f, fontsize=z)
    tw.write_text(page, color=None, overlay=True)
    return tw


def _base_y(g, em):
    """水平文字: 给定中心 y -> 基线 y (cm)"""
    t = g['t']
    up = any(c in UPCH for c in t)
    dn = any(c in DNCH for c in t)
    top = CAP if up else XH
    bot = -DESC if dn else 0.0
    return g['y'] + (top + bot) * em / 2


def write_spans(page, items):
    """批量写文字。items: [{t,x,y,z,bold?,c?,rot?}]  x,y = 墨迹中心 (cm);
    rot=90 自下向上读。按 (颜色,字重) 分组, 每组一个 TextWriter 一次提交
    (PyMuPDF 的 TextWriter 颜色在构造时给定; 逐条 insert_text 在大页面上极慢)。"""
    hor = [g for g in items if not g.get('rot')]
    ver = [g for g in items if g.get('rot')]
    groups = {}
    OFFS = ((0, -1), (0, 1), (-1, 0), (1, 0), (-.72, -.72), (.72, -.72), (-.72, .72), (.72, .72))
    for g in hor + ver:
        if not g.get('halo'):
            continue
        for dx, dy in OFFS:
            h = dict(g, x=g['x'] + dx * 0.018, y=g['y'] + dy * 0.018,
                     c=(1, 1, 1), halo=False)
            (ver if g.get('rot') else hor).append(h)
    for g in hor:
        key = (tuple(round(c, 3) for c in g.get('c', (0, 0, 0))), bool(g.get('bold')))
        groups.setdefault(key, []).append(g)
    okeys = sorted(groups, key=lambda k: 0 if min(k[0]) > 0.98 else 1)
    for col, bold in okeys:
        gs = groups[(col, bold)]
        f = FONT['b' if bold else 'r']
        tw = fitz.TextWriter(page.rect, color=col)
        for g in gs:
            em_cm = g['z'] / 72 * 2.54
            L = fitz.Font.text_length(f, g['t'], fontsize=g['z'])
            _tw_words(tw, f, g['t'], g['x'] * CM - L / 2, _base_y(g, em_cm) * CM, g['z'])
        tw.write_text(page, color=col)
    for g in ver:
        f = FONT['b' if g.get('bold') else 'r']
        em_cm = g['z'] / 72 * 2.54
        L = fitz.Font.text_length(f, g['t'], fontsize=g['z'])
        col = g.get('c', (0, 0, 0))
        tw = fitz.TextWriter(page.rect, color=col)
        tw.append(fitz.Point(g['x'] * CM - L / 2, _base_y(g, em_cm) * CM), g['t'],
                  font=f, fontsize=g['z'])
        tw.write_text(page, color=col,
                      morph=(fitz.Point(g['x'] * CM, g['y'] * CM),
                             fitz.Matrix(90 if g['rot'] == 90 else -90)))


def _tw_words(tw, f, t, x, y, z):
    """含空格的串**按词逐个 append**（空格用 text_length 的推进量补）。
    原因: TextWriter 整串提交时给这个 TTF 子集建的 ToUnicode 把空格映射成 U+00A0,
    于是 `get_text()` 读出 'Ependymal\\xa0cells' —— 图上看不出, 但文字层检索/校对/
    “有没有改成 Ependymal cells” 的自动检查全部漏判。拆词后回读是正常 U+0020,
    实测渲染等价：总墨量到小数第 4 位不变，串尾坐标差 0.028 pt（内容流 3 位小数取整），
    1200 dpi 下 0.029% 的像素因反锯齿翻转 —— 肉眼与 RIP 都看不出。竖排分支不走这里
    （旋转矩阵按整串算，拆词会散架）。
    """
    if ' ' not in t:
        tw.append(fitz.Point(x, y), t, font=f, fontsize=z)
        return
    sp = fitz.Font.text_length(f, ' ', fontsize=z)
    for i, w in enumerate(t.split(' ')):
        if i:
            x += sp
        tw.append(fitz.Point(x, y), w, font=f, fontsize=z)
        x += fitz.Font.text_length(f, w, fontsize=z)


def cover(page, rects, color=(1, 1, 1)):
    """一次 Shape 提交多个白色覆盖框 (rects: fitz.Rect 列表, pt)"""
    sh = page.new_shape()
    for r in rects:
        sh.draw_rect(r)
    sh.finish(color=None, fill=color, width=0)
    sh.commit()


_PDFFONT = {}


def pdf_font(kind='r'):
    """按 TTF 建一个 fitz.Font (度量用); kind in r/b/i/bi"""
    if kind not in _PDFFONT:
        _PDFFONT[kind] = fitz.Font(fontfile=FFILES[kind])
    return _PDFFONT[kind]


def _fk(bold=False, italic=False):
    return 'bi' if (bold and italic) else ('b' if bold else ('i' if italic else 'r'))


def put_text(page, pt, text, size=6.5, bold=False, italic=False, color=(0, 0, 0), fname=None):
    """insert_text 的正字 Arial 版。

    旧脚本一律写 fontname='helv'/'hebo' (base-14 Helvetica): 那种字**不嵌进 PDF**,
    期刊 preflight 直接判不合格, 而且 Helvetica 与 Arial 度量不同 -> 拼版宽度会算歪。
    这里强制带 fontfile, 字体随文件嵌入。"""
    k = _fk(bold, italic)
    return page.insert_text(pt, text, fontsize=size, fontfile=FFILES[k],
                            fontname=(fname or FN[k]), color=color)


def text_width(text, size=6.5, bold=False, italic=False):
    """与 put_text 同字体的文本宽度 (cm) —— 替代 fitz.get_text_length(fontname='helv')"""
    return pdf_font(_fk(bold, italic)).text_length(text, fontsize=size) * 2.54 / 72.0


def letter(page, x, y, ch, size=8.0, fname=None):
    """面板字母: (x,y) = 墨迹左上角 (cm)
    fname: 字体资源名。拼版页上如果已经 show_pdf_page 进来了同名资源 (子集字体),
    PyMuPDF 会复用那个旧子集 -> 新字母缺字 (变 \x00), 所以拼版时请传唯一名。"""
    em = size / 72 * 2.54
    page.insert_text(fitz.Point(x * CM, (y + CAP * em) * CM), ch, fontsize=size,
                     fontfile=F_BOLD, fontname=fname or FN['b'], color=(0, 0, 0))


def save_doc(doc, path, clean=False):
    import tempfile
    fd, tmp = tempfile.mkstemp(suffix='.pdf', dir=os.path.dirname(os.path.abspath(path)))
    os.close(fd)
    doc.save(tmp, garbage=4, deflate=True, clean=clean)
    doc.close()
    os.replace(tmp, path)
    return path


def audit(path):
    """字号分布 + NUL 计数"""
    sizes = collections.Counter(); nul = 0
    for pg in fitz.open(path):
        for b in pg.get_text('dict')['blocks']:
            for l in b.get('lines', []):
                for s in l['spans']:
                    t = s['text']
                    nul += t.count(chr(0))
                    t = t.replace(chr(0), '').strip()
                    if t:
                        sizes[round(s['size'], 1)] += len(t)
    return dict(sorted(sizes.items())), nul


def _strip_bt_et(data):
    """从 PDF 内容流里去掉 BT...ET 文本块。**字符串面量 / hex 串 / inline image 里的
    "BT"/"ET" 不算关键字** —— 旧版正则 \bBT\b.*?\bET\b 在本项目踩过雷: FigS2 的细胞类型
    标签就叫 "ET" (写成 (ET) 在字符串里), 正则把文本块提前终止, 内容流被截坏, MuPDF 报
    syntax error, 该段之后的文字整段消失。这里按 PDF 词法扫一遍。
    返回 (新流, 是否有改动)"""
    WS = b' \t\r\n\f\x00'
    DELIM = WS + b'()<>[]{};/\x7f'
    n = len(data)
    out = bytearray()
    i = 0
    in_text = False

    def read_string(j):                       # data[j] == '(' -> 返回右括号之后位置
        depth = 1
        j += 1
        while j < n and depth:
            c = data[j]
            if c == 0x5C:                     # 反斜杠转义, 跳两字节
                j += 2
                continue
            if c == 0x28:
                depth += 1
            elif c == 0x29:
                depth -= 1
            j += 1
        return j

    def read_dict(j):                         # data[j:j+2] == '<<' -> 返回 '>>' 之后
        depth = 1
        j += 2
        while j < n and depth:
            two = data[j:j + 2]
            if two == b'<<':
                depth += 1
                j += 2
                continue
            if two == b'>>':
                depth -= 1
                j += 2
                continue
            if data[j:j + 1] == b'(':
                j = read_string(j)
                continue
            j += 1
        return j

    while i < n:
        ch = data[i]
        if ch == 0x28:                        # '('
            j = read_string(i)
            if not in_text:
                out += data[i:j]
            i = j
            continue
        if ch == 0x3C:                        # '<'
            if data[i:i + 2] == b'<<':
                j = read_dict(i)
            else:
                k = data.find(b'>', i)
                j = n if k < 0 else k + 1
            if not in_text:
                out += data[i:j]
            i = j
            continue
        if ch in WS:
            if not in_text:
                out.append(ch)
            i += 1
            continue
        j = i
        while j < n and data[j] not in DELIM:
            j += 1
        tok = data[i:j]
        i = j
        if tok == b'BT':
            in_text = True
            continue
        if tok == b'ET':
            in_text = False
            continue
        if tok == b'BI' and not in_text:       # inline image: 原样搬到 EI
            k = i
            while True:
                k = data.find(b'EI', k)
                if k < 0:
                    break
                if k + 2 >= n or data[k + 2] in DELIM:
                    break
                k += 2
            k = n if k < 0 else k + 2
            out += data[i - 2:k]
            i = k
            continue
        if not in_text:
            out += tok
            if i < n and data[i] not in WS:   # 保留 token 之间的分隔符, 免得粘连
                pass
    return bytes(out), len(out) != n


def strip_text_stream(doc, page):
    """只删文字, 不碰图形/图像: 直接从页面内容流里去掉 BT...ET 文本块。
    (apply_redactions 会破坏密集微矢量点, 见 README_Fig3.md 的说明)"""
    import re
    pat = re.compile(rb'\bBT\b.*?\bET\b', re.S)
    n = 0
    for xref in page.get_contents():
        cont = doc.xref_stream(xref)
        if cont is None:
            continue
        new = pat.sub(b'', cont)
        if new != cont:
            doc.update_stream(xref, new)
            n += 1
    return n


def text_block_ranges(b):
    """PDF 内容流字节 -> [(start, end)] 每个 BT..ET 文本对象区间 (含头尾)。

    为什么不用正则 \bBT\b.*?\bET\b (figkit.strip_text_stream 的老做法, 2026-10-14 发现会吃图形):
    字符串 (…)/<…>/名字 /xx 以及 % 注释 里出现的 "BT"/"ET" 也会被正则当作定界符 —— 贴进来的
    参考图用 8-bit 编码, 串里碰到 'ET' 后非贪婪匹配一路吞到下一个真 ET, 中间的矢量线直接消失
    (Fig1_v3 实测丢 4 条 1px 横线, 位于面板 C y 9.03-10.11 cm)。
    这里按 PDF 词法扫一遍: 字符串/名字/注释 整体跳过, 只对真正的操作符做 BT/ET 配对。"""
    DELIM = b' \t\r\n/()<>[]{}%'
    out, stack, i, n = [], [], 0, len(b)
    while i < n:
        c = b[i]
        if c == 0x25:                                        # % 注释到行尾
            j = b.find(b'\n', i)
            i = n if j < 0 else j + 1
            continue
        if c == 0x28:                                        # ( 字符串, 可嵌套, 可转义
            i += 1
            depth = 1
            while i < n and depth:
                if b[i] == 0x5c:                             # backslash 转义
                    i += 2
                    continue
                if b[i] == 0x28:
                    depth += 1
                elif b[i] == 0x29:
                    depth -= 1
                i += 1
            continue
        if c == 0x3c:                                        # < 十六进制串
            j = b.find(b'>', i)
            i = n if j < 0 else j + 1
            continue
        if c == 0x2f:                                        # / 名字
            i += 1
            while i < n and b[i] not in DELIM:
                i += 1
            continue
        if c == 0x5b:                                      # [ 数组开始
            out.append((i, i + 1, 'open', b'['))
            i += 1
            continue
        if c == 0x5d:                                      # ] 数组结束
            out.append((i, i + 1, 'close', b']'))
            i += 1
            continue
        if c in DELIM:
            i += 1
            continue
        j = i
        while j < n and b[j] not in DELIM:
            j += 1
        w = b[i:j]
        if w == b'BT':
            stack.append(i)
        elif w == b'ET' and stack:
            out.append((stack.pop(), j))
        i = j
    out.extend((s, n) for s in stack)                        # 未闭合 -> 到结尾
    return sorted(out)


SHOW_OPS = (b'Tj', b'TJ', b"'", b'"')
TEXT_STATE = (b'BT', b'ET', b'Tf', b'Td', b'TD', b'Tm', b'TL', b'Ts', b'Tr', b'TC',
              b'TW', b'TZ', b'gs', b'q', b'Q', b'cm', b'w', b'J', b'j', b'M', b'd', b'ri',
              b'i', b'CS', b'scn', b'sc', b'S', b'f', b'f*', b'F', b'B', b'B*', b'b', b'b*',
              b'W', b'W*', b'm', b'l', b'c', b'v', b'y', b'h', b're', b'BI', b'ID', b'EI',
              b'D0', b'D1', b'MP', b'DP', b'BMC', b'BDC', b'EMC', b'BX', b'EX', b'n', b'X0',
              b'X1', b'DL', b'RL', b'K', b'k', b'RG', b'rg', b'G', b'g')


def _tokens(b):
    """PDF 内容流词法扫描 -> [(start, end, kind, text)]
    kind: 'str' ( (…) ) | 'hex' ( <…> ) | 'name' ( /X ) | 'num' | 'word' | 'delim'"""
    DELIM = b' \t\r\n/()<>[]{}%'
    out, i, n = [], 0, len(b)
    while i < n:
        c = b[i]
        if c == 0x25:                                        # % 注释
            j = b.find(b'\n', i)
            i = n if j < 0 else j + 1
            continue
        if c == 0x28:                                        # ( 字符串
            k = i + 1
            depth = 1
            while k < n and depth:
                if b[k] == 0x5c:
                    k += 2
                    continue
                if b[k] == 0x28:
                    depth += 1
                elif b[k] == 0x29:
                    depth -= 1
                k += 1
            out.append((i, k, 'str', b[i:k]))
            i = k
            continue
        if c == 0x3c:                                        # < 十六进制串
            j = b.find(b'>', i)
            j = n if j < 0 else j + 1
            out.append((i, j, 'hex', b[i:j]))
            i = j
            continue
        if c == 0x2f:                                        # / 名字
            j = i + 1
            while j < n and b[j] not in DELIM:
                j += 1
            out.append((i, j, 'name', b[i:j]))
            i = j
            continue
        if c == 0x5b:                                      # [ 数组开始
            out.append((i, i + 1, 'open', b'['))
            i += 1
            continue
        if c == 0x5d:                                      # ] 数组结束
            out.append((i, i + 1, 'close', b']'))
            i += 1
            continue
        if c in DELIM:
            i += 1
            continue
        j = i
        while j < n and b[j] not in DELIM:
            j += 1
        w = b[i:j]
        try:
            float(w)
            kind = 'num'
        except ValueError:
            kind = 'word'
        out.append((i, j, kind, w))
        i = j
    return out


OPERAND = ('num', 'str', 'hex', 'name', 'arr')


def scrub_text_ops(b):
    """只删 BT..ET 里"显示文字"的算子 (Tj / TJ / ' / ") 连同它的操作数;
    文本对象里的其它算子 (路径、颜色、gs、q/Q…) 一律保留。

    为什么不能整段删 BT..ET: Illustrator/cairo 导出的成品里, 文本对象内部偶尔夹着画图算子
    (Fig1_v3 页面级有 3 个绘图对象就在 BT..ET 里, 整段删掉 -> 图例细线消失, 逐像素差 32 px)。"""
    _t = _tokens(b)                      # 只扫一次 (之前写成 while k < len(_tokens(b)) 是 O(n^2) 重扫, 直接卡死)
    toks, k = [], 0
    while k < len(_t):
        st, en, kind, txt = _t[k]
        if kind == 'open':                                 # [ … ] 整体当一个操作数
            depth = 1
            k += 1
            while k < len(_t) and depth:
                if _t[k][2] == 'open':
                    depth += 1
                elif _t[k][2] == 'close':
                    depth -= 1
                k += 1
            toks.append((st, _t[k - 1][1], 'arr', b'[]'))
            continue
        toks.append((st, en, kind, txt))
        k += 1
    cuts, in_bt, op_start = [], False, None
    for st, en, kind, txt in toks:
        if kind == 'word':
            if txt == b'BT':
                in_bt, op_start = True, None
            elif txt == b'ET':
                in_bt, op_start = False, None
            elif in_bt and txt in SHOW_OPS:
                cuts.append((op_start if op_start is not None else st, en))
                op_start = None
            else:
                op_start = None                     # 普通算子: 操作数栈清空
        elif kind in OPERAND:
            if op_start is None:
                op_start = st
    if not cuts:
        return b, 0
    out, p = bytearray(), 0
    for s, e in sorted(cuts):
        if s < p:
            s = p
        if s >= e:
            continue
        out += b[p:s]
        p = e
    out += b[p:]
    return bytes(out), len(cuts)


def strip_text_blocks(b):
    """内容流 -> (删掉全部文字后的内容流, 删了几处显示算子)"""
    if b'BT' not in b:
        return b, 0
    return scrub_text_ops(b)


def strip_text_stream_v2(doc, page, xrefs=None):
    """同 strip_text_stream, 但用词法安全的 strip_text_blocks。xrefs=None -> 页面内容流。"""
    n = 0
    for xref in (page.get_contents() if xrefs is None else xrefs):
        cont = doc.xref_stream(xref)
        if cont is None:
            continue
        new, k = strip_text_blocks(cont)
        if k:
            doc.update_stream(xref, new)
            n += k
    return n


def gradient_rect(page, rect, pix, n=60, vertical=True, flip=False):
    """在 rect 内用 n 个细条画渐变。pix: fitz.Pixmap (1xH 竖条 或 Wx1 横条)。
    vertical=True -> 第 0 号像素画在顶部; flip=True 反向。"""
    sh = page.new_shape()
    W, H = rect.width, rect.height
    for i in range(n):
        k = i if not flip else (n - 1 - i)
        if pix.width == 1:
            c = pix.pixel(0, int(k * (pix.height - 1) / max(n - 1, 1)))
        else:
            c = pix.pixel(int(k * (pix.width - 1) / max(n - 1, 1)), 0)
        col = (c[0] / 255, c[1] / 255, c[2] / 255)
        if vertical:
            y0 = rect.y0 + H * i / n
            r = fitz.Rect(rect.x0, y0, rect.x1, rect.y0 + H * (i + 1) / n + 0.2)
        else:
            x0 = rect.x0 + W * i / n
            r = fitz.Rect(x0, rect.y0, rect.x0 + W * (i + 1) / n + 0.2, rect.y1)
        sh.draw_rect(r)
        sh.finish(color=None, fill=col, width=0)
    sh.commit()
