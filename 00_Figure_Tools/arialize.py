import os
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
arialize.py — 把成品 PDF 的文字层整体换成真 Arial (只动文字, 不碰任何图形) v3
================================================================================
为什么需要: 成品里有些文字不是我们的脚本画的, 是**贴进来的外部源 PDF** (R/ggplot 早年在
fontconfig alias 下产出的 Liberation Sans、base-14 Helvetica、matplotlib 默认 DejaVu)。
无法从数据重跑, 但字体必须统一成 Arial -> 在成品页上"删字重绘"。

为什么可以原位换 (2026-10-13 fontTools 实测):
    Arial vs Liberation Sans: hmtx 推进宽度逐字符相同 (总宽比值 1.0000);
                              hhea 1854/-434, OS/2 typo 1491/-431, win 1854/434 全同
    => 同一串文字从同一基线原点写出, 左右边界与上下伸出量都落在原位 (竖排也一样)。
    DejaVu Sans 比 Arial 宽 7~10% -> 按 bbox 中心补偿重锚 (报告会列出)。

两个坑 (都踩过, 已在本版本修掉):
 1) 不能用 apply_redactions: figkit.strip_text_stream 的注释/README_Fig3.md 记录它会破坏
    密集微矢量点 (UMAP 点云); 而且本项目 Fig1_v3/Fig2 上它直接抛
    "Cannot create unique resource name"。 -> 改用 "删 BT..ET 文本块" (含 Form XObject 递归)。
 2) 不能用 PyMuPDF 的 TextWriter/insert_text 写 Arial: 它给嵌入字体生成的 /ToUnicode 会把
    glyph **反查成 cmap 里的另一个码位** (实测 Arial: 空格 U+0020 -> U+00A0, · U+00B7 -> U+A78F),
    字形/位置全对, 但文本层取出的字符错 -> Illustrator 里编辑/搜索会拿到 nbsp 和怪字符。
    matplotlib 的 PDF 后端自己写 CMap+ToUnicode, 文本层是干净的 (FigS1 实测无 nbsp)。
    => 重绘一律走 matplotlib 透明叠加层 (use_matplotlib_arial 注册真 Arial, pdf.fonttype=42),
       再用 show_pdf_page 贴回原页 (与拼版用的同一套机制, 文字保持可编辑)。

前导空格: 原 PDF 里空格常是用 TJ 核距摆出来的, 推进与 Arial 的空格宽度不同 (实测差 1.5 pt)。
    => 用"第一个可见字符的原始 x"反锚: origin_new = inkx0 - Arial(前导空格串)。

v4 (2026-10-14) 两处必修:
 1) **隐藏文字**: get_text() 不看可见性 —— 贴外部源 PDF 带进来的旧文字层 (render mode 3 /
    被不透明图元盖住 / 被 clip 掉) 一样会被取出。若照原样重绘, 这些**本来看不见的字会变成
    看得见** (noc2 实测 5063 个隐藏字符, 含旧 A-F 面板字母)。故 full 模式先做"有字渲染 vs
    删字渲染"逐像素差分, 墨量为 0 的 span 不重绘 = 顺手把旧文字层清掉 (画面零变化,
    可用 --keep-hidden 关掉)。
 2) **未嵌入字体残留**: 删掉 BT..ET 后, 页面/XObject 的 /Resources /Font 里那些 base-14
    (Helvetica 等, 无 FontFile) 仍被声明, 期刊 preflight 照样报。故存盘前递归剪掉
    内容流里不再被 Tf 引用的字体名, 再 garbage=4 回收对象。

安全阀 (每次都跑): 逐页矢量绘图对象数必须相同; 字符数不减(隐藏层除外); 逐词位置残差;
                  只余 Arial; 无 nbsp; 文本往返一致; 画面像素逐点比对 (只允许字形边缘抖动)。

用法:
    python3 arialize.py Fig3_v10.pdf [FigS2.pdf ...]   # 就地改, 备份 *.pre_arial.pdf
    python3 arialize.py --list a.pdf                   # 只报告字体现状
    python3 arialize.py --restore Fig3_v10.pdf         # 从备份还原
"""
import io, math, os, re, sys, shutil, collections

sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use('Agg')
import fitz
import numpy as np
import figkit as K
from matplotlib.figure import Figure

FIG = os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures'))
CM = 72 / 2.54
K.use_matplotlib_arial()
IS_ARIAL = lambda n: n.split('+')[-1].lower().startswith('arial')
# PyMuPDF 自己嵌入时用的名字 (带 "Regular/Bold/Italic" 字样) —— 是我们贴的字, 且 ToUnicode 被写坏
# (空格变 U+00A0), 所以也算"要重绘的目标"; 期刊原图的名字是 ArialMT / Arial-Bold / Arial-ItalicMT
OURS_PYMUPDF = lambda n: re.match(r'^Arial (Regular|Bold|Italic|Bold Italic)$', n) is not None
IS_TARGET = lambda n: (not IS_ARIAL(n)) or OURS_PYMUPDF(n) or n == 'Arial'
BOLDISH = lambda n, fl: bool(fl & 16) or 'bold' in n.lower() or 'black' in n.lower()
ITALISH = lambda n, fl: bool(fl & 2) or 'italic' in n.lower() or 'oblique' in n.lower()


def fixtext(t):
    """文本层里的 ToUnicode 伪影: PyMuPDF 自己嵌 Arial 时写的反查表把空格反查成 U+00A0,
    把连字符反查成 U+00AD(软连字符)。Arial 里两者都有字形, 照原样重绘会把伪影带进新文件
    (实测 Fig1_v3 96 个 / Fig5_v2 35 个 nbsp), 回读时又变成 "非 Arial 痕迹"。"""
    return t.replace('\xa0', ' ').replace('\xad', '-')


def spans(page):
    out = []
    for b in page.get_text('rawdict')['blocks']:
        if b['type'] != 0:
            continue
        for l in b['lines']:
            d = tuple(l.get('dir') or (1.0, 0.0))
            for s in l['spans']:
                chars = [fixtext(c['c']) for c in s['chars']]
                t = ''.join(chars)
                if not t.strip():
                    continue
                # 第一个可见字符的原始墨位 x (用于前导空格反锚)
                i0 = next((i for i, c in enumerate(chars) if c.strip()), 0)
                i1 = max((i for i, c in enumerate(chars) if c.strip()), default=0)
                inkx0 = s['chars'][i0]['bbox'][0]
                inkx1 = s['chars'][i1]['bbox'][2]      # 末个可见字符右缘 (不受尾空格影响)
                out.append(dict(t=t, font=s['font'], size=s['size'], color=s['color'],
                                origin=fitz.Point(s['origin']), bbox=fitz.Rect(s['bbox']),
                                dir=d, flags=s.get('flags', 0), i0=i0, inkx0=inkx0, inkx1=inkx1,
                                wsum=sum(c['bbox'][2] - c['bbox'][0] for c in s['chars'])))
    return out


def text_ink_mask(doc, page, dpi=150):
    """"这个 span 到底看不看得见"的判据: 整页有字渲染 vs 删掉全部文字后渲染 的逐像素差。
    返回 (bool 数组, dpi)。

    坑 (2026-10-14 实测, 两次三番把判定搞翻): **同一页, 先 get_text() 再 get_pixmap() 与先
    get_pixmap() 再 get_text() 渲染结果不一样** (MuPDF 的页面缓存被文本设备预热后走了另一条
    路径), FigS3_v3 上文字墨量差 4% (108015 vs 112278 px)。所以这里两次渲染都从磁盘原件重开
    文档, 与调用方 doc 的"是否已被读过/改过"完全无关 -> 判据可复现。
    前提: 调用时该文件在磁盘上就是要看的样子 (arialize 里保险/隐藏层判定都在改动之前跑)。"""
    import numpy as np
    path = doc.name
    if path and os.path.exists(path):
        d1 = fitz.open(path)
        pm = d1[page.number].get_pixmap(dpi=dpi)
        d1.close()
        d2 = fitz.open(path)
        strip_everywhere(d2, d2[page.number])
        pm2 = d2[page.number].get_pixmap(dpi=dpi)
        d2.close()
    else:
        pm = page.get_pixmap(dpi=dpi)
        tmp = fitz.open('pdf', doc.tobytes())
        strip_everywhere(tmp, tmp[page.number])
        pm2 = tmp[page.number].get_pixmap(dpi=dpi)
        tmp.close()
    A_ = np.frombuffer(pm.samples, np.uint8).reshape(pm.height, pm.width, pm.n)[:, :, :3].astype(np.int16)
    B_ = np.frombuffer(pm2.samples, np.uint8).reshape(pm2.height, pm2.width, pm2.n)[:, :, :3].astype(np.int16)
    return (np.abs(A_ - B_).max(axis=2) > 10), dpi


def split_hidden(mask_ddpi, S):
    """(文字墨迹图, dpi) + S -> (可见 spans, 隐藏 spans); 隐藏 = 该 span 框内文字墨迹为 0"""
    mask, dd = mask_ddpi
    sc = dd / 72.0
    vis, hid = [], []
    for s in S:
        x0, y0, x1, y1 = [int(v) for v in (s['bbox'].x0 * sc, s['bbox'].y0 * sc,
                                           s['bbox'].x1 * sc, s['bbox'].y1 * sc)]
        sub = mask[max(0, y0):max(y0 + 1, y1), max(0, x0):max(x0 + 1, x1)]
        (hid if sub.size == 0 or sub.sum() == 0 else vis).append(s)
    return vis, hid


FONTNAME = re.compile(rb'/([A-Za-z0-9_.\-]+)\s+[0-9.]+\s+Tf')


def prune_fonts(doc):
    """递归删掉 /Resources /Font 里内容流已不再引用的字体名 (清 base-14 未嵌入残留)。"""
    def used_names(stream):
        if not stream:
            return set()
        return {m.decode() for m in FONTNAME.findall(stream)}

    def clean_obj(xref, is_page):
        try:
            xo = doc.xref_object(xref, compressed=True)
        except Exception:
            return 0
        if not xo or '/Font' not in xo:
            return 0
        m = re.search(r'/Font\s*<<(.*?)>>', xo, re.S)
        if not m:
            return 0
        entries = re.findall(r'/([A-Za-z0-9_.\-]+)\s+(\d+)\s+0\s+R', m.group(1))
        if not entries:
            return 0
        conts = set()
        if is_page:
            for cx in re.findall(r'/Contents\s+(\d+)\s+0\s+R', xo) + \
                    [x for grp in re.findall(r'/Contents\s+\[(.*?)\]', xo, re.S)
                     for x in re.findall(r'(\d+)\s+0\s+R', grp)]:
                conts.add(int(cx))
        else:
            conts.add(xref)
        used = set()
        for cx in conts:
            try:
                used |= used_names(doc.xref_stream(cx))
            except Exception:
                pass
        keep = [(nm, xr) for nm, xr in entries if nm in used]
        if len(keep) == len(entries):
            return 0
        newdict = ''.join('/%s %s 0 R ' % (nm, xr) for nm, xr in keep).strip()
        xo2 = xo[:m.start()] + '/Font <<%s>>' % newdict + xo[m.end():]
        if not keep:
            xo2 = xo2.replace('/Font <<>> ', '').replace('/Font <<>>', '')
        doc.update_object(xref, xo2)
        return len(entries) - len(keep)

    n = 0
    for pg in doc:
        n += clean_obj(pg.xref, True)
    for x in {xx[0] for pg in doc for xx in pg.get_xobjects()}:
        try:
            s = doc.xref_object(x, compressed=True)
        except Exception:
            continue
        if '/XObject' in s:
            for gx in re.findall(r'/XObject\s*<<(.*?)>>', s, re.S):
                for y in re.findall(r'(\d+)\s+0\s+R', gx):
                    n += clean_obj(int(y), False)
        n += clean_obj(x, False)
    return n


def _fig(W, H):
    fig = Figure(figsize=(W / 72.0, H / 72.0), dpi=72)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis('off')
    ax.set_facecolor('none'); fig.patch.set_alpha(0.0)
    return fig, ax


def _pdf_of(fig):
    buf = io.BytesIO(); fig.savefig(buf, format='pdf', transparent=True, pad_inches=0)
    buf.seek(0); return fitz.open('pdf', buf.read())


def fackey(s):
    f = face_of(s)
    return ('bi' if f['weight'] == 'bold' and f['style'] == 'italic' else
            'b' if f['weight'] == 'bold' else 'i' if f['style'] == 'italic' else 'r')



def face_of(s):
    return dict(weight='bold' if BOLDISH(s['font'], s['flags']) else 'normal',
                style='italic' if ITALISH(s['font'], s['flags']) else 'normal')


def rgb(v):
    return ((v >> 16 & 255) / 255, (v >> 8 & 255) / 255, (v & 255) / 255)


def page_level_only_spans(doc, page):
    """判断哪些文字是"贴进来的 XObject 里带的" (如期刊参考图): 在副本上只删页面级 BT..ET,
    剩下的就是 XObject 里的文字。返回 (页面级 spans, XObject 内 spans)。选择性模式据此
    做到**完全不碰参考图原文字**。"""
    import copy as _c
    tmp = fitz.open('pdf', doc.tobytes())
    allS = spans(tmp[page.number])
    K.strip_text_stream_v2(tmp, tmp[page.number])
    inx = spans(tmp[page.number])
    tmp.close()
    key = lambda s: (s['t'], round(s['inkx0'], 2), round(s['bbox'].y0, 2), round(s['size'], 2))
    left = collections.Counter(map(key, inx))
    ours = []
    for s in allS:
        k = key(s)
        if left.get(k, 0) > 0:
            left[k] -= 1
        else:
            ours.append(s)
    return ours, inx


def strip_everywhere(doc, page):
    """删掉页面**全部**文字: 页面内容流 + 递归所有 Form XObject (含流直接挂本体那种)。"""
    n = K.strip_text_stream_v2(doc, page)
    seen = set()

    def clean(xref):
        if xref in seen or xref <= 0:
            return 0
        seen.add(xref)
        c = 0
        try:
            st = doc.xref_stream(xref)
        except Exception:
            st = None
        if st is not None:
            new, _nk = K.strip_text_blocks(st)
            if new != st:
                doc.update_stream(xref, new)
                c += 1
        try:
            xo = doc.xref_object(xref, compressed=True)
        except Exception:
            return c
        cm = re.search(r'/Contents\s+(\d+)\s+0\s+R', xo)
        subs = re.findall(r'/Contents\s+\[(\d+)\s+0\s+R[^\]]*\]', xo)
        for cx in ([cm.group(1)] if cm else []) + subs:
            try:
                st = doc.xref_stream(int(cx))
            except Exception:
                continue
            if st is None:
                continue
            new, _nk = K.strip_text_blocks(st)
            if new != st:
                doc.update_stream(int(cx), new)
                c += 1
        mx = re.search(r'/XObject\s*<<(.*?)>>', xo, re.S)
        if mx:
            for gx in re.finditer(r'(\d+)\s+0\s+R', mx.group(1)):
                c += clean(int(gx.group(1)))
        return c

    for x in page.get_xobjects():
        n += clean(x[0])
    return n


def _overlay_pdf(S, anchors, W, H):
    """把 spans 按给定锚点画到一张与页面同尺寸的透明 matplotlib PDF (纯文字层, 便于回读)。
    rotation_mode='anchor' 是关键: 'default' 是绕"旋转前的盒角"转, 回读原点恒比锚点偏
    (-4.3, -1.2) pt @6pt 字, 竖排标签怎么迭代都对不上; 'anchor' 绕基线锚点转, 回读原点 ==
    锚点 (实测差 0.000 pt), 且 bbox 与原 PDF 完全重合 (MOL2: 98.490 vs 98.491)。"""
    fig, ax = _fig(W, H)
    for s, (x, y) in zip(S, anchors):
        vx, vy = s['dir']
        rot = -math.degrees(math.atan2(vy, vx))     # PyMuPDF 的 dir 在屏幕坐标(y 向下); mpl 逆时针为正
        ax.text(x, H - y, s['t'], fontsize=s['size'], color=rgb(s['color']),
                family='Arial', rotation=rot, rotation_mode='anchor',
                ha='left', va='baseline', zorder=3, **face_of(s))
    return _pdf_of(fig)


def fit_anchors(S, iters=5, tol=0.015, flat_only_guess=True):
    """闭环求落笔锚点: 画 -> 回读 -> 按"期望落笔点 - 回读原点"修正 -> 再画, 直到 <tol cm。
    返回 (anchors, worst)。期望落笔点已含前导空格反锚 (空格在源 PDF 里常是 TJ 核距, 宽度与
    Arial 不同, 用起笔墨位 x 反推原点才对得齐)。"""
    targets = [(s['origin'].x, s['origin'].y) for s in S]
    anchors = list(targets)
    for i, s in enumerate(S):
        if abs(s['dir'][1]) >= 0.01:
            continue
        face = K.FONT[fackey(s)]
        if s['i0']:
            anchors[i] = targets[i] = (s['inkx0'] - face.text_length(s['t'][:s['i0']], fontsize=s['size']),
                                       s['origin'].y)
        # 注: 不对 DejaVu->Arial 的宽度差做\"居中补偿\"。补偿会把整串字平移 (实测质心移 1.6mm),
        # 那就是动内容; 换字体的正确语义是\"原点不动, 从原点重新排版\"(与 Illustrator 一致)。
    norm = lambda t: ' '.join(t.split())
    worst = 0.0
    CLAMP, RAD = 6.0, 14.0
    for it in range(iters):
        src = _overlay_pdf(S, anchors, CUR['W'], CUR['H'])
        R = spans(src[0])
        byindex = (len(R) == len(S) and
                   all(norm(a['t']) == norm(b['t']) for a, b in zip(S, R)))
        used, new, worst = set(), [], 0.0
        for i, s in enumerate(S):
            t = norm(s['t'])
            if byindex:
                r = R[i]
            else:
                cand = [(abs(R[j]['origin'].x - targets[i][0]) + abs(R[j]['origin'].y - targets[i][1]), j)
                        for j in range(len(R))
                        if j not in used and (norm(R[j]['t']) == t or norm(R[j]['t']).startswith(t))]
                cand = [c for c in cand if c[0] <= RAD]
                if not cand:
                    new.append(anchors[i]); continue
                r = R[min(cand)[1]]; used.add(min(cand)[1])
            dx = targets[i][0] - r['origin'].x
            dy = targets[i][1] - r['origin'].y
            if abs(dx) > CLAMP or abs(dy) > CLAMP:
                new.append(anchors[i]); continue
            worst = max(worst, max(abs(dx), abs(dy)) / CM)
            new.append((anchors[i][0] + dx, anchors[i][1] + dy))
        anchors = new
        src.close()
        if worst < tol:
            break
    return anchors, worst


CUR = {'W': 0.0, 'H': 0.0}


def place_overlay(doc, page, S, anchors, W, H):
    """把 Arial 文字层贴到当前页 (返回新加的内容流 xref 列表, 供容器模式并回容器流)"""
    before = set(page.get_contents())
    src = _overlay_pdf(S, anchors, W, H)
    page.show_pdf_page(fitz.Rect(0, 0, W, H), src, 0)
    src.close()
    after = [x for x in page.get_contents()]
    return [x for x in after if x not in before]


DELIM = b' \t\r\n/()<>[]{}%'


def stream_state_at(b, off):
    """扫内容流 [0,off) 求: (q/Q 深度, 该处的 CTM, 是否有生效中的裁剪 W n)。
    插入重绘层时必须知道这三件事, 否则:
      * 深度>0 处插入 -> 继承外面的 cm, 整层文字被平移/缩放 (实测 noc2 挪了 3cm);
      * 有裁剪时插入 -> 被 clip 路径切掉一半;
      所以要么在深度 0 处插, 要么用逆矩阵补偿。"""
    import numpy as np
    depth = 0
    M = fitz.Matrix(1, 0, 0, 1, 0, 0)
    stack = []
    clip = 0
    nums = []
    i, n = 0, min(off, len(b))
    while i < n:
        c = b[i]
        if c == 0x25:
            j = b.find(b'\n', i); i = n if j < 0 else j + 1; continue
        if c == 0x28:
            i += 1; dep = 1
            while i < n and dep:
                if b[i] == 0x5c: i += 2; continue
                if b[i] == 0x28: dep += 1
                elif b[i] == 0x29: dep -= 1
                i += 1
            continue
        if c == 0x3c:
            j = b.find(b'>', i); i = n if j < 0 else j + 1; continue
        if c == 0x2f:
            i += 1
            while i < n and b[i] not in DELIM: i += 1
            continue
        if c in DELIM:
            i += 1; continue
        j = i
        while j < n and b[j] not in DELIM: j += 1
        w = b[i:j]
        if w in (b'q', b'Q'):
            if w == b'q':
                stack.append((M, clip))
            elif stack:
                M, clip = stack.pop()
            nums = []
        elif w == b'cm':
            if len(nums) >= 6:
                a, bb, c_, dd, e, f = nums[-6:]
                M = fitz.Matrix(a, bb, c_, dd, e, f) * M
            nums = []
        elif w in (b'W', b'W*', b'BI'):
            clip += 1 if w != b'BI' else 0
            nums = []
        elif w in (b'f', b'f*', b'F', b'B', b'B*', b'S', b's'):
            clip = max(0, clip - 1)
            nums = []
        elif w in (b'RG', b'rg', b'K', b'k', b'CS', b'scn', b'gs', b'Tf', b'Td', b'TD',
                   b'Tm', b'TJ', b'Tj', b"'", b'"', b'D1', b'D0', b'MP', b'DP', b'BMC', b'BDC', b'EMC'):
            nums = []
        else:
            try:
                nums.append(float(w))
            except ValueError:
                nums = []
        i = j
    return len(stack), M, clip


def insert_point_for(st, first_off, prefer_top=True):
    """返回 (插入偏移, 补偿矩阵字节串或 None)

    **[2026-10-15 修正] 默认插到流末尾, 不是插到第一个文本块处。**
    原来插在首文本块是为了"保住容器内的 z 序", 但那是错的: 容器里凡是画在文字**之后**的
    不透明图元 (matplotlib 的白底 patch、相邻色块) 都会把重绘层整个盖住 ——
    实测 Fig5_v2 面板 A 的竖排刻度 'COP_NFOL' 墨迹 661px -> 199px, 只剩尾巴可见,
    而文本层完全正确 (所以只查文本层的审计查不出来)。
    插到末尾 = 文字在**本容器内**永远最顶; 跨容器的遮挡顺序由各容器自身的调用顺序决定, 不受影响。
    副作用 (被本容器内图元压住的字会重新露出来) 由 audit 的"揭开"计数兜底, >0 即判失败。
    末尾若不干净 (q 未闭合 / 有活动裁剪) 才退回旧策略。"""
    if prefer_top:
        dE, ME, cE = stream_state_at(st, len(st))
        if dE == 0 and cE == 0:
            return len(st), None
    depth, M, clip = stream_state_at(st, first_off)
    if depth == 0 and clip == 0:
        return first_off, None                      # 顶层且无裁剪: 原位插, CTM 就是容器基空间
    if clip == 0:
        try:
            mi = ~M
            return first_off, (b'%g %g %g %g %g %g cm' % (mi.a, mi.b, mi.c, mi.d, mi.e, mi.f))
        except Exception:
            pass
    return len(st), None                            # 有裁剪 -> 只能追加到末尾 (那里深度 0)


def _gray(page, dpi=150):
    import numpy as np
    pm = page.get_pixmap(dpi=dpi)
    return np.frombuffer(pm.samples, np.uint8).reshape(pm.height, pm.width, pm.n)[:, :, :3].astype(np.int16)


def container_fonts(doc, page, cx):
    """容器流里 Tf 用到的字体资源名 -> BaseFont 名集合。
    解析不了 (资源字典是间接引用且拿不到 / 名字没映射 / BaseFont 缺失) 就返回 None = "不确定"。
    只用来做**便宜的预筛**: 一个容器里的字全是 Arial 时, 没必要跑 _own_spans
    (那玩意儿每容器要整页渲染两次 + 整文档序列化一次, Fig1_v3 有 300+ 容器, 直接跑到天荒地老)。"""
    st = doc.xref_stream(cx) or b''
    names = set(re.findall(rb'/([A-Za-z0-9._+-]+)\s+[\d.]+\s+Tf', st))
    if not names:
        return set()
    res = doc.xref_get_key(cx, 'Resources/Font')[1]
    if res in ('null', 'None'):
        res = doc.xref_get_key(page.xref, 'Resources/Font')[1]
    if not res.startswith('<<'):
        return None
    d = dict((n, int(x)) for n, x in re.findall(r'/([A-Za-z0-9._+-]+)\s+(\d+)\s+\d+\s+R', res))
    out = set()
    for n in names:
        xr = d.get(n.decode())
        if xr is None:
            return None
        bf = doc.xref_get_key(xr, 'BaseFont')[1]
        if bf in ('null', 'None'):
            return None
        out.add(bf.strip())
    return out


def form_xrefs(doc):
    """页面上所有含文字的 Form XObject xref (含嵌套的); 按 xref 升序 = 外层先处理"""
    out = []
    for x in range(1, doc.xref_length()):
        so = doc.xref_object(x, compressed=True) or ''
        if not so.startswith('<<') or '/Subtype/Form' not in so.replace(' ', ''):
            continue
        try:
            st = doc.xref_stream(x)
        except Exception:
            continue
        if st and K.text_block_ranges(st):
            out.append(x)
    return out


def _reload(page):
    """改过内容流后强制 MuPDF 重建这页的显示列表 (PyMuPDF 1.24 没有 Page.clean())。
    不刷新的话删字前后的两张渲染图完全一样, 可见性判定会全错。"""
    d = page.parent
    k, v = d.xref_get_key(page.xref, 'Contents')
    d.xref_set_key(page.xref, 'Contents', v)


def _own_spans(doc, page, cxref):
    """本容器自己(不含嵌套)的文字 span, 按当前页面坐标系给出; 同时返回删掉它们后的页面墨迹差。
    做法: 记下容器流原文 -> 删容器文字 -> 页面 span 表做多重集差 -> 复原容器流。"""
    st0 = doc.xref_stream(cxref) or b''
    allS = spans(page)
    stripped, nk = K.strip_text_blocks(st0)
    if not nk:
        return [], None, st0
    pixA = _gray(page)                              # 当前状态 (含本容器文字)
    doc.update_stream(cxref, stripped)
    tmp = fitz.open('pdf', doc.tobytes())           # 另开一份渲染: 同一 doc 改流后 MuPDF 会用缓存
    pixB = _gray(tmp[0])
    tmp.close()
    left = spans(page)
    key = lambda s: (s['t'], round(s['inkx0'], 2), round(s['bbox'].y0, 2), round(s['size'], 2))
    cnt = collections.Counter(map(key, left))
    own = []
    for t in allS:
        k = key(t)
        if cnt.get(k, 0) > 0:
            cnt[k] -= 1
        else:
            own.append(t)
    doc.update_stream(cxref, st0)                      # 复原, 真正重绘时再删一次
    return own, (np.abs(pixA - pixB) > 10), st0


def fix_container(doc, page, cxref, drop_hidden=True, iters=5, tol=0.015):
    """**原位**重画一个 Form XObject 自己的文字 (容器在页面里的 z 序、外部遮挡、内部遮挡全不变)。

    可见性一律在**整页坐标系**里量 (_own_spans): 只看容器自身渲染会误删"页面上其实看得见"的
    标签 (noc2 的 'CDS UTR'/'8.88' 就是这么丢的), 只看整页又清不掉贴进来的旧文字层
    (那层在自己容器里是画得完完整整的, 只是在页面上被后面的面板盖住 = 零墨迹)。
    画的时候再把同一批 span 换回容器局部坐标 (两边同一内容顺序, 按序号配对)。"""
    own_page, vmask, st0 = _own_spans(doc, page, cxref)
    if not own_page:
        return 0, 0, 0.0
    if vmask is not None:                     # 逐容器保险: 本容器有"看得见但提不出文本"的字就整块不碰
        sc = 150 / 72.0
        m = np.zeros(vmask.shape, bool)
        for t in own_page:
            x0, y0, x1, y1 = [int(v * sc) for v in (t['bbox'].x0 - 1, t['bbox'].y0 - 1,
                                                    t['bbox'].x1 + 1, t['bbox'].y1 + 1)]
            m[max(0, y0):max(y0 + 1, y1), max(0, x0):max(x0 + 1, x1)] = True
        tot, lost = int(vmask.sum()), int((vmask & ~m).sum())
        if tot and lost > 0.002 * tot:
            print('   xref%-5d 跳过: 本容器 %d/%d 文字墨 px 提不出文本 (Type3 无 ToUnicode), '
                  '重绘会静默丢字' % (cxref, lost, tot))
            return 0, 0, 0.0
    keep_idx = []
    if vmask is not None:
        sc = 150 / 72.0
        for i, t in enumerate(own_page):
            x0, y0, x1, y1 = [int(v * sc) for v in t['bbox']]
            sub = vmask[max(0, y0):max(y0 + 1, y1), max(0, x0):max(x0 + 1, x1)]
            if (not drop_hidden) or sub.sum() > 0:
                keep_idx.append(i)
    todo = [i for i in keep_idx if IS_TARGET(own_page[i]['font'])]
    if not todo:
        return 0, len(own_page) - len(keep_idx), 0.0
    rect_str = doc.xref_get_key(cxref, 'BBox')[1]
    if not rect_str:
        return 0, 0, 0.0
    rect = fitz.Rect([float(v) for v in rect_str.strip('[]').split()])
    if rect.is_empty or rect.width <= 0 or rect.height <= 0:
        return 0, 0, 0.0
    save = dict(contents=doc.xref_get_key(page.xref, 'Contents')[1],
                res=doc.xref_get_key(page.xref, 'Resources')[1],
                media=doc.xref_get_key(page.xref, 'MediaBox')[1],
                crop=doc.xref_get_key(page.xref, 'CropBox')[1])
    if save['media'] in (None, 'null'):
        r0 = page.rect
        save['media'] = '[%g %g %g %g]' % (r0.x0, r0.y0, r0.x1, r0.y1)
    kres, res = doc.xref_get_key(cxref, 'Resources')
    inline_res = (kres != 'null' and res and res.startswith('<<'))
    doc.xref_set_key(page.xref, 'Contents', '%d 0 R' % cxref)
    doc.xref_set_key(page.xref, 'Resources', res if (kres != 'null' and res) else 'null')
    page.set_mediabox(rect)
    worst = 0.0
    nred = 0
    try:
        allS = spans(page)
        K.strip_text_blocks(st0)                       # 触发一次编译期检查, 无副作用
        stripped, _nk = K.strip_text_blocks(doc.xref_stream(cxref) or b'')
        doc.update_stream(cxref, stripped)
        left = spans(page)
        key = lambda t: (t['t'], round(t['inkx0'], 2), round(t['bbox'].y0, 2), round(t['size'], 2))
        cnt = collections.Counter(map(key, left))
        own_local = []
        for t in allS:
            k = key(t)
            if cnt.get(k, 0) > 0:
                cnt[k] -= 1
            else:
                own_local.append(t)
        if len(own_local) != len(own_page):
            print('   xref%-5d 容器 span 数不一致 (页面 %d / 局部 %d) -> 复原并跳过' %
                  (cxref, len(own_page), len(own_local)))
            doc.update_stream(cxref, st0)            # 必须复原: 已经删过字了, 不写回就是丢内容
            return 0, 0, 0.0
        S = [own_local[i] for i in keep_idx]           # 页面坐标 -> 局部坐标, 同序配对
        W, H = rect.width, rect.height
        CUR['W'], CUR['H'] = W, H
        anchors, worst = fit_anchors(S, iters=iters, tol=tol)
        added = place_overlay(doc, page, S, anchors, W, H)
        extra = b'\n'.join(doc.xref_stream(x) or b'' for x in added)
        cur = doc.xref_stream(cxref) or b''
        rg0 = K.text_block_ranges(st0)
        cut, mat = insert_point_for(cur, rg0[0][0] if rg0 else len(cur))
        doc.update_stream(cxref, cur[:cut] + b'\nq\n' + (mat + b'\n' if mat else b'') +
                          extra + b'\nQ\n' + cur[cut:])
        if inline_res:
            k2, r2 = doc.xref_get_key(page.xref, 'Resources')
            if r2:
                doc.xref_set_key(cxref, 'Resources', r2)
        for x in added:
            try:
                doc.update_stream(x, b' ')
            except Exception:
                pass
        nred = len(todo)
    finally:
        doc.xref_set_key(page.xref, 'Contents', save['contents'])
        doc.xref_set_key(page.xref, 'Resources', save['res'])
        page.set_mediabox(fitz.Rect([float(v) for v in save['media'].strip('[]').split()]))
        if save['crop'] not in (None, 'null'):
            doc.xref_set_key(page.xref, 'CropBox', save['crop'])
    return nred, len(own_page) - len(keep_idx), worst


def redraw_page(doc, page, S, keep_xobject_text=False, iters=5, tol=0.015):
    """整页拍平重绘: 删掉页面级(或全部)文字 -> 在页面最上层重画。
    会把"被别的图元盖住的字"也揭开, 所以只适用于无遮挡的图; 有遮挡用 fix_container()。"""
    if keep_xobject_text:
        K.strip_text_stream_v2(doc, page)
    else:
        strip_everywhere(doc, page)
    W, H = page.rect.width, page.rect.height
    CUR['W'], CUR['H'] = W, H
    anchors, worst = fit_anchors(S, iters=iters, tol=tol)
    place_overlay(doc, page, S, anchors, W, H)
    return worst


def _graphics_only(doc, dpi=110):
    """删掉全部文字后的图形渲染 = "内容"的指纹 (字体改动不该碰它)"""
    import numpy as np
    tmp = fitz.open(doc.name) if (doc.name and os.path.exists(doc.name)) else fitz.open('pdf', doc.tobytes())
    strip_everywhere(tmp, tmp[0])
    pm = tmp[0].get_pixmap(dpi=dpi)
    tmp.close()
    return np.frombuffer(pm.samples, np.uint8).reshape(pm.height, pm.width, pm.n).copy()


def audit(before, after):
    """改字前后对比: 图形没被碰 / 字没被挪 / 没揭开遮挡 / 只剩 Arial / 文本层干净"""
    d0, d1 = fitz.open(before), fitz.open(after)
    r = {}
    g0, g1 = _graphics_only(d0), _graphics_only(d1)
    r['graphics_identical'] = (g0.shape == g1.shape and bool((g0 == g1).all()))
    r['drawings'] = ([len(p.get_drawings()) for p in d0], [len(p.get_drawings()) for p in d1])
    A_all = spans(d0[0])
    A, A_hid = split_hidden(text_ink_mask(d0, d0[0]), A_all)
    B = spans(d1[0])
    r['hidden_dropped'] = (len(A_hid), sum(len(x['t'].strip()) for x in A_hid))
    r['chars'] = (sum(len(x['t'].strip()) for x in A), sum(len(x['t'].strip()) for x in B))
    norm = lambda t: ' '.join(t.split())
    # "文字有没有消失" 按 (文本, 位置) 去重比, 不按实例数: 拼版图里同一段标签常有 3-5 份重叠
    # 拷贝, 只有最上面那份在页面上有墨, 其余在各自容器里是被盖住的 -> 原位重绘时正确丢弃,
    # 但按实例计数会被误判成"丢了字" (noc2 实测 95 条误报)。
    # 每个"改前可见"的文字实例, 改后附近 (6pt 内) 是否还有同一段文字 —— 允许被合并成更长的
    # 一段 (PyMuPDF 会把共线相邻的两段读成一段, 如 'CDS UTR'+'intron line' -> 'CDS UTR intron line')
    lost = []
    for a in A:
        ta = norm(a['t'])
        # 位置容差 6pt 只适用于"同一段字还在原位"; 合并读取时 b 比 a 长, 起点会左移
        # (实测 noc2: 'CDS UTR'@1.27 + 'intron line'@2.36 -> 读成一段 'CDS UTR intron line'@1.27),
        # 所以再给一条: 同一条基线上, a 的水平中心落在 b 框内 = 被并进更长的结果, 不算丢。
        cx = (a['bbox'].x0 + a['bbox'].x1) / 2.0
        ok = any(((abs(b['bbox'].x0 - a['bbox'].x0) <= 6 or
                   (b['bbox'].x0 - 1 <= cx <= b['bbox'].x1 + 1)) and
                  abs(b['bbox'].y0 - a['bbox'].y0) <= 6 and
                  (norm(b['t']) == ta or ta in norm(b['t']) or norm(b['t']) in ta)) for b in B)
        if not ok:
            lost.append((ta[:30], round(a['bbox'].x0 / CM, 2), round(a['bbox'].y0 / CM, 2)))
    r['words_lost'] = lost[:20]
    r['n_lost'] = len(lost)
    r['fonts_after'] = sorted({x['font'] for x in B})
    r['sizes_after'] = sorted({round(x['size'], 2) for x in B})
    r['nbsp'] = sum(x['t'].count('\xa0') for x in B) + sum(x['t'].count('\xad') for x in B)
    mA, _ = text_ink_mask(d0, d0[0])
    mB, _ = text_ink_mask(d1, d1[0])
    inter = int((mA & mB).sum()); union = int((mA | mB).sum())
    r['ink_iou'] = inter / union if union else 1.0
    r['ink_before'] = int(mA.sum()); r['ink_after'] = int(mB.sum())

    def _g(dd, dpi=150):
        pm = dd[0].get_pixmap(dpi=dpi)
        return np.frombuffer(pm.samples, np.uint8).reshape(pm.height, pm.width, pm.n)[:, :, :3].astype(np.int16)
    dif = np.abs(_g(d0) - _g(d1)).max(axis=2)
    r['pixdiff_pct'] = 100.0 * float((dif > 40).sum()) / float(dif.size)

    # 逐 span: 改前零墨改后有墨 = 揭开了被遮住的字 (改内容); 起笔/基线动了 = 字被挪了
    mA3, _ = text_ink_mask(d0, d0[0], dpi=300)
    mB3, _ = text_ink_mask(d1, d1[0], dpi=300)
    s3 = 300 / 72.0

    def stat(m, bb):
        x0, y0, x1, y1 = [int(v * s3) for v in bb]
        sub = m[max(0, y0):max(y0 + 1, y1), max(0, x0):max(x0 + 1, x1)]
        n = int(sub.sum())
        if not n:
            return 0, None, None, None
        ys, xs = np.nonzero(sub)
        rows = np.bincount(ys, minlength=sub.shape[0]); cols = np.bincount(xs, minlength=sub.shape[1])
        rr = np.nonzero(rows >= 3)[0]; cc = np.nonzero(cols >= 3)[0]
        if len(rr) == 0 or len(cc) == 0:
            return n, xs.min() / s3, ys.min() / s3, ys.max() / s3
        return (n, cc.min() / s3, rr.min() / s3, rr.max() / s3)

    exposed, moved = [], []
    usedB = set()
    for a in A:
        cand = [(abs(b['bbox'].x0 - a['bbox'].x0) + abs(b['bbox'].y0 - a['bbox'].y0), j)
                for j, b in enumerate(B) if b['t'] == a['t'] and j not in usedB]
        cand = [c for c in cand if c[0] <= 6.0]
        if not cand:
            continue
        j = min(cand)[1]; usedB.add(j); b = B[j]
        win = fitz.Rect(min(a['bbox'].x0, b['bbox'].x0) - 2, min(a['bbox'].y0, b['bbox'].y0) - 2,
                        max(a['bbox'].x1, b['bbox'].x1) + 2, max(a['bbox'].y1, b['bbox'].y1) + 2)
        na, la, ta, ba = stat(mA3, win)
        nb, lb, tb, bb_ = stat(mB3, win)
        if nb > 1.6 * na + 40:            # 改前被盖住 (全隐藏或部分隐藏), 改后墨量暴涨 = 揭开了
            exposed.append((a['t'][:16], round(a['bbox'].x0 / CM, 2), round(a['bbox'].y0 / CM, 2),
                            na, nb))
        # 只有"窗口内墨量确实变了 (>15%)"才判定为字被挪/被盖:
        # 换字体的正确语义是原点与文本层逐字段不动, 只重排字形 -> 同一 span 框里的墨量只会差几个百分点
        # (Libe->Arial 的 '1': 1584 vs 1583 px), 而真出事故时墨量是断崖式的
        # (Fig5_v2 竖排 'COP_NFOL' 被后画的白底盖住: 661 -> 199 px)。
        # 绝对阈值 60px@300dpi 是必须的: 单个 6pt 数字的窗口只有几十 px,
        # 光抗锯齿就能差出 20%, 只有"相对+绝对"双超阈才是真被盖/真被挪。
        elif na > 40 and nb > 40 and abs(nb - na) > max(0.15 * na, 60) and \
                None not in (la, lb, ta, tb, ba, bb_) and \
                (max(abs(lb - la), abs(bb_ - ba)) > 1.2 or abs(tb - ta) > 2.0):
            moved.append((a['t'][:16], round(a['bbox'].x0 / CM, 2), round(a['bbox'].y0 / CM, 2),
                          round(lb - la, 2), round(tb - ta, 2), round(bb_ - ba, 2)))
    r['exposed'] = exposed
    r['moved'] = moved
    d0.close(); d1.close()
    return r


def bakname(path):
    return (path[:-4] if path.lower().endswith('.pdf') else path) + '.pre_arial.pdf'


def arialize_inplace(path, drop_hidden=True):
    """**逐容器原位 Arial 化** (多图拼版/有遮挡的成品用这个, 别用整页拍平):
      1) 页面级文字 (拼版脚本自己写的字母/图例) -> 页面层删字重绘;
      2) 每个含文字的 Form XObject (贴进来的面板源 PDF) -> 在它自己的坐标系里删字重绘,
         画完并回该容器流。容器的调用位置、矩阵、被谁盖住 一概不动。
    返回备份路径 (无改动则 None)。"""
    doc = fitz.open(path)
    page = doc[0]
    nchg = ntot = nhid = 0
    resid = 0.0
    # --- 1) 页面级 ---
    ours, inx = page_level_only_spans(doc, page)
    S = [s for s in ours if IS_TARGET(s['font'])]
    if S:
        ntot += len(ours); nchg += len(S)
        CUR['W'], CUR['H'] = page.rect.width, page.rect.height
        anchors, worst = fit_anchors(ours)
        resid = max(resid, worst)
        K.strip_text_stream_v2(doc, page)
        place_overlay(doc, page, ours, anchors, page.rect.width, page.rect.height)
        print('   [页面级] %d 段重绘 (XObject 内 %d 段交给容器模式)' % (len(ours), len(inx)))
    # --- 2) 容器级 ---
    forms = form_xrefs(doc)
    done = 0
    nskip = 0
    for cx in forms:
        cf = container_fonts(doc, page, cx)
        if cf is not None and cf and not any(IS_TARGET(f) for f in cf):
            nskip += 1                      # 整容器都是 Arial, 不必碰
            continue
        try:
            n, h, w = fix_container(doc, page, cx, drop_hidden=drop_hidden)
        except Exception as e:
            print('   xref%-5d 跳过 (%s)' % (cx, e))
            continue
        if n:
            done += 1; nchg += n; resid = max(resid, w)
            print('   xref%-5d 原位重绘 %d 段%s, 残余 %.4f cm' % (cx, n,
                  ' (删零墨 %d 段)' % h if h else '', w))
        nhid += h
    if not nchg:
        doc.close()
        print('   无需改动 (文字层已是 Arial)')
        return None
    print('   共 %d 个容器需重绘 (预筛跳过 %d 个纯 Arial 容器); 重绘 %d 段, 删零墨 %d 段; 最大残余误差 %.4f cm'
          % (done, nskip, nchg, nhid, resid))
    npf = prune_fonts(doc)
    if npf:
        print('   剪掉 %d 个不再被引用的字体声明 (base-14 未嵌入残留)' % npf)
    bak = bakname(path)
    if not os.path.exists(bak):
        shutil.copyfile(path, bak)
    tmp = path + '.tmp'
    doc.save(tmp, garbage=4, deflate=True, clean=True)
    doc.close()
    os.replace(tmp, path)
    return bak


def _guard(mask_ddpi, S):
    """整页拍平重绘前的硬保险: 每一处文字墨迹都要落在某个"可提取" span 的框里。
    框外还有墨 = 存在看得见却提不出文本的字 (Type3 无 ToUnicode 等), 删字重绘会静默丢字。"""
    mask, dd = mask_ddpi
    sc = dd / 72.0
    m = np.zeros(mask.shape, bool)
    for s in S:
        x0, y0, x1, y1 = [int(v * sc) for v in (s['bbox'].x0 - 1, s['bbox'].y0 - 1,
                                                s['bbox'].x1 + 1, s['bbox'].y1 + 1)]
        m[max(0, y0):max(y0 + 1, y1), max(0, x0):max(x0 + 1, x1)] = True
    tot, lost = int(mask.sum()), int((mask & ~m).sum())
    if tot and lost > 0.002 * tot:
        print('   ❌ 保险拦下: %d/%d 文字墨像素 (%.1f%%) 不在任何可提取 span 框内'
              % (lost, tot, 100.0 * lost / tot))
        print('      = 有提不出文本的可见字 (Type3 无 ToUnicode 等), 重绘会静默丢字; 未改动文件')
        return False
    print('   保险: 文字墨 %d px 全部由可提取 span 覆盖 ✅' % tot)
    return True


def arialize(path, mode='auto', drop_hidden=True):
    """mode: 'full' = 整页重绘 (含贴进来的我们自己面板里的字);
             'selective' = 只重绘页面级文字 (参考图/XObject 里的原文字完全不动);
             'auto' = 页面里有"已经是 Arial 的 XObject 文字"就用 selective, 否则 full
       drop_hidden: full 模式下不重绘"零墨迹" span = 顺手删掉外部源 PDF 带进来的隐藏旧文字层"""
    doc = fitz.open(path)
    ntot = nchg = 0
    nhidden = hchars = 0
    resid = 0.0
    for page in doc:
        S = spans(page)
        keep = False
        if mode in ('selective', 'auto'):
            ours, inx = page_level_only_spans(doc, page)
            ref_in_xobj = [s for s in inx if IS_ARIAL(s['font']) and not OURS_PYMUPDF(s['font'])]
            if mode == 'selective' or (mode == 'auto' and ref_in_xobj):
                keep = True
                S = ours
                print('   [selective] 页面级 %d 段重绘; XObject 内 %d 段 (其中参考 Arial %d 段) 原样保留'
                      % (len(S), len(inx), len(ref_in_xobj)))
        md = None
        if not keep:                                  # 整页拍平: 先验覆盖, 再摘隐藏层
            md = text_ink_mask(doc, page)
            if not _guard(md, S):
                doc.close(); return None
            if drop_hidden:
                S, H_ = split_hidden(md, S)
                if H_:
                    nhidden += len(H_)
                    hchars += sum(len(s['t'].strip()) for s in H_)
                    print('   [隐藏层] %d 段 / %d 字符零墨迹 (如 %s) -> 不重绘 = 删除'
                          % (len(H_), hchars, ', '.join(sorted({repr(s['t'].strip()[:10]) for s in H_[:5]}))))
        todo = [s for s in S if IS_TARGET(s['font'])]
        ntot += len(S); nchg += len(todo)
        if not todo:
            continue
        if keep:
            S2 = S
        resid = max(resid, redraw_page(doc, page, S, keep_xobject_text=keep))
    if not nchg:
        doc.close()
        print('   无需改动 (文字层已是 Arial)')
        return None
    print('   全页文字 %d 段, 需重绘 %d 段; 闭环收敛后残余锚点误差 %.4f cm' % (ntot, nchg, resid))
    npf = prune_fonts(doc)
    if npf:
        print('   剪掉 %d 个不再被引用的字体声明 (base-14 未嵌入残留)' % npf)
    bak = bakname(path)
    if not os.path.exists(bak):
        shutil.copyfile(path, bak)
    tmp = path + '.tmp'
    doc.save(tmp, garbage=4, deflate=True, clean=True)
    doc.close()
    os.replace(tmp, path)
    return bak


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    if not args or '-h' in sys.argv or '--help' in sys.argv:
        print(__doc__); sys.exit(0)
    allok = True
    for p in args:
        p = p if os.path.exists(p) else os.path.join(FIG, p)
        print('===', os.path.basename(p))
        if '--restore' in sys.argv:
            shutil.copyfile(bakname(p), p)
            print('   已从备份还原'); continue
        if '--list' in sys.argv:
            c = collections.Counter(s['font'] for s in spans(fitz.open(p)[0]))
            print('   现状:', dict(c)); continue
        bak = (arialize_inplace(p, drop_hidden=('--keep-hidden' not in sys.argv))
               if '--inplace' in sys.argv else
               arialize(p, mode=('selective' if '--selective' in sys.argv else
                                 'full' if '--full' in sys.argv else 'auto'),
                        drop_hidden=('--keep-hidden' not in sys.argv)))
        if bak is None:
            continue
        r = audit(bak, p)
        same = r['drawings'][0] == r['drawings'][1]
        bad = [f for f in r['fonts_after'] if not IS_ARIAL(f)]
        ok = same and not bad and not r['nbsp'] \
            and r['graphics_identical'] and not r['exposed'] and not r['moved'] \
            and r['pixdiff_pct'] < 1.5 and not r['n_lost']
        print('   校验: 绘图对象 %s %s | 可见字符 %d->%d (删隐藏层 %d 段/%d 字) | 丢词 %d | nbsp %d'
              % (r['drawings'][0], '相同✅' if same else '变了❌', r['chars'][0], r['chars'][1],
                 r['hidden_dropped'][0], r['hidden_dropped'][1], r['n_lost'], r['nbsp']))
        print('   内容指纹: 删字后图形渲染 %s' % ('与改字前逐像素相同✅ (只动了文字层)' if r['graphics_identical']
                                            else '变了❌ -> 碰到了图形, 必须 --restore'))
        print('   逐 span 墨迹边缘: 被揭开的隐藏字 %d 处 | 起笔/基线漂移>1.2pt %d 处 %s'
              % (len(r['exposed']), len(r['moved']),
                 '✅ 画面内容未变' if not r['exposed'] and not r['moved'] else '❌'))
        for e in r['exposed'][:4]:
            print('        揭开: %r @%.2f,%.2f 墨 %d -> %d px' % e)
        for m in r['moved'][:4]:
            print('        漂移: %r @%.2f,%.2f  d左=%.2f d上=%.2f d下=%.2f pt' % m)
        print('   遮蔽检查: %s' % ('无文字被揭开 ✅' if not r['exposed'] else
              '❌ %d 处原本被遮住的字重绘后露出来了' % len(r['exposed'])))
        print('   字体: %s\n   字号: %s  ->  %s' % (', '.join(r['fonts_after']), r['sizes_after'],
              'PASS ✅ 全 Arial・原位・文本层干净' if ok else
              'FAIL ❌ 残留 %s / 丢字 %s' % (bad, r['words_lost'][:3])))
        allok &= ok
    print('总判定:', 'PASS' if allok else '有需要复核项')
