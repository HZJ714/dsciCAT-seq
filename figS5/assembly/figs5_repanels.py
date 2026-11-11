#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
figs5_repanels.py — 重出 FigS5 的 b/e/f 三块密集面板 (槽位真实尺寸, 6-8pt, Arial)
================================================================================
b  RSS top10 网格 (9.0 x 10.7 cm): 从 rss_scores_filtered_96.csv 复刻官方
   plot_rss (2 列 x 3 行, 蓝色 eRSS 曲线, top10 红点 + adjustText 标注 TF 名);
   原图 25.8x30.8cm 缩 0.35 倍 -> 文字 3.2pt, 必须按槽位重出
e  Nrxn3 locus ATAC track (18.0 x 10.3 cm): 复用 grn_mef2c_nrxn3_track.py,
   figsize 47x27cm -> 槽位, 字号 clamp 6-8pt (脚本内自动改写)
f  GRN 网络 x3 (5.5 x 5.5 cm): 从 cytoscape CSV 重建 (networkx spring layout,
   seed 固定); 节点 = 全部, 标签 = TF + marker + degree>=3 (mol2 21/58,
   mol6 按 marker degree 取前 20, mol3_6 31/56), 6pt + adjustText
   —— 原 49/56/58 节点全标注在 5.5cm 内物理不可行
输出: <work>/figures/replot/S5{b,e,f}_*.pdf
用法: python figs5_repanels.py          (e 约 2-5 min, 读 12 个 h5)
"""
import os, re, glob, pickle, sys
import numpy as np
# LEGACY：作者机上那棵分析树的根（本文件只对那棵树里的脚本做文本替换，搜索串必须逐字一致）
LEGACY = os.environ.get('DSCI_LEGACY_ROOT', '')
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(LEGACY, 'pylibs')))
import fitz
CM = 72 / 2.54
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

# ---------- 字体: Liberation Sans (Arial 度量), TrueType 嵌入 ----------
import sys, os as _os
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from figkit import use_matplotlib_arial as _arial
_FAM = _arial()          # 真 Arial (无则 Liberation Sans) + pdf.fonttype=42
plt.rcParams.update({
    'axes.linewidth': 0.4, 'xtick.major.width': 0.4, 'ytick.major.width': 0.4,
    'xtick.major.size': 1.6, 'ytick.major.size': 1.6,
})
assert font_manager.findfont('Liberation Sans').endswith('.ttf')
OUT = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(LEGACY, 'figures')), 'replot')
os.makedirs(OUT, exist_ok=True)
PSC = os.environ.get('DSCI_SCENIC_ROOT', os.path.join(LEGACY, 'pi/plots-scenic+'))
BLUE, RED, GREY = '#3A6DB5', '#C03028', '#9A9A9A'


# Liberation Sans / Arial 的 em 比例 (与 figkit 同一套, 用于算真实墨迹盒)
CAP, XH, DESC = 0.716, 0.519, 0.212
UPCH = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789bdfhklt')
DNCH = set('gjpqy_')


def resolve_ink(fig, ax, texts, pts=(), dot_ms=2.4, gap=2.5, iters=1500,
                margin_pt=4.0, verbose=True):
    """按**真实墨迹盒** (display px) 确定性消解标签重叠: 文字-文字推开 + 文字-红点避开。
    每轮**先夹紧再判碰** -> 退出时既在轴内又无重叠 (旧版用 matplotlib 布局盒 +
    数据单位间隙, 与成品墨迹不守恒, 9.0cm 版残留 16 对、17.6cm 版残留 17 对重叠)。
    返回残余重叠对数 (adjust_text 会改 ha/va, 故锚点偏移逐文字预存)。
    gap 取 2.5px(≈0.06cm) 而不是 0: Agg 的文本行盒与 PDF 后端真正落笔的行盒在竖直方向
    差 ~1.3pt (va=center/top 时更大), 零间隙消解后成品仍会残留 2mm² 级别的斜角擦碰,
    留出这个余量后实测 0 重叠。"""
    ren = fig.canvas.get_renderer()
    dpi = fig.dpi

    def metrics(t):
        """(kx, ky, dn): 锚点相对墨迹盒左下角的偏移 (px) + 下伸部高度 (px).
        adjust_text 会把 ha/va 改成 center 等 (optimally_align_text), 所以**不能**假设
        锚点=基线左端; 基线一律从渲染器拿: base = e.y0 + descent."""
        e = t.get_window_extent(renderer=ren)
        _w, _h, d = ren.get_text_width_height_descent(t.get_text(),
                                                      t.get_fontproperties(), False)
        base = e.y0 + d
        ax_, ay_ = ax.transData.transform(t.get_position())
        s_ = t.get_text()
        dn = DESC * t.get_fontsize() * dpi / 72.0 if any(c in DNCH for c in s_) else 0.0
        return ax_ - e.x0, ay_ - base, dn

    m0 = [metrics(t) for t in texts]

    def one(t, i):
        e = t.get_window_extent(renderer=ren)
        _w, _h, d = ren.get_text_width_height_descent(t.get_text(),
                                                      t.get_fontproperties(), False)
        base = e.y0 + d
        s_ = t.get_text()
        em = t.get_fontsize() * dpi / 72.0
        top = (CAP if any(c in UPCH for c in s_) else XH) * em
        dn = DESC * em if any(c in DNCH for c in s_) else 0.0
        return [e.x0, base - dn, e.x1, base + top, dn, m0[i][0], m0[i][1]]

    def allb():
        return [one(t, i) for i, t in enumerate(texts)]

    r = dot_ms * dpi / 144.0 + gap
    dots = []
    for px_, py_ in pts:
        x, y = ax.transData.transform((px_, py_))
        dots.append([x - r, y - r, x + r, y + r, 0.0])
    ab = ax.get_window_extent()
    bx0, by0, bx1, by1 = ab.x0 - margin_pt, ab.y0 - margin_pt, ab.x1 + margin_pt, ab.y1 + margin_pt

    def clamp(bs):
        for b in bs:
            w, h = b[2] - b[0], b[3] - b[1]
            if b[0] < bx0:
                b[0], b[2] = bx0, bx0 + w
            if b[2] > bx1:
                b[2], b[0] = bx1, bx1 - w
            if b[1] < by0:
                b[1], b[3] = by0, by0 + h
            if b[3] > by1:
                b[3], b[1] = by1, by1 - h

    inv = ax.transData.inverted()
    left, it = -1, 0
    for it in range(iters):
        bs = allb()
        clamp(bs)
        pool = bs + dots
        nlab = len(bs)
        bad = []
        for i in range(len(pool)):
            for j in range(i + 1, len(pool)):
                a, c = pool[i], pool[j]
                ox = min(a[2], c[2]) - max(a[0], c[0])
                oy = min(a[3], c[3]) - max(a[1], c[1])
                if ox > 0 and oy > 0:
                    bad.append((i, j, ox, oy))
        if not bad:
            left = 0
            break
        for i, j, ox, oy in bad:
            a, c = pool[i], pool[j]
            ai, ci = i < nlab, j < nlab
            if not ai and not ci:
                continue
            d = max(min(ox, oy) * 0.6, 0.4)
            if ox <= oy:
                sg = 1 if (a[0] + a[2]) >= (c[0] + c[2]) else -1
                if ai and ci:
                    a[0] += sg * d / 2; a[2] += sg * d / 2
                    c[0] -= sg * d / 2; c[2] -= sg * d / 2
                elif ai:
                    a[0] += sg * d; a[2] += sg * d
                else:
                    c[0] -= sg * d; c[2] -= sg * d
            else:
                sg = 1 if (a[1] + a[3]) >= (c[1] + c[3]) else -1
                if ai and ci:
                    a[1] += sg * d / 2; a[3] += sg * d / 2
                    c[1] -= sg * d / 2; c[3] -= sg * d / 2
                elif ai:
                    a[1] += sg * d; a[3] += sg * d
                else:
                    c[1] -= sg * d; c[3] -= sg * d
        for t, b in zip(texts, bs):
            x, y = inv.transform((b[0] + b[5], b[1] + b[4] + b[6]))
            t.set_position((x, y))
    if verbose:
        print(f'    resolve_ink[{ax.get_title()}]: {len(texts)} 标签 {it} 轮, 残余 {left}')
    return left



# ================================================================
# b: RSS top10 网格  (默认 9.0 x 10.7 cm; FigS5_noc2 用 17.6 x 11.2 cm)
# ================================================================
def panel_b(w=9.0, h=10.7, out=None):
    """RSS top10 网格 (复刻官方 plot_rss): 2 列 x 3 行, 蓝曲线 + top10 红点 + 红标签。
    标签位置: adjust_text 初排 -> resolve_ink 按真实墨迹盒消解 (含红点避让/轴内夹紧)。
    导线在**第二遍**画: 第一遍之后还有 title/tick/label 等设置会改变几何, 先画完导线
    再改几何会让导线与标签错位。"""
    out = out or f'{OUT}/S5b_rss_grid.pdf'
    rss = pd.read_csv(f'{PSC}/rss_scores_filtered_96.csv', index_col=0)
    order = ['COP_NFOL', 'MFOL', 'MOL2', 'MOL3', 'MOL3_6', 'MOL6']   # 原图顺序
    fig, axes = plt.subplots(3, 2, figsize=(w / 2.54, h / 2.54))
    # 左边距: 旋转轴标题(两行 6.5pt ≈0.5cm) + y 刻度(0.4cm) + 缝 ≈ 1.30cm。
    # **必须在消解之前**定好版面: 之后再 subplots_adjust 会改变 axes 盒, 标签按数据坐标
    # 存的位移会整体错位, 之前的夹紧/消解全部作废。
    fig.subplots_adjust(left=min(0.13, 1.30 / w), right=0.985, top=0.955, bottom=0.065,
                        wspace=0.10, hspace=0.42)
    from adjustText import adjust_text
    lab = {}
    for ax, ct in zip(axes.flat, order):
        s = rss.loc[ct].sort_values(ascending=False)
        x = np.arange(len(s))
        ax.plot(x, s.values, '-', color=BLUE, lw=0.8, zorder=1)
        ax.plot(x, s.values, 'o', color=BLUE, ms=1.7, zorder=2)
        top = s.nlargest(10)
        xi = [s.index.get_loc(n) for n in top.index]
        ax.plot(xi, top.values, 'o', color=RED, ms=2.4, zorder=3)
        ax.set_xlim(-2, 97)                      # 先固定坐标范围
        texts = [ax.text(xx, yy, n.split('_+')[0].split('_(')[0]
                         .replace('_extended', ' ext').replace('_direct', ' dir'),
                         fontsize=6, color=RED, zorder=4)
                 for xx, yy, n in zip(xi, top.values, top.index)]
        adjust_text(texts, x=xi, y=list(top.values), ax=ax,
                    expand_text=(1.1, 1.25), force_text=(0.15, 0.25),
                    align=False, autoalign='')
        # adjust_text 会把 ha/va 改成 center; 统一显式设成 center/center 并把锚点设为盒中心,
        # 使**所有标签共用同一个 Agg->PDF 竖直系统差** (~1.7pt), 相对几何因此是精确的。
        fig.canvas.draw()
        ren0 = fig.canvas.get_renderer()
        inv0 = ax.transData.inverted()
        for t in texts:
            e = t.get_window_extent(renderer=ren0)
            t.set_ha('center'); t.set_va('center')
            t.set_position(inv0.transform(((e.x0 + e.x1) / 2, (e.y0 + e.y1) / 2)))
        fig.canvas.draw()
        resolve_ink(fig, ax, texts, pts=list(zip(xi, top.values)), dot_ms=2.4)
        ax.autoscale(False)                      # 冻结 limits: 之后加导线/标题不再改几何
        lab[ct] = (ax, texts, xi, list(top.values))
        ax.set_title(ct, fontsize=7)
        ax.tick_params(labelsize=6)
        col = 0 if ax in axes[:, 0] else 1
        if col == 0:
            ax.set_ylabel('eRegulon specificity\nscore (eRSS)', fontsize=6.5)
            ax.yaxis.set_label_coords(-0.105, 0.5)   # 让开 y 刻度 (默认 labelpad 会与刻度擦碰)
        else:
            ax.set_yticklabels([])
        if ax in axes[2, :]:
            ax.set_xlabel('eRegulon rank', fontsize=6.5)
        else:
            ax.set_xticklabels([])
    # ---- 第二遍: 存盘 -> 在**成品 PDF 的墨迹空间**里复核消解 + 画导线 ----
    tmp = f'{OUT}/.S5b_raw.pdf'
    fig.savefig(tmp)
    plt.close(fig)
    postfix_b(tmp, out, verbose=True)
    return out


def postfix_b(src, out, verbose=True):
    """面板 b 的成品级消解: Agg 行盒与 PDF 落笔差 1.3-1.7pt, matplotlib 里"0 重叠"
    存出来仍有 2mm² 级擦碰; 这里在 PDF 自己的墨迹坐标系里把红色标签推开 (整行当一个
    单元), 再按原基线精确重写文字层, 并补画 标签->红点 的导线。"""
    import pdfix, fitz
    from collections import defaultdict
    doc = fitz.open(src)
    page = doc[0]
    rects, dots = [], []
    for d in page.get_drawings():
        r = d['rect']
        if d['type'] in ('s', 'fs') and 5 < r.width / CM < 12 and 1 < r.height / CM < 4:
            rects.append((r.x0 / CM, r.y0 / CM, r.x1 / CM, r.y1 / CM))
        elif (d['type'] in ('f', 'fs') and r.width / CM < 0.2 and r.height / CM < 0.2
              and d.get('fill') and max(abs(a - b) for a, b in
                                        zip(d['fill'], (0.753, 0.188, 0.157))) < 0.02):
            dots.append(((r.x0 + r.x1) / 2 / CM, (r.y0 + r.y1) / 2 / CM))
    ws = pdfix.words(page)
    doc.close()
    RED = (0.753, 0.188, 0.157)
    us = pdfix.units_of(ws)

    def is_red(u):
        return any(max(abs(a - b) for a, b in zip(ws[i]['c'], RED)) < 0.02
                   and ws[i]['z'] < 6.2 for i in u['idx'])

    def rect_of(u):
        cx = (u['box'][0] + u['box'][2]) / 2
        cy = (u['box'][1] + u['box'][3]) / 2
        for r in rects:
            if r[0] - 0.35 <= cx <= r[2] + 0.35 and r[1] - 0.35 <= cy <= r[3] + 0.35:
                return r
        return None

    groups = defaultdict(list)
    for k, u in enumerate(us):
        r = rect_of(u)
        if r and is_red(u):
            groups[r].append(k)

    def soft_bounds(r):
        """上沿只留 0.15cm (再上面是本图标题); 下沿: 非末行可以溢到行间距里
        (1-2 行没有 x 刻度, 间距 ah*hspace=1.10cm, 减去下一图标题 0.25cm)。"""
        ymax = max(x[3] for x in rects)
        down = 0.15 if r[3] > ymax - 0.5 else 0.85
        return (r[0] - 0.1, r[1] - 0.15, r[2] + 0.1, r[3] + down)

    tot = 0
    for r, mv in groups.items():
        obs = [(x - 0.05, y - 0.05, x + 0.05, y + 0.05) for x, y in dots
               if r[0] <= x <= r[2] and r[1] <= y <= r[3]]
        left = pdfix.resolve_units(us, mv, bounds=soft_bounds(r), obstacles=obs,
                                   gap=0.012, iters=3000)
        if left:                       # 没收敛 -> 放宽边界再来一次
            left = pdfix.resolve_units(us, mv, bounds=(r[0] - 0.6, r[1] - 0.4,
                                                       r[2] + 0.6, r[3] + 1.3),
                                       obstacles=obs, gap=0.012, iters=4000)
        if left:                       # 还不行 -> 去掉软边界, 只保证不压别的东西
            left = pdfix.resolve_units(us, mv, bounds=None, obstacles=obs,
                                       gap=0.012, iters=6000)
        if left and verbose:
            print(f'    [b] rect({r[0]:.1f},{r[1]:.1f}) 仍未收敛: {len(mv)} 标签')
        tot += left
    pdfix.apply_units(ws, us)
    leaders = []
    for r, mv in groups.items():
        cand = [(x, y) for x, y in dots
                if r[0] - 0.5 <= x <= r[2] + 0.5 and r[1] - 0.5 <= y <= r[3] + 0.5]
        for k in mv:
            u = us[k]
            cx = (u['box'][0] + u['box'][2]) / 2
            cy = (u['box'][1] + u['box'][3]) / 2
            if not cand:
                continue
            ax_, ay_ = min(cand, key=lambda p: (p[0] - cx) ** 2 + (p[1] - cy) ** 2)
            leaders.append((ax_, ay_, u['idx'][0]))
    pdfix.rewrite(src, out, ws, leaders=leaders, lcolor=(0.753, 0.188, 0.157),
                  lw=0.25, lthresh=0.06)
    if verbose:
        n, bad = pdfix.check(out)
        print(f'  postfix_b: {len(us)} 行单元/{len(ws)} 词, 可动 {sum(len(v) for v in groups.values())} '
              f'标签, 导线 {len(leaders)} 条, 残余重叠 {tot} | 成品 {n} 词, 自检重叠 {len(bad)}')
        for a, b, ar in bad[:6]:
            print(f'      残留: {a["t"]!r} <> {b["t"]!r} {ar*100:.2f}mm2')
    return out


# ================================================================
# e: Nrxn3 locus track  18.0 x 10.3 cm  (改写字号后重跑原脚本)
# ================================================================
def panel_e():
    src = open(os.path.join(os.environ.get('DSCI_SCENIC_SCRIPTS', os.path.join(LEGACY, 'pi/scripts-scenic+')),
                  'grn_mef2c_nrxn3_track.py')).read()
    src = src.replace('figsize=(18.5, 11)', 'figsize=(18.0/2.54, 10.3/2.54)')
    # 字号 clamp 到 6-8 (原 5.5-10.5, 在 0.38 倍画布上只剩 2-4pt)
    src = re.sub(r'fontsize=([\d.]+)',
                 lambda m: f'fontsize={min(max(float(m.group(1)), 6.0), 8.0):.1f}',
                 src)
    src = src.replace('labelpad=28', 'labelpad=12')

    # ---- 布局修正 (18cm 槽位上的碰撞) ----
    # 标题: 原两行 >18cm 且 4 行 8pt 越界 -> 3 行 6.5pt (正则替换, 不依赖原文精确字符)
    _TITLE_NEW = """    axes[0].set_title(
        f'Nrxn3 locus (chr12:{WIN_S:,}-{WIN_E:,})  |  ATAC insertions per cell state (RPM-normalized, 200bp bins)\\n'
        'Mef2c eRegulon → Nrxn3 enhancers: 4 regions (2000hvr)  ·  MEF2 motif metacluster_179.4 (JASPAR MA0052.4), NES/hits shown\\n'
        'Right: Nrxn3 expression ('
        + ('boxplot + jitter, ' if style == 'box' else '') + 'linearized 2^X−1) per state — MOL3-specific '
        '(median ≈7.0 vs ≈3.1–3.8, ≈2×; MWU p<1e-22)',
        fontsize=6.5, loc='left', pad=16)"""
    src, nt = re.subn(
        r"    axes\[0\]\.set_title\(.*?loc='left'\)",
        lambda m: _TITLE_NEW, src, flags=re.S)
    assert nt == 1, 'title 替换未命中'

    reps = [
        # 上边距: 标题 3 行需要 ~12mm
        ('top=0.95, bottom=0.065', 'top=0.86, bottom=0.065'),
        # 边距: 无 tight 输出时所有元素须在画布内
        ("left=0.05,",
         "left=0.098,"),
        ("right=0.985, top=0.86, bottom=0.065",
         "right=0.99, top=0.84, bottom=0.075"),
        # violin 主体裁剪到轴内 (MOL3 长尾越界)
        ("body.set_linewidth(0.6)",
         "body.set_linewidth(0.6); body.set_clip_path(ax_vio.patch)"),
        # 保存: 不用 tight -> 画布恰为 figsize (18x10.3), 字号 1:1
        ("fig.savefig(out_pdf, bbox_inches='tight')", "fig.savefig(out_pdf)"),
        # 全局刻度字号 (原默认 ~10pt)
        ("fig = plt.figure(figsize=(18.0/2.54, 10.3/2.54))",
         "fig = plt.figure(figsize=(18.0/2.54, 10.3/2.54))\n"
         "    plt.rcParams['xtick.labelsize'] = 7\n"
         "    plt.rcParams['ytick.labelsize'] = 7\n"
         "    plt.rcParams['axes.labelsize'] = 7\n"
         "    plt.rcParams['axes.titlesize'] = 7\n"
         "    plt.rcParams['font.family'] = 'sans-serif'\n"
         "    plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'DejaVu Sans']\n"
         "    plt.rcParams['pdf.fonttype'] = 42\n"
         "    plt.rcParams['xtick.major.width'] = 0.4\n"
         "    plt.rcParams['ytick.major.width'] = 0.4\n"
         "    plt.rcParams['axes.linewidth'] = 0.4"),
        # MEF2 motif 注释: 0.66 -> 0.04 (原压在转录本行上)
        ("ax.text(WIN_S + (WIN_E - WIN_S) * 0.005, 0.66,",
         "ax.text(WIN_S + (WIN_E - WIN_S) * 0.36, 0.04,"),
        ("'MEF2 motif hit (metacluster_179.4 / JASPAR MA0052.4)',",
         "'MEF2 motif hits · metacluster_179.4',"),

        # 位置注释分层: E1/E2 在 0.30, E3/E4 在 0.10
        ("ax.text(mid, 0.25, loc_note[lab]",
         "ax.text(mid, {'E1': 0, 'E2': 0.75, 'E3': 0.10, 'E4': 0.10}[lab], loc_note[lab]"),
        # 距离标注分层: E1/E2 在 -0.08, E3/E4 在 -0.24
        ("dists = {88_707_581: '-14.8 kb (TSS)'",
         "dists = {88_707_581: 'intergenic, −14.8 kb (TSS)'"),
        ("ax.text((s + e) / 2, -0.10, dists[s]",
         "ax.text((s + e) / 2, {'E1': -0.08, 'E2': -0.08, 'E3': -0.24, 'E4': -0.40}[lab], dists[s]"),
        # E 标签改单行 (两行版 NES 文本与相邻标签相碰); E4 右移避开 E3
        (r"ax.text(mid, 1.00, f'{lab}\nNES={nes}', fontsize=7.5,",
         r"ax.text(mid + {'E4': 25_000}.get(lab, 0), 1.00, f'{lab}  NES={nes}', fontsize=7.5, "
         r"bbox=dict(boxstyle='square,pad=0.15', fc='white', ec='none', alpha=0.85),"),
        # 图例移到 gene model 轴下方 (原压在底部注释上)
        ("ax.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.0, 0.03),",
         "ax_vio.legend(handles=handles, loc='upper right', bbox_to_anchor=(1.0, 1.0),"),
        # xlabel 缩短 (2^X-1 说明移图注; 两处 box/violin 分支)
        ("ax_vio.set_xlabel('Nrxn3 expr. (2$^X$−1)', fontsize=8.0, labelpad=2)",
         "ax_vio.set_xlabel('Nrxn3 expr.', fontsize=6.5, labelpad=2)"),
        # 图例竖排 (3 列横排越页)
        ("frameon=False, fontsize=7.0, ncol=3, handlelength=1.5,",
         "frameon=False, fontsize=6.5, ncol=1, handlelength=1.2,"),
        # E1 位置注释并入距离标注, 删除原文本
        ("'intergenic · TSS −14.8 kb'", "''"),
        # violin 绿框移入轴内 + 4 行短行 + 6.5pt (原越页 2cm)
        ("ax_vio.text(hi * 0.48, 5.05,",
         "if False: ax_vio.text(hi * 0.48, 5.05,"),
        (r"'MOL3-specific\nmedian ≈7.0 vs ≈3.1–3.8\n≈2×; MWU p<1e-22',",
         r"'MOL3-specific\nmedian ≈7.0\nvs ≈3.1–3.8\n≈2×; MWU p<1e-22',"),
        (r"fontsize=7.0, color='#1b7a24', ha='left', va='center',",
         r"fontsize=6.5, color='#1b7a24', ha='left', va='center',"),
        # violin 刻度只留 0/5/10 (15 越页)
        ("ax_vio.set_xlim(0, hi)", "ax_vio.set_xlim(0, hi); ax_vio.set_xticks([0, 5, 10])"),
        # 输出到 replot
        ("f'{OUT}/grn_mef2c_nrxn3_track.pdf'", "os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(LEGACY, 'figures')), 'replot/S5e_nrxn3_track.pdf')"),
    ]
    for old, new in reps:
        assert old in src, f'panel_e patch 未命中: {old[:60]!r}'
        src = src.replace(old, new)
    # load_nrxn3_expr: 原需 8GB pkl+scenicplus 模块, 改读预提取缓存
    # (cache/nrxn3_expr.npz, 由 scplus_obj_final.pkl stub-unpickle 提取,
    #  中位数与原图标注 'median 7.0 vs 3.1 3.8' 经 2^X-1 换算完全一致)
    src = re.sub(r'def load_nrxn3_expr\(\):.*?(?=# ===== CELL: load =====)',
                 '''def load_nrxn3_expr():\n    z = np.load(os.path.join(LEGACY, 'figures/scripts/cache/nrxn3_expr.npz'))\n    return {k: z[k] for k in z.files}\n\n\n''',
                 src, flags=re.S)
    src = src.replace("out_pdf = '", "out_pdf = '")     # 占位
    src = src.replace("grn_mef2c_nrxn3_track.pdf",
                      "S5e_nrxn3_track.pdf")
    # 输出目录: 保持原脚本逻辑, 但把写盘路径改到 replot
    src = src.replace("OUTDIR = '" + LEGACY + "/pi/plots-scenic+'",
                      f"OUTDIR = '{OUT}'") if "OUTDIR = '" in src else src
    path = '/tmp/track_small.py'
    open(path, 'w').write(src)
    g = {'__name__': '__main__', '__file__': path}
    old_argv = sys.argv
    sys.argv = [path]                     # 原脚本自带 argparse, 隔离 sys.argv
    try:
        exec(compile(open(path).read(), path, 'exec'), g)
    finally:
        sys.argv = old_argv


# ================================================================
# f: GRN 网络 x3  5.5 x 5.5 cm
# ================================================================
def panel_f():
    import networkx as nx
    from adjustText import adjust_text
    tags = ['mol2', 'mol6', 'mol3_6']
    for tag in tags:
        nodes = pd.read_csv(f'{PSC}/grn_{tag}_cytoscape_nodes.csv')
        edges = pd.read_csv(f'{PSC}/grn_{tag}_cytoscape_edges.csv')
        G = nx.Graph()
        for _, r in nodes.iterrows():
            G.add_node(r.Node_ID, type=r.Type, marker=str(r.is_marker) == 'yes',
                       deg=int(r.degree))
        for _, r in edges.iterrows():
            G.add_edge(r.source, r.target, w=float(r.weight))
        pos = nx.spring_layout(G, weight='w', seed=11, iterations=500)
        # 布局: median 质心居中 + 按 75 分位半径缩放主簇, 离散节点钳在 r=1.08 圈上
        xs = sorted(p[0] for p in pos.values()); ys = sorted(p[1] for p in pos.values())
        cx_, cy_ = xs[len(xs) // 2], ys[len(ys) // 2]
        rads = sorted(((x - cx_) ** 2 + (y - cy_) ** 2) ** 0.5
                      for x, y in pos.values())
        sc = 0.78 / max(rads[int(0.75 * len(rads))], 1e-9)
        pos = {n: ((x - cx_) * sc, (y - cy_) * sc) for n, (x, y) in pos.items()}
        pos = {n: ((lambda r: (p[0] / r * 1.08, p[1] / r * 1.08) if r > 1.08 else p)(\
                    (p[0] ** 2 + p[1] ** 2) ** 0.5))
               for n, p in pos.items()}

        fig, ax = plt.subplots(figsize=(5.5 / 2.54, 5.5 / 2.54))
        ax.set_axis_off()
        # 边 (权重 -> 线宽)
        wmax = max(d['w'] for *_, d in G.edges(data=True))
        for u, v, d in G.edges(data=True):
            p, q = pos[u], pos[v]
            ax.plot([p[0], q[0]], [p[1], q[1]], '-', color='#B5B5B5',
                    lw=0.25 + 0.9 * d['w'] / wmax, zorder=1, solid_capstyle='round')
        # 节点: TF 大红, marker 蓝, 其余灰; 尺寸随 degree
        for n, d in G.nodes(data=True):
            p = pos[n]
            if d['type'] == 'TF':
                ms, c, z = 26, '#D2413B', 4
            elif d['marker']:
                ms, c, z = 10, '#4DAECB', 3
            else:
                ms, c, z = 6, '#C9C9C9', 2
            ax.plot(*p, 'o', ms=ms * 0.12, mfc=c, mec='white', mew=0.2, zorder=z)
        # 标签: TF + marker + degree>=3, cap 22 (按 marker degree 优先)
        cand = [n for n, d in G.nodes(data=True)
                if d['type'] == 'TF' or d['marker'] or d['deg'] >= 3]
        cand.sort(key=lambda n: (G.nodes[n]['type'] != 'TF',
                                 -(G.nodes[n]['deg'] if G.nodes[n]['marker'] else 0),
                                 G.nodes[n]['marker'] is False))
        cand = cand[:22]
        t2node = {}
        texts = []
        for n in cand:
            t = ax.text(*pos[n], n, fontsize=6, zorder=5,
                        color='#1A1A1A' if G.nodes[n]['type'] == 'TF' else '#333333')
            texts.append(t); t2node[id(t)] = n
        adjust_text(texts, ax=ax, expand_text=(1.15, 1.3),
                    expand_points=(1.15, 1.3),
                    force_text=(0.12, 0.2), force_points=(0.05, 0.1))
        # TF 名加粗
        for t in texts:
            if t.get_text() in HUB_TFS:
                t.set_fontweight('bold'); t.set_color('#B02020')
        # 确定性消解残余重叠 (adjustText 后处理): 文字-文字推开 + 文字-节点避开
        fig.canvas.draw()
        ren = fig.canvas.get_renderer()
        inv = ax.transData.inverted()
        node_r = {'TF': 0.075, 'marker': 0.045}.get  # 数据坐标近似半径
        GAP = 0.007   # 最小字间距 (数据单位, ~0.015cm)
        def dbb(t):
            e = t.get_window_extent(renderer=ren)
            (x0, y0), (x1, y1) = inv.transform([(e.x0, e.y0), (e.x1, e.y1)])
            return (x0 - GAP, y0 - GAP, x1 + GAP, y1 + GAP)
        for it in range(400):
            hit = False
            bs = [(t,) + dbb(t) for t in texts]
            for i in range(len(bs)):
                for j in range(i + 1, len(bs)):
                    t, a0, a1, a2, a3 = bs[i]
                    u, b0, b1, b2, b3 = bs[j]
                    ox = min(a2, b2) - max(a0, b0)
                    oy = min(a3, b3) - max(a1, b1)
                    if ox > 0.0003 and oy > 0.0003:
                        hit = True
                        dx = (a0 + a2) / 2 - (b0 + b2) / 2
                        dy = (a1 + a3) / 2 - (b1 + b3) / 2
                        if abs(dx) < 1e-9 and abs(dy) < 1e-9:
                            dx, dy = 0.02, -0.012
                        n = (dx * dx + dy * dy) ** 0.5
                        push = max(min(ox, oy) * 0.6, 0.014)
                        t.set_position((t.get_position()[0] + dx / n * push,
                                        t.get_position()[1] + dy / n * push))
                        bs[i] = (t,) + dbb(t)
            # 后 100 轮只做文字-文字 (避免与节点回避相互振荡)
            for t, a0, a1, a2, a3 in bs if it < 300 else ():
                for n, d in G.nodes(data=True):
                    px, py = pos[n]
                    r = node_r(d['type'], 0.03)
                    if a0 - r < px < a2 + r and a1 - r < py < a3 + r:
                        hit = True
                        dx = (a0 + a2) / 2 - px
                        dy = (a1 + a3) / 2 - py
                        if abs(dx) < 1e-9 and abs(dy) < 1e-9:
                            dx, dy = 0.02, 0.02
                        m = (dx * dx + dy * dy) ** 0.5
                        t.set_position((t.get_position()[0] + dx / m * 0.05,
                                        t.get_position()[1] + dy / m * 0.05))
            if not hit:
                break
        if tag == 'mol6':
            bs = [(t,) + dbb(t) for t in texts]
            for i in range(len(bs)):
                for j in range(i + 1, len(bs)):
                    t, a0, a1, a2, a3 = bs[i]
                    u, b0, b1, b2, b3 = bs[j]
                    ox = min(a2, b2) - max(a0, b0)
                    oy = min(a3, b3) - max(a1, b1)
                    if ox > 0 and oy > 0:
                        print(f'  [dbg] {t.get_text()} x {u.get_text()} '
                              f'ox={ox:.4f} oy={oy:.4f}')
        # 导线: 按最终文字位置, 从文本框边缘连到对应节点
        for t in texts:
            p = pos[t2node[id(t)]]
            e = t.get_window_extent(renderer=ren)
            (x0, y0), (x1, y1) = inv.transform([(e.x0, e.y0), (e.x1, e.y1)])
            cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
            dxs, dys = p[0] - cx, p[1] - cy
            if (dxs * dxs + dys * dys) ** 0.5 < 0.15:
                continue
            cands = [v for v in (((x1 - cx) / dxs) if dxs else 1e9,
                                 ((x0 - cx) / dxs) if dxs else 1e9,
                                 ((y1 - cy) / dys) if dys else 1e9,
                                 ((y0 - cy) / dys) if dys else 1e9) if v > 0]
            tp = min(cands)
            ax.plot([cx + dxs * tp, p[0]], [cy + dys * tp, p[1]], '-',
                    color='#999999', lw=0.25, zorder=2.5)
        ax.set_xlim(-1.25, 1.25); ax.set_ylim(-1.25, 1.25)
        fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
        fig.savefig(f'{OUT}/S5f_{tag}_network.pdf')
        plt.close(fig)
        print(f'saved S5f_{tag} ({len(G)} nodes, {len(texts)} labels)')


HUB_TFS = set()


if __name__ == '__main__':
    import sys
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    # TF 名单 (用于 f 加粗): 各网络 TF
    for tag in ['mol2', 'mol6', 'mol3_6']:
        n = pd.read_csv(f'{PSC}/grn_{tag}_cytoscape_nodes.csv')
        HUB_TFS.update(n[n.Type == 'TF'].Node_ID)
    if which in ('all', 'b'): panel_b()
    if which in ('all', 'f'): panel_f()
    if which in ('all', 'e'): panel_e()
    print('ALL DONE')
