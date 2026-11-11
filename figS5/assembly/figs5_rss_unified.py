#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
figs5_rss_unified.py — Figure S5 面板 A（eRSS 排序）按 Codex #90 重出：**六格统一 y 量程**
=========================================================================================
为什么要改（Codex 批注 #45/#90 的根因，实测定位，不是它"眼神不好"就算了）：
  旧版 `figs5_repanels.panel_b()` 里
      else: ax.set_yticklabels([])
  —— 右列三格（MFOL / MOL3 / MOL6）**一个 y 刻度都没有**，而且六格各自 autoscale：
      COP_NFOL 0.167–0.621 / MOL2 0.167–0.378 / MOL3_6 0.167–0.331 ...
  读者于是把左列的刻度搬到右列上读：
      · MOL6 真值 max 0.513，按 MOL3_6 的 0.20–0.30 轴去读 → "约 0.33"（正是它报的数）
      · MFOL 真值 max 0.343，按 COP_NFOL 的 0.2–0.6 轴去读 → "接近 0.6"（也正是它报的数）
  ⇒ 数值本身没错（交付表 S5a_rss_scores_96.csv 复核：MOL6 恰有 5 个 >0.5），
     **错的是版式允许跨格串读**。本版三条修法：
       1) 六格共用同一 y 量程 (0.15, 0.65) 与同一组刻度 (0.2/0.4/0.6)，**每格都画刻度**；
       2) 画一条 eRSS = 0.5 的虚线基准（正文那句 "five MOL6 regulons above 0.5" 变成可数的）；
       3) 红点与标签由 top-10 收到 **top-5**（MOL6 的 top-5 恰就是那 5 个 >0.5 的），
          标签密度减半，只给 rank-1 附数值锚定。
  完整 top-10（含 rank 与 eRSS）另出 **Table S5**（`figS5/assembly/make_table_s5_rss.py`）。

数值一个都没动：数据仍读同一份 `rss_scores_filtered_96.csv`（与交付表逐字符相同）。

⚠ 为什么本文件自带一份 `postfix_b2` 而不复用 `figs5_repanels.postfix_b`：
  旧 postfix 用 `page.get_drawings()` 的**尺寸窗口**去认六个坐标框，而 MuPDF 会把
  样式相同的相邻路径**合并成一个 drawing**（本版面下六个框只认出 3 个，其中一个 bbox
  还是两框并集 7.12 × 1.81 cm），于是 30 个红标签只有 13 个进入消解、成品仍有碰字。
  这里改成**由 matplotlib 自己报坐标框**（存 sidecar JSON），红点仍从矢量里认（实测 30 个）。

用法: PYTHONPATH=<work>/pylibs python3 scripts/figs5_rss_unified.py [--top 5]
输出: <work>/figures/replot/S5a_rss_grid_v2.pdf  (17.60 × 11.20 cm，与旧槽位同尺寸)
"""
import os, sys
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import figs5_repanels as R          # 复用 resolve_ink / 字体与配色
from figkit import CM
import pdfix
from collections import defaultdict

OUTDIR = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'replot')
SRC = os.path.join(os.environ.get('DSCI_SCENIC_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pi/plots-scenic+')), 'rss_scores_filtered_96.csv')
ORDER = ['COP_NFOL', 'MFOL', 'MOL2', 'MOL3', 'MOL3_6', 'MOL6']      # 与交付版同序
YLO, YHI = 0.15, 0.65
YTICKS = [0.2, 0.4, 0.6]
REF = 0.50                          # 正文用的阈值线
W, H = 17.60, 11.20                 # 槽位尺寸（与交付版一致，页面不动）


def short(n):
    return (n.split('_+')[0].split('_(')[0]
             .replace('_extended', ' ext').replace('_direct', ' dir'))


def build(out=None, ntop=4, spread=6.0, verbose=True):
    out = out or f'{OUTDIR}/S5a_rss_grid_v2.pdf'
    rss = pd.read_csv(SRC, index_col=0)
    fig, axes = plt.subplots(3, 2, figsize=(W / 2.54, H / 2.54))
    fig.subplots_adjust(left=min(0.13, 1.30 / W), right=0.985, top=0.955, bottom=0.065,
                        wspace=0.16, hspace=0.40)
    for ax, ct in zip(axes.flat, ORDER):
        s = rss.loc[ct].sort_values(ascending=False)
        x = np.arange(len(s))
        ax.plot(x, s.values, '-', color=R.BLUE, lw=0.8, zorder=1)
        ax.plot(x, s.values, 'o', color=R.BLUE, ms=1.7, zorder=2)
        # 基准线：先画，压在曲线下面
        ax.axhline(REF, color=R.GREY, lw=0.5, linestyle=(0, (2, 2)), zorder=1)
        top = s.nlargest(ntop)
        # 红点集合 = top-n ∪ 全部 >0.5（MOL6 有 5 个 >0.5，只标 top-4 的话线上会有
        # 一个蓝点，"数红点" 就不成立）；标签仍只给 top-n。MOL6 的 top-4 恰是正文
        # 点名的 Creb5/Tbx2/Nfia/Mbnl2，第 5 个 Nfix 只红不标。
        above = {n for n in s.index if s[n] > REF}
        redset = list(top.index) + [n for n in s.sort_values(ascending=False).index
                                    if n in above and n not in top.index]
        xi_all = [s.index.get_loc(n) for n in redset]
        ax.plot(xi_all, [s[n] for n in redset], 'o', color=R.RED, ms=2.4, zorder=3)
        xi = [s.index.get_loc(n) for n in top.index]
        ax.set_xlim(-2, 97)
        ax.set_ylim(YLO, YHI)                     # ← 六格统一量程
        ax.set_yticks(YTICKS)                     # ← 每格都画刻度
        ax.tick_params(axis='y', labelsize=6)
        labels = [short(n) + (f' {v:.2f}' if i == 0 else '')     # rank-1 附数值
                  for i, (n, v) in enumerate(zip(top.index, top.values))]
        # 初始位置就**沿 x 轴错开**（每个标签比前一个右移 spread 个 rank）：
        # 前 4–5 名的点几乎同高、又都挤在 rank 0–5，从点位置起步的消解器推不开（实测
        # COP_NFOL 框内 4000 轮不收敛、成品上 "Rxrg/Lcorl/Tcf12" 三个标签互碰）。
        texts = [ax.text(xx + spread * (i + 1), yy + 0.018, t, fontsize=6, color=R.RED, zorder=4)
                 for i, (xx, yy, t) in enumerate(zip(xi, top.values, labels))]
        from adjustText import adjust_text
        adjust_text(texts, x=[xx for xx in xi], y=list(top.values), ax=ax,
                    expand_text=(1.1, 1.25), force_text=(0.15, 0.25),
                    align=False, autoalign='')
        fig.canvas.draw()
        ren0 = fig.canvas.get_renderer()
        inv0 = ax.transData.inverted()
        for t in texts:
            e = t.get_window_extent(renderer=ren0)
            t.set_ha('center'); t.set_va('center')
            t.set_position(inv0.transform(((e.x0 + e.x1) / 2, (e.y0 + e.y1) / 2)))
        fig.canvas.draw()
        R.resolve_ink(fig, ax, texts, pts=list(zip(xi_all, [s[n] for n in redset])), dot_ms=2.4)
        ax.autoscale(False)
        ax.set_title(ct, fontsize=7)
        ax.tick_params(axis='x', labelsize=6)
        if ax in axes[:, 0]:
            ax.set_ylabel('eRegulon specificity\nscore (eRSS)', fontsize=6.5)
            ax.yaxis.set_label_coords(-0.105, 0.5)
        if ax in axes[2, :]:
            ax.set_xlabel('eRegulon rank', fontsize=6.5)
        else:
            ax.set_xticklabels([])
    # 基准线只在左上格命名一次
    axes[0, 0].text(48, REF + 0.012, 'eRSS = 0.5', fontsize=6, color=R.GREY, ha='left',
                    va='bottom')
    tmp = f'{OUTDIR}/.S5a_unified_raw.pdf'
    DPI = 100
    fig.savefig(tmp, dpi=DPI)
    # 坐标框由 matplotlib 自己报（PDF 与 Agg 同为左下原点，翻成 MuPDF 的自上而下 cm）
    rend = fig.canvas.get_renderer()
    rects = []
    for ax in axes.flat:
        bb = ax.get_window_extent(renderer=rend)
        x0, x1 = bb.x0 / DPI * 2.54, bb.x1 / DPI * 2.54
        y0, y1 = bb.y0 / DPI * 2.54, bb.y1 / DPI * 2.54
        rects.append((x0, H - y1, x1, H - y0))
    plt.close(fig)
    postfix_b2(tmp, out, rects, verbose=verbose)   # 成品级墨迹消解 + 导线
    return out


def postfix_b2(src, out, rects, verbose=True):
    """同 postfix_b，但坐标框来自调用方（见模块头 ⚠ 说明）。"""
    doc = fitz.open(src)
    page = doc[0]
    dots = []
    for d in page.get_drawings():
        r = d['rect']
        if (d['type'] in ('f', 'fs') and r.width / CM < 0.2 and r.height / CM < 0.2
                and d.get('fill')
                and max(abs(a - b) for a, b in zip(d['fill'], (0.753, 0.188, 0.157))) < 0.02):
            dots.append(((r.x0 + r.x1) / 2 / CM, (r.y0 + r.y1) / 2 / CM))
    ws = pdfix.words(page)
    doc.close()
    us = pdfix.units_of(ws)
    RED = (0.753, 0.188, 0.157)

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
    ymax = max(r[3] for r in rects)
    tot = 0
    for r, mv in groups.items():
        obs = [(x - 0.05, y - 0.05, x + 0.05, y + 0.05) for x, y in dots
               if r[0] <= x <= r[2] and r[1] <= y <= r[3]]
        down = 0.15 if r[3] > ymax - 0.5 else 0.85        # 末行不许溢出，其余可进行距
        # 左沿内收 0.03cm：否则 rank-1 标签会被推到正好压在坐标框左边线上（成品上
        # "Pax5 0.34"/"Usf1 0.38" 首字符与竖线相接）。
        bounds = (r[0] + 0.10, r[1] - 0.15, r[2] + 0.1, r[3] + down)
        left = pdfix.resolve_units(us, mv, bounds=bounds, obstacles=obs, gap=0.012, iters=4000)
        if left:
            left = pdfix.resolve_units(us, mv, bounds=(r[0] + 0.10, r[1] - 0.4,
                                                       r[2] + 0.6, r[3] + 1.3),
                                       obstacles=obs, gap=0.012, iters=6000)
        if left and verbose:
            print(f'    [b2] rect({r[0]:.1f},{r[1]:.1f}) 仍未收敛: {len(mv)} 标签')
        tot += left
    if verbose:
        print('    每框可动标签数: %s / 共 %d' % (
            sorted(len(v) for v in groups.values()), sum(len(v) for v in groups.values())))
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
        print(f'  postfix_b2: {len(us)} 行单元/{len(ws)} 词, 导线 {len(leaders)} 条, '
              f'残余重叠 {tot} | 成品 {n} 词, 自检重叠 {len(bad)}')
        for a, b, ar in bad[:8]:
            print(f'      残留: {a["t"]!r} <> {b["t"]!r} {ar*100:.2f}mm2')
    return out


if __name__ == '__main__':
    a = sys.argv
    def opt(k, d, cast=float):
        p = [s for s in a if s.startswith('--%s=' % k)]
        return cast(p[0].split('=', 1)[1]) if p else d
    ntop = opt('top', 4, int)
    spread = opt('spread', 6.0)
    out = opt('out', None, str)
    p = build(out=out, ntop=ntop, spread=spread)
    import fitz
    d = fitz.open(p)[0]
    print('%s  %.2f × %.2f cm' % (os.path.basename(p), d.rect.width / CM, d.rect.height / CM))
