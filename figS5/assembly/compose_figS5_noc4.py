#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compose_figS5_noc4.py — Figure S5 第四版：Codex #90/#91 两条落地（G3 + G4）
==========================================================================
相对交付版 `FigS5_noc3.pdf` 只动两处，其余**逐像素不变**：

  A 面板  `replot/S5b_rss_grid_wide.pdf` → **`replot/S5a_rss_grid_v2.pdf`**
          （六格统一 y 量程 0.15–0.65 + 每格都画刻度 + eRSS=0.5 虚线 + 标签收到 top-4）
          —— 旧版右列三格**一个 y 刻度都没有**且六格各自 autoscale，
             读者把左列刻度搬到右列上读，于是"MOL6 只有 0.33 / MFOL 接近 0.6"
             （Codex #45/#90 的判错正是这两个串位值；数值本身没错）。
  B 面板  第二行左侧空白区 (x 0.15, y 21.45, 4.6 cm 宽) 加**网络图例块**
          （`scripts/network_key.py`：方向 / 边宽 / 边色 / 符号 / 节点色五条，全部实测自成品矢量）
          三张网络的几何、标签、消解结果一字未动（复用 noc3 的 `network_layer()`）。

页面尺寸、纵向常数、C/D 两块与 noc3 完全相同（18.00 × 46.75 cm）。

自检（脚本自带）：
  1) 与交付版 noc3 逐块像素比对 —— A 块**允许**不同，B/C/D 三块 + 字母必须 0 可见差；
  2) 图例块落笔前实测目标空白盒墨占比 = 0；
  3) 成品跑 `final_audit` + `inkcheck`（见 README/CHANGELOG 记录）。

用法: python3 figS5/assembly/compose_figS5_noc4.py    （从仓库根目录或任意位置均可）
输出: <work>/figures/FigS5_noc4.pdf
"""
import os, sys
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
for _d in ('00_Figure_Tools', os.path.join('fig5', 'assembly')):
    _p = os.path.join(_ROOT, _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)
import fitz
import numpy as np
from figkit import CM, letter, save_doc
import compose_figS5_noc3 as C3          # 复用 network_layer / paste / 纵向常数
from network_key import draw_key

B = os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures'))
REP = f'{B}/replot'
OUT = f'{B}/FigS5_noc4.pdf'
A_SRC = f'{REP}/S5a_rss_grid_v2.pdf'          # ← G3 新 A 面板
KEY_AT = (0.15, 21.45)                          # ← G4 图例块左上角 (cm)
KEY_MAXW = 4.60
KEY_BOX = (0.0, 21.0, 4.9, 29.8)                # 落笔前实测必须全空的盒子


def main():
    # ---- 0) 落笔前：目标盒子确实全空 ----
    ref = fitz.open(f'{B}/FigS5_noc3.pdf')[0]
    x0, y0, x1, y1 = KEY_BOX
    pm = ref.get_pixmap(dpi=200, clip=fitz.Rect(x0 * CM, y0 * CM, x1 * CM, y1 * CM), alpha=False)
    a = np.frombuffer(pm.samples, np.uint8).reshape(pm.height, pm.width, pm.n)[:, :, :3]
    frac = (a.min(2) < 235).mean()
    print('图例目标盒 x[%.1f,%.1f] y[%.1f,%.1f] 在交付版上的墨占比 = %.4f%%' % (
        x0, x1, y0, y1, 100 * frac))
    assert frac < 0.0005, '目标盒不空，换位置'
    assert os.path.exists(A_SRC), '先跑 scripts/figs5_rss_unified.py'

    # ---- 1) 组版（除 A 源件外全部沿用 noc3 的函数与常数）----
    out = fitz.open()
    pg = out.new_page(width=C3.W * CM, height=46.75 * CM)
    C3.paste(pg, A_SRC, C3.A_Y0, label='A')
    for f, ttl, x0n, ty, ny in C3.NETS:
        C3.network_layer(pg, f'{C3.PSC}/{f}', x0n, ny, C3.NET_W, ttl, ty)
    letter(pg, 0.02, C3.B_T1 - 0.36, 'B', C3.LETTER_SZ)
    draw_key(pg, KEY_AT[0], KEY_AT[1], KEY_MAXW, imp=(0.5, 5.0))
    C3.paste(pg, f'{REP}/S5c_nrxn3_track_box.pdf', C3.C_Y0, x0=0.0, label='C')
    C3.paste(pg, f'{REP}/S5_stability_wide.pdf', C3.D_Y0, x0=0.0, label='D')
    save_doc(out, OUT, clean=False)

    # ---- 2) 自检：逐块像素比对 noc3 ----
    new = fitz.open(OUT)[0]
    print('成品 %.2f × %.2f cm' % (new.rect.width / CM, new.rect.height / CM))
    DPI = 150
    blocks = [('A 面板(允许变)', 0.0, 11.6), ('B 行1标题带', 11.6, 12.2),
              ('B 行1网络', 12.2, 20.6), ('B 行2标题带', 20.6, 21.12),
              ('B 行2网络(图例区除外)', 24.0, 29.5), ('C 面板', 30.0, 41.4),
              ('D 面板', 41.4, 46.75)]
    for name, by0, by1 in blocks:
        r = fitz.Rect(0, by0 * CM, new.rect.width, by1 * CM)
        p_new = new.get_pixmap(dpi=DPI, clip=r, alpha=False)
        p_old = ref.get_pixmap(dpi=DPI, clip=r, alpha=False)
        A = np.frombuffer(p_new.samples, np.uint8).reshape(p_new.height, p_new.width, p_new.n)[:, :, :3].astype(int)
        Bm = np.frombuffer(p_old.samples, np.uint8).reshape(p_old.height, p_old.width, p_old.n)[:, :, :3].astype(int)
        d = np.abs(A - Bm).max(2)
        tot = d.size
        print('  %-24s y[%.1f,%.1f] 可见差 %d / %d px (%.4f%%)' % (
            name, by0, by1, int((d > 24).sum()), tot, 100 * (d > 24).mean()))
    return OUT


def finish(path=OUT, verbose=True):
    """noc3 配方里的后三步（少了任何一步城市门就会红）：
         1) arialize --inplace —— 贴进来的 C 带源件是 Liberation，必须换成真 Arial
            （不做这一步：C 带与交付版 noc3 差 1.38%，且审计报“非 Arial 痕迹”）；
         2) prune_dead_fonts —— 三张网络去字后留下的死 Type3/DejaVu 声明；
         3) fix_letters_regular —— 面板字母粗体→常规（2026-10-21 全稿口径，交付版已是常规）；
         3) doc.subset_fonts() —— 把 526 KB 整字库 Arial 裁到 30 KB，
            实测逐像素 0 差（150/300 dpi）、文字层逐字符相同
            ⇒ 交付版 noc3 上挂的【豁免: 未子集=Arial,Liberation】在本版**自动取消**。"""
    import subprocess, fitz
    S = os.path.dirname(os.path.abspath(__file__))
    env = dict(os.environ, PYTHONPATH=os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
    for cmd in ([sys.executable, f'{S}/arialize.py', '--inplace', path],
                [sys.executable, f'{S}/prune_dead_fonts.py', path]):
        r = subprocess.run(cmd, env=env, capture_output=True, text=True)
        if verbose:
            print('  $ %s\n%s' % (' '.join(cmd[1:]), (r.stdout or r.stderr).strip()[-400:]))
    # 3) 面板字母 A–D 粗体→常规（2026-10-21 全稿口径；交付版 Figure_S5.pdf 已是常规字重，
    #    少了这步就会把字母退回粗体 —— 与交付件逐块比会多出 125 px）
    step = path + '.letters.tmp'
    subprocess.run([sys.executable, f'{S}/fix_letters_regular.py', path, step, '--prune'],
                   env=env, check=True, capture_output=True, text=True)
    os.replace(step, path)
    d = fitz.open(path)
    d.subset_fonts(verbose=False)
    tmp = path + '.sub.tmp'
    d.save(tmp, garbage=4, deflate=True)
    d.close()
    os.replace(tmp, path)
    return path


if __name__ == '__main__':
    main()
    finish()
