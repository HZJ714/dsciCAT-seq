#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig1c_v4.py — Figure 1 panel C **整块重画**（不再从旧参考图矢量裁贴）
=====================================================================
为什么重画（Codex #237 + 作者 2026-10-28 要求）：
  旧 C 是从 `Figure 1.pdf` 裁下来的 Signac CoveragePlot/TilePlot 拼贴，缺
    ① 坐标轴（窗口内没有任何刻度，读者无法把峰对到基因的哪个部位）；
    ② 基因组版本（hg38 / mm10 只在图注里，图内没有）；
    ③ 归一化口径（CoveragePlot 的 norm 在"每文库单独一组"时整个抵消，图上其实是**原始 cut 数**，
       三行之间不可比 —— 读 Signac 1.9.0 源码得到，见 fig1c_compute.py 头注）；
    ④ 基因名写法：人 `Mir17hg` → 应为 `MIR17HG`；鼠 `Hex7` → 该坐标上的基因是 `Hes7`
       （mm10 chr11:69,120,452-69,123,259，UCSC refGene NM_033041；注释库里没有 Hex7）；
    ⑤ 单细胞 tile 那一块旧图只标 "300 single cells"，**没说是哪个文库**（三套文库都画了才说得通）。
  作者另外问：RNA 轨当年是 BAM→bw→IGV 手工截的，能不能整个用代码重做 → 能：
  RNA 覆盖直接用 bw（pyBigWig）/ bedgraph（awk），ATAC 用 fragments.tsv.gz + 该文库峰集，
  单细胞 tile 用同一批 filtered cells 的 cut 位点，基因模型用当初出图同一套 EnsDb
  （EnsDb.Hsapiens.v86 = GRCh38 / EnsDb.Mmusculus.v79 = GRCm38，dump_genes.R 导出）。

bw 的身份（2026-10-28 xcheck_bw.py，1-Mb 大窗 2×2 同/跨模态对照）：
  r(bw, HT RNA bedgraph)  = 0.95 / 0.92 / (0.40 低信号区) / 0.94
  r(bw, 本文库 ATAC cuts) = 0.21 / 0.21 / -0.01 / 0.23
  正对照 r(ATAC50k, ATAC_HT) = 0.85 / 0.79 / 0.57 / 0.89；跨模态基线 = 0.21 / 0.22 / 0.00 / 0.20
  峰内/峰外深度：bw 1.11–1.37x，ATAC cuts 2.9–29.7x
  ⇒ human.bw / mouse.bw 是 **RNA 覆盖**，不是 ATAC；panel C 的 RNA 行可以放心用它。

单位（图内写明，图注照抄）：
  ATAC = 落在该文库自己 MACS2 峰集内的 cut 位点，每 10⁶ 个 cut（该文库 filtered cells）每 100-bp 窗；
  RNA  = 同一 track 的比对碱基，每 10⁶ 个比对碱基、每 100-bp 窗；
  tile = 该 cell 在该 100-bp 窗内的 cut 数（原始值，色标 0→3+ 截断），每文库取窗内 cut 最多的 120 个细胞。
  两个细胞系各自独立量程（HEK 与 NIH 的绝对信号差 ~7x，共用刻度会把一边压平），
  所以**每一半的每一组都画自己的 y 刻度 + 横向网格线**，并在底部写明。

用法: <conda-env>/bin/python3 fig1c_v4.py
输出: <work>/figures/replot/C_v4_tracks.pdf  （尺寸 = fig1_layout.C_SLOT，1:1 不再缩放）
"""
import os, sys
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))          # 持久盘 PyMuPDF（/opt/conda 每轮重置）
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle, Polygon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig1_layout as FL
from figkit import use_matplotlib_arial

CCH = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'scripts/cache/region')
OUT = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'replot/C_v4_tracks.pdf')
CM = 72 / 2.54

# ---------------------------------------------------------------- 几何（cm）----
W, H = FL.CW, FL.ROWH                       # 11.85 x 6.53 = C 槽位，出图 1:1
LBL_X = 1.45                                # 左标签列宽（组名 + 方法名 + y 刻度）
GRP_W = 0.40                                # 组名（旋转）+ 竖括号 占 0.02–0.42
MTH_R = 1.00                                # 方法名右对齐位置（y 刻度留给 1.00–LBL_X）
RGT, GAP = 0.12, 0.50                       # 右边距 / 两半之间（右半的 y 刻度要放这里）
HALF = (W - LBL_X - RGT - GAP) / 2          # 每一半的绘图宽 = 4.89
X0 = {'HEK293T': LBL_X, 'NIH3T3': LBL_X + HALF + GAP}

ROWS = [('legend', 0.26), ('title', 0.26),
        ('atac', 0.44), ('atac', 0.44), ('atac', 0.44),
        ('tile', 0.52), ('tile', 0.52), ('tile', 0.52),
        ('rna', 0.44), ('rna', 0.44), ('rna', 0.44),
        ('gene', 0.52), ('xax', 0.40), ('note', 0.22)]
GAP_IN, GAP_GRP = 0.02, 0.085               # 组内 / 组间行距
TILE_CELLS = 120                            # 每文库画多少个细胞（3 块 ×120 才放得下）

FS_MTH, FS_TK, FS_GRP, FS_TITLE, FS_NOTE = 4.4, 3.6, 4.6, 5.6, 3.2


def row_y():
    """把 ROWS 从槽位顶部往下排，返回 [(kind, y_top, h)] 与缩放因子 k（超高就等比压）"""
    ngrp = sum(1 for i, (k, _) in enumerate(ROWS) if i and k != ROWS[i - 1][0])
    nin = len(ROWS) - 1 - ngrp
    need = sum(h for _, h in ROWS) + ngrp * GAP_GRP + nin * GAP_IN
    k = min(1.0, (H - 0.02) / need)
    out, y = [], 0.0
    for i, (kind, h) in enumerate(ROWS):
        h *= k
        out.append((kind, y, h))
        nxt = ROWS[i + 1][0] if i + 1 < len(ROWS) else None
        y += h + (GAP_GRP if (nxt and nxt != kind) else GAP_IN) * k
    return out, k


# ---------------------------------------------------------------- 数据 ----------
LIB = ['dsciCAT-seq-50k', 'dsciCAT-seq-100k', 'HT-scCAT-seq']
SHORT = {'dsciCAT-seq-50k': '50k', 'dsciCAT-seq-100k': '100k', 'HT-scCAT-seq': 'HT'}
# 三色 = 旧 C 里逐像素取样核对的同一套（红/绿/蓝），保证作者"色没变"的直觉成立
COL = {'dsciCAT-seq-50k': '#EE2E2B', 'dsciCAT-seq-100k': '#6EBF4B', 'HT-scCAT-seq': '#273B93'}
WIN_INFO = {'HEK293T': dict(chrom='chr13', lo=91340000, hi=91370000, build='hg38'),
            'NIH3T3': dict(chrom='chr11', lo=69090000, hi=69130000, build='mm10')}
TILE_CMAP = LinearSegmentedColormap.from_list('wtblu', ['#FFFFFF', '#5B84C4', '#0B2A6B'])
PERMIL = 'per 10$^6$'


def load():
    z = np.load(f'{CCH}/fig1c_tracks.npz')
    d = {}
    for spn in WIN_INFO:
        for lib in LIB:
            prof = z[f'{spn}|{lib}|atac'].astype(float)
            rna = z[f'{spn}|{lib}|rna'].astype(float)
            tile = z[f'{spn}|{lib}|tile'].astype(float)
            mat_tot, rna_tot, ncell, ncell_reg, npk = z[f'{spn}|{lib}|tot']
            nb = len(prof)
            win = (WIN_INFO[spn]['hi'] - WIN_INFO[spn]['lo']) / nb
            d[(spn, lib)] = dict(atac=prof / mat_tot * 1e6, rna=rna / rna_tot * 1e6,
                                 tile=tile[:TILE_CELLS], nb=nb, win=win,
                                 ncell=int(ncell), ncell_reg=int(ncell_reg), npk=int(npk),
                                 raw_cut=prof.sum())
    genes = {}
    with open(f'{CCH}/genes.tsv') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        for ln in f:
            g = dict(zip(hdr, ln.rstrip('\n').split('\t')))
            genes.setdefault('HEK293T' if g['sp'] == 'hs' else 'NIH3T3', []).append(g)
    return d, genes


def kill_axis(a):
    """彻底关掉一个 axes 的刻度/刻度标签/偏移量文本（gene 行用数据坐标但不许它长轴）"""
    a.set_xticks([]); a.set_yticks([])
    a.tick_params(bottom=False, top=False, left=False, right=False,
                  labelbottom=False, labelleft=False)
    a.xaxis.get_offset_text().set_text('')
    a.yaxis.get_offset_text().set_text('')
    for s in a.spines.values():
        s.set_visible(False)


# ---------------------------------------------------------------- 画图 ----------
def main():
    D, GENES = load()
    fam = use_matplotlib_arial()
    matplotlib.rcParams['mathtext.rm'] = fam
    matplotlib.rcParams['mathtext.default'] = 'regular'
    plt.rcParams.update({'savefig.dpi': 600})

    fig = plt.figure(figsize=(W / 2.54, H / 2.54))     # cm → inch
    layout, kk = row_y()
    Y = {}
    for kind in ('legend', 'title', 'atac', 'tile', 'rna', 'gene', 'xax', 'note'):
        Y[kind] = [(y, h) for (kd, y, h) in layout if kd == kind]

    def axes(x, ytop, w, h):
        return fig.add_axes([x / W, 1 - (ytop + h) / H, w / W, h / H])

    # ---- 顶部图例：方法色线（左） + tile 色标（右）------------------------------
    ly, lh = Y['legend'][0]
    lax = axes(0, ly, W, lh); lax.axis('off')
    lax.set_xlim(0, W); lax.set_ylim(0, 1)              # 这一行用**厘米**当数据坐标
    x = LBL_X
    for lib in LIB:
        lax.plot([x, x + 0.40], [0.5, 0.5], color=COL[lib], lw=1.2, clip_on=False, solid_capstyle='butt')
        lax.text(x + 0.48, 0.5, lib, fontsize=4.6, va='center', ha='left', color='#222222')
        x += 0.48 + 0.098 * len(lib) + 0.42
    cbw = 0.85
    cbx = W - 0.12 - cbw - 1.05              # 右端：色标 + 两端文字
    cax = fig.add_axes([cbx / W, 1 - (ly + lh * 0.86) / H, cbw / W, lh * 0.58 / H])
    cax.imshow(np.linspace(0, 3, 64).reshape(1, 64), aspect='auto', cmap=TILE_CMAP, vmin=0, vmax=3)
    cax.set_xticks([]); cax.set_yticks([])
    for s in cax.spines.values():
        s.set_linewidth(0.3); s.set_color('#666666')
    lax.text(cbx - 0.05, 0.5, 'tile 0', fontsize=FS_TK + 0.2, va='center', ha='right')
    lax.text(cbx + cbw + 0.05, 0.5, '3+ cuts/cell', fontsize=FS_TK + 0.2, va='center', ha='left')

    for spn, info in WIN_INFO.items():
        x0 = X0[spn]
        lo, hi, chrom, build = info['lo'], info['hi'], info['chrom'], info['build']
        ext, step = [lo, hi], (5000 if (hi - lo) <= 30000 else 10000)
        ticks = list(range(int(np.ceil(lo / step) * step), hi + 1, step))
        ymax_a = max(D[(spn, l)]['atac'].max() for l in LIB) * 1.10
        ymax_r = max(D[(spn, l)]['rna'].max() for l in LIB) * 1.10

        # 标题：左=细胞系，右=版本 + 坐标 + 窗宽
        ty, th = Y['title'][0]
        a = axes(x0, ty, HALF, th); a.axis('off')
        a.text(0, 0.30, spn, transform=a.transAxes, fontsize=FS_TITLE, fontweight='bold',
               va='bottom', ha='left', color='#111111')
        a.text(1, 0.30, '%s · %s:%s–%s · %.0f kb' % (build, chrom, format(lo, ','),
               format(hi, ','), (hi - lo) / 1000), transform=a.transAxes, fontsize=4.2,
               va='bottom', ha='right', color='#333333')

        # 覆盖度行（ATAC / RNA 共用）
        def cov_row(i, lib, key, ymax, bg=None):
            ry, rh = Y[key][i]
            a = axes(x0, ry, HALF, rh)
            a.set_xlim(*ext); a.set_ylim(0, ymax)
            if bg:
                a.set_facecolor(bg)
            d = D[(spn, lib)][key]
            xc = (np.arange(len(d)) + 0.5) * D[(spn, lib)]['win'] + lo
            a.fill_between(xc, 0, d, color=COL[lib], lw=0, zorder=3)
            a.plot(xc, d, color=COL[lib], lw=0.25, zorder=4)
            a.set_xticks([])                            # x 轴只在最下面那条画（旧图整块没有 x 轴）
            if i == 0:
                a.set_yticks([0, ymax / 2.2, ymax / 1.1])
                a.set_yticklabels(['0', '%.2f' % (ymax / 2.2), '%.2f' % (ymax / 1.1)])
            else:
                a.set_yticks([0]); a.set_yticklabels(['0'])
            a.tick_params(axis='y', labelsize=FS_TK, length=1.0, width=0.3, pad=1.0)
            # v4.1（作者 2026-10-28：“把 c 和 d 的网格都去掉，有网格感觉不好看”）：
            #   删掉 v4 加的横向网格（0/中/顶三根）+ 竖向定位虚线。
            #   删了不丢信息：ymax_a / ymax_r 是**整组（3 个文库）共用**的量程（第 176/177 行
            #   取的是 max over LIB），所以第一行的三个刻度数字就是这三行的尺；
            #   对坐标靠底部 x 轴 + 基因模对齐。
            for s in a.spines.values():
                s.set_linewidth(0.35); s.set_color('#444444')
            lab = axes(0, ry, MTH_R - 0.02, rh); lab.axis('off')
            lab.text(1.0, 0.5, SHORT[lib], transform=lab.transAxes, fontsize=FS_MTH,
                     va='center', ha='right', color=COL[lib], fontweight='bold')
            return a

        for i, lib in enumerate(LIB):
            cov_row(i, lib, 'atac', ymax_a)
        for i, lib in enumerate(LIB):
            cov_row(i, lib, 'rna', ymax_r, bg='#F6F6F6')

        # 组名（旋转 90°，跨该组三行）+ 竖括号
        for key, txt in (('atac', 'ATAC\n' + PERMIL + ' cuts'),
                         ('tile', 'single-cell tiles\n%d cells each' % TILE_CELLS),
                         ('rna', 'RNA\n' + PERMIL + ' bases')):
            ys = Y[key][0][0]; ye = Y[key][-1][0] + Y[key][-1][1]
            a = axes(0.02, ys, GRP_W, ye - ys); a.axis('off')
            a.text(0.55, 0.5, txt, transform=a.transAxes, fontsize=FS_GRP, va='center',
                   ha='center', rotation=90, color='#111111', linespacing=1.15)
            a.plot([0.98, 0.98], [0.02, 0.98], transform=a.transAxes, lw=0.5,
                   color='#222222', clip_on=False)

        # 单细胞 tile：每文库一块，块左标细胞数
        for i, lib in enumerate(LIB):
            ry, rh = Y['tile'][i]
            a = axes(x0, ry, HALF, rh)
            t = D[(spn, lib)]['tile']
            a.imshow(t, aspect='auto', cmap=TILE_CMAP, vmin=0, vmax=3,
                     interpolation='nearest', extent=[lo, hi, t.shape[0], 0])
            a.set_xlim(*ext); a.set_xticks([])
            a.set_yticks([])                     # 细胞数已经写在左列，不再加 y 刻度
            for s in a.spines.values():
                s.set_linewidth(0.35); s.set_color('#444444')
            lab = axes(0, ry, MTH_R - 0.02, rh); lab.axis('off')
            lab.text(1.0, 0.5, SHORT[lib], transform=lab.transAxes, fontsize=FS_MTH,
                     va='center', ha='right', color=COL[lib], fontweight='bold')

        # 基因模型（EnsDb：人 v86/GRCh38，鼠 v79/GRCm38；模型 = 外显子总长最长的转录本）
        gy, gh = Y['gene'][0]
        a = axes(x0, gy, HALF, gh); kill_axis(a)
        a.set_xlim(*ext); a.set_ylim(0, 1)
        last_end = -1
        for g in sorted(GENES.get(spn, []), key=lambda g: (int(g['start']) + int(g['end'])) / 2):
            s, e = max(int(g['start']), lo), min(int(g['end']), hi)
            ex = []
            for piece in g['exons'].split(';'):
                if not piece:
                    continue
                p, q = piece.split('-')
                p, q = max(int(p), lo), min(int(q), hi)
                if q > p:
                    ex.append((p, q))
            c = (s + e) / 2
            lev = 1 if (last_end > 0 and c - 2500 < last_end) else 0     # 名字打架就分两层
            yb = 0.66 if lev == 0 else 0.30
            last_end = max(last_end, e) if lev == 0 else max(last_end, e)
            a.plot([s, e], [yb + 0.11, yb + 0.11], lw=0.4, color='#222222',
                   solid_capstyle='butt', zorder=2)
            for p, q in ex:
                a.add_patch(Rectangle((p, yb + 0.03), q - p, 0.17, facecolor='#1F3B73',
                                      edgecolor='none', zorder=3))
            ar = 1 if g['strand'] == '+' else -1
            px = e if ar > 0 else s
            a.add_patch(Polygon([[px, yb + 0.115], [px + ar * 600, yb + 0.175],
                                 [px + ar * 600, yb + 0.055]], closed=True,
                                facecolor='#1F3B73', edgecolor='none', zorder=4))
            a.text(c, yb - 0.01, g['symbol'], fontsize=4.2, style='italic', ha='center',
                   va='top', color='#111111', zorder=6,
                   bbox=dict(fc='white', ec='none', pad=0.12, alpha=0.80))

        # x 轴（每一半各自一条：两个细胞系的坐标系不同）
        #   ⚠ 轴脊放在这一行的**最上沿**（axes 高 0.015 cm），刻度标签往下挂进行内，
        #     否则标签会画到行外、和下面那行口径说明叠在一起（v4 第一版就是这么撞的）。
        xy, xh = Y['xax'][0]
        a = axes(x0, xy, HALF, 0.015)
        a.set_xlim(*ext); a.set_ylim(0, 1)
        a.set_xticks(ticks); a.set_xticklabels(['%.3f' % (t / 1e6) for t in ticks])
        # 首/末刻度标签靠内对齐，否则右半最后一个刻度（=窗右界）会戳到页边被裁掉
        _lb = a.get_xticklabels()
        if _lb:
            _lb[0].set_ha('left'); _lb[-1].set_ha('right')
        a.set_yticks([])
        a.tick_params(axis='x', length=1.2, width=0.35, pad=1.2, labelsize=3.8)
        for s_ in ('top', 'right', 'left'):
            a.spines[s_].set_visible(False)
        a.spines['bottom'].set_linewidth(0.4); a.spines['bottom'].set_color('#333333')
        a = axes(x0, xy + 0.20, HALF, 0.18); a.axis('off')
        a.text(0.5, 0.5, '%s (%s) — genomic coordinate (Mb)' % (chrom, build),
               transform=a.transAxes, fontsize=3.8, ha='center', va='center', color='#222222')

    # ---- 底部口径说明（图注照抄）----------------------------------------------
    ny, nh = Y['note'][0]
    a = axes(0.02, ny, W - 0.04, nh); a.axis('off')
    # ⚠ 这一行带 mathtext（10$^6$），matplotlib 对**含公式的字符串**里的非 ASCII 标点会退回
    #   DejaVuSansDisplay（实测 ’ 和 – 两个字符就退，同一字符在纯文本行里又正常），
    #   整页因此多嵌一套字体 → 这里一律用 ASCII 标点（’ → 改写避开，– → -）。
    a.text(0, 0.9, 'ATAC: cuts inside per-library MACS2 peaks / 10$^6$ cuts · '
                   'RNA: aligned read-bases / 10$^6$ read-bases · 100-bp windows · '
                   'tiles: 120 top cells, cuts/100 bp (0-3+) · y-scale shared within each block',
           transform=a.transAxes, fontsize=FS_NOTE, va='top', ha='left', color='#444444')

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, facecolor='white')
    plt.close(fig)
    import fitz
    doc = fitz.open(OUT)
    p = doc[0]
    Wp, Hp = p.rect.width / CM, p.rect.height / CM
    for nm, r in (('Cv4_full', (0, 0, Wp, Hp)), ('Cv4_left', (0, 0, Wp / 2, Hp)),
                  ('Cv4_right', (Wp / 2, 0, Wp, Hp))):
        p.get_pixmap(dpi=800 if nm == 'Cv4_full' else 1100,
                     clip=fitz.Rect(r[0] * CM, r[1] * CM, r[2] * CM, r[3] * CM)).save('/tmp/%s.png' % nm)
    doc.close()
    print('[C] saved %s  页 %.2f x %.2f cm（槽位 %.2f x %.2f）' % (OUT, Wp, Hp, FL.CW, FL.ROWH))
    for spn in WIN_INFO:
        for lib in LIB:
            d = D[(spn, lib)]
            print('   %-8s %-16s ATAC max %.2f (窗内 %d cuts / %d cells) | RNA max %.2f | tile %dx%d'
                  % (spn, lib, d['atac'].max(), d['raw_cut'], d['ncell_reg'], d['rna'].max(),
                     d['tile'].shape[0], d['nb']))


if __name__ == '__main__':
    main()
