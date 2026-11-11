#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig1c_compute.py (v2) — Figure 1 panel C 重画取数：ATAC cuts / 单细胞 tile / RNA depth
=====================================================================================
为什么不用 Signac CoveragePlot（本环境 2026-10-28 读 Signac 1.9.0 源码后的结论）：
  norm.value = count / group.scale.factors × scale.factor，
  group.scale.factors = AverageCounts(group) × n_cells(group) = 该组总计数，
  scale.factor 缺省 = median(该组总计数)。
  → 三个文库分开画（每组只 1 个 group）时 scale.factor 恰等于该组总值，归一化整个抵消，
    图上显示的是**原始 cut 数**，行间不可比。要可比必须显式给同一个分母口径。
  → 这里自己取数并统一成「每 10⁶ 个该文库分子的信号 / 100-bp 窗」，ATAC 与 RNA 同一族单位。

口径（会写进图注）：
  ATAC 分子数分母 = 该文库 ATAC peak 矩阵的计数总和（= 分配给 filtered cells 的 cut 数）。
      分子只计**落在该文库自己的 MACS2 峰内**的 cut 位点（片段两端各 1 个 cut，与 Signac
      CutMatrix 同义），这样分子/分母同一个 universe；峰外的弥散信号不进图。
  tile   = 同一区间内每个 filtered cell 的 cut 数（原始值，不归一化），按区间内 cut 数
           降序取前 300 个 cell（= Signac TilePlot(order.by='top', tile.cells=300)）。
  RNA    = 比对深度（bigWig 每碱基 depth / bedgraph 分箱 depth）在 100-bp 窗内的碱基和；
           分母 = 该 track 全基因组的比对碱基和（bigWig header sumData；bedgraph 全文件 Σ宽×深）。
           read 长在分子分母同时出现、自动抵消，所以不需要假设 read 长度。
  单位统一：signal per 10^6 molecules of that library, per 100-bp window。

数据源核对（2026-10-28）：
  * human.bw / mouse.bw（HG12、HG14 的 sub/ 目录）经三项检查判定为 **RNA 覆盖**，不是 ATAC：
      (1) 与同文库 ATAC 峰集无富集关系（峰内/峰外深度仅 1.35–1.55x，峰区占 20–37% 面积，
          top20 窗只有 10–11/20 落在峰内 = 随机期望）；
      (2) 与 fragments 文件的 1-kb 剖面相关系数仅 0.32–0.45；
      (3) MIR17HG 基因体/侧翼深度比 3.2x（转录本特征）。
  * 细胞数与图例一致：HEK293T ATAC 16,165 / 27,760 / 1,970（dsciCAT-seq-50k/100k/HT-scCAT-seq）。
输出：cache/region/fig1c_tracks.npz
"""
import os, subprocess, json
import numpy as np

D    = os.environ.get('FIG1_DATA', '')
CCH  = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'scripts/cache/region')
OUT  = f'{CCH}/fig1c_tracks.npz'
WIN  = 100
NCELL_TILE = 300

LIB  = ['dsciCAT-seq-50k', 'dsciCAT-seq-100k', 'HT-scCAT-seq']
TAG  = {'dsciCAT-seq-50k': 'HG12', 'dsciCAT-seq-100k': 'HG14', 'HT-scCAT-seq': 'HT'}
SPS  = {
 'HEK293T': dict(build='hg38', chrom='chr13', lo=91340000, hi=91370000, sp='Human',
   frag={'HG12':'HG/HG12/sub/human_fragments.tsv.gz',
         'HG14':'HG/HG14/sub/human_fragments.tsv.gz',
         'HT'  :'scCAT-seq-HT/IL9/ATAC-sub/human_fragments.tsv.gz'},
   mtx ={'HG12':'HG/HG12/sub/ATAC/Human','HG14':'HG/HG14/sub/ATAC/Human',
         'HT' :'scCAT-seq-HT/IL9/ATAC-sub/Human'},
   rna ={'HG12':'HG/HG12/sub/human.bw','HG14':'HG/HG14/sub/human.bw',
         'HT' :'scCAT-seq-HT/IL9/RNA-sub/outs/human_cleaned.bedgraph'}),
 'NIH3T3': dict(build='mm10', chrom='chr11', lo=69090000, hi=69130000, sp='Mouse',
   frag={'HG12':'HG/HG12/sub/mouse_fragments.tsv.gz',
         'HG14':'HG/HG14/sub/mouse_fragments.tsv.gz',
         'HT'  :'scCAT-seq-HT/IL9/ATAC-sub/mouse_fragments.tsv.gz'},
   mtx ={'HG12':'HG/HG12/sub/ATAC/Mouse','HG14':'HG/HG14/sub/ATAC/Mouse',
         'HT' :'scCAT-seq-HT/IL9/ATAC-sub/Mouse'},
   rna ={'HG12':'HG/HG12/sub/mouse.bw','HG14':'HG/HG14/sub/mouse.bw',
         'HT' :'scCAT-seq-HT/IL9/RNA-sub/outs/mouse_cleaned.bedgraph'}),
}

def read_totals(path, key=None):
    d = {}
    if not os.path.exists(path): return d
    for ln in open(path):
        p = ln.rstrip('\n').split('\t')
        if len(p) < 2: continue
        d[p[0]] = float(p[1])
    return d

def peaks_in_region(features, chrom, lo, hi):
    out = []
    for ln in open(features):
        loc = ln.split('\t')[0].strip()
        try:
            c, rz = loc.split(':'); a, z = rz.split('-')
        except ValueError:
            continue
        if c == chrom and int(z) > lo and int(a) < hi:
            out.append((int(a), int(z)))
    return np.array(sorted(out), dtype=np.int64) if out else np.zeros((0, 2), np.int64)

def cut_bins(frag_gz, barcodes, chrom, lo, hi, peaks):
    """tabix 只读目标区间的片段；返回 (每 100bp 窗的峰内 cut 数, {cell: 该 cell 的 cut 向量}, filtered cell 数)"""
    nb = (hi - lo) // WIN
    cells = set(x.strip() for x in open(barcodes))
    prof = np.zeros(nb, np.int64); per = {}
    starts, ends = peaks[:, 0], peaks[:, 1]
    def inside(pos):
        if len(peaks) == 0: return False
        j = np.searchsorted(ends, pos, side='left')       # 第一个 end>pos 的峰
        return j < len(peaks) and starts[j] <= pos
    p = subprocess.Popen(['tabix', frag_gz, f'{chrom}:{lo+1}-{hi}'],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for ln in p.stdout:
        f = ln.split('\t')
        c = f[3]
        if c not in cells: continue
        a, b = int(f[1]), int(f[2])
        d = per.get(c)
        if d is None:
            d = per[c] = np.zeros(nb, np.int8)
        for cut in (a, b):
            if lo <= cut < hi and inside(cut):
                k = (cut - lo) // WIN
                prof[k] += 1
                if d[k] < 120: d[k] += 1
    p.wait()
    return prof, per, len(cells)

def rna_bw(path, chrom, lo, hi):
    import pyBigWig
    b = pyBigWig.open(path)
    v = np.nan_to_num(np.asarray(b.values(chrom, lo, hi), dtype=float))
    tot = float(b.header()['sumData']); b.close()
    nb = (hi - lo) // WIN
    return v[:nb * WIN].reshape(nb, WIN).sum(1), tot

def rna_bedgraph(path, chrom, lo, hi, totfile):
    nb = (hi - lo) // WIN
    acc = np.zeros(nb)
    prog = '$1==c && $3>s && $2<e {lo=($2>s?$2:s); hi=($3<e?$3:e); if(hi>lo) printf "%d %d %s\\n", lo, hi-lo, $4}'
    p = subprocess.run(['awk', '-v', f'c={chrom}', '-v', f's={lo}', '-v', f'e={hi}', prog, path],
                       capture_output=True, text=True, check=True)
    for ln in p.stdout.split('\n'):
        if not ln: continue
        s0, w, d = ln.split(' '); s0, w, d = int(s0), int(w), float(d)
        k0 = (s0 - lo) // WIN; k1 = min((s0 + w - 1 - lo) // WIN, nb - 1)
        if k0 == k1: acc[k0] += d * w
        else:
            ks = np.arange(k0, k1 + 1)
            acc[k0:k1 + 1] += d * (np.minimum(s0 + w, lo + (ks + 1) * WIN) - np.maximum(s0, lo + ks * WIN))
    key = 'human' if chrom == 'chr13' else 'mouse'
    tot = read_totals(totfile).get(key + '_cleaned', np.nan)
    return acc, tot

def main():
    atot = read_totals(f'{CCH}/totals_atac.txt')          # key = 矩阵目录相对路径
    btot = f'{CCH}/totals_rna_bedgraph.txt'
    store = {}
    for spn, s in SPS.items():
        lo, hi, chrom = s['lo'], s['hi'], s['chrom']
        for lib in LIB:
            t = TAG[lib]; mdir = s['mtx'][t]
            peaks = peaks_in_region(f"{D}/{mdir}/features.tsv", chrom, lo, hi)
            prof, per, ncell = cut_bins(f"{D}/{s['frag'][t]}", f"{D}/{mdir}/barcodes.tsv",
                                        chrom, lo, hi, peaks)
            mat_tot = atot.get(mdir, np.nan)
            order = sorted(per.items(), key=lambda kv: -int(kv[1].sum()))[:NCELL_TILE]
            tile = np.vstack([kv[1] for kv in order]).astype(np.int8) if order else np.zeros((0, (hi-lo)//WIN), np.int8)
            rp = s['rna'][t]
            if rp.endswith('.bw'):
                rprof, rtot = rna_bw(f'{D}/{rp}', chrom, lo, hi)
            else:
                rprof, rtot = rna_bedgraph(f'{D}/{rp}', chrom, lo, hi, btot)
            store[f'{spn}|{lib}|atac'] = prof
            store[f'{spn}|{lib}|rna']  = rprof
            store[f'{spn}|{lib}|tile'] = tile
            store[f'{spn}|{lib}|tot']  = np.array([mat_tot, rtot, ncell, len(per), len(peaks)], float)
            print('%-8s %-16s ATAC 峰内cut=%7d/%-12.3g (%.2f per1e6/max窗) max窗=%5d | '
                  'RNA bases=%11.4g/%-11.4g max窗=%8.0f | cells %d/%d peaks %d' %
                  (spn, lib, prof.sum(), mat_tot, prof.max()/mat_tot*1e6, prof.max(),
                   rprof.sum(), rtot, rprof.max(), len(per), ncell, len(peaks)), flush=True)
    np.savez_compressed(OUT, **store)
    print('wrote', OUT, os.path.getsize(OUT) // 1024, 'KB')

if __name__ == '__main__':
    main()
