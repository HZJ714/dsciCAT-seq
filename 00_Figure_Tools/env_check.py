#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
env_check.py — 报告本仓库需要的环境变量、是否已设置、解析到哪里
================================================================
用法:
    python3 00_Figure_Tools/env_check.py          # 人看的表格
    python3 00_Figure_Tools/env_check.py --json   # 机器读
    python3 00_Figure_Tools/env_check.py --require  # 必需项缺失时退出码 1（CI 用）

必需（不设置就读不到数据）：DSCI_WORK_ROOT、FIG1_DATA、DSCI_REFDATA_ROOT
其余都有 `<DSCI_WORK_ROOT>/…` 的派生默认值；DSCI_REPO_ROOT 默认 '.'（在仓库根目录运行）。
契约的权威说明在根 README《Data layout and environment variables》。
"""
import json, os, sys

WORK = 'DSCI_WORK_ROOT'
REQUIRED = [
    (WORK, '分析树根：FY/ figures/ pi/ scripts/ fonts/ …（原作者机的工作盘）'),
    ('FIG1_DATA', '原始 count matrix（ResultData/…），仓库不分发，从 CNGBdb CNP0010263 取'),
    ('DSCI_REFDATA_ROOT', '参考注释与索引：GTF / TSS BED / chrom.sizes / bead whitelist'),
]
DERIVED = [
    ('DSCI_FY_ROOT', 'FY', 'oligo/dsciCAT-seq 处理树（ArchR arrows、scDRS 输出）'),
    ('DSCI_FIG_ROOT', 'figures', '交付图与分析输出（routeA_out/…）'),
    ('DSCI_PYLIBS', 'pylibs', '可选：额外的 site-packages（固定版本 PyMuPDF 等）'),
    ('DSCI_FONTS', 'fonts', 'Arial TTF 目录（成图要求 Arial only）'),
    ('DSCI_FONTS_NAMED', 'fonts_named', '按字重命名的字体目录'),
    ('DSCI_SEURAT_ROOT', 'pi/plots-seurat', 'Seurat replot 产物'),
    ('DSCI_SCENIC_ROOT', 'pi/plots-scenic+', 'SCENIC+ 产物（eRSS、stability、增强子序列）'),
    ('DSCI_SCENIC_SCRIPTS', 'pi/scripts-scenic+', 'SCENIC+ 脚本（Nrxn3 位点轨道等）'),
    ('DSCI_SCRIPTS_ROOT', 'scripts', '分析脚本树（MS GWAS sweep 等）'),
    ('DSCI_PROBE_ROOT', 'routeA_probe', 'fig6 dualomics 探针产物'),
    ('DSCI_REF_ROOT', 'referrence', '参考件（liftOver chain 等，见 figS6/README.md 的 md5）'),
    ('DSCI_MOTIF_ROOT', 'aertslab_motif_colleciton', 'Aerts lab motif 集合'),
    ('DSCI_PI_ROOT', 'pi', 'pi/ 本身（少数脚本直接用）'),
]
OTHER = [
    ('DSCI_REPO_ROOT', '.', '本仓检出；默认当前目录，要求在仓库根目录运行 notebook/脚本'),
    ('DSCI_RESULTDATA_ROOT', '', 'ResultData 根，FIG1_DATA 的上一级（两者给一个即可）'),
    ('DSCI_GS_PROBE_DIR', '', 'fig6 dualomics 探针目录（gene-score 输入）'),
    ('DSCI_GS_RDS', '', 'fig6 用的 gene-score RDS 全路径（人鼠共享那批文库）'),
    ('DSCI_OWN_RNA_H', '', 'own-data RNA h5ad 全路径（fig6 同核对照）'),
    ('DSCI_SCDRS_SUBDIR', '', 'scDRS 表达通道输出子目录名（默认相对 FY 树）'),
    ('DSCI_SCDRS_SUBDIR_ATAC', '', 'scDRS 可及性通道输出子目录名'),
    ('DSCI_TAB_SUBDIR', '', '交付台账 CSV 的子目录名（fig6/figS6 读表用）'),
    ('DSCI_LEGACY_ROOT', '', '只被 figS5/assembly/figs5_repanels.py 用：它按文本改写作者树上那批脚本'),
    ('DSCI_HOMOLOGENE', '', 'homologene.data 路径（fig6 面板 D 的人鼠同源映射）'),
]


def rows():
    out = []
    w = os.environ.get(WORK, '')
    for var, why in REQUIRED:
        v = os.environ.get(var, '')
        out.append((var, 'required', v or '(unset)', why, bool(v)))
    for var, sub, why in DERIVED:
        v = os.environ.get(var, '')
        resolved = v or ((os.path.join(w, sub) if w else '(unset: %s)' % WORK))
        out.append((var, 'derived', resolved, why, bool(v) or bool(w)))
    for var, dflt, why in OTHER:
        v = os.environ.get(var, '')
        out.append((var, 'optional', v or (dflt or '(unset)'), why, True))
    return out


def main():
    js = '--json' in sys.argv
    rws = rows()
    if js:
        print(json.dumps([{'var': a, 'kind': b, 'resolved': c, 'why': d, 'ok': e} for a, b, c, d, e in rws],
                         ensure_ascii=False, indent=1))
    else:
        print('%-22s %-9s %-56s %s' % ('VARIABLE', 'KIND', 'RESOLVED TO', 'SET?'))
        print('-' * 100)
        for a, b, c, d, e in rws:
            print('%-22s %-9s %-56s %s' % (a, b, c[:56], 'yes' if e else 'NO'))
        print()
        for a, b, c, d, e in rws:
            if not e:
                print('  待设置 %s — %s' % (a, d))
    miss = [a for a, b, c, d, e in rws if not e and b == 'required']
    if '--require' in sys.argv and miss:
        print('缺少必需变量: %s' % ', '.join(miss), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
