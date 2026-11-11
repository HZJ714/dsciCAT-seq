# -*- coding: utf-8 -*-
"""c_scdrs_genelevel.py — ② 面板的正规做法 + 交付分数可复现性校验。

做法：把 MS 基因集里的每个基因写成 gs 文件的一行（TRAIT -> 单个基因），用**生成交付分数时的
同一 scdrs compute-score 参数**（n-ctrl 1000, cov-file, weight-opt vs, flag-raw-count False,
flag-filter-data True, gs-species hsapiens, h5ad-species mmusculus）一次跑完 → 得到 890 个基因各自的
逐核 scDRS 分数；再按亚型求成熟−前体差。这是 scDRS 框架内的 gene-level 分数，不是线性重构。

阶段 1（本脚本 --check）：只跑 MS 全集 + 3 个单基因，校验能否复现交付的 full_score（Pearson≈1）。
阶段 2（--full）：跑全部 890 个单基因（分 batch 并行）。
用法：conda run -n <你的分析环境> python c_scdrs_genelevel.py --check|--full   （即 $DSCI_PY；python 3.9 + scdrs 1.0.2）
"""
import os, sys, glob, time, subprocess, shutil
import pandas as pd, numpy as np
from scipy.stats import pearsonr, spearmanr

WORK = os.environ.get('DSCI_WORK_ROOT', '')
SCDRS = os.path.expanduser('"~/.local/bin/scdrs')
H5AD = f"{WORK}/FY/oligo/3k_dsci_oligo_subset.h5ad"
COV = f"{WORK}/FY/oligo/scdrs_ms_out/cov.tsv"
OUT = f"{WORK}/figures/routeA_out/GENELEVEL_v1"
GSDIR, RESDIR = OUT + "/gs", OUT + "/out"
DELIV = f"{WORK}/FY/oligo/scdrs_traits_draft/full74"
os.makedirs(GSDIR, exist_ok=True); os.makedirs(RESDIR, exist_ok=True)
log = lambda *a: print("[%6.0fs]" % (time.time() - T0), *a, flush=True)
T0 = time.time()

# MS 基因集（带权重的那份：magma_10kb_top1000_zscore.74_traits.rv1.gs 的 MS 行）
line = [l.rstrip("\n") for l in open(f"{WORK}/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs")
        if l.startswith("PASS_Multiple_sclerosis\t")][0]
W = {}
for tok in line.split("\t")[1].split(","):
    if ":" in tok:
        g, z = tok.rsplit(":", 1)
        try: W[g.strip().upper()] = float(z)
        except ValueError: pass
import anndata as ad
vn = set(str(x).upper() for x in ad.read_h5ad(H5AD, backed="r").var_names)
GENES = [g for g in sorted(W) if g in vn]
log("MS 加权集 %d 基因，矩阵内命中 %d" % (len(W), len(GENES)))


def write_gs(path, rows):
    """rows: [(trait_name, "G1,G2,...")]；scDRS gs 文件为 tab 分隔 TRAIT<TAB>GENE,GENE"""
    with open(path, "w") as f:
        f.write("TRAIT\tGENESET\n")      # scDRS 读 gs 时要 header，否则首行被当表头吃掉
        for tr, gl in rows:
            f.write("%s\t%s\n" % (tr, ",".join(gl)))


def run(gspath, tag, n_ctrl=1000):
    d = f"{RESDIR}/{tag}"
    os.makedirs(d, exist_ok=True)
    cmd = [SCDRS, "compute-score",
           "--h5ad-file", H5AD, "--h5ad-species", "mmusculus",
           "--cov-file", COV, "--gs-file", gspath, "--gs-species", "hsapiens",
           "--ctrl-match-opt", "mean_var", "--weight-opt", "vs", "--adj-prop", "None",
           "--flag-filter-data", "True", "--flag-raw-count", "False",
           "--n-ctrl", str(n_ctrl), "--flag-return-ctrl-raw-score", "False",
           "--flag-return-ctrl-norm-score", "True", "--out-folder", d]
    lp = open(d + "/compute.log", "w")
    p = subprocess.Popen(cmd, stdout=lp, stderr=subprocess.STDOUT)
    return p, d, lp


def scores_from(d):
    out = {}
    for f in sorted(glob.glob(d + "/*/*.score.gz") + glob.glob(d + "/*/*.full_score.gz")):
        tr = os.path.basename(f).split(".full_score")[0].split(".score")[0]
        if tr in out:
            continue
        out[tr] = pd.read_csv(f, sep="\t", index_col=0)
    return out


if "--check" in sys.argv:
    msl = [l.rstrip("\n") for l in open(f"{WORK}/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs") if l.startswith("PASS_Multiple_sclerosis\t")][0]
    ms_plain = [g.strip().upper() for g in msl.split("\t")[1].split(",") if g.strip().upper() in vn]
    probes = []      # 单基因集会被 scDRS 以 small size 跳过，这里只校验全集
    gs = f"{GSDIR}/check.gs"
    write_gs(gs, [("PASS_Multiple_sclerosis", ms_plain)] + [("G_" + g, [g]) for g in probes])
    log("校验 gs：%d 行（MS 全集 %d 基因 + %d 个单基因）" % (1 + len(probes), len(ms_plain), len(probes)))
    p, d, lp = run(gs, "check")
    log("已启动 PID %d → %s" % (p.pid, d))
    while p.poll() is None:
        time.sleep(20); log("   运行中 %.1f min（尾日志：%s）" % ((time.time() - T0) / 60, open(d + "/compute.log").read().strip().split("\n")[-1][:110]))
    lp.close()
    log("结束 rc=%d" % p.returncode)
    S = scores_from(d)
    log("输出 %d 个性状：%s" % (len(S), list(S)[:6]))
    ref = pd.read_csv(glob.glob(DELIV + "/*/*Multiple_sclerosis*full_score.gz")[0], sep="\t", index_col=0)
    for tr, df in sorted(S.items()):
        if tr in ref.index.name or tr == "PASS_Multiple_sclerosis":
            i = df.index.intersection(ref.index)
            if len(i) > 100:
                for c in ["raw_score", "norm_score"]:
                    log("  %-22s vs 交付 %-10s n=%d Pearson=%.4f Spearman=%.4f 中位 %.3f vs %.3f"
                        % (tr, c, len(i), pearsonr(df.loc[i, c], ref.loc[i, c])[0],
                           spearmanr(df.loc[i, c], ref.loc[i, c]).correlation,
                           df.loc[i, c].median(), ref.loc[i, c].median()))
        else:
            c = df["norm_score"]
            log("  %-22s 逐核 norm_score 中位=%.3f IQR=%.3f  有值核=%d" % (tr, c.median(), c.quantile(.75) - c.quantile(.25), int(c.notna().sum())))
elif "--full" in sys.argv:
    NB = 30
    batches = [GENES[i::NB] for i in range(NB)]
    procs = []
    for i, b in enumerate(batches):
        gs = f"{GSDIR}/gl_b%02d.gs" % i
        write_gs(gs, [("G_" + g, [g]) for g in b])
        p, d, lp = run(gs, "gl_b%02d" % i)
        procs.append((i, p, d, lp, len(b)))
        log("  启动 batch %02d（%d 基因）" % (i, len(b)))
    log("共 %d 个并行作业，%d 基因（并行度受 n-ctrl 与内存约束）" % (len(procs), sum(x[4] for x in procs)))
    pend = list(procs)
    while pend:
        time.sleep
        time.sleep(30)
        still = []
        for i, p, d, lp, n in pend:
            if p.poll() is None:
                still.append((i, p, d, lp, n))
            else:
                lp.close(); log("  batch %02d 完成 rc=%d" % (i, p.returncode))
        pend = still
        if pend:
            log("   进行中 %d/%d 批 ..." % (len(pend), len(procs)))
    log("全部完成，开始汇总")
    rows = []
    for f in sorted(glob.glob(RESDIR + "/gl_b*/*/*.score.gz")):
        tr = os.path.basename(f).split(".full_score")[0].split(".score")[0]
        g = tr[2:]
        df = pd.read_csv(f, sep="\t", index_col=0)
        if "norm_score" not in df:
            continue
        ct = ad.read_h5ad(H5AD, backed="r").obs["cell_type"].astype(str)
        cells = [c for c in df.index if c in ct.index]
        s = df.loc[cells, "norm_score"]; lab = ct.loc[cells]
        mat, cop = s[lab == "MOL6"], s[lab == "COP/NFOL"]
        rows.append(dict(gene=g, magma_z=W.get(g, np.nan), scdrs_gene_score_norm_median=float(s.median()),
                         n_mol6=int(len(mat)), n_cop=int(len(cop)),
                         delta1_gene=float(mat.median() - cop.median()),
                         rna_expr_mean_log=float(np.nan), pct_detect=float("nan")))
    t = pd.DataFrame(rows).sort_values("delta1_gene", ascending=False)
    t.to_csv(OUT + "/genelevel_scdrs_delta1.csv", index=False)
    log("交付 %s/genelevel_scdrs_delta1.csv  %d 基因  Δ1 中位 %.3f，正向 %d 个"
        % (OUT, len(t), t.delta1_gene.median(), int((t.delta1_gene > 0).sum())))
    log("最大 5 个：%s" % t.head(5)[["gene", "magma_z", "delta1_gene"]].to_string(index=False))
