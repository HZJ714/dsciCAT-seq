# -*- coding: utf-8 -*-
"""c2_logo.py — ② 面板的正规算法：leave-one-gene-out（LOGO）对 scDRS 成熟−前体差的贡献分解。

为什么不用线性重构：交付的 scDRS 分数是"分层内经验零分布归一化"的结果，实测 w_g·z_g 重构与它只有
r=0.334（见 recon_check），所以任何"Σ基因贡献=Δscore"的加法分解都会与正文报告的分数不一致。
LOGO 直接在交付统计量上做消融：把 MS 基因集里第 i 个基因去掉，用**同一 scdrs compute-score 参数**
重打一次分，Δ1 的变化就是该基因的边际贡献。附带好处：基因集里那些在数据中无同源命中的基因，
其 LOGO 结果必须等于全集聚，于是它们给出本分析自身的噪声尺度（内部负对照）。

内存/并行：每个作业自己读 h5ad（约 2–4 GB 峰值），并行上限 MAXP；启动前与运行中检查空闲内存，
低于 FLOOR_GB 就暂停启动新作业（防止被 OOM killer 干掉）。
用法：conda run -n <你的分析环境> python c2_logo.py [--maxp 20] [--nctrl 1000]   （即 $DSCI_PY）
"""
import os, sys, glob, time, subprocess, json
import numpy as np, pandas as pd

WORK = os.environ.get('DSCI_WORK_ROOT', '')
SCDRS = os.path.expanduser('"~/.local/bin/scdrs')
H5AD = f"{WORK}/FY/oligo/3k_dsci_oligo_subset.h5ad"
COV = f"{WORK}/FY/oligo/scdrs_ms_out/cov.tsv"
GSFULL = f"{WORK}/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs"
OUT = f"{WORK}/figures/routeA_out/GENELEVEL_v1"
GSD, RESD = OUT + "/gs_logo", OUT + "/out_logo"
for d in (GSD, RESD):
    os.makedirs(d, exist_ok=True)
T0 = time.time()
log = lambda *a: print("[%6.0fs]" % (time.time() - T0), *a, flush=True)
MAXP = 20
NCTRL = 1000
FLOOR_GB = 70
MAT, COP = "MOL6", "COP/NFOL"

args = sys.argv[1:]
if "--maxp" in args:
    MAXP = int(args[args.index("--maxp") + 1])
if "--nctrl" in args:
    NCTRL = int(args[args.index("--nctrl") + 1])

# ---- 基因集（与交付同一输入文件、同一行）
row = [l.rstrip("\n") for l in open(GSFULL) if l.startswith("PASS_Multiple_sclerosis\t")][0]
GENES = [g.strip().upper() for g in row.split("\t")[1].split(",") if g.strip()]
log("MS 基因集 %d 个基因（人名义，交给 scDRS 自己做人->鼠映射）" % len(GENES))

# ---- 细胞类型标签（一次读 obs 即可，backed 模式很轻）
import anndata as ad
ct = ad.read_h5ad(H5AD, backed="r").obs["cell_type"].astype(str)


def free_gb():
    return float(os.popen("free -g | awk 'NR==2{print $7}'").read().strip() or 999)


def write_gs(path, rows):
    with open(path, "w") as f:
        f.write("TRAIT\tGENESET\n")
        for tr, gl in rows:
            f.write("%s\t%s\n" % (tr, ",".join(gl)))


def launch(gsfile, tag):
    d = f"{RESD}/{tag}"
    os.makedirs(d, exist_ok=True)
    cmd = [SCDRS, "compute-score", "--h5ad-file", H5AD, "--h5ad-species", "mmusculus",
           "--cov-file", COV, "--gs-file", gsfile, "--gs-species", "hsapiens",
           "--ctrl-match-opt", "mean_var", "--weight-opt", "vs", "--adj-prop", "None",
           "--flag-filter-data", "True", "--flag-raw-count", "False",
           "--n-ctrl", str(NCTRL), "--flag-return-ctrl-raw-score", "False",
           "--flag-return-ctrl-norm-score", "True", "--out-folder", d]
    lp = open(d + "/compute.log", "w")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    return subprocess.Popen(cmd, stdout=lp, stderr=subprocess.STDOUT, env=env), d, lp


def delta1_from(d):
    """返回 {trait: Δ1}（median(MOL6) − median(COP/NFOL)，逐核 norm_score）"""
    out = {}
    for f in sorted(glob.glob(d + "/*.score.gz") + glob.glob(d + "/*/*.score.gz")):
        tr = os.path.basename(f).split(".full_score")[0].split(".score")[0]
        if tr in out:
            continue
        s = pd.read_csv(f, sep="\t", index_col=0)["norm_score"]
        cells = [c for c in s.index if c in ct.index]
        v, lab = s.loc[cells], ct.loc[cells]
        m, p = v[lab == MAT], v[lab == COP]
        if len(m) > 5 and len(p) > 5:
            out[tr] = float(m.median() - p.median())
    return out


# ---- 打包：每批 BATCH 个 LOGO 作业合成一个 gs（scDRS 一次调用可算多个 trait，省掉重复 preprocess）
BATCH = int(os.environ.get("LOGO_BATCH", "50"))
LIMIT = int(os.environ.get("LOGO_LIMIT", "0"))
jobs = [(i, g) for i, g in enumerate(GENES)][:LIMIT or None]
if LIMIT:
    log("校准模式：只做前 %d 个基因（LOGO_BATCH=%s）" % (len(jobs), BATCH))
batches = [jobs[k:k + BATCH] for k in range(0, len(jobs), BATCH)]
log("共 %d 个 LOGO 基因，分 %d 批（每批最多 %d 个性状），并行上限 %d" % (len(jobs), len(batches), BATCH, MAXP))

running, done = [], []
for bi, b in enumerate(batches):
    # 每批内同时跑一次完整基因集 FULL_<tag>：LOGO 与它配对相减，把"去掉一个基因引起的
    # 分层/经验归一化漂移"从基因信号里分离出来（scDRS 非线性，加法不闭合，必须同批对照）。
    rows = [("FULL_b%03d" % bi, list(GENES))] + \
           [("LOO_%04d_%s" % (gi, g), [x for j, x in enumerate(GENES) if j != gi]) for gi, g in b]
    gs = f"{GSD}/b{bi:03d}.gs"
    write_gs(gs, rows)
    while len(running) >= MAXP:
        time.sleep(20)
        wr = [r for r in running if r[1].poll() is not None]
        for tag, p, d, lp in wr:
            lp.close(); done.append((tag, d, p.returncode))
            log("  完成 %s rc=%d（累计 %d/%d 批，用时 %.1f min）" % (tag, p.returncode, len(done), len(batches), (time.time() - T0) / 60))
        running = [r for r in running if r[1].poll() is None]
        if free_gb() < FLOOR_GB:
            log("  ! 空闲内存 %.0f GB < %d GB → 等待" % (free_gb(), FLOOR_GB)); time.sleep(60)
    while free_gb() < FLOOR_GB:
        log("  ! 空闲内存不足（%.0f GB），暂缓启动" % free_gb()); time.sleep(45)
    p, d, lp = launch(gs, "b%03d" % bi)
    running.append(("b%03d" % bi, p, d, lp))
    log("  启动 b%03d（%d 个性状）| 运行中 %d | 空闲内存 %.0f GB" % (bi, len(rows), len(running), free_gb()))
while running:
    time.sleep(25)
    wr = [r for r in running if r[1].poll() is not None]
    for tag, p, d, lp in wr:
        lp.close(); done.append((tag, d, p.returncode))
        log("  完成 %s rc=%d（%d/%d）" % (tag, p.returncode, len(done), len(batches)))
    running = [r for r in running if r[1].poll() is None]

# ---- 参照：全集 Δ1（用已交付分数，口径一致）
ref = pd.read_csv(glob.glob(f"{WORK}/FY/oligo/scdrs_traits_draft/full74/*/*Multiple_sclerosis*full_score.gz")[0], sep="\t", index_col=0)
cells = [c for c in ref.index if c in ct.index]
v, lab = ref["norm_score"].loc[cells], ct.loc[cells]
DELTA_FULL = float(v[lab == MAT].median() - v[lab == COP].median())
log("参照 Δ1(全集, 交付分数) = %.5f（MOL6 n=%d, COP/NFOL n=%d）" % (DELTA_FULL, (lab == MAT).sum(), (lab == COP).sum()))

rows, perbatch = [], {}
for tag, d, rc in sorted(done):
    if rc != 0:
        log("  ! %s 非零退出，跳过" % tag); continue
    got = delta1_from(d)
    full_here = [v for k, v in got.items() if k.startswith("FULL_b")]
    bloo = {k: v for k, v in got.items() if k.startswith("LOO_")}
    if not full_here or not bloo:
        log("  ! %s 缺同批 FULL 对照（full=%d loo=%d），跳过" % (tag, len(full_here), len(bloo))); continue
    bf = float(np.mean(full_here))
    perbatch[tag] = bf
    for tr, dv in bloo.items():
        parts = tr.split("_"); gi, g = int(parts[1]), parts[2]
        rows.append(dict(gene=g, gene_index=gi, delta1_LOO=round(dv, 6),
                         delta1_full_same_batch=round(bf, 6),
                         contribution=round(bf - dv, 6), batch=tag))
L = pd.DataFrame(rows).drop_duplicates("gene_index")
wv = {}
wl = [l.rstrip("\n") for l in open(f"{WORK}/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs") if l.startswith("PASS_Multiple_sclerosis\t")][0]
for tok in wl.split("\t")[1].split(","):
    if ":" in tok:
        gg, zz = tok.rsplit(":", 1)
        try: wv[gg.strip().upper()] = float(zz)
        except ValueError: pass
L["magma_z"] = L.gene.map(lambda g: wv.get(g, np.nan))
L["effect_on_own_delta"] = L.contribution / DELTA_FULL if DELTA_FULL else np.nan
L = L.sort_values("contribution", ascending=False)
L["rank_desc"] = np.arange(1, len(L) + 1)
L.to_csv(OUT + "/logo_delta1_contributions.csv", index=False)

# 内部负对照：本应无影响的基因（名单里在数据无命中者）给出噪声尺度
noise = L[L.gene_index >= len([g for g in GENES])] if False else None
pb = pd.Series(perbatch)
log("同批 FULL 对照 Δ1：%d 批 均值 %.5f SD %.5f（= 批次间归一化漂移尺度）；交付分数参照 %.5f"
    % (len(pb), pb.mean(), (pb.std(ddof=1) if len(pb) > 1 else float("nan")), DELTA_FULL))
log("LOGO 完成：%d 个基因；Σ贡献=%.5f，Δ1(全集)=%.5f → 相加不闭合（scDRS 非线性，不能用加法分解）"
    % (len(L), L.contribution.sum(), DELTA_FULL))
log("Top10（提升成熟−前体差最多）:\n%s" % L.head(10)[["gene", "magma_z", "contribution", "delta1_LOO", "effect_on_own_delta"]].to_string(index=False))
log("Bottom5（去掉后反而增大）:\n%s" % L.tail(5)[["gene", "magma_z", "contribution"]].to_string(index=False))
zero = L[L.contribution.abs() < 1e-6]
log("贡献≈0 的基因 %d 个（这些即在数据中无同源命中的名单项，可作本分析噪声尺度：|c| 的 p95=%.2e）"
    % (len(zero), float(zero.contribution.abs().quantile(.95)) if len(zero) else np.nan))
json.dump(dict(delta1_full=DELTA_FULL, n_logo=len(L), nctrl=NCTRL, sum_contrib=float(L.contribution.sum()),
               noise_genes=len(zero)), open(OUT + "/logo_meta.json", "w"), indent=1)
