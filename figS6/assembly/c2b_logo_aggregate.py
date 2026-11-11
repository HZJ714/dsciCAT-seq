#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""c2b_logo_aggregate.py — 只聚合已完成的 LOGO 结果（不重跑 scDRS）。

c2_logo.py 是"生成 gs → 启动作业 → 聚合"一体且**无 resume**，直接重跑会把 1,050 个作业
再算一遍（≈2 h）。本脚本从 out_logo/b*/ 读现成结果，走与 c2_logo 完全相同的聚合口径：

    contribution(g) = Δ1_full(same batch) − Δ1_LOO(g)

同批 FULL 对照是必需的：scDRS 的经验归一化与控制核抽样会随基因集组成漂移，
不同批的 Δ1 不能直接相减。

输出：logo_delta1_contributions.csv、logo_meta.json（覆盖校准批的 5 行版本）
"""
import os, sys, glob, json
import numpy as np, pandas as pd

WORK = os.environ.get('DSCI_WORK_ROOT', '')
V1 = WORK + "/figures/routeA_out/GENELEVEL_v1"
OUT = V1 + "/out_logo"
H5AD = WORK + "/FY/oligo/3k_dsci_oligo_subset.h5ad"
GSFULL = WORK + "/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs"
MAT, COP = "MOL6", "COP/NFOL"
NCTRL = int(os.environ.get("NCTRL", "1000"))


def log(*a):
    print(*a, flush=True)


row = [l.rstrip("\n") for l in open(GSFULL) if l.startswith("PASS_Multiple_sclerosis\t")][0]
GENES = [g.strip().upper() for g in row.split("\t")[1].split(",") if g.strip()]
import anndata as ad
ct = ad.read_h5ad(H5AD, backed="r").obs["cell_type"].astype(str)


BAD = []


def delta1_from(d):
    out = {}
    global BAD
    for f in sorted(glob.glob(d + "/*.score.gz") + glob.glob(d + "/*/*.score.gz")):
        tr = os.path.basename(f).split(".full_score")[0].split(".score")[0]
        if tr in out:
            continue
        # ⚠ 不要用 usecols=["norm_score"]：这些文件表头第一列是**空名**（index 列），
        #    pandas 会静默返回 0 列（实测 n=5962 但 columns=[]），下游 KeyError 看起来像坏文件。
        #    实测全读单个文件 0.0 s，1050 个 ≈2 min，直接按 c2_logo 原口径全读。
        try:
            s = pd.read_csv(f, sep="\t", index_col=0)["norm_score"]
        except Exception as e:
            BAD.append((os.path.relpath(f, OUT), type(e).__name__)); continue
        cells = [c for c in s.index if c in ct.index]
        v, lab = s.loc[cells], ct.loc[cells]
        m, p = v[lab == MAT], v[lab == COP]
        if len(m) > 5 and len(p) > 5:
            out[tr] = float(m.median() - p.median())
    return out


# ---- 参照 Δ1（交付分数，与主图同口径）
ref = pd.read_csv(glob.glob(WORK + "/FY/oligo/scdrs_traits_draft/full74/*/*Multiple_sclerosis*full_score.gz")[0],
                  sep="\t", index_col=0)
cells = [c for c in ref.index if c in ct.index]
v, lab = ref["norm_score"].loc[cells], ct.loc[cells]
DELTA_FULL = float(v[lab == MAT].median() - v[lab == COP].median())
log("参照 Δ1(全集, 交付分数) = %.5f" % DELTA_FULL)

batches = sorted(os.path.basename(p) for p in glob.glob(OUT + "/b*"))
rows, perbatch, skipped = [], {}, []
for tag in batches:
    got = delta1_from("%s/%s" % (OUT, tag))
    full_here = [x for k, x in got.items() if k.startswith("FULL_b")]
    bloo = {k: x for k, x in got.items() if k.startswith("LOO_")}
    if not full_here or not bloo:
        skipped.append(tag); continue
    bf = float(np.mean(full_here))
    perbatch[tag] = bf
    for tr, dv in bloo.items():
        parts = tr.split("_")
        gi, g = int(parts[1]), parts[2]
        rows.append(dict(gene=g, gene_index=gi, delta1_LOO=round(dv, 6),
                         delta1_full_same_batch=round(bf, 6),
                         contribution=round(bf - dv, 6), batch=tag))
L = pd.DataFrame(rows)
assert len(L), "聚合结果为空——先确认 out_logo/b*/ 里有可读的 score.gz"
L = L.drop_duplicates("gene_index")
log("批次数 %d（缺同批 FULL 被跳过 %d）；基因行 %d（名单 %d）；坏文件 %d"
    % (len(perbatch), len(skipped), len(L), len(GENES), len(BAD)))
if BAD:
    log("  坏文件前 5：%s" % BAD[:5])
    pd.DataFrame(BAD, columns=["file", "error"]).to_csv(V1 + "/logo_bad_files.csv", index=False)

wv = {}
for tok in row.split("\t")[1].split(","):
    if ":" in tok:
        gg, zz = tok.rsplit(":", 1)
        try:
            wv[gg.strip().upper()] = float(zz)
        except ValueError:
            pass
L["magma_z"] = L.gene.map(lambda g: wv.get(g, np.nan))
L["effect_on_own_delta"] = L.contribution / DELTA_FULL
L = L.sort_values("contribution", ascending=False).reset_index(drop=True)
L["rank_desc"] = np.arange(1, len(L) + 1)
L.to_csv(V1 + "/logo_delta1_contributions.csv", index=False)

pb = pd.Series(perbatch)
zero = L[L.contribution.abs() < 1e-6]
log("同批 FULL 对照 Δ1：%d 批 均值 %.5f SD %.5f（= 批次间归一化漂移尺度）" % (len(pb), pb.mean(), pb.std(ddof=1)))
log("Σ贡献 %.5f vs Δ1(全集) %.5f → 相加不闭合（scDRS 非线性，禁止写成加法分解）" % (L.contribution.sum(), DELTA_FULL))
log("Top8:\n%s" % L.head(8)[["gene", "magma_z", "contribution", "delta1_LOO", "effect_on_own_delta"]].to_string(index=False))
log("Bottom4:\n%s" % L.tail(4)[["gene", "magma_z", "contribution"]].to_string(index=False))
log("贡献≈0 的基因 %d 个（数据中无同源命中的名单项 = 本分析噪声尺度；|c| p95=%.2e）"
    % (len(zero), float(zero.contribution.abs().quantile(.95)) if len(zero) else float("nan")))
json.dump(dict(delta1_full=DELTA_FULL, n_logo=int(len(L)), nctrl=NCTRL, n_batches=int(len(pb)),
               full_batch_sd=float(pb.std(ddof=1)), sum_contrib=float(L.contribution.sum()),
               noise_genes=int(len(zero))),
          open(V1 + "/logo_meta.json", "w"), indent=1)
log("已写 logo_delta1_contributions.csv / logo_meta.json")
