# -*- coding: utf-8 -*-
"""a1b_splithalf_74.py — 关闭 Table S6 行 199：74 个性状的 split-half Δ1（两个通道）。

口径（必须照读，这是 199 的全部含义）：
  * 逐核 scDRS norm_score 直接取盘上已交付的两通道打分（RNA full74/、ATAC out_gs74/），
    同核交集 = O1 的 5,962 核；Δ1 = median(MOL6) − median(COP/NFOL)（与 O1 的 MAT/COP 常量一致）。
  * 分半是**按细胞类型分层的随机二分**（保证每半里稀有的 COP/NFOL 都不为空），R 次重复。
  * 因为 norm_score 是在全体 5,962 核上做的经验归一化，两半共享该归一化 —— 本分析量化的是
    **抽样变异**（n 减半后 Δ1 还稳不稳、性状间排序还认不认），不是"两半各自独立重跑 scDRS"的
    归一化变异。后者需要重打分，另见 _missing.json / 回执。
用法：/opt/conda/bin/python3 a1b_splithalf_74.py
"""
import os, glob, time, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr

WORK = os.environ.get('DSCI_WORK_ROOT', '')
D_RNA = f"{WORK}/FY/oligo/scdrs_traits_draft/full74"
D_ATAC = f"{WORK}/FY/oligo/fig6_rerun_20261019/dualomics_probe/out_gs74"
H5_R = f"{WORK}/FY/oligo/3k_dsci_oligo_subset.h5ad"
OUT = f"{WORK}/figures/routeA_out/SUPPORT_v2/tables_v12"
os.makedirs(OUT, exist_ok=True)
MAT, COP = "MOL6", "COP/NFOL"
R, SEED = 1000, 20261112
t0 = time.time()
log = lambda *a: print("[%5ds]" % (time.time() - t0), *a, flush=True)


def load_scores(root):
    out = {}
    for f in sorted(glob.glob(f"{root}/*/*.score.gz")):
        tr = os.path.basename(f).split(".full_score")[0].split(".score")[0]
        if tr in out:
            continue
        out[tr] = pd.read_csv(f, sep="\t", index_col=0)["norm_score"]
    return pd.DataFrame(out)


log("读两通道 74 性状逐核分数 ...")
Sr, Sa = load_scores(D_RNA), load_scores(D_ATAC)
cells = Sr.index.intersection(Sa.index)
log("RNA %s | ATAC %s | 同核交集 %d" % (Sr.shape, Sa.shape, len(cells)))
Sr, Sa = Sr.loc[cells], Sa.loc[cells]

import anndata as ad
ct = ad.read_h5ad(H5_R, backed="r").obs.loc[cells, "cell_type"].astype(str)
traits = sorted(Sr.columns)
im = {c: np.where(ct.values == c)[0] for c in (MAT, COP)}
log("核数：%s=%d  %s=%d（分半后每半约 %d / %d）"
    % (MAT, len(im[MAT]), COP, len(im[COP]), len(im[MAT]) // 2, len(im[COP]) // 2))

rng = np.random.default_rng(SEED)
# 分半方式与 O1 的 T5 完全一致：MOL6、COP/NFOL 各自 permutation 后取前/后半（精确对半），
# 两个通道共用同一组半（这样 xchan_same_half / xchan_half1_vs_half2 才有定义）。
idx_mat = im[MAT]; idx_cop = im[COP]
Xr = Sr[traits].to_numpy(dtype=float); Xa = Sa[traits].to_numpy(dtype=float)
dAr = np.zeros((R, len(traits))); dBr = np.zeros((R, len(traits)))
dAa = np.zeros((R, len(traits))); dBa = np.zeros((R, len(traits)))
draws = []
for k in range(R):
    pa = rng.permutation(len(idx_mat)); ha1, ha2 = idx_mat[pa[:len(pa) // 2]], idx_mat[pa[len(pa) // 2:]]
    pb = rng.permutation(len(idx_cop)); ka1, ka2 = idx_cop[pb[:len(pb) // 2]], idx_cop[pb[len(pb) // 2:]]
    dAr[k] = np.median(Xr[ha1], 0) - np.median(Xr[ka1], 0)
    dBr[k] = np.median(Xr[ha2], 0) - np.median(Xr[ka2], 0)
    dAa[k] = np.median(Xa[ha1], 0) - np.median(Xa[ka1], 0)
    dBa[k] = np.median(Xa[ha2], 0) - np.median(Xa[ka2], 0)
    draws.append(dict(draw=k, seed=SEED,
                      rna_within_channel=spearmanr(dAr[k], dBr[k]).correlation,
                      atac_within_channel=spearmanr(dAa[k], dBa[k]).correlation,
                      xchan_same_half=spearmanr(dAr[k], dAa[k]).correlation,
                      xchan_half1_vs_half2=spearmanr(dAr[k], dBa[k]).correlation))
DR = pd.DataFrame(draws)
DR.to_csv(f"{OUT}/splithalf_T5_reliability_draws.csv", index=False)
T5REF = pd.read_csv(f"{WORK}/figures/routeA_out/O1_own_multiome_paired/tables/T5_split_half_reliability.csv").set_index("metric")
summ = []
for m in ["rna_within_channel", "atac_within_channel", "xchan_same_half", "xchan_half1_vs_half2"]:
    v = DR[m].to_numpy(dtype=float)
    ref = float(T5REF.loc[m, "median_split_half_rho"]) if m in T5REF.index else float("nan")
    summ.append(dict(metric=m, median_split_half_rho=float(np.median(v)),
                     q25=float(np.quantile(v, .25)), q75=float(np.quantile(v, .75)),
                     n_draw=len(v), seed=SEED,
                     delivered_T5_median_20draw=round(ref, 5),
                     delta_vs_delivered=round(float(np.median(v)) - ref, 5),
                     note="subsample split-half on the delivered per-nucleus norm_score; "
                          "shared empirical normalisation, quantifies sampling variability only"))
SM = pd.DataFrame(summ)
SM.to_csv(f"{OUT}/splithalf_T5_reliability_summary.csv", index=False)
log("T5 复算（%d draw，逐值已存盘）:" % R)
log(SM[["metric", "median_split_half_rho", "q25", "q75", "delivered_T5_median_20draw", "delta_vs_delivered"]].round(4).to_string(index=False))

rows, per_trait = [], []
for name, dA_, dB_ in [("RNA", dAr, dBr), ("ATAC", dAa, dBa)]:
    X_ = Xr if name == "RNA" else Xa
    full = np.median(X_[idx_mat], 0) - np.median(X_[idx_cop], 0)
    diff = dA_ - dB_
    rho_ab = spearmanr(dA_.mean(0), dB_.mean(0))
    for i, tr in enumerate(traits):
        sd = float(np.std(diff[:, i], ddof=1))
        rows.append(dict(channel=name, trait=tr, n_full=int(len(cells)),
                         n_mature=int(len(idx_mat)), n_precursor=int(len(idx_cop)),
                         delta1_full=round(float(full[i]), 5),
                         delta1_halfA_mean=round(float(dA_[:, i].mean()), 5),
                         delta1_halfB_mean=round(float(dB_[:, i].mean()), 5),
                         delta1_halfA_sd=round(float(dA_[:, i].std(ddof=1)), 5),
                         half_diff_mean=round(float(diff[:, i].mean()), 6),
                         half_diff_se=round(sd / np.sqrt(R), 6),
                         ci95_lo=round(float(np.percentile(diff[:, i], 2.5)), 5),
                         ci95_hi=round(float(np.percentile(diff[:, i], 97.5)), 5),
                         sign_concordance=round(float((np.sign(dA_[:, i]) == np.sign(dB_[:, i])).mean()), 4),
                         both_halves_same_sign_as_full=bool(np.sign(dA_[:, i].mean()) == np.sign(dB_[:, i].mean()) == np.sign(full[i])),
                         reps=R, seed=SEED,
                         split_half_type="subsample split-half (shared empirical normalisation)",
                         what_it_quantifies="sampling variability of Delta1 at n/2, NOT re-normalisation variability"))
    t = pd.DataFrame(rows).query("channel == @name")
    log("%s: 74 性状 Δ1 全样本中位 %.3f；两半均值之差中位 %.4f；符号一致率中位 %.3f；两半排序 Spearman=%.3f；同号性状 %d/74"
        % (name, np.median(full), np.median(t.half_diff_mean.astype(float)),
           t.sign_concordance.median(), rho_ab.correlation, int(t.both_halves_same_sign_as_full.sum())))
    per_trait.append(dict(channel=name, rho_halfA_vs_halfB=round(float(rho_ab.correlation), 4),
                          rho_p=float(rho_ab.pvalue), n_traits=len(traits),
                          median_sign_concordance=round(float(t.sign_concordance.median()), 4),
                          median_half_diff_abs=round(float(t.half_diff_mean.abs().median()), 5),
                          traits_all_concordant=int(t.both_halves_same_sign_as_full.sum()),
                          seed=SEED, reps=R))

res = pd.DataFrame(rows).sort_values(["channel", "trait"])
res.to_csv(f"{OUT}/splithalf_delta1_all74.csv", index=False)
pd.DataFrame(per_trait).to_csv(f"{OUT}/splithalf_across74_summary.csv", index=False)
log("交付 %s  %d 行 → splithalf_delta1_all74.csv（199 由 3 性状补全到 74×2）" % (OUT, len(res)))
worst = res.assign(x=res.half_diff_mean.abs()).sort_values("x", ascending=False).head(3)
log("最不稳定的 3 个性状：\n%s" % worst[["channel", "trait", "delta1_full", "half_diff_mean", "ci95_lo", "ci95_hi", "sign_concordance"]].to_string(index=False))
# 更新 missing 台账（199 已闭环）
mp = f"{OUT}/_missing.json"
m = json.load(open(mp)) if os.path.exists(mp) else {"missing": []}
m["missing"] = [x for x in m["missing"] if "split" not in x.lower()]
m["closed_199"] = dict(rows=len(res), seed=SEED, reps=R, caveat="subsample split-half; shared empirical normalisation")
json.dump(m, open(mp, "w"), indent=1)
log("done")
