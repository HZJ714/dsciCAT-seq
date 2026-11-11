# -*- coding: utf-8 -*-
"""S6 v11 F 格新增梯级数据源：trait_O1_npeak（2026-11-06，去混杂 Tier A 执行）。

口径与交付 S6v10_bootstrap_units.py 完全同语法（B=1,000、seed=20261104、按 74 性状重抽、
spearmanr），唯一区别：两侧 Δ1 先各自对 ln n_peak 做 OLS 残差化（**带截距**——
DECONFOUND_TIERA 产物④§5 的教训），再进入同一 bootstrap。协变量列 = 人侧公共 trait 级
ln n_peak（seven4_confound_ALL_cohorts.csv，GSE231584/GRCh38 来源，见产物⑤）。

输出（全部新文件，不改任何 v8/v10 交付表）：
  SUPPORT_v2/tables_v11/figS6v11_bootstrap_npeak.csv   tidy: unit,b,rho
  SUPPORT_v2/tables_v11/figS6v11_unit_summary_ext.csv  交付 summary + 本行
对照锚：D10（B=10,000, seed 20261102）给 0.708（CI 0.550–0.815, perm P 0.0001）。
"""
import os, sys

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
sys.path.insert(0, WORK + "/scripts")
import numpy as np, pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr

R = WORK + "/figures/routeA_out"
V8 = R + "/SUPPORT_v2/tables_v8"
V11 = R + "/SUPPORT_v2/tables_v11"; os.makedirs(V11, exist_ok=True)
B, SEED = 1000, 20261104

T4 = pd.read_csv(R + "/O1_own_multiome_paired/tables/T4_delta1_74traits_both_channels.csv")
S = pd.read_csv(R + "/FIG7_v1/tables/seven4_confound_ALL_cohorts.csv")[["trait", "lnp"]]
T = T4.merge(S, on="trait", how="left", validate="one_to_one")
assert T.lnp.notna().all(), "74 性状 ln n_peak 未对齐"

def resid(y, x):                      # OLS 残差化，带截距（add_constant）
    return np.asarray(y, float) - sm.OLS(np.asarray(y, float), sm.add_constant(np.asarray(x, float))).fit().fittedvalues

x = resid(T.d1_rna, T.lnp); z = resid(T.d1_atac, T.lnp)
point = float(spearmanr(x, z).statistic)
rng = np.random.default_rng(SEED)
vals = np.array([spearmanr(x[i], z[i]).statistic
                 for i in (rng.integers(0, len(x), len(x)) for _ in range(B))])
lo, hi = np.percentile(vals, [2.5, 97.5])

# 残余混杂自检（产物④§1 的 +0.055 / −0.015）：调整后各通道 Δ1 与 lnp 的 Spearman
rr = (spearmanr(T.lnp, x).statistic, spearmanr(T.lnp, z).statistic)

# 2026-11-06 作者复审追加：n_hit 调整级也要上面（raw→n_hit→n_peak 三连读，让读者自己看到
# n_hit ≈ 无操作 +0.007、n_peak 才动真格 −0.111）。协变量 = 各通道自身 ln n_hit（交付 D6 规则），
# 同样带截距残差化；bootstrap 用**独立随机流**（同 seed 重播种），不影响上面 n_peak 列的逐位数值。
xh = resid(T.d1_rna, np.log(T.n_hit_rna)); zh = resid(T.d1_atac, np.log(T.n_hit_atac))
nh = float(spearmanr(xh, zh).statistic)
rng2 = np.random.default_rng(SEED)
valsh = np.array([spearmanr(xh[i], zh[i]).statistic
                  for i in (rng2.integers(0, len(xh), len(xh)) for _ in range(B))])
loh, hih = np.percentile(valsh, [2.5, 97.5])

pd.DataFrame([dict(unit="trait_O1_npeak", b=b, rho=float(v)) for b, v in enumerate(vals)]
             ).to_csv(V11 + "/figS6v11_bootstrap_npeak.csv", index=False)
pd.DataFrame([dict(unit="trait_O1_nhit", b=b, rho=float(v)) for b, v in enumerate(valsh)]
             ).to_csv(V11 + "/figS6v11_bootstrap_nhit.csv", index=False)
SUM = pd.read_csv(V8 + "/figS6v10_unit_summary.csv")
new = [dict(unit="trait_O1_npeak", n=len(x), rho_point=round(point, 4),
            boot_median=round(float(np.median(vals)), 4), ci_lo=round(float(lo), 4),
            ci_hi=round(float(hi), 4), B=B, seed=SEED,
            pairs="O1 same-nucleus Delta1, ln n_peak-adjusted (n = 74)"),
       dict(unit="trait_O1_nhit", n=len(xh), rho_point=round(nh, 4),
            boot_median=round(float(np.median(valsh)), 4), ci_lo=round(float(loh), 4),
            ci_hi=round(float(hih), 4), B=B, seed=SEED,
            pairs="O1 same-nucleus Delta1, ln n_hit-adjusted (n = 74)")]
pd.concat([SUM, pd.DataFrame(new)], ignore_index=True).to_csv(
    V11 + "/figS6v11_unit_summary_ext.csv", index=False)

print("point rho (n_peak adj) = %.3f   [D10 锚: 0.708]" % point)
print("boot median %.3f  95%% CI [%.3f, %.3f]  (B=1,000 seed 20261104)  [D10 锚: 0.550–0.815, B=10,000]"
      % (np.median(vals), lo, hi))
print("残余混杂 ρ(ln n_peak, 残差) RNA %+.3f / ATAC %+.3f  [回执锚: +0.055 / −0.015]" % rr)
print("n_hit 级: point %.3f [回执锚 0.826]  boot median %.3f  CI [%.3f, %.3f]"
      % (nh, np.median(valsh), loh, hih))
print("written -> %s/figS6v11_bootstrap_{npeak,nhit}.csv, figS6v11_unit_summary_ext.csv" % V11)
print("written -> %s/figS6v11_bootstrap_npeak.csv, figS6v11_unit_summary_ext.csv" % V11)
