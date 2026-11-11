# -*- coding: utf-8 -*-
"""S6 v10 F 格数据源：单元梯各行的 bootstrap rho 分布（B=1,000, seed=20261104）。

每行 = 一个独立抽样单元上的跨通道 Spearman rho；bootstrap 在该单元的"单元层"重抽：
  nucleus      按 5,962 个核重抽（descriptive，核非独立单元——保留灰显）
  library      按 74 性状重抽（FY13/FY14 各自的库内 Δ1 对，T4 表现成列）
  mouse        按 18 只鼠重抽（M2b 逐鼠 rna/atac 轴对比）
  trait O1     按 74 性状重抽（同核 Δ1 对，MOL6-COP/NFOL）
  trait O1 wide 同上，宽对比 Δ1b
  trait D-own  按 74 性状重抽（own 臂 Oligo-OPC，rho_raw）
输出：figS6v10_bootstrap_units.csv（tidy: unit,b,rho）+ figS6v10_unit_summary.csv
"""
import sys, os

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
import os
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
sys.path.insert(0, WORK + "/scripts")
import numpy as np, pandas as pd
from scipy.stats import spearmanr

V8 = WORK + "/figures/routeA_out/SUPPORT_v2/tables_v8"
R_OUT = WORK + "/figures/routeA_out"
B, SEED = 1000, 20261104
rng = np.random.default_rng(SEED)

def boot(x, y, n_unit):
    x = np.asarray(x, float); y = np.asarray(y, float)
    return np.array([spearmanr(x[i], y[i]).statistic
                     for i in (rng.integers(0, len(x), len(x)) for _ in range(n_unit))])

units = {}
T4 = pd.read_csv(R_OUT + "/O1_own_multiome_paired/tables/T4_delta1_74traits_both_channels.csv")
CELL = pd.read_csv(R_OUT + "/O1_own_multiome_paired/tables/cell_level_table.csv", low_memory=False)
M2B = pd.read_csv(R_OUT + "/M2b_authoranno_v1/tables/M2b_per_mouse.csv")
PERC = pd.read_csv(R_OUT + "/D_crosschemistry/tables_log1p_20261104/fig6D_delta1_74traits_per_tech.csv")

# nucleus (MS) —— 按核重抽
xr = CELL["PASS_Multiple_sclerosis_rna"].astype(float).values
xa = CELL["PASS_Multiple_sclerosis_atac"].astype(float).values
units["nucleus"] = (xr, xa, "per-nucleus pairs (n = %d)" % len(xr))
# library FY13 / FY14 —— 按 74 性状重抽（库内成对 Δ1）
for L in ("FY13", "FY14"):
    units["library_" + L] = (T4["d1_rna_" + L].values, T4["d1_atac_" + L].values,
                             "per-trait Delta1 within library %s (n = 74)" % L)
# mouse —— 按 18 只鼠重抽
mm = M2B.dropna(subset=["atac_axis_contrast", "rna_axis_contrast"])
units["mouse"] = (mm["rna_axis_contrast"].values, mm["atac_axis_contrast"].values,
                  "per-mouse axis contrasts (n = %d)" % len(mm))
# trait O1 raw / wide —— 按 74 性状重抽
units["trait_O1"] = (T4["d1_rna"].values, T4["d1_atac"].values, "O1 same-nucleus Delta1 (n = 74)")
units["trait_O1_wide"] = (T4["d1b_rna"].values, T4["d1b_atac"].values, "O1 wide contrast Delta1b (n = 74)")
# trait D-own raw —— own 臂 Oligo-OPC
own = PERC[PERC.tech == "dsciCAT-seq (own)"]
units["trait_D_own"] = (own["d1_rna"].values, own["d1_atac"].values,
                        "D-arm own Delta1, rho_raw (n = 74)")

tidy, summ = [], []
for u, (x, y, note) in units.items():
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    point = float(spearmanr(x, y).statistic)
    vals = boot(x, y, B)
    lo, hi = np.percentile(vals, [2.5, 97.5])
    for b, v in enumerate(vals):
        tidy.append(dict(unit=u, b=b, rho=float(v)))
    summ.append(dict(unit=u, n=len(x), rho_point=round(point, 4), boot_median=round(float(np.median(vals)), 4),
                     ci_lo=round(float(lo), 4), ci_hi=round(float(hi), 4), B=B, seed=SEED, pairs=note))
    print("%-16s n=%-6d rho %.3f  boot 95%% CI [%.3f, %.3f]" % (u, len(x), point, lo, hi))

pd.DataFrame(tidy).to_csv(V8 + "/figS6v10_bootstrap_units.csv", index=False)
S = pd.DataFrame(summ); S.to_csv(V8 + "/figS6v10_unit_summary.csv", index=False)
print("written: figS6v10_bootstrap_units.csv (%d rows), figS6v10_unit_summary.csv" % len(tidy))
