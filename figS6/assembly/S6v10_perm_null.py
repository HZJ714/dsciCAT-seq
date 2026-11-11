# -*- coding: utf-8 -*-
"""S6 v10：Fig S6-E 的置换零分布——74 性状 Δ1 RNA vs ATAC 配对 rho 的性状标签置换。

已交付 T4b（O1_own_multiome_paired/tables/T4b_cross_channel_rho_variants.csv）只存了
P_perm 与 null_975 两个汇总数，零分布向量没有落盘。本脚本按同一判据重算：
  统计量 = spearmanr(d1_rna, d1_atac)（74 个性状，O1 同核，dsciCAT-seq own 臂）
  零分布 = 固定 ATAC 顺序打乱（trait-pairing permutation），N=10,000，seed 20261104
对账：观测 rho、单侧 P、null 97.5% 必须与 T4b 交付值在 MC 误差内一致，否则报错退出。
产出 tables_v8/figS6v10_permnull_o1.csv（向量）+ figS6v10_permnull_summary.csv（汇总+对账行）。
"""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

ROUTE = WORK + "/figures/routeA_out"
OUT = ROUTE + "/SUPPORT_v2/tables_v8"
T4b = f"{ROUTE}/O1_own_multiome_paired/tables/T4b_cross_channel_rho_variants.csv"
T4 = f"{ROUTE}/O1_own_multiome_paired/tables/T4_delta1_74traits_both_channels.csv"
N, SEED = 10000, 20261104

D = pd.read_csv(T4)
d_r = D["d1_rna"].values.astype(float)
d_a = D["d1_atac"].values.astype(float)
assert len(d_r) == len(d_a) == 74, (len(d_r), len(d_a))
obs, p_sp = spearmanr(d_r, d_a)

rng = np.random.default_rng(SEED)
null = np.empty(N)
for i in range(N):
    null[i] = spearmanr(d_r, d_a[rng.permutation(74)]).statistic
p_perm = (1.0 + np.sum(null >= obs)) / (N + 1.0)
nl975 = float(np.percentile(null, 97.5))

# ---- 与已交付 T4b（variant=raw）汇总数对账（MC 误差容限）----
TB = pd.read_csv(T4b); r0 = TB[TB["variant"] == "raw"].iloc[0]
rho_del, P_del, n975_del = float(r0["rho_cross_channel"]), float(r0["P_perm"]), float(r0["null_975"])
ok = (abs(obs - rho_del) < 1e-6) and (p_perm < 0.01 and P_del < 0.01 or abs(p_perm - P_del) < 0.05) \
     and abs(nl975 - n975_del) < 0.05
print("[对账] obs %.4f vs 交付 %.4f | P %.4f vs %.4f | null975 %.3f vs %.3f -> %s"
      % (obs, rho_del, p_perm, P_del, nl975, n975_del, "OK" if ok else "MISMATCH"))
if not ok:
    raise SystemExit("[对账失败] 零分布判据与已交付 T4b 不一致，图不得使用")

pd.DataFrame({"rho_perm": null}).to_csv(OUT + "/figS6v10_permnull_o1.csv", index=False)
pd.DataFrame([
    dict(stat="observed cross-channel rho", value=round(obs, 4), source=T4b, source_column="rho_cross_channel"),
    dict(stat="permutation P (one-sided)", value=round(p_perm, 6), source=T4b, source_column="P_perm (delivered %.4f)" % P_del),
    dict(stat="null 97.5 pct", value=round(nl975, 4), source=T4b, source_column="null_975 (delivered %.4f)" % n975_del),
    dict(stat="null mean", value=round(float(np.mean(null)), 6), source=T4b, source_column="null_mean (delivered %.6f)" % float(r0["null_mean"])),
    dict(stat="permutation draws", value=N, source="computed here"),
    dict(stat="seed", value=SEED, source="computed here"),
]).to_csv(OUT + "/figS6v10_permnull_summary.csv", index=False)
print("written figS6v10_permnull_o1.csv (%d 行) + summary" % N)
