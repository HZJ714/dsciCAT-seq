# -*- coding: utf-8 -*-
"""p0_8c_p2_robust.py —— 38 vs 同面板 293 的 P≈0.015 是真是假？

三重体检：
  R1 留一基因（LOO over the 38）：任意剔掉 1 个基因后 P 的最坏值 —— 防"一个基因撑起来"
  R2 检测率门槛：只保留 max 检测率 ≥30% 的基因，防低检测基因放大效应
  R3 效应量与方向：中位数差、稳健置换（对极端值不敏感）
"""
import os
import numpy as np, pandas as pd
TAB = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'routeA_out/P0_8C/tables')
PLAT = os.path.join(os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')), 'oligo/scdrs_ms_out/genes_of_interest.txt')
SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
NPERM, SEED = 5000, 20261116
M = pd.read_csv(TAB + "/p08c_dual_landing_genomewide.csv")
C38 = set(str(x).upper() for x in pd.read_csv(TAB + "/p08c_c1_rna_atac_by_subtype.csv").mouse.astype(str))
P331 = set(x.strip().upper() for x in open(PLAT) if x.strip())
M["gu"] = M.gene.astype(str).str.upper()

def strat_perm(d38, d293, pct38, pct293, mn38, mn293, nperm=NPERM, seed=SEED, stat=np.mean):
    d = np.concatenate([d38, d293]); pct = np.concatenate([pct38, pct293]); mn = np.concatenate([mn38, mn293])
    lab = np.concatenate([np.ones(len(d38)), np.zeros(len(d293))]).astype(bool)
    qb = np.quantile(pct, [.25, .5, .75]); qm = np.quantile(mn, [.25, .5, .75])
    cellb = np.digitize(pct, qb) * 4 + np.digitize(mn, qm)
    cells = {c: np.flatnonzero(cellb == c) for c in np.unique(cellb)}
    rng = np.random.default_rng(seed)
    obs = stat(d[lab]) - stat(d[~lab]); null = np.empty(nperm)
    for k in range(nperm):
        pl = lab.copy()
        for c, idx in cells.items():
            if len(idx) > 1:
                pl[idx] = rng.permutation(pl[idx])
        null[k] = stat(d[pl]) - stat(d[~pl])
    p = (np.sum(null >= obs) + 1) / (nperm + 1)
    return obs, null.mean(), null.std(), (obs - null.mean()) / null.std(), p

for tag in ("rna", "atac"):
    pctM = M[tag + "_pct"].values      # 6 亚型平均检测率（与主检验同一口径）
    meanM = M[[tag + "_" + s for s in SUB]].mean(1).values
    d = M[tag + "_d"].values
    if os.environ.get("P2_CONTRAST") == "mol6":   # 与主图同口径：Δ = MOL6 − COP/NFOL
        d = M[tag + "_" + SUB[5]].values - M[tag + "_" + SUB[0]].values
    base = np.isfinite(d) & np.isfinite(meanM) & (pctM >= 10)
    i38 = np.flatnonzero(base & M.gu.isin(C38).values)
    i293 = np.flatnonzero(base & M.gu.isin(P331 - C38).values)
    print("\n=== %s ===  38→%d, 293→%d" % (tag.upper(), len(i38), len(i293)))
    o, nm, ns, z, p = strat_perm(d[i38], d[i293], pctM[i38], pctM[i293], meanM[i38], meanM[i293])
    print("R0 基线：差 %+0.4f  零%+.4f±%.4f  z=%+.2f  P=%.4f" % (o, nm, ns, z, p))
    # R1 留一基因
    worst, wname = (p, "")
    for j in i38:
        k = np.array([x for x in i38 if x != j])
        _, _, _, _, pp = strat_perm(d[k], d[i293], pctM[k], pctM[i293], meanM[k], meanM[i293], nperm=1500)
        if pp > worst:
            worst, wname = pp, M.gu.iloc[j]
    print("R1 留一基因：最坏 P=%.4f（去掉 %s 后失去显著？%s）" % (worst, wname, "是" if worst > 0.05 else "否"))
    # R2 检测率门槛
    for thr in (30, 50):
        a = i38[pctM[i38] >= thr]; b = i293[pctM[i293] >= thr]
        if len(a) >= 8:
            o2, _, _, z2, p2 = strat_perm(d[a], d[b], pctM[a], pctM[b], meanM[a], meanM[b])
            print("R2 检测≥%d%%：n=%d vs %d  差%+.4f  z=%+.2f  P=%.4f" % (thr, len(a), len(b), o2, z2, p2))
    # R3 中位数（稳健）
    o3, _, _, z3, p3 = strat_perm(d[i38], d[i293], pctM[i38], pctM[i293], meanM[i38], meanM[i293], stat=np.median)
    print("R3 中位数差（稳健）：%+.4f  z=%+.2f  P=%.4f" % (o3, z3, p3))
