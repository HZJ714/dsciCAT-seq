import os
# -*- coding: utf-8 -*-
"""p0_8c_power_reconcile.py —— 为什么 Fig6A/B 能出显著分数，而 top1000 基因层是阴性？

同一个效应量，两种抽样误差：
  · 细胞层（A/B 的做法）：把**基因集当作固定**，只让细胞抽样 ⇒ n=5,962，SE 极小
  · 基因层（我们的检验）：把**细胞集当作固定**，让基因抽样 ⇒ n=223，SE 大一个量级
本脚本用同一批细胞级数据把两个 SE 都算出来，直接对账。
"""
import numpy as np, pandas as pd
TAB = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'routeA_out/P0_8C/tables')
SUB = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
W = np.array([203, 490, 1016, 688, 1006, 2574], float)
fA = W[4] / (W[4] + W[5])          # 0.2810 = C1 富集池里 MOL3/6 的细胞占比

z = np.load(TAB + "/p08c_dual_celllevel223.npz", allow_pickle=True)
XG, AG, key, genes = z["XG"], z["AG"], z["key"], z["genes"]
rows = {s: np.flatnonzero(key == s) for s in SUB}
print("细胞 %d × 基因 %d | 富集权重 f = %.4f" % (XG.shape[0], XG.shape[1], fA))

def landing_cell(M, cols):                      # 细胞层：先对基因平均 → 每细胞一个分数
    S = np.nanmean(M[:, cols], 1)
    mA, mB = S[rows["MOL3/6"]].mean(), S[rows["MOL6"]].mean()
    m0 = S[rows["COP/NFOL"]].mean()
    vA = S[rows["MOL3/6"]].var(ddof=1) / len(rows["MOL3/6"])
    vB = S[rows["MOL6"]].var(ddof=1) / len(rows["MOL6"])
    v0 = S[rows["COP/NFOL"]].var(ddof=1) / len(rows["COP/NFOL"])
    return fA * (mA - m0) + (1 - fA) * (mB - m0), np.sqrt(fA**2 * (vA + v0) + (1 - fA)**2 * (vB + v0))

def landing_gene(M, cols):                      # 基因层：先对细胞类型求均值差 → 每基因一个 Δ
    with np.errstate(invalid="ignore"):
        Mn = np.column_stack([np.nanmean(M[rows[s]][:, cols], 0) for s in SUB])
    L = Mn - Mn[:, :1]
    return fA * L[:, 4] + (1 - fA) * L[:, 5]

rng = np.random.default_rng(20261116)
out = []
for tag, M in (("RNA", XG), ("ATAC", AG)):
    all_c = np.arange(M.shape[1])
    eff, se_cell = landing_cell(M, all_c)
    d = landing_gene(M, all_c)
    dg = d[np.isfinite(d)]
    se_gene = dg.std(ddof=1) / np.sqrt(len(dg))  # 解析式（基因独立近似）
    bs = [np.nanmean(landing_gene(M, rng.choice(all_c, len(all_c), replace=True))) for _ in range(2000)]
    se_boot = float(np.std(bs, ddof=1))          # 基因自助法（含基因间相关）
    # 二维匹配背景给出的"超额"（来自 matched2d 表，与图面同源）
    T = pd.read_csv(TAB + "/p08c_dual_setlevel_matched2d.csv")
    r = T[(T.channel == tag.lower()) & (T["set"] == "top1000可测223")].iloc[0]
    excess = r.mean_d - r.null_mean
    print("\n=== %s ===" % tag)
    print("  风险集着陆量（同一个数）  Δ̄ = %+.4f log2" % eff)
    print("  对匹配背景的超额（被检验效应）  E = %+.4f log2  = %.1f%% 相对差异" % (excess, (2**abs(excess) - 1) * 100))
    print("  ① 细胞抽样误差（A/B 口径，基因集固定）SE = %.5f → z = %.1f" % (se_cell, abs(excess) / se_cell))
    print("  ② 基因抽样误差（基因层口径）        SE = %.5f（解析） / %.5f（bootstrap） → z = %.2f / %.2f"
          % (se_gene, se_boot, abs(excess) / se_gene, abs(excess) / se_boot))
    print("  误差比 ②/① = %.1f 倍   ⇒ 同一个效应：细胞层 %+.1f\u03c3，基因层 %+.1f\u03c3"
          % (se_gene / se_cell, excess / se_cell, excess / se_gene))
    sd_g = dg.std(ddof=1)
    print("  每基因 Δ 的基因间 SD = %.4f；要在基因层达到 P<0.05 需 n ≈ %d 个可测基因（现有 %d）"
          % (sd_g, int(np.ceil((1.96 * sd_g / abs(excess))**2)), len(dg)))
    out.append(dict(channel=tag, effect_log2=eff, excess=excess, se_cell=se_cell, se_gene=se_gene,
                    se_gene_boot=se_boot, ratio=se_gene / se_cell, sd_per_gene=d.std(ddof=1),
                    n_have=len(dg), n_needed=int(np.ceil((1.96 * sd_g / abs(excess))**2))))
pd.DataFrame(out).to_csv(TAB + "/p08c_power_reconcile.csv", index=False)
print("\n[写出] p08c_power_reconcile.csv")
