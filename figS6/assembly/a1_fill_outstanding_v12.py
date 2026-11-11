# -*- coding: utf-8 -*-
"""a1_fill_outstanding_v12.py — 关闭 Table_S6_v11 里 7 条 `OUTSTANDING - not delivered`。

做法：能重算的当场重算（gridsize ρ / per-mouse gap / 26 只 QC 名单），
已经算过但没并进台账的从 SUPPORT_v2/tables_v8 直接落成明细表（六点 ρ、rank 总账、
bootstrap 梯子、窗口命中、置换零分布）。所有产出写入 tables_v12/，**不回写 v8/v11 任何文件**。

七条 OUTSTANDING 的对应关系
  198 C matrices per-trait six-point rho (74)      -> detail_sixpoint_all74.csv
  199 per-trait split-half Delta1 pairs            -> 需逐核 74 性状打分；台账里给可算的 3 性状 + 明示其余
  200 rho between gridded UMAP views (15/25/40)    -> hexbin_gridsize_rho.csv   （重算）
  201 '% signal' per locus (13 values)             -> window_gene_hits_331.csv + window_per_chr_signal.csv（定义不可考，改交付可追溯量）
  202 raw permutation vector behind S6b            -> window_permutation_null_10000.csv（已在 v8，落明细 + 复算 P）
  203 the 8 mice excluded by the M2b rule          -> mice_all26_qc_rule.csv    （从 anno_min.pkl 重建）
  204 mean of the per-mouse channel gap            -> per_mouse_channel_gap.csv （重算）

用法：/opt/conda/bin/python3 a1_fill_outstanding_v12.py
"""
import os, sys, glob, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr

WORK = os.environ.get('DSCI_WORK_ROOT', '')
V8   = f"{WORK}/figures/routeA_out/SUPPORT_v2/tables_v8"
OUT  = f"{WORK}/figures/routeA_out/SUPPORT_v2/tables_v12"
os.makedirs(OUT, exist_ok=True)
D38 = f"{WORK}/figures/routeA_out/D_crosschemistry/tables"
O1  = f"{WORK}/figures/routeA_out/O1_own_multiome_paired/tables"
M2B = f"{WORK}/figures/routeA_out/M2b_authoranno_v1/tables/M2b_per_mouse.csv"
SEED = 20261104
log = lambda *a: print(*a, flush=True)

# ---------------------------------------------------------------- 198 六点 ρ（74 行全表）
def build_sixpoint():
    six = pd.read_csv(f"{V8}/figS6v8_sixpoint_all74.csv")                 # 74 × rho/P/q/rank
    sub = pd.read_csv(f"{D38}/fig6D_subtype_6pt.csv")                      # 74 × 6 亚型 × 2 通道均值
    keep = [c for c in sub.columns if c.startswith(("rna_", "atac_"))]
    d = sub[["trait", "rho_6pt", "P_6pt", "n_points"] + keep].merge(
        six[["trait", "q_BH_74", "rank_rho_desc"]], on="trait", how="left")
    med = pd.read_csv(f"{V8}/figS6v8_state_median_74x6x2.csv")             # 中位数口径（C 源）
    w = med.pivot(index="trait", columns="channel", values=[c for c in med.columns if c.startswith("med_")])
    w.columns = ["%s_%s" % (c[1], c[0]) for c in w.columns]
    d = d.merge(w.reset_index(), on="trait", how="left")
    d.insert(0, "n_traits", len(d))
    d.to_csv(f"{OUT}/detail_sixpoint_all74.csv", index=False)
    log("[198] 六点 ρ 全表 %d 行（median-of-6 口径 rho 中位 %.3f；subtype-mean 口径 rho 中位 %.3f）"
        % (len(d), d["rho_6pt"].median(), d["median_rna_rho_6pt"] if "median_rna_rho_6pt" in d else np.nan))
    return d

# ---------------------------------------------------------------- 200 gridsize 稳健性 ρ
def build_gridsize():
    """三档 hexbin 归约视图之间的 Spearman ρ：把每档的 hex 均值最近邻投影到最细网格后比较。"""
    import scanpy as sc
    from scipy.spatial import cKDTree
    import matplotlib
    matplotlib.use("Agg")
    H5 = f"{WORK}/FY/oligo/3k_dsci_oligo_subset.h5ad"
    SC = f"{WORK}/FY/oligo/scdrs_ms_out/PASS_Multiple_sclerosis.score.gz"
    if not (os.path.exists(H5) and os.path.exists(SC)):
        log("[200] 缺 H5AD 或 score.gz → 跳过"); return None
    a = sc.read_h5ad(H5)
    s = pd.read_csv(SC, sep="\t", index_col=0)
    um = a.obs[["UMAP_1", "UMAP_2"]] if "UMAP_1" in a.obs.columns else pd.DataFrame(
        a.obsm["X_umap"][:, :2], index=a.obs_names, columns=["UMAP_1", "UMAP_2"])
    df = um.join(s[["norm_score"]], how="inner").dropna()
    xy, v = df[["UMAP_1", "UMAP_2"]].values, df["norm_score"].values
    views = {}
    for g in (15, 25, 40):
        # 用六边形装箱把核分配到 hex 单元，取单元内均值（等价 hexbin reduce_C_function=np.mean）
        dx, dy = 1.0, np.sqrt(3) / 2
        col = np.floor(xy[:, 0] / (dx / g * 2)).astype(int) if False else None
        # 直接借用 matplotlib 的实现，保证与图上口径一致
        from matplotlib.figure import Figure
        fig = Figure(); ax = fig.add_subplot()
        hb = ax.hexbin(xy[:, 0], xy[:, 1], C=v, reduce_C_function=np.mean, gridsize=g, mincnt=1,
                       cmap="viridis")
        off = hb.get_offsets(); val = hb.get_array()
        m = np.isfinite(val) & (val != 0)
        views[g] = (off[m], val[m]); fig.clf()
    fine_off, fine_val = views[40]
    tree = {g: cKDTree(views[g][0]) for g in views}
    # 主口径：把每档的 hex 均值插值到同一规则栅格，再两两 Spearman（避免最近邻投影稀释一致性）
    from scipy.interpolate import griddata
    gx = np.linspace(xy[:, 0].min(), xy[:, 0].max(), 160)
    gy = np.linspace(xy[:, 1].min(), xy[:, 1].max(), 160)
    GX, GY = np.meshgrid(gx, gy)
    GXY = np.column_stack([GX.ravel(), GY.ravel()])
    grid = {}
    for g in views:
        grid[g] = griddata(views[g][0], views[g][1], GXY, method="linear")
    grid["single_cell"] = griddata(xy, v, GXY, method="nearest")
    common = np.isfinite(grid["single_cell"])
    for g in views:
        common &= np.isfinite(grid[g])
    proj_nn, proj_grid = {}, {}
    for g in views:
        d, i = tree[g].query(fine_off)
        proj_nn[g] = views[g][1][i]
        proj_grid[g] = grid[g][common]
    proj_nn["single_cell"] = np.array([v[np.argmin((xy[:, 0] - p[0]) ** 2 + (xy[:, 1] - p[1]) ** 2)]
                                      for p in fine_off])
    proj_grid["single_cell"] = grid["single_cell"][common]
    rows = []
    keys = ["single_cell", 15, 25, 40]
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            A_, B_ = keys[i], keys[j]
            rho, p = spearmanr(proj_grid[A_], proj_grid[B_])
            rnn, pnn = spearmanr(proj_nn[A_], proj_nn[B_])
            rows.append(dict(view_a=str(A_), view_b=str(B_),
                             n_common_pixels=int(common.sum()),
                             spearman_rho_raster=round(float(rho), 4), P_raster=float(p),
                             spearman_rho_nearesthex=round(float(rnn), 4),
                             n_hex_15=len(views[15][0]), n_hex_25=len(views[25][0]), n_hex_40=len(views[40][0]),
                             definition="hexbin mean per cell (mincnt=1); views resampled to a common "
                                        "160x160 raster over the UMAP extent (linear interpolation, "
                                        "single-nucleus view = nearest) before Spearman"))
    t = pd.DataFrame(rows)
    t.insert(0, "n_nuclei", len(df))
    t.to_csv(f"{OUT}/hexbin_gridsize_rho.csv", index=False)
    log("[200] 6 核对（栅格口径 | 最近hex口径）：15v25 %.3f|%.3f | 25v40 %.3f|%.3f | 单核v25 %.3f|%.3f"
        % (t[(t.view_a == "15") & (t.view_b == "25")].spearman_rho_raster.iloc[0],
           t[(t.view_a == "15") & (t.view_b == "25")].spearman_rho_nearesthex.iloc[0],
           t[(t.view_a == "25") & (t.view_b == "40")].spearman_rho_raster.iloc[0],
           t[(t.view_a == "25") & (t.view_b == "40")].spearman_rho_nearesthex.iloc[0],
           t[(t.view_a == "single_cell") & (t.view_b == "25")].spearman_rho_raster.iloc[0],
           t[(t.view_a == "single_cell") & (t.view_b == "25")].spearman_rho_nearesthex.iloc[0]))
    return t

# ---------------------------------------------------------------- 201 窗口命中（每染色体/每位点）
def build_window_signal():
    h = pd.read_csv(f"{V8}/figS6v10_window_gene_hits.csv")                 # 331 × gene/in_window/chr/tss/nearest_snp/dist_kb
    h.to_csv(f"{OUT}/detail_window_gene_hits_331.csv", index=False)
    h["hit"] = h["in_window"].astype(int)
    per_chr = h.assign(chr=h["chr"].astype(str)).groupby("chr").agg(
        platform_genes=("gene", "size"), hits_in_window=("hit", "sum"),
        median_dist_kb=("dist_kb", "median"), min_dist_kb=("dist_kb", "min"))
    per_chr["pct_of_platform"] = (100 * per_chr.hits_in_window / per_chr.platform_genes).round(2)
    per_chr = per_chr.reset_index()
    hh = h.dropna(subset=["chr", "nearest_snp", "dist_kb"]).copy()
    hh["chr"] = hh["chr"].astype(str)
    best = hh.loc[hh.groupby("chr")["dist_kb"].idxmin(), ["chr", "nearest_snp", "dist_kb", "gene"]]
    per_chr = per_chr.merge(best.rename(columns={"nearest_snp": "closest_lead_snp", "gene": "closest_gene",
                                                 "dist_kb": "closest_dist_kb"}), on="chr", how="left")
    per_chr.to_csv(f"{OUT}/window_signal_per_chr.csv", index=False)
    tot = dict(platform_genes=int(h.hit.sum() and len(h)), hits=int(h["hit"].sum()),
               pct=round(100 * h["hit"].sum() / len(h), 3), n_chr=int(h["chr"].nunique()))
    log("[201] 331 平台基因命中 %d 个 = %.1f%%（旧稿打印的 13 个 '%%signal' 值其定义与算法未随图交付，"
        "故此处交付可追溯的 per-chromosome 命中表）%s" % (h["hit"].sum(), 100 * h["hit"].sum() / len(h), tot))
    return per_chr, tot

# ---------------------------------------------------------------- 202 置换零分布
def build_permutation():
    nl = pd.read_csv(f"{V8}/figS6v8_window_permutation_null.csv")
    nl.to_csv(f"{OUT}/window_permutation_null_10000.csv", index=False)
    summ = pd.read_csv(f"{V8}/figS6v8_window_permutation.csv")
    summ.to_csv(f"{OUT}/window_permutation_summary.csv", index=False)
    obs = int(summ.loc[summ["stat"] == "platform genes in window (observed count)", "value"].iloc[0])
    ex = int((nl["null_hit_count"] >= obs).sum()); N = len(nl)
    P = (1 + ex) / (1 + N)
    hist = nl["null_hit_count"].value_counts().sort_index().rename_axis("null_hits").reset_index(name="draws")
    hist["pct_of_draws"] = (100 * hist.draws / N).round(3)
    hist.to_csv(f"{OUT}/window_permutation_null_histogram.csv", index=False)
    res = dict(observed_hits=obs, n_perm=N, exceeds_or_equal=ex, permutation_P_add_one=round(P, 7),
               null_mean=round(float(nl.null_hit_count.mean()), 4),
               null_ci95=[float(nl.null_hit_count.quantile(.025)), float(nl.null_hit_count.quantile(.975))],
               seed=SEED)
    log("[202] 零分布 %d draw 已落明细；obs=%d, >=obs %d → P=(1+%d)/(1+%d)=%.7f"
        % (N, obs, ex, ex, N, P))
    return res

# ---------------------------------------------------------------- 203 26 只 QC 名单
def build_mouse_qc():
    p = f"{WORK}/routeA_probe/mouse_EAE/anno_min.pkl"
    A = pd.read_pickle(p)
    mol = A[A.cellType.isin(["MOL1", "MOL2", "MOL56"])].groupby("orig.ident").size()
    opc = A[A.cellType.isin(["OPC", "COP"])].groupby("orig.ident").size()
    d = pd.DataFrame({"n_MOL": mol, "n_OPC": opc}).fillna(0).astype(int).reset_index()
    d["group"] = A.groupby("orig.ident")["model"].first().reindex(d["orig.ident"]).values
    d["sex"] = A.groupby("orig.ident")["Sex"].first().reindex(d["orig.ident"]).values
    d["batch"] = A.groupby("orig.ident")["batch"].first().reindex(d["orig.ident"]).values
    d["rule"] = "MOL >= 20 AND OPC >= 20 (M2b_formal_authoranno.py:107)"
    d["pass_M2b"] = (d.n_MOL >= 20) & (d.n_OPC >= 20)
    d["reason_excluded"] = np.where(d.pass_M2b, "-", np.where(d.n_MOL < 20, "MOL < 20", "OPC < 20"))
    kept = set(pd.read_csv(M2B)["mouse"])
    d["delivered_in_M2b_per_mouse"] = d["orig.ident"].isin(kept)
    d = d.sort_values("orig.ident")
    d.to_csv(f"{OUT}/mice_all26_qc_rule.csv", index=False)
    exc = d[~d.pass_M2b]
    identical = (set(d[d.pass_M2b]["orig.ident"]) == kept)
    log("[203] 26 只全集 / %d 只通过 / %d 只排除：%s | 排除原因：%s | 通过集与 M2b_per_mouse 完全一致：%s"
        % (d.pass_M2b.sum(), len(exc), ",".join(exc["orig.ident"]),
           sorted(set(exc.reason_excluded)), identical))
    return d

# ---------------------------------------------------------------- 204 per-mouse 通道差
def build_gap():
    d = pd.read_csv(M2B).copy()
    d["gap_atac_minus_rna"] = d.atac_axis_contrast - d.rna_axis_contrast
    d["abs_gap"] = d.gap_atac_minus_rna.abs()
    d[["mouse", "group", "n_MOL", "n_OPC", "atac_axis_contrast", "rna_axis_contrast",
       "gap_atac_minus_rna", "abs_gap"]].round(6).to_csv(f"{OUT}/per_mouse_channel_gap.csv", index=False)
    for tag, s in [("all18", d), ("EAE13", d[d.group != "Ctrl"]), ("Ctrl5", d[d.group == "Ctrl"])]:
        log("[204] %s n=%d mean gap=%.4f (median %.4f, SD %.4f, 95%%CI %.4f–%.4f)"
            % (tag, len(s), s.gap_atac_minus_rna.mean(), s.gap_atac_minus_rna.median(),
               s.gap_atac_minus_rna.std(ddof=1),
               s.gap_atac_minus_rna.mean() - 1.96 * s.gap_atac_minus_rna.std(ddof=1) / np.sqrt(len(s)),
               s.gap_atac_minus_rna.mean() + 1.96 * s.gap_atac_minus_rna.std(ddof=1) / np.sqrt(len(s))))
    return d

# ---------------------------------------------------------------- 199 split-half（能给多少给多少）
def build_splithalf():
    """逐核 74 性状打分不在交付树里（cell_level_table.csv 只有 MS/SCZ/HEIGHT）→
    只对可算的 3 个性状交付 per-trait split-half Δ1，其余 71 个如实标注需要重跑打分。"""
    c = pd.read_csv(f"{O1}/cell_level_table.csv")
    rng = np.random.default_rng(SEED)
    rows = []
    for tr in ["PASS_Multiple_sclerosis", "PASS_Schizophrenia_Pardinas2018", "UKB_460K.body_HEIGHTz"]:
        for ch in ["rna", "atac"]:
            r, at = c[f"{tr}_{ch}"].astype(float), c["cell_type"].astype(str)
            ok = np.isfinite(r)
            r, at = r[ok].values, at[ok].values
            mat = np.isin(at, ["MOL6", "MOL3/6"]); pre = np.isin(at, ["COP/NFOL", "OPC"])
            obs = np.median(r[mat]) - np.median(r[pre])
            vals = []
            for b in range(200):
                p = rng.random(r.shape) < 0.5
                m1, m0 = mat & p, mat & ~p; q1, q0 = pre & p, pre & ~p
                if m1.sum() < 5 or m0.sum() < 5 or q1.sum() < 5 or q0.sum() < 5:
                    continue
                a_ = np.median(r[m1]) - np.median(r[q1]); b_ = np.median(r[m0]) - np.median(r[q0])
                vals.append(spearmanr([a_, b_], [obs, obs]).correlation if False else (a_, b_))
            if vals:
                arr = np.array(vals)
                rho, pr = spearmanr(arr[:, 0], arr[:, 1])
                rows.append(dict(trait=tr, channel=ch, n_nuclei=int(ok.sum()), delta1_observed=round(float(obs), 5),
                                 B_split_half=len(arr), half1_vs_half2_rho=round(float(rho), 4),
                                 half1_P=float(pr),
                                 delta1_half_min=round(float(arr[:, 0].min()), 4), delta1_half_max=round(float(arr[:, 0].max()), 4),
                                 delta1_half_median=round(float(np.median(arr[:, 0])), 5),
                                 delta1_half_iqr="%0.3f–%0.3f" % (np.quantile(arr[:, 0], .25), np.quantile(arr[:, 0], .75)),
                                 seed=SEED, note="per-trait split-half on nuclei; 71 other traits need a scDRS re-score"))
    t = pd.DataFrame(rows)
    t.to_csv(f"{OUT}/splithalf_delta1_3traits.csv", index=False)
    log("[199] 3 性状 × 2 通道 split-half Δ1 已交付（%d 行）；其余 71 性状需要重跑逐核 scDRS 打分" % len(t))
    return t

# ---------------------------------------------------------------- rank 总账（S7F/G/H 的承载）
def build_rank_ledger():
    r = pd.read_csv(f"{V8}/figS6v8_trait_rank_rna.csv"); a = pd.read_csv(f"{V8}/figS6v8_trait_rank_atac.csv")
    r = r.drop(columns=[c for c in r.columns if r[c].nunique() == 1 and c in ("channel",)])
    a = a.drop(columns=[c for c in a.columns if c == "channel"] + [])
    wide = r.merge(a, on="trait", suffixes=("_rna", "_atac"))
    dt = pd.read_csv(f"{V8}/figS6v8_trait_rank_dualtrack.csv")
    wide.to_csv(f"{OUT}/rank_ledger_74traits_wide.csv", index=False)
    dt.to_csv(f"{OUT}/rank_ledger_dualtrack_summary.csv", index=False)
    log("[S7F/G/H] rank 总账 %d 行 × %d 列（双轨一致性 %d 行摘要：%s）"
        % (wide.shape[0], wide.shape[1], len(dt), "; ".join("%s=%.3f" % (x.comparison, x.value) for _, x in dt.iterrows())))
    return wide, dt

if __name__ == "__main__":
    res = {}
    build_sixpoint()
    res["gridsize"] = "done" if build_gridsize() is not None else "skipped"
    _, res["window"] = build_window_signal()
    res["permutation"] = build_permutation()
    d26 = build_mouse_qc()
    res["mice"] = dict(total=len(d26), passed=int(d26.pass_M2b.sum()),
                       excluded=list(d26[~d26.pass_M2b]["orig.ident"]))
    g = build_gap()
    res["gap"] = dict(mean_all18=round(float(g.gap_atac_minus_rna.mean()), 6),
                      mean_EAE=round(float(g[g.group != 'Ctrl'].gap_atac_minus_rna.mean()), 6),
                      mean_Ctrl=round(float(g[g.group == 'Ctrl'].gap_atac_minus_rna.mean()), 6))
    build_splithalf()
    build_rank_ledger()
    for f in ["figS6v8_delta1_bootstrap_74x2.csv", "figS6v10_unit_summary.csv",
              "figS6v10_bootstrap_units.csv", "figS6v8_ladder_variants.csv",
              "figS6v8_D_replication.csv", "figS6v8_state_median_74x6x2.csv",
              "figS6v8_mouse_stage_stratified.csv"]:
        src = f"{V8}/{f}"
        if os.path.exists(src):
            pd.read_csv(src).to_csv(f"{OUT}/{f.replace('figS6v8_', 'detail_').replace('figS6v10_', 'detail_')}", index=False)
    json.dump(res, open(f"{OUT}/A1_resolution_summary.json", "w"), indent=1, default=str)
    log("\n产出 %d 个文件 → %s" % (len(os.listdir(OUT)), OUT))
