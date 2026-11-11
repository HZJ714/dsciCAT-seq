# -*- coding: utf-8 -*-
"""b_new_panels.py — 新面板所需的全部计算（① ② ③ + Fig7D 靶点矩阵 + S7C 留一法 + S7D 组成）。

原则：只用盘上已有数据；每个输出一律 long-format CSV + 一句 provenance；
算不出的（如 71 个性状的 split-half）不硬造，写进 MISSING 清单。
产出 → figures/routeA_out/NEWPANELS_v1/tables/
用法：/opt/conda/bin/python3 b_new_panels.py
"""
import os, re, gzip, glob, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr, fisher_exact, mannwhitneyu

WORK = os.environ.get('DSCI_WORK_ROOT', '')
OUT = f"{WORK}/figures/routeA_out/NEWPANELS_v1/tables"
os.makedirs(OUT, exist_ok=True)
H5 = f"{WORK}/FY/oligo/3k_dsci_oligo_subset.h5ad"
GS = f"{WORK}/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs"
GSPLAIN = f"{WORK}/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs"
PLAT = f"{WORK}/FY/oligo/scdrs_ms_out/genes_of_interest.txt"
GTF = f"{WORK}/FY/oligo/scdrs_ms_out/gencode.v45.basic.annotation.gtf.gz"
BUB = f"{WORK}/figures/S8_bubble/tables"
BED = f"{WORK}/FY/oligo/multivelo/peaks/by_celltype"
OL_FC = f"{WORK}/figures/routeA_out/A2_GSE277430_v1/tables/OL_lineage_log2FC_mature_vs_precursor.csv"
log = lambda *a: print(*a, flush=True)
MISS = []

SUBTYPES = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]   # 与 Fig6A 亚型轴一致
BETA = r"PASS_Multiple_sclerosis"

# ---------------------------------------------------------------- 载入自有 multiome
def load_own():
    import scanpy as sc
    a = sc.read_h5ad(H5)
    X = a.X
    mx = float(X.max()) if hasattr(X, "max") else np.nanmax(X)
    log("[own] shape=%s type=%s X.max=%.2f var_names[:4]=%s obs=%s"
        % (X.shape, type(X).__name__, mx, list(a.var_names[:4]), list(a.obs.columns)[:6]))
    return a, X, ("log1p" if mx < 12 else "counts")


def _densify(X, rows, cols):
    """只在需要的子块上转稠密，避免 5,977 x 31,263 全矩阵进内存。"""
    sub = X[np.ix_(rows, cols)]
    return np.asarray(sub.todense()) if hasattr(sub, "todense") else np.asarray(sub)

# ---------------------------------------------------------------- ① 风险基因 × 6 亚型
def riskgene_subtype_heatmap(a, X, scale):
    """行 = MS 基因集里在 OL 谱系内变动的基因（按 |MOL6-COP| 表达差排序取 48）；
       RNA 块 = 各亚型平均表达；ATAC 块 = 各亚型 peak bed 落在 TSS±100 kb 的 peak 数。"""
    genes = [g.strip().upper() for g in open(PLAT) if g.strip()]
    plat = set(genes)
    # 亚型标签（obs 里可能叫 cell_type，取值可能是 SUBPAL 名或近似）
    ct = a.obs["cell_type"].astype(str)
    mapping = {}
    for s in SUBTYPES:
        hit = ct[ct.str.upper().str.replace("-", "", regex=False).str.replace("/", "", regex=False) == s.upper().replace("/", "").replace("-", "")]
        mapping[s] = list(hit.index)
    if sum(len(v) for v in mapping.values()) == 0:
        log("[①] obs cell_type 取值=%s → 需映射表" % dict(ct.value_counts().head(8))); MISS.append("① subtype mapping")
        return None
    var = {g.upper(): i for i, g in enumerate(a.var_names)}
    idx = [var[g] for g in genes if g.upper() in var]
    gsub = [genes[i] for i in range(len(genes)) if genes[i].upper() in var]
    idx = [var[g.upper()] for g in gsub]
    means, counts = {}, {}
    for s, cells in mapping.items():
        rows = [i for i, n in enumerate(a.obs_names) if n in set(cells)]
        blk = _densify(X, rows, idx)
        means[s] = (np.log1p(blk).mean(0) if scale == "counts" else blk.mean(0))
        counts[s] = len(rows)
    M = pd.DataFrame(means, index=gsub)
    M["n_in_MS_set"] = 1
    M["abs_gap"] = (M["MOL6"] - M["COP/NFOL"]).abs() if "MOL6" in M and "COP/NFOL" in M else M.iloc[:, :6].max(1) - M.iloc[:, :6].min(1)
    top = M.sort_values("abs_gap", ascending=False).head(48)
    rna = top[[c for c in SUBTYPES if c in top.columns]].round(4)
    rna.insert(0, "gene", rna.index)
    rna["mol6_minus_cop"] = (top["MOL6"] - top["COP/NFOL"]).round(4).values
    rna["nuclei_per_subtype"] = ";".join("%s=%d" % (s, counts[s]) for s in SUBTYPES if s in counts)
    rna.to_csv(f"{OUT}/panel6D_riskgene_subtype_RNA.csv", index=False)
    # ATAC 块：peak bed ∩ TSS±100 kb
    tss = {}
    with gzip.open(GTF, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            nm = re.search(r'gene_name "([^"]+)"', c[8])
            if nm and nm.group(1).upper() in set(g.upper() for g in rna.gene):
                tss.setdefault(nm.group(1).upper(), (c[0].replace("chr", ""), int(c[3]) if c[6] == "+" else int(c[4])))
    rows = []
    if os.path.isdir(BED):
        peaks = {}
        for bf in sorted(glob.glob(BED + "/*.bed")):
            st = os.path.basename(bf)[:-4].replace("-6", "/6").replace("MOL3/6", "MOL3/6")
            st = {"COP-NFOL": "COP/NFOL", "MFOL": "MFOL", "MOL2": "MOL2", "MOL3": "MOL3",
                  "MOL3-6": "MOL3/6", "MOL6": "MOL6"}.get(os.path.basename(bf)[:-4], st)
            arr = []
            for ln in open(bf):
                p = ln.split("\t")
                if len(p) >= 3:
                    arr.append((p[0].replace("chr", ""), int(p[1]), int(p[2])))
            peaks[st] = arr
        for g in rna.gene:
            if g.upper() not in tss:
                rows.append(dict(gene=g, **{s: np.nan for s in peaks}, note="no TSS in GTF basic")); continue
            ch, ts = tss[g.upper()]
            d = {}
            for s, arr in peaks.items():
                d[s] = sum(1 for (c, lo, hi) in arr if c == ch and lo - 100000 <= ts <= hi + 100000)
            rows.append(dict(gene=g, **d, note=""))
        at = pd.DataFrame(rows)
        at.to_csv(f"{OUT}/panel6D_riskgene_subtype_ATAC_peakcount.csv", index=False)
        log("[①] RNA 块 %d 基因 × %d 亚型；ATAC 块 %d 亚型 peak bed（核数 %s）；TSS 命中 %d/%d"
            % (len(rna), len(SUBTYPES), len(peaks), rna.nuclei_per_subtype.iloc[0], at[list(peaks)].notna().sum().sum(), len(at) * len(peaks)))
    else:
        log("[①] 无 by_celltype peak bed → 只交付 RNA 块"); MISS.append("① ATAC block")
    return rna

# ---------------------------------------------------------------- ② 基因级贡献分解
def contribution(a, X, scale):
    """w_g = MAGMA 基因水平 z（scDRS 输入基因集原样）；
       分解口径 = w_g × 该基因逐核 z 分数在 MOL6 与前体之间的均值差 / Σ|w|。
       恒等式只对"重构 score"成立，故同表交付重构 score 与交付 norm_score 的 Spearman。"""
    line = None
    with open(GS) as f:
        for ln in f:
            if ln.startswith(BETA + "\t"):
                line = ln.rstrip("\n").split("\t")[1]; break
    if line is None:
        log("[②] gs 文件里找不到 MS 行"); MISS.append("② weights"); return None
    w = {}
    for tok in line.split(","):
        if ":" in tok:
            g, z = tok.rsplit(":", 1)
            try: w[g.strip().upper()] = float(z)
            except ValueError: pass
    log("[②] MS 加权基因集 %d 个基因（|z| 中位 %.2f，最大 %.2f）" % (len(w), np.median(list(map(abs, w.values()))), max(w.values())))
    ct = a.obs["cell_type"].astype(str)
    names = list(a.obs_names)
    pos = {g.upper(): i for i, g in enumerate(a.var_names)}
    use = [(g, pos[g]) for g in w if g in pos]
    log("[②] 与矩阵交集 %d/%d" % (len(use), len(w)))
    E = _densify(X, list(range(X.shape[0])), [i for _, i in use])
    E = np.log1p(E) if scale == "counts" else E
    zE = (E - E.mean(0)) / (E.std(0, ddof=0) + 1e-9)
    grp = {}
    for s in SUBTYPES:
        rows = [k for k, n in enumerate(names) if str(ct.iloc[k]).upper().replace("-", "").replace("/", "") == s.upper().replace("-", "").replace("/", "")]
        if rows: grp[s] = rows
    if not grp: MISS.append("② subtypes"); return None
    mat = "MOL6" if "MOL6" in grp else max(grp, key=lambda k: len(grp[k]))
    pre = "COP/NFOL" if "COP/NFOL" in grp else min(grp, key=lambda k: len(grp[k]))
    dz = zE[grp[mat]].mean(0) - zE[grp[pre]].mean(0)
    wv = np.array([w[g] for g, _ in use])
    contrib = wv * dz / np.abs(wv).sum()
    ordi = np.argsort(-np.abs(contrib))[:20]
    recon = zE @ (wv / np.abs(wv).sum())
    rows = []
    for k in ordi:
        g = use[k][0]
        rows.append(dict(gene=g, magma_z=wv[k], expr_z_diff_mature_minus_precursor=round(float(dz[k]), 5),
                         contribution_to_score_difference=round(float(contrib[k]), 6),
                         pct_of_total_gap=round(100 * contrib[k] / contrib.sum(), 2) if contrib.sum() != 0 else np.nan,
                         direction="higher in %s" % mat if dz[k] > 0 else "higher in %s" % pre))
    t = pd.DataFrame(rows)
    t.insert(0, "mature_subtype", mat); t.insert(1, "precursor_subtype", pre)
    t.insert(2, "sum_top20_contrib", round(float(contrib[ordi].sum()), 6))
    t.insert(3, "sum_all_contrib", round(float(contrib.sum()), 6))
    t.insert(4, "reconstructed_vs_delivered_rho", np.nan)
    # 与交付的逐核 scDRS 分数对照（若可对齐 cell_id）
    cl = f"{WORK}/figures/routeA_out/O1_own_multiome_paired/tables/cell_level_table.csv"
    if os.path.exists(cl):
        c = pd.read_csv(cl).set_index("cell_id")
        inter = [n for n in names if n in c.index]
        if inter:
            rr = spearmanr(recon[[names.index(x) for x in inter]],
                           c.loc[inter, "%s_rna" % BETA].astype(float).fillna(0)).correlation
            t.loc[0, "reconstructed_vs_delivered_rho"] = round(float(rr), 4)
            log("[②] 重构 score 与交付 scDRS norm_score 的 Spearman = %.3f (n=%d) → 分解口径%s"
                % (rr, len(inter), "可用" if rr > 0.5 else "只能作示意，不可称 driver"))
    t.to_csv(f"{OUT}/panel6E_gene_level_contribution.csv", index=False)
    log("[②] top20 交付；Σ贡献=%.4f，全部基因 Σ=%.4f（恒等校验）" % (contrib[ordi].sum(), contrib.sum()))
    return t

# ---------------------------------------------------------------- ③ eRegulon × MS 风险基因
def eregulon(a):
    f = f"{WORK}/FY/oligo/scenic_tables/TableSX_eregulon_target_regions.summary.csv"
    d = pd.read_csv(f)
    plat = set(g.strip().upper() for g in open(PLAT) if g.strip())
    uni = set()
    # 全集：矩阵里出现的基因符号
    uni = set(g.upper() for g in a.var_names)
    rows = []
    for _, r in d.iterrows():
        gl = str(r.get("genes", ""))
        tg = set(x.strip().upper() for x in gl.replace(";", ",").split(",") if x.strip() and x.strip() != "nan")
        if not tg:
            tg = set()
        hit = len(tg & plat); n = len(tg & uni)
        if n < 5:
            continue
        a_ = hit; b_ = n - hit; c_ = len(plat & uni) - hit; dd = len(uni) - n - c_
        try:
            OR, p = fisher_exact([[a_, b_], [c_, dd]], alternative="greater")
        except Exception:
            continue
        rows.append(dict(eregulon=r["eRegulon"], TF=r["TF"], peak_state=r["peak_state"],
                         n_target_genes_in_matrix=n, n_MS_platform_targets=a_, OR=round(float(OR), 3), P=float(p),
                         pct_targets_that_are_MS=round(100 * a_ / n, 2)))
    t = pd.DataFrame(rows)
    if len(t) == 0:
        log("[③] 汇总表里没有 target gene 列 → 需从 SCENIC regulon 源文件取（列=%s）" % list(d.columns))
        MISS.append("③ target-gene lists absent in summary table")
        return None
    from statsmodels.stats.multitest import multipletests
    t["q_BH"] = multipletests(t.P, method="fdr_bh")[1]
    t = t.sort_values("P").round(6)
    t.to_csv(f"{OUT}/panelS6C_eregulon_MS_enrichment.csv", index=False)
    log("[③] %d 个 eRegulon；最小 P=%.3g，q<0.05 的个数=%d" % (len(t), t.P.min(), (t.q_BH < 0.05).sum()))
    return t

# ---------------------------------------------------------------- Fig7D 四套靶点矩阵
STATUS = {  # 只在 DailyMed/label 实测过的标 verified，其余一律 needs-label-verification（不凭记忆写获批）
    "siponimod": "verified (DailyMed MAYZENT setid 44492772-5aed-4627-bd85-e8e89f308bb3 v17)",
}
def target_matrix():
    frames = []
    f = f"{BUB}/A3_GSE118257_target_stats.csv"
    d = pd.read_csv(f); d.columns = [c.strip().lower() for c in d.columns]
    frames.append(d.rename(columns={"group": "celltype", "pct": "pct_nuclei_with_read", "mean_exp": "mean_log_detecting"})
                   .assign(dataset="GSE118257"))
    d = pd.read_csv(f"{BUB}/A2_GSE277430_target_stats.csv"); d.columns = [c.strip().lower() for c in d.columns]
    frames.append(d.rename(columns={"group": "celltype", "pct": "pct_nuclei_with_read", "mean_log": "mean_log_detecting"})
                   .assign(dataset="GSE277430"))
    d = pd.read_csv(f"{BUB}/MO_dsciCAT_target_stats.csv"); d.columns = [c.strip().lower() for c in d.columns]
    frames.append(d.rename(columns={"group": "celltype", "pct": "pct_nuclei_with_read", "mean_log": "mean_log_detecting"})
                   .assign(dataset="own_dsciCAT"))
    d = pd.read_csv(f"{BUB}/A4_GSE301908_target_stats.csv")
    gcol = [c for c in d.columns if c.lower() == "gene"][0]
    long = []
    for c in d.columns:
        if c.lower() == "gene" or c.endswith("_me"):
            continue
        long.append(d[[gcol, c]].rename(columns={c: "pct_nuclei_with_read"}).assign(celltype=c, dataset="GSE301908"))
    m = d[[gcol] + [c for c in d.columns if c.endswith("_me")]].melt(id_vars=gcol)
    m["celltype"] = m["variable"].str.replace("_me", "", regex=False)
    m = m.rename(columns={"value": "mean_log_detecting"})[[gcol, "celltype", "mean_log_detecting"]]
    a4 = pd.concat(long, ignore_index=True).merge(m, on=[gcol, "celltype"], how="left")
    a4 = a4.rename(columns={gcol: "gene"})
    frames.append(a4)
    all_ = pd.concat(frames, ignore_index=True, sort=False)
    all_["gene"] = all_["gene"].astype(str).str.upper()
    all_["celltype"] = all_["celltype"].astype(str)
    all_.to_csv(f"{OUT}/panel7D_target_detection_all.csv", index=False)
    common = all_.groupby("gene")["dataset"].nunique()
    four = sorted(common[common == 4].index)
    piv = all_.pivot_table(index=["dataset", "gene"], columns="celltype", values="pct_nuclei_with_read")
    piv.round(3).to_csv(f"{OUT}/panel7D_detection_freq_wide.csv")
    log("[Fig7D] 4 套合并 %d 行；数据集×靶点：%s；四套全覆盖靶点 %d 个：%s"
        % (len(all_), dict(all_.groupby("dataset")["gene"].nunique()), len(four), four))
    log("[Fig7D] 各套细胞类型数：%s" % {k: int(v) for k, v in all_.groupby("dataset")["celltype"].nunique().items()})
    return all_, four

# ---------------------------------------------------------------- S7C 逐供体留一法
def leave_one_out():
    f = f"{WORK}/figures/routeA_out/A2_GSE277430_v1/tables/celltype_donors_2.csv"
    d = pd.read_csv(f)
    rows = []
    for tr, g in d.groupby("trait"):
        piv = g.pivot_table(index="patient", columns="celltype", values="median_score")
        matc = [c for c in piv.columns if "mature" in c.lower()]
        prec = [c for c in piv.columns if "OPC" in c.upper() or "COP" in c.upper()]
        if not matc or not prec:
            continue
        delta = (piv[matc].mean(1) - piv[prec].mean(1))
        full = float(np.sign(np.median(delta)))
        for p in delta.index:
            loo = delta.drop(p)
            rows.append(dict(trait=tr, donor_left_out=p, n_donors=len(loo),
                             median_delta_full=round(float(np.median(delta)), 4),
                             median_delta_LOO=round(float(np.median(loo)), 4),
                             sign_flip=bool(np.sign(np.median(loo)) != np.sign(full)),
                             frac_delta_sign_changes=round(float((np.sign(loo) != np.sign(full)).mean()), 3)))
    t = pd.DataFrame(rows)
    t.to_csv(f"{OUT}/panelS7C_leave_one_donor_out_GSE277430.csv", index=False)
    if len(t):
        log("[S7C] %d 行（%d 性状 × %d 供体）；留一后变号的性状数=%d；至少一个供体剔除即变号的性状=%d"
            % (len(t), t.trait.nunique(), t.donor_left_out.nunique(),
               int(t.sign_flip.sum()), t[t.sign_flip].trait.nunique()))
        s = t.groupby("trait").sign_flip.sum().sort_values(ascending=False)
        log("[S7C] 最依赖单个供体的前 5 个性状：%s" % {k: int(v) for k, v in s.head(5).items()})
    else:
        MISS.append("S7C no mature/precursor columns")
    return t

# ---------------------------------------------------------------- S7D 疾病状态 × OL 亚型组成
def composition():
    rows = []
    # 人源自有：每供体亚型占比（A7 有 sample + celltype）
    f = f"{WORK}/figures/routeA_out/A7_GSE231585_rna_v1/tables/obs_qc_annot.csv"
    d = pd.read_csv(f)
    ct = d[[c for c in d.columns if c in ("sample", "celltype")]].dropna()
    tab = ct.pivot_table(index="sample", columns="celltype", values=None, aggfunc="size", fill_value=0)
    if tab.empty:
        g = ct.groupby(["sample", "celltype"]).size().unstack(fill_value=0)
    else:
        g = tab
    frac = (g.T / g.sum(1)).T
    for s in frac.index:
        for c in frac.columns:
            rows.append(dict(dataset="GSE231585", unit="sample", label=str(s), celltype=str(c),
                             n_nuclei=int(g.loc[s, c]), fraction=round(float(frac.loc[s, c]), 4)))
    # 自有小鼠：每鼠 MOL/OPC 比
    m = pd.read_csv(f"{WORK}/figures/routeA_out/M2b_authoranno_v1/tables/M2b_per_mouse.csv")
    for _, r in m.iterrows():
        tot = r.n_MOL + r.n_OPC
        rows.append(dict(dataset="own_EAE_mouse", unit="mouse", label="%s(%s)" % (r.mouse, r.group), celltype="MOL",
                         n_nuclei=int(r.n_MOL), fraction=round(float(r.n_MOL / tot), 4)))
        rows.append(dict(dataset="own_EAE_mouse", unit="mouse", label="%s(%s)" % (r.mouse, r.group), celltype="OPC/COP",
                         n_nuclei=int(r.n_OPC), fraction=round(float(r.n_OPC / tot), 4)))
    t = pd.DataFrame(rows)
    t.to_csv(f"{OUT}/panelS7D_composition_by_donor.csv", index=False)
    # MS vs Ctrl 的人源组成检验（GSE231585 若有分组列则分列，否则按 sample）
    ol = t[(t.dataset == "own_EAE_mouse")]
    eae = ol[ol.label.str.contains("Early|Peak|Late") & ~ol.label.str.contains("Ctrl")]
    ctrl = ol[ol.label.str.contains("Ctrl")]
    for ct_ in ["MOL", "OPC/COP"]:
        x = eae[eae.celltype == ct_].fraction.values; y = ctrl[ctrl.celltype == ct_].fraction.values
        if len(x) > 1 and len(y) > 1:
            U, p = mannwhitneyu(x, y, alternative="two-sided")
            log("[S7D] 鼠 %s 比例 EAE n=%d %.3f vs Ctrl n=%d %.3f  MWU P=%.4f" % (ct_, len(x), np.mean(x), len(y), np.mean(y), p))
    log("[S7D] %d 行（%d 个数据集单元 × 细胞类型）" % (len(t), t.label.nunique()))
    return t

if __name__ == "__main__":
    a, X, scale = load_own()
    riskgene_subtype_heatmap(a, X, scale)
    contribution(a, X, scale)
    eregulon(a)
    target_matrix()
    leave_one_out()
    composition()
    json.dump(dict(missing=MISS), open(f"{OUT}/_missing.json", "w"), indent=1)
    log("\n产出 %d 表 → %s\nMISS: %s" % (len(glob.glob(OUT + "/*.csv")), OUT, MISS))
