#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S6v8_derived_stats.py — fig S6 v8 的双轨统计（先在派生表里落数，画图脚本一律只读表）

产出（全部写 SUPPORT_v2/tables_v8/，每行带 source 列供反向溯源）：
  1) figS6v8_group_assoc_74x6x2.csv   scDRS **原生**组级 assoc_mcp / assoc_mcz（q95 对照法，n_ctrl=1000），
     74 性状 × 6 状态 × 2 通道。公式与 scripts/MS_GWAS/u2_assoc_matrix.py 逐字相同（该实现已对 12/12
     交付 group 文件校验过），此处只是把通道从 RNA 扩到染色质，并对每个 (通道,状态) 家族做 BH。
  2) figS6v8_delta1_bootstrap_74x2.csv  Δ1 = median(MOL6) − median(COP/NFOL) 的点估计 + 亚型内细胞
     bootstrap（B=1000，seed 固定）的 95% CI 与符号一致率；同时给宽对比 Δ1b。
  3) figS6v8_trait_rank_dualtrack.csv   同一通道内三种性状排名的对照：
     Δ1（我们的对比）、assoc_mcz(MOL6)（原文口径的组级富集）、Δz = z(MOL6) − z(COP/NFOL)（把原文统计量放进对比形式的我们自己构造的版本）。
  4) figS6v8_sixpoint_all74.csv         6 点状态剖面 ρ（O1 同核，两通道），补 BH q 与名次 —— 关掉旧图的 "other 71 traits NOT DELIVERED"。
  5) figS6v8_ladder_variants.csv        梯子 ρ 的全部口径（raw / depth-resid / n_hit-adj / 宽对比 / 分库 / D-own 另一套定义），只从交付表读，不重算。
  6) figS6v8_D_replication.csv          修正后 D 车道：跨臂 ρ（RNA/ATAC 各两两）+ 梯子 rho_raw/rho_adj/CI/天花板 + offdiag z。
  7) figS6v8_mouse.csv                  18 只鼠逐鼠两通道对比与 split-half AUC（**不算跨模态幅度比**）。

红线（来自 METHODS_HARD_AUDIT_v1 / SCDRS_DELTA_REASONING_v1）：
  * 染色质通道不用 scDRS 的 p 做推断；此处仍算出 p/z，但用途是**排名对照**，不是显著性断言。
  * 不比较跨模态的 Δ 幅度；只报秩/方向。
  * Δ1 的显著性来自我们自建的 bootstrap/置换，不借 scDRS 细胞级 p。
用法：python3 S6v8_derived_stats.py [--quick]     --quick = 只做 6)7) 与 5)（不读 full_score）
"""
import os, sys, glob, gzip
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

FY = WORK + "/FY/oligo"
ROUTE = WORK + "/figures/routeA_out"
CH_DIRS = {'RNA': f'{FY}/scdrs_traits_draft/full74',
           'ATAC(gene score)': f'{FY}/fig6_rerun_20261019/dualomics_probe/out_gs74'}
OUTD = f'{ROUTE}/SUPPORT_v2/tables_v8'
os.makedirs(OUTD, exist_ok=True)
STATES = ['COP/NFOL', 'MFOL', 'MOL2', 'MOL3', 'MOL3/6', 'MOL6']
MAT, PRE = 'MOL6', 'COP/NFOL'          # O1 主口径
WBROAD_M, WBROAD_P = ['MOL6', 'MOL2'], ['MFOL', 'COP/NFOL']   # 宽对比 Δ1b
B = int(os.environ.get('BOOT', 1000)); SEED = 20261104
QC = '--quick' in sys.argv


def SAY(*a):
    print(a[0] % a[1:] if len(a) > 1 else a[0], flush=True)


# ---------------- 标签 ----------------
ctab = pd.read_csv(f'{FY}/scdrs_ms_out/cell_to_type.csv')
lab = ctab.set_index('cell_id')['cell_type']
SAY('细胞标签（scDRS 打分集）: %s，共 %d 核', lab.value_counts().to_dict(), len(lab))


def load_full(p):
    with gzip.open(p, 'rt') as f:
        fs = pd.read_csv(f, sep='\t', index_col=0)
    ctrl = fs[[c for c in fs.columns if c.startswith('ctrl_norm_score')]].values.astype(np.float64)
    return fs['norm_score'].values.astype(np.float64), np.asarray(fs.index), ctrl


def assoc_q95(norm, ctrl, cell_lab):
    """scDRS 组级 assoc（q95 对照法）：与 u2_assoc_matrix.py 逐字同式"""
    out = {}
    for st in STATES:
        m = (cell_lab == st)
        if m.sum() < 5:
            out[st] = (np.nan, np.nan, int(m.sum()))
            continue
        q95 = np.quantile(norm[m], 0.95)
        cq95 = np.quantile(ctrl[m], 0.95, axis=0)
        p = ((cq95 >= q95).sum() + 1) / (cq95.shape[0] + 1)
        z = (q95 - cq95.mean()) / cq95.std()
        out[st] = (p, z, int(m.sum()))
    return out


def boot_delta(v_mat, v_pre, rng):
    """亚型内细胞 bootstrap：两侧各自重抽核，再取中位差。返回 (点, lo, hi, frac>0, B)"""
    obs = np.median(v_mat) - np.median(v_pre)
    bm = rng.integers(0, len(v_mat), size=(B, len(v_mat)))
    bp = rng.integers(0, len(v_pre), size=(B, len(v_pre)))
    d = np.median(v_mat[bm], axis=1) - np.median(v_pre[bp], axis=1)
    return obs, float(np.quantile(d, .025)), float(np.quantile(d, .975)), float((d > 0).mean())


# ---------------- 5) 梯子口径（只读交付表） ----------------
T4 = pd.read_csv(f'{ROUTE}/O1_own_multiome_paired/tables/T4_delta1_74traits_both_channels.csv')
T4b = pd.read_csv(f'{ROUTE}/O1_own_multiome_paired/tables/T4b_cross_channel_rho_variants.csv')
LIB = pd.read_csv(f'{ROUTE}/O1_own_multiome_paired/tables/L1_by_library.csv')
DLAD = pd.read_csv(f'{ROUTE}/D_crosschemistry/tables_log1p_20261104/fig6D_crosschem_ladder.csv')
DTECH = pd.read_csv(f'{ROUTE}/D_crosschemistry/tables_log1p_20261104/fig6D_delta1_74traits_per_tech.csv')
rows = []
for _, r in T4b.iterrows():
    rows.append(dict(family='O1 同核跨通道梯子', variant=r['variant'], rho=round(float(r.rho_cross_channel), 4),
                     P_spearman=float(r.P), perm_null_975=round(float(r.null_975), 4), perm_P=float(r.P_perm),
                     unit='74 traits', n=74, definition='Δ1 = MOL6 − COP/NFOL' if r['variant'] != 'sens: (MOL6+MOL2)-(MFOL+COP)' else 'Δ1b = (MOL6+MOL2) − (MFOL+COP/NFOL)',
                     source='O1/tables/T4b_cross_channel_rho_variants.csv'))
for _, r in LIB.iterrows():
    rows.append(dict(family='O1 分库复现', variant=r['lib'], rho=round(float(r.rho_xchan), 4), P_spearman=np.nan,
                     perm_null_975=round(float(r.null975), 4), perm_P=float(r.P_perm), unit='74 traits', n=74,
                     definition='Δ1 = MOL6 − COP/NFOL (within one library)',
                     source='O1/tables/L1_by_library.csv'))
for _, r in DLAD.iterrows():
    rows.append(dict(family='D 同臂跨通道（另一套定义）', variant=r['tech'] + ' rho_raw', rho=round(float(r.rho_raw), 4),
                     P_spearman=np.nan, perm_null_975=np.nan, perm_P=np.nan, unit='74 traits', n=74,
                     definition='Δ1 = Oligo − OPC', source='D/tables_log1p_20261104/fig6D_crosschem_ladder.csv'))
    rows.append(dict(family='D 同臂跨通道（另一套定义）', variant=r['tech'] + ' rho_adj', rho=round(float(r.rho_adj), 4),
                     P_spearman=np.nan, perm_null_975=np.nan, perm_P=float(r.perm_P_adj) if 'perm_P_adj' in r else np.nan,
                     unit='74 traits', n=74, definition='Δ1 = Oligo − OPC, n_hit 残差化',
                     source='D/tables_log1p_20261104/fig6D_crosschem_ladder.csv'))
pd.DataFrame(rows).to_csv(f'{OUTD}/figS6v8_ladder_variants.csv', index=False)
SAY('[5] ladder_variants: %d 行', len(rows))


# ---------------- 6) D 复制数（修正表） ----------------
TECHS = DTECH.tech.unique().tolist()
rows = []
for ch in ['rna', 'atac']:
    for i in range(len(TECHS)):
        for j in range(i + 1, len(TECHS)):
            a, b = TECHS[i], TECHS[j]
            pa = DTECH[DTECH.tech == a].set_index('trait')['d1_' + ch]
            pb = DTECH[DTECH.tech == b].set_index('trait')['d1_' + ch]
            pa, pb = pa.align(pb, join='inner')
            rho, pv = spearmanr(pa.values, pb.values)
            rows.append(dict(family='跨制库化学 同通道（同一性状 Δ1 两两一致）', channel=ch.upper(), pair=f'{a} vs {b}',
                             stat='Spearman rho of Δ1', value=round(rho, 4), P=pv, n=len(pa),
                             definition='Δ1 = Oligo − OPC', source='D/tables_log1p_20261104/fig6D_delta1_74traits_per_tech.csv'))
lad = pd.read_csv(f'{ROUTE}/D_crosschemistry/tables_log1p_20261104/fig6D_crosschem_ladder.csv')
for _, r in lad.iterrows():
    for k, v in r.items():
        rows.append(dict(family='同臂跨通道梯子（预注册判据）', channel='RNA+ATAC', pair=r['tech'], stat=k,
                         value=v, P=np.nan, n=int(r.get('n_traits', 74)) if 'n_traits' in r else 74,
                         definition='Δ1 = Oligo − OPC', source='D/tables_log1p_20261104/fig6D_crosschem_ladder.csv'))
MS = pd.read_csv(f'{ROUTE}/D_crosschemistry/tables_log1p_20261104/fig6D_MS_row.csv')
for _, r in MS.iterrows():
    for k in ['d1_rna', 'd1_atac', 'rank_rna_raw', 'rank_atac_raw', 'rank_rna_adj', 'rank_atac_adj', 'offdiag_z']:
        if k in r:
            rows.append(dict(family='D MS 行（性状级）', channel='RNA+ATAC', pair=r['tech'], stat=k,
                             value=r[k], P=np.nan, n=74, definition='Δ1 = Oligo − OPC',
                             source='D/tables_log1p_20261104/fig6D_MS_row.csv'))
pd.DataFrame(rows).to_csv(f'{OUTD}/figS6v8_D_replication.csv', index=False)
SAY('[6] D_replication: %d 行 | 跨臂 RNA/ATAC rho = %s', len(rows), [x['value'] for x in rows if x['stat'] == 'Spearman rho of Δ1'])


# ---------------- 7) 小鼠（不算幅度比） ----------------
MP = f'{ROUTE}/M2b_authoranno_v1/tables/M2b_per_mouse.csv'
if os.path.exists(MP):
    M = pd.read_csv(MP)[['mouse', 'group', 'n_MOL', 'n_OPC', 'auc_test', 'p_perm',
                         'atac_axis_contrast', 'rna_axis_contrast']].rename(columns={
        'auc_test': 'split_half_AUC', 'p_perm': 'perm_P',
        'atac_axis_contrast': 'ATAC_axis_contrast', 'rna_axis_contrast': 'RNA_axis_contrast'})
    MSRC = 'M2b_authoranno_v1/tables/M2b_per_mouse.csv'
else:
    M = pd.read_csv(f'{ROUTE}/SUPPORT_v2/tables/Table_S7.csv'); MSRC = 'SUPPORT_v2/tables/Table_S7.csv'
rows = []
for _, r in M.iterrows():
    for k in ['ATAC_axis_contrast', 'RNA_axis_contrast', 'split_half_AUC', 'perm_P']:
        if k in r:
            rows.append(dict(mouse=r['mouse'], group=r['group'], stat=k, value=r[k], source=MSRC))
sign = (np.sign(M.RNA_axis_contrast) == np.sign(M.ATAC_axis_contrast)).sum()
rho_mb, p_mb = spearmanr(M.RNA_axis_contrast, M.ATAC_axis_contrast)
rows.append(dict(mouse='ALL', group='ALL', stat='n_mice', value=len(M), source='count'))
rows.append(dict(mouse='ALL', group='ALL', stat='same_sign_fraction', value=round(sign / len(M), 4), source='computed here, direction only'))
rows.append(dict(mouse='ALL', group='ALL', stat='cross_channel_rho_per_mouse', value=round(rho_mb, 4), source='computed here'))
rows.append(dict(mouse='ALL', group='ALL', stat='cross_channel_perm_P_delivered', value=0.049, source='M2b delivered (perm of mouse labels, 10000)'))
pd.DataFrame(rows).to_csv(f'{OUTD}/figS6v8_mouse.csv', index=False)
SAY('[7] mouse: n=%d, 同号 %d/%d, 跨通道 rho=%.3f （**未计算任何跨模态幅度比**）', len(M), sign, len(M), rho_mb)

# ---------------- 8b) 推断单位梯（回应审计 P0-1 伪重复：图上每一级都写明单位与 n） -------------
rows = []
RHO_NUC = float(pd.read_csv(f'{ROUTE}/O1_own_multiome_paired/tables/T2_MS_per_nucleus_rna_vs_atac.csv')
                .query("trait.str.contains('Multiple_sclerosis')").rho_raw.iloc[0])
rows.append(dict(level='nucleus', unit='1 nucleus', n=int(os.environ.get('N_NUC', 5962)), rho=round(RHO_NUC, 4),
    P_or_perm=np.nan, inferential='NO - descriptive only (nuclei are not independent sampling units)',
    quantity='per-nucleus RNA vs ATAC rho, MS',
    source='O1_own_multiome_paired/tables/T2_MS_per_nucleus_rna_vs_atac.csv'))
for _, r in LIB.iterrows():
    rows.append(dict(level='library', unit='1 library (independent replicate)', n=int(r.n_cells),
                     rho=round(float(r.rho_xchan), 4), P_or_perm=float(r.P_perm),
                     inferential='YES - permutation over traits, one library = one independent replicate',
                     quantity='cross-channel ladder rho, within one library',
                     source='O1_own_multiome_paired/tables/L1_by_library.csv'))
MX = pd.read_csv(f'{ROUTE}/M2b_authoranno_v1/tables/M2b_crosschannel.csv').iloc[0]
rows.append(dict(level='mouse', unit='1 mouse', n=int(len(M)), rho=round(float(MX.spearman_rho), 4),
                 P_or_perm=float(MX.perm_P), inferential='YES - permutation over mouse labels',
                 quantity='cross-channel rho of per-mouse axis contrasts',
                 source='M2b_authoranno_v1/tables/M2b_crosschannel.csv'))
raw = T4b[T4b.variant == 'raw'].iloc[0]
broadv = T4b[T4b.variant.str.startswith('sens')].iloc[0]
rows.append(dict(level='trait', unit='1 trait (dsciCAT-seq nuclei)', n=74, rho=round(float(raw.rho_cross_channel), 4),
                 P_or_perm=float(raw.P_perm), inferential='YES - permutation over trait labels',
                 quantity='ladder rho of Delta1 (MOL6 - COP/NFOL)',
                 source='O1_own_multiome_paired/tables/T4b_cross_channel_rho_variants.csv'))
rows.append(dict(level='trait', unit='1 trait (dsciCAT-seq nuclei)', n=74, rho=round(float(broadv.rho_cross_channel), 4),
                 P_or_perm=float(broadv.P_perm), inferential='YES - permutation over trait labels',
                 quantity='ladder rho of broad contrast Delta1b (MOL6+MOL2 - MFOL-COP/NFOL)',
                 source='O1_own_multiome_paired/tables/T4b_cross_channel_rho_variants.csv'))
DLAD_OWN = DLAD[DLAD.tech.astype(str).str.contains('dsciCAT')].iloc[0]
rows.append(dict(level='trait', unit='1 trait (dsciCAT-seq, Oligo-OPC definition)', n=int(DLAD_OWN.get('n_traits', 74)),
                 rho=round(float(DLAD_OWN.rho_adj), 4), P_or_perm=float(DLAD_OWN.perm_P_adj),
                 inferential='YES - permutation over trait labels',
                 quantity='ladder rho_adj, D lane same-arm cross-channel',
                 source='D_crosschemistry/tables_log1p_20261104/fig6D_crosschem_ladder.csv'))
pd.DataFrame(rows).to_csv(f'{OUTD}/figS6v8_unit_ladder.csv', index=False)
SAY('[8b] unit ladder: %s', [(r['level'], round(r['rho'], 3)) for r in rows])

# ---------------- 8c) 小鼠疾病对比的阶段分层置换（回应审计问题 4） ----------------
FA = pd.read_csv(f'{ROUTE}/FIG7_v11/tables/fig7A_donors_delta1.csv')
mm = FA[FA.study == 'M GSE250589'].copy()
mm['stage'] = mm.id.str.extract(r'^(Early|Peak|Late)')[0]
mm['arm'] = mm.id.str.extract(r'_(Ctrl|EAE)_')[0]


def cliff(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    return (np.sum(x[:, None] > y[None, :]) - np.sum(x[:, None] < y[None, :])) / (len(x) * len(y))


d_obs = cliff(mm[mm.arm == 'EAE'].d1.values, mm[mm.arm == 'Ctrl'].d1.values)
rng2 = np.random.default_rng(SEED)
null = []
for _ in range(10000):
    lab2 = []
    for st, g in mm.groupby('stage'):
        n1 = int((g.arm == 'EAE').sum()); vv = np.array([1] * n1 + [0] * (len(g) - n1))
        rng2.shuffle(vv); lab2.append(pd.Series(vv, index=g.index))
    l = pd.concat(lab2).sort_index()
    null.append(cliff(mm.d1.values[l == 1], mm.d1.values[l == 0]))
null = np.array(null)
P_strat = (1 + np.sum(np.abs(null) >= abs(d_obs))) / (1 + len(null))
rows8 = [dict(stat='Cliff delta pooled (as delivered in Fig 7A)', value=round(d_obs, 4), P_delivered_MWU=0.0092,
              n_EAE=int((mm.arm == 'EAE').sum()), n_Ctrl=int((mm.arm == 'Ctrl').sum()),
              note='stage-stratified permutation P below', source='computed here from FIG7_v11/tables/fig7A_donors_delta1.csv'),
         dict(stat='stage-stratified permutation P', value=P_strat, P_delivered_MWU=np.nan, n_EAE=10000, n_Ctrl=np.nan,
              note='labels permuted within Early/Peak/Late; ctrl per stage = 3/1/2', source='computed here')]
for st, g in mm.groupby('stage'):
    e, c = g[g.arm == 'EAE'].d1.values, g[g.arm == 'Ctrl'].d1.values
    rows8.append(dict(stat=f'within-stage median diff ({st})', value=float(np.median(e) - np.median(c)),
                      P_delivered_MWU=float(mannwhitneyu(e, c, alternative='two-sided').pvalue) if len(c) > 1 and len(e) > 1 else np.nan,
                      n_EAE=len(e), n_Ctrl=len(c), note='no test when a stage has <2 libraries', source='computed here'))
pd.DataFrame(rows8).to_csv(f'{OUTD}/figS6v8_mouse_stage_stratified.csv', index=False)
SAY('[8] 小鼠疾病对比：池化 δ=%+.3f → 阶段分层置换 P=%.4f | 各阶段库数 %s', d_obs, P_strat,
    {st: (int((g.arm == 'EAE').sum()), int((g.arm == 'Ctrl').sum())) for st, g in mm.groupby('stage')})
if QC:
    sys.exit(0)


# ---------------- 1)+2)+3)+4) 需要 full_score 的重计算 ----------------
TR = {}
for chn, root in CH_DIRS.items():
    files = {}
    for f in glob.glob(f'{root}/magma_10kb_1000.batch*/*.full_score.gz'):
        files[os.path.basename(f).replace('.full_score.gz', '')] = f
    TR[chn] = files
    SAY('%s: %d 个 full_score.gz', chn, len(files))
traits = sorted(set(TR['RNA']) & set(TR['ATAC(gene score)']))
SAY('两通道交集性状 = %d', len(traits))

assoc_rows, boot_rows, med_rows = [], [], []
rng = np.random.default_rng(SEED)
for chn in ['RNA', 'ATAC(gene score)']:
    for k, tr in enumerate(traits):
        v, ids, ctrl = load_full(TR[chn][tr])
        L = lab.reindex(ids).values
        a = assoc_q95(v, ctrl, L)
        for st in STATES:
            p_, z_, n_ = a[st]
            assoc_rows.append(dict(channel=chn, trait=tr, state=st, assoc_mcp=p_, assoc_mcz=z_, n_cells=n_))
        vm, vp = v[L == MAT], v[L == PRE]
        obs, lo, hi, fr = boot_delta(vm, vp, rng)
        broad = np.median(np.concatenate([v[L == s] for s in WBROAD_M])) - \
            np.median(np.concatenate([v[L == s] for s in WBROAD_P]))
        obsb, lob, hib, frb = boot_delta(np.concatenate([v[L == s] for s in WBROAD_M]),
                                         np.concatenate([v[L == s] for s in WBROAD_P]), rng)
        boot_rows.append(dict(channel=chn, trait=tr, definition='Δ1 = MOL6 − COP/NFOL', point=obs, ci_lo=lo, ci_hi=hi,
                              frac_pos=fr, n_mat=len(vm), n_pre=len(vp), B=B, seed=SEED))
        boot_rows.append(dict(channel=chn, trait=tr, definition='Δ1b = (MOL6+MOL2) − (MFOL+COP/NFOL)', point=obsb,
                              ci_lo=lob, ci_hi=hib, frac_pos=frb, n_mat=int(np.isin(L, WBROAD_M).sum()),
                              n_pre=int(np.isin(L, WBROAD_P).sum()), B=B, seed=SEED))
        med_rows.append(dict(channel=chn, trait=tr, **{('med_' + s): float(np.median(v[L == s])) if (L == s).sum() else np.nan for s in STATES}))
        if (k + 1) % 10 == 0:
            SAY('  %s: %d/%d traits done', chn, k + 1, len(traits))

A = pd.DataFrame(assoc_rows)
# BH 逐 (通道,状态) 家族（74 个 p）挂回表
A['q_BH_within_state'] = np.nan
for chn, g in A.groupby('channel'):
    for st, gg in g.groupby('state'):
        pv = gg.assoc_mcp.values.astype(float)
        _, qq, _, _ = multipletests(pv, method='fdr_bh')
        A.loc[gg.index, 'q_BH_within_state'] = qq
A['source'] = 'scDRS q95 group assoc, recomputed here from *.full_score.gz (公式同 scripts/MS_GWAS/u2_assoc_matrix.py，该实现已对 12/12 交付 group 文件校验)'
A.to_csv(f'{OUTD}/figS6v8_group_assoc_74x6x2.csv', index=False)
SAY('[1] group_assoc: %d 行', len(A))

BO = pd.DataFrame(boot_rows)
BO['source'] = 'cell bootstrap within subtype, B=%d seed=%d (computed here)' % (B, SEED)
BO.to_csv(f'{OUTD}/figS6v8_delta1_bootstrap_74x2.csv', index=False)
SAY('[2] bootstrap: %d 行', len(BO))

MED = pd.DataFrame(med_rows)
MED['source'] = 'subtype medians of scDRS norm_score, from *.full_score.gz (computed here)'
MED.to_csv(f'{OUTD}/figS6v8_state_median_74x6x2.csv', index=False)

# --- RNA 轨校验：必须复现已交付 assoc_matrix_74x6.tsv ---
REF = pd.read_csv(f'{FY}/fig6_rerun_20261019/assoc_matrix_74x6.tsv', sep='\t')
chk = A[A.channel == 'RNA'].pivot(index='trait', columns='state', values='assoc_mcz')
chk = chk.rename(columns={s: f'{s}_z' for s in STATES})
ref = REF.set_index('trait')[[f'{s}_z' for s in STATES]]
ref.index.name = 'trait'
d = np.nanmax(np.abs(chk.reindex(ref.index)[ref.columns].values - ref.values))
SAY('   校验 RNA assoc_mcz vs 已交付 assoc_matrix_74x6.tsv: 最大绝对差 = %.3e  %s', d, 'PASS' if d < 1e-9 else 'FAIL')

# --- 3) 双轨排名 ---
rk = []
for chn in ['RNA', 'ATAC(gene score)']:
    bo = BO[(BO.channel == chn) & (BO.definition.str.startswith('Δ1 ='))].set_index('trait')
    zM = A[(A.channel == chn) & (A.state == MAT)].set_index('trait').assoc_mcz
    zC = A[(A.channel == chn) & (A.state == PRE)].set_index('trait').assoc_mcz
    t = pd.DataFrame({'delta1': bo.point, 'd1_lo': bo.ci_lo, 'd1_hi': bo.ci_hi, 'frac_pos': bo.frac_pos,
                      'z_MOL6': zM, 'z_COPNFOL': zC, 'dz': zM - zC, 'q_MOL6': A[(A.channel == chn) & (A.state == MAT)].set_index('trait').q_BH_within_state})
    t['rank_delta1'] = t.delta1.rank(ascending=False).astype(int)
    t['rank_zMOL6'] = t.z_MOL6.rank(ascending=False).astype(int)
    t['rank_dz'] = t.dz.rank(ascending=False).astype(int)
    r1 = spearmanr(t.delta1, t.z_MOL6); r2 = spearmanr(t.delta1, t.dz); r3 = spearmanr(t.z_MOL6, t.dz)
    SAY('[3] %s: Δ1 vs z(MOL6) rho=%.3f P=%.2g | Δ1 vs Δz rho=%.3f P=%.2g | z vs Δz rho=%.3f',
        chn, r1.correlation, r1.pvalue, r2.correlation, r2.pvalue, r3.correlation)
    for nm in ['PASS_Multiple_sclerosis', 'PASS_Schizophrenia_Pardinas2018', 'UKB_460K.body_HEIGHTz']:
        SAY('     %-32s Δ1 #%d (CI %+.3f..%+.3f, frac>0 %.2f) | z(MOL6) #%d q=%.3f | Δz #%d', nm,
            t.loc[nm, 'rank_delta1'], t.loc[nm, 'd1_lo'], t.loc[nm, 'd1_hi'], t.loc[nm, 'frac_pos'],
            t.loc[nm, 'rank_zMOL6'], t.loc[nm, 'q_MOL6'], t.loc[nm, 'rank_dz'])
    for pair, rr in [('rank_delta1 vs rank_zMOL6', r1), ('rank_delta1 vs rank_dz', r2), ('rank_zMOL6 vs rank_dz', r3)]:
        rk.append(dict(channel=chn, comparison=pair, stat='Spearman rho over 74 trait ranks', value=round(rr.correlation, 4),
                       P=rr.pvalue, n=74, source='computed here'))
    t.assign(channel=chn).to_csv(f'{OUTD}/figS6v8_trait_rank_%s.csv' % ('rna' if chn == 'RNA' else 'atac'), index_label='trait')
pd.DataFrame(rk).to_csv(f'{OUTD}/figS6v8_trait_rank_dualtrack.csv', index=False)

# --- 4) 六点剖面 ρ（用本脚本自算的状态中位，两通道对齐 T4） ---
pm = MED[MED.channel == 'RNA'].set_index('trait')[['med_' + s for s in STATES]]
pa = MED[MED.channel == 'ATAC(gene score)'].set_index('trait')[['med_' + s for s in STATES]]
sp = []
for tr in pm.index:
    rho, pv = spearmanr(pm.loc[tr].values.astype(float), pa.loc[tr].values.astype(float))
    sp.append(dict(trait=tr, rho_6pt=rho, P_6pt=pv, n=6))
SP = pd.DataFrame(sp)
_, SP['q_BH_74'], _, _ = multipletests(SP.P_6pt.values, method='fdr_bh')
SP['rank_rho_desc'] = SP.rho_6pt.rank(ascending=False).astype(int)
SP['source'] = '6 state medians per channel (computed here); 与 O1/tables/T4 同源同定义'
SP.to_csv(f'{OUTD}/figS6v8_sixpoint_all74.csv', index=False)
SAY('[4] sixpoint: 中位 rho=%+.3f, rho>0 %.0f%% , P<0.05 未校正 %d/74, BH q<0.05 %d/74',
    SP.rho_6pt.median(), 100 * (SP.rho_6pt > 0).mean(), (SP.P_6pt < .05).sum(), (SP.q_BH_74 < .05).sum())

# 与已交付 T4 的 Δ1 对齐校验（同一批核同一式子应当逐位相同）
T4i = T4.set_index('trait')
bR = BO[(BO.channel == 'RNA') & (BO.definition.str.startswith('Δ1 ='))].set_index('trait').point
bA = BO[(BO.channel == 'ATAC(gene score)') & (BO.definition.str.startswith('Δ1 ='))].set_index('trait').point
SAY('校验 vs 交付 T4：RNA 最大差 %.3e | ATAC 最大差 %.3e',
    np.nanmax(np.abs(bR.reindex(T4i.index).values - T4i.d1_rna.values)),
    np.nanmax(np.abs(bA.reindex(T4i.index).values - T4i.d1_atac.values)))
rr, pp = spearmanr(bR.reindex(T4i.index).values, bA.reindex(T4i.index).values)
SAY('梯子复现：本脚本自算 Δ1 跨通道 rho = %.4f（交付 0.8191）', rr)
SAY('全部产出写入 %s', OUTD)
