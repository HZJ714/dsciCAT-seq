# Figure S3c（100k COP/NFOL「两簇」四格验证图）· 结果回执 2026-10-29

判据冻结文件：`S3C_COP2_FROZEN_CRITERIA_20261029.md`（**先看规则后看数据**，本轮没有改判据）
对象：`<work>/pi/plots-seurat/fy14_cop_nfol_subset_reclustered.RDS`（133 核，res 0.5 → 81 / 52）

---

## 0. 一句话结论

**三条冻结判据里 S1 过、S2 与 S3 不过 ⇒ 按冻结时的约定，正文那句"100k 把 COP/NFOL 分成两群"
必须降级为描述性报告，不能命名成两个亚状态。**
但本轮最有价值的一条不是"没过"，而是：**把两组一起下采样到同一深度重跑同一条判据，
251 / 306 个独立差异基因仍然成立、方向 100% 不变** ⇒ 这个分裂**不是**纯深度假象；
过不了的是"它到底是状态还是捕获深度"这个问题在**单库单动物**里无法区分。

---

## 1. 三条判据的结果（数字全部现算，脚本见 §5）

| 判据 | 冻结要求 | 实测 | 判定 |
|---|---|---|---|
| **S1 独立基因** | 排除 4 个判据基因后，≥3 个基因满足 \|簇均值差\|≥0.5 且 **BH-FDR<0.05** | **260 个**（31,263 基因宇宙）；按 Seurat `p_val_adj`（其实是 Bonferroni 口径）是 146 个 | ✅ **PASS**（两种口径都远超 3） |
| **S2 分裂稳定性** | 扰动网格（5 分辨率 × 5 种子 × {全量, 80%×20 draws} = 525 次）里"恰好 2 簇且各 ≥10 核"复现率 ≥70% 且 ARI 均值 ≥0.70 | **41.9% / ARI 0.624** | ❌ **FAIL** |
| ↳ 同网格按分辨率拆开 | — | res 0.1：0/105；**res 0.3：105/105，ARI 0.94**；**res 0.5：105/105，ARI 0.93**；res 0.8：10/105；res 1.0：0/105 | 见 §3 的口径说明 |
| **S3 无技术混杂** | 9 个 QC 量每一个都要满足 BH≥0.05 **或** \|Cliff δ\|≤0.2 | 可评 6 项里 **2 项不过**：`nCount_RNA`（2,682 vs 4,292，BH 6.0e-4，δ=−0.40）、`nFeature_RNA`（1,740 vs 2,322，BH 2.1e-3，δ=−0.35）。ATAC 侧 `nCount_ATAC`/`nFeature_ATAC`、`percent.mt`、`RNA.weight` 全过（BH≥0.92，\|δ\|≤0.07） | ❌ **FAIL** |

**深度匹配重检**（看到 S3 挂了之后加的敏感性检验，方向写死：只能削弱结论）：
把两组**一起**多项式下采样到 2,682 UMI，再跑同一条 S1 判据（8,081 个 ≥10 核表达的基因宇宙）：
**306 → 251 个基因仍合格（保留 82%），方向一致率 100%，\|Δ\| 中位只缩到 0.95 倍。**
⇒ 差异不是深度的影子；但两组的深度差本身仍然存在，且在少突 lineage 里"越成熟 RNA 越多"是真的生物学，
所以 S3 这条挂得**并不代表分裂是假的**，只代表**在这份数据里分不清"状态"和"捕获深度"**。

---

## 2. C1ql1 那条批注（Codex #87：15-cell 候选被判失败，规则是不是事后挑的）

新脚本 `scripts/s3c_c1ql1_sensitivity.py` → `u13_cop_4k/s3c/s3c_c1ql1_sensitivity.csv`。
把冻结规则里的 `C1ql1 ≤ 0.2` 一项**整条拿掉**重算：

| 运行 | 簇数 | 带 C1ql1（冻结） | 去掉 C1ql1 |
|---|---|---|---|
| FY13_full（50k） | 12 | 1 | 1 |
| FY14_full（100k） | 13 | 1 | 1 |
| HTp_full（4 库合并） | 8 | 1 | 1 |
| **EL14_full（HT 单库）** | 6 | **0** | **1**（簇 4，n=15，C1ql1=0.432） |
| 全脑抽 4,029 核 × 20 draws | 7–10 | 7/20 | **7/20**（逐 draw 完全一致） |

⇒ **C1ql1 这一项只会"取消候选"，从来不会"制造候选"**：24 个运行里只有 1 个判定被它改变，
20 次下采样 draw 一次都没改变。这正是"规则预先冻结"该有的样子，
所以**不删基因、不删那个 15-cell 候选**，把这张对照表放进图注/回复即可。

---

## 3. 图（方案 A：独立小图，两版都给作者看效果）

| 文件 | 尺寸 | 内容 |
|---|---|---|
| `replot/S3c/FigS3c_full.pdf` | **18.00 × 9.75 cm** | A 26 行（21 独立基因 + 4 灰色判据基因）+ B/C/D + 底部判据条 |
| `replot/S3c/FigS3c_compact.pdf` | **18.00 × 7.75 cm** | 去掉 4 行灰色判据基因 + 判据条进图注，只剩四格 |

体检（`final_audit.py` / `inkcheck.py`）：两版**硬门全过** —— 最小可见字号 5.5 pt、
只有 ArialMT、无未嵌入、无 Type3、未子集化=无、隐藏字符 0、面板字母 ABCD 8 pt 常规（house rule）、
紧墨盒重叠 **0 对**。

四格读法（图注草稿见 §4）：
- **A** 点面积 = 表达百分比，颜色 = 该基因在 5 组间的 z 分数（coolwarm ±2.2）。
  sub-1 一侧 `Mobp/Cnp/Trf/Sox10/Gpc5` 偏高（更成熟），sub-0 一侧 `Sox6/Cspg4/Pdgfra` 偏高（更 OPC 样）；
  灰色 = 当初用来定义 COP/NFOL 的四个基因（不参与 S1）；淡灰 = 五组里都测不到（`Cdk1/Dcn/Fn1/Mki67/Myc/Pdgfrb/Top2a`，
  增殖/纤维化信号在两簇都为零 ⇒ **不是循环细胞假象**）。
- **B** 空心蓝 = 原样，实心红 = 两组一起下采样到同一深度；右沿数字 = 合格基因数 306 → 251。
- **C** 柱 = 该分辨率下"恰好 2 簇且各 ≥10 核"的运行百分比，菱形 = 与出图那次分裂的 ARI（×100）。
- **D** Cliff δ 的哑铃，灰带 = \|δ\|≤0.2，红 = BH<0.05 且 \|δ\|>0.2；右侧两栏是两组各自的中位数。

---

## 4. 图注草稿 —— 作者看过 v8 后重排为**三格 H/I/J**（`FigS3_v9.pdf`）

作者意见与落地：
| 意见 | 处理 |
|---|---|
| "至少得放 100k COP 再聚类的 UMAP，不然不明所以" | **新增 H**：50k 与 100k 的 COP/NFOL 用**同一套参数**再聚类的 UMAP 对，标 sub-0/sub-1 |
| "K 图没啥信息量" | **撤 K**（6 行 lollipop）。但两条数不能丢 ⇒ 中位 UMI 2,682 vs 4,292 直接写进 H 的注释，统计量进图注 |
| "I 图是啥意思？" | **撤 I**（深度匹配重检散点）。结论压成图注一句：同深度重检后 306 → 251 个基因仍合格、方向 100% 不变 |
| "气泡图可以去掉 oligo 那列？" | 照做，z 分数按余下 4 组重算。**代价见下面 ⚠** |
| "核心就是气泡图和 J 图" | I（气泡图）与 J（稳定性）保留，宽度重新分配 |

**H 这块的前提是实测过的，不是叙述出来的**（`scripts/s3c_cop_umap_pair.R`）：
50k 60 核 → **1 簇**（res 0.1–0.5 都不分，0.8 才分 32/28）；100k 133 核 → **81/52**，
与交付分区**完全一致（ARI = 1.000，两库都是）**。出图脚本里写成了断言：
同参数再聚类若复现不出交付分区，直接 `RuntimeError`，图不会生成。

> ⚠ **去掉 Oligo 列丢掉的正是"sub-1 像成熟少突"这个锚点**：`Mobp` sub-1 = 1.52，同库 Oligo = 1.50，
> sub-0 只有 0.17；`Mog` 0.48 vs 1.34、`Mal` 0.15 vs 1.80、`Trf` 0.45 vs 1.81。
> ⇒ 图注里补了一句给数（见下），或作者想留一列时把 `GRP` 加回 `'100k Oligo (same lib)'` 即可（z 会自动按 5 组重算）。

图注（接在 Figure S3 图注末尾，三格 = **H–J**；已按作者第二轮意见改：J 用图例不用文字、
H 的参数行放大到 6.2 pt、I 收成 18 行并分组、中位 UMI 与 pipeline 细节进图注）：

> **H–J**, the two sub-clusters seen inside the 100k COP/NFOL subset (the 133 nuclei circled in Figure 3C), tested against rules fixed before the data were seen. **H**, the COP/NFOL nuclei of both dsciCAT-seq libraries re-clustered with identical parameters — 3,000 highly variable genes, scaling, PCA, graph neighbours and Leiden clustering on dimensions 1–10, resolution 0.5, seed 1, i.e. the pipeline that produced the delivered partition, reproduced here at ARI = 1.00 for both libraries. The 50k library (60 nuclei, median 3,558 UMI) yields one cluster at these parameters and splits only when the resolution is raised to 0.8 (32/28); the 100k library (133 nuclei) yields sub-0 (n = 81, median 2,682 UMI) and sub-1 (n = 52, median 4,292 UMI). **I**, the frozen panel of markers that took no part in defining COP/NFOL, shown for the two sub-clusters, the 50k COP/NFOL pool and the OPCs of the same 100k library. Dot area gives the percentage of nuclei expressing, colour the mean log-normalised expression z-scored across these four groups (−2.2 to +2.2). The fourteen markers detected in at least one group are listed above the dashed line in order of decreasing expression; the four below it, in grey type, are undetectable in all four groups (mean < 0.1) and are kept as negative controls — proliferation (`Cdk1`, `Mki67`, `Top2a`) and matrix/fibroblast (`Dcn`) markers, i.e. neither sub-cluster is a cycling or stromal contaminant. Three further undetectable panel members (`Fn1`, `Myc`, `Pdgfrb`, mean < 0.05) are not plotted. The oligodendrocyte reference of the same library (n = 3,614) is likewise not shown: sub-1 reaches it for the myelin genes (`Mobp` 1.52 vs 1.50; `Mog` 0.48 vs 1.34; `Cnp` 0.71 vs 0.57) while sub-0 does not (`Mobp` 0.17, `Mog` 0.44), and sub-1 stays below it for `Mal` (0.15 vs 1.80). **J**, recovery of the two-cluster split across the frozen perturbation grid (dark bars, five whole-library runs per resolution; pale bars, one hundred runs at 80% downsampling), with adjusted rand index against the delivered partition (green diamonds, ×100) and the frozen acceptance line (red dashes). Three criteria were fixed in advance (thresholds derived from the FY13 cluster 29 and FY14 cluster 31 reference clusters): **S1** at least three markers excluded from the call differ by ≥0.5 at BH-FDR < 0.05 — **met**, 260 genes over the whole transcriptome (146 under the Bonferroni convention), of which three are in the panel above (`Gpc5`, `Sox6`, `Mobp`); **S2** the split recovered in ≥70% of the grid with mean ARI ≥ 0.70 — **not met** (41.9% / ARI 0.62 over the whole grid, but 210/210 runs with ARI 0.93 at resolution 0.3–0.5, the grid also counting resolutions where a two-cluster answer is not expected); **S3** no pre-specified quality metric differing between the sub-clusters at both BH < 0.05 and |Cliff's δ| > 0.2 — **not met** for `nCount_RNA` (2,682 vs 4,292, BH 6.0e-4, δ = −0.40) and `nFeature_RNA` (1,740 vs 2,322, BH 2.1e-3, δ = −0.35), while both ATAC metrics, `percent.mt` and the WNN RNA weight are clean (BH ≥ 0.92, |δ| ≤ 0.07). Re-sampling both sub-clusters to the same depth of 2,682 UMI leaves 251 of the 306 qualifying genes, all in the same direction, so the signature is not simply the depth difference; still, within one library from one animal the two groups differ in captured RNA, so they are reported descriptively and are not named as two substates.

---

## 5. 本轮修掉的两个**真错**（不是版式问题）

1. **BH-FDR 的秩与 p 错配**（`scripts/s3c_depth_matched_recheck.py`）。
   旧写法 `rk[order] = arange(...)` 得到的是"按原始位置索引的秩"，却又拿它去乘**已排序**的 `p[m][order]`
   ⇒ 相当于把 FDR 随机打乱。实测 `Prom1` 真 p = 4.3e-21 却被算成 q = 0.099。
   修好后 as-is 从 **62 → 306**，深度匹配后从 **36 → 251**（结论方向不变，但之前那对数是错的）。
   已加 `bh()` 函数 + 与 statsmodels `multipletests(method='fdr_bh')` 逐值对账的自检 assert。
2. **S1 的口径**：Seurat `FindMarkers` 返回的 `p_val_adj` 是 `p_val × 检验基因数`（Bonferroni），
   而冻结判据写的是 **BH-FDR**。两个口径都算并写进图/表（260 vs 146），**不**只报好看的那个。
3. （上一轮直接崩掉的原因，记一笔免得再踩）`merge` 之后列名带后缀 ⇒ 是 `mean_diff_a/mean_diff_m`，
   不是 `diff_a/diff_m`。

---

## 6. 脚本与产物

| 路径 | 作用 |
|---|---|
| `scripts/s3c_cop2_analysis.R` | A 格 marker 表 + B 格 Seurat DE（S1 主口径） |
| `scripts/s3c_fix_stability_qc.R` | C 格扰动网格（与出图脚本同一套多分辨率调用） |
| `scripts/s3c_fix_qc.R` | D 格技术量（FetchData 取数值 + 正确的 Cliff δ / AUC / BH） |
| `scripts/s3c_export_counts.R` | 导出 133 核 raw counts + 元信息，供深度匹配 |
| `scripts/s3c_depth_matched_recheck.py` | 深度匹配重检（**已修 BH**） |
| `scripts/s3c_score_stability.py` | S2 打分（主口径=冻结原文的整个网格，res 匹配档只作次要参考值） |
| `scripts/s3c_c1ql1_sensitivity.py` | §2 的 C1ql1 拿掉重算 |
| `scripts/s3c_cop2_figure.py` | 出图（`--compact` 出 7.75 cm 版） |
| `u13_cop_4k/s3c/*.csv` | 全部中间表（含负结果） |

---

## 7. 需要作者点头的（**2026-10-29 已点头，1/2/3/4 已落地，见 §8**）

1. **正文 [72] 那句**现在写 "the higher-throughput library resolved internal substructure within COP/NFOL
   (two subclusters of 81 and 52 cells)"。按冻结约定应降级为描述性；我建议改成
   "…resolved two reproducible sub-clusters of 81 and 52 nuclei that differ in a myelination-associated
   gene programme which survives depth matching, but that also differ in nuclear RNA capture depth within
   a single library and animal, so they are reported descriptively rather than as two substates."
   —— **等作者点头再动 docx**。
2. **Figure 3 图注 [1303]** 现在那句"they are defined by the embedding and by Gpr17, not by an
   independent marker set"是本轮分析**之前**写的，现在已被 S1 推翻（有 260 个独立基因）。
   建议改成引用 Figure S3c 的口径：有独立基因支持、但深度混杂 ⇒ 仍只作描述。
3. **版式**：S3c 是独立小图（推荐，两版尺寸都在上面），还是硬塞进 Figure S3 底部（页高 34.56 → 42.3/44.3 cm）？
4. **点径回退那条（上一份回执 §11）已按"外科补丁"路线解决**：
   `FigS3_v7.pdf` = 以已放大的交付件为基准，只换四行字（F 标题 / 两处 separable 措词 / `fails C1ql1`→`C1ql1 not met`），
   逐像素差异 0.166% 且**全部落在四个标签框内**，挖掉标签区后各区块墨量与母版逐块相同 ⇒ 点径不可能变。
   细节与三条机制坑见 `FIGURES_CODEX_ROUND2_20261029.md` §12.2。**是否用它替换定稿，等作者点头。**
5. `FigS3c_full.pdf` 要不要顺带把 S3 页里的 G 面板灰色判据基因行去掉（避免同一批基因画两遍）？
   现在两版都保留，去留由作者定。

---

## 8. 定稿落地记录（作者："可以定稿了目前，记得一并修改好图注"）

| 上面第几条 | 结果 |
|---|---|
| 1 正文 [72] | 已改，v6。措辞按实测收紧：H 图注里那句"50k 只有把分辨率提到 0.8 才分开"、"260 个独立基因（面板内 3 个）"、"41.9% 网格复现"、"2,682 vs 4,292 UMI"全部进正文一句，结论=**只作描述、不命名为两个亚态** |
| 2 图注 [1303] | 已改：删掉"not by an independent marker set"（被 S1 推翻），改成"有独立 marker 支持（S3I），但 S2/S3 不过 ⇒ 仍只作描述" |
| 3 版式 | 按作者意见并入 Figure S3 底部 = **H–J 三格**（页高 42.88 cm），独立小图版留在 `replot/S3c/` 不再对外 |
| 4 替换定稿 | 已替换 `FINAL_SET_20261021/Figure_S3.pdf`（旧件进 `_superseded/`），S3b 同步 E–J，PNG/md5/GENERATED_TABLE/块内 figures.docx 一起同步 |
| 5 G 面板灰色判据行 | 未动（H–J 是 marker 面板 I，与 S3 的 G 面板不重复同一批基因的画法：G 是 7 个谱系 marker × 12 簇，I 是 21 个独立 marker × 4 组） |

**还没改、需要作者单独点头的一条**：正文 [72] 前半句 "Gpr17 feature plots **confirmed** …
**bona fide** … **not annotation artifacts**" 与 Figure 3 图注（v5 起写 "does not by itself exclude OPC"）
口径不一致。本轮只按 S3c 结论改了后半句，前半句未动。
