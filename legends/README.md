# legends/ — 定稿图注（2026-11-05 同步）

投稿可直接复制的图注正文。来源：另一环境 `<work>/figures/`，
2026-11-05 由本仓库统一收录。**原则（作者口径）：图注只写图面看不到的细节；结果解读一律在正文。**

| 文件 | 内容 | 权威源 |
|---|---|---|
| `figS6_v10_legend_final.md` | Supplementary Figure S6（v10 合成页 A–G，12 段 journal-ready）；草稿在 `SUPPORT_v2/figS6_v10_legend_draft.md`（英文段与定稿 diff=0） | `SUPPORT_v2/figS6_v10_LEGEND_FINAL.md` |
| `fig6_v21_legend_final.md` | Figure 6 主图 + Fig S6g（支持 6C）+ Fig S6h（支持 6D）+ Fig S6m（逐性状散点） | `routeA_out/FIG6_v21/fig6_v21_legend_final.md` |
| `fig7_v13_legend_final.md` | Figure 7 主图（v13，A/B/C 三格，7.93×16.90 cm）+ Fig S7 支持说明 | `routeA_out/FIG7_v13/fig7_v13_legend_final.md` |
| `RESULTS_TEXT_FIG6_FIG7_20261105.md` | Fig 6 / Fig 7 对应正文解读段（替换稿件相应段落用） | 同上目录 |
| `figure_legends.md`（master） | 全图图注 master 文件（Fig1–S5 + Fig6/S6 历史段） | `online/figures/figure_legends.md` |

✅ **硬伤 #1–#3 已修**（2026-11-05，6 文件 26 处替换，见 `GITHUB_UNIFY_RECEIPT_20261105.md` §3a）：

1. **编号同步**：v13 图注 6 处、正文解读段 7 处（3 处 "S7d"→Table S8、3 处 "Fig S7b"→S7a、
   1 处 →S7c;Table S7-8）、fig7/README 2 处、Figure_7.ipynb 4 处、build 脚本 4 处、
   assembly 注释 1 处——S7 v3 字母 = **S7a 去混杂 / S7b 三轴账本 / S7c 双组学**。
2. **"chromosome-position control (36→38)" 列名误读已更正**：图注 2 处 + 正文 1 处 →
   "the snATAC chromatin-channel rank (MS 36 → 38)"（= GSE231584 snATAC 通道 Δ1 名次；
   全仓无染色体分布性状表，此前的"染色体位置对照"说法作废）。
3. **小鼠 EAE 组成校正成对披露已补**：图注 (A) 总体句改为
   "Δ = −0.419, δ = −0.744, P = 0.0092; composition-adjusted Δ = −0.338, P = 0.387, n.s."。

⚠ 注意：**那边环境目录只读**，其原始文件仍是旧文本；若那边再动图注/notebook，须先照搬本仓
`legends/` 与 `fig7/` 的 diff（见回执 §3a 清单），否则下次同步会把这些修复覆盖回去。

✅ **S6 编号已定（2026-11-11，随稿件 v16 重排）**：表号统一为 **Table S1–S8**，
其中 **eRSS 表 = Table S5**（旧 S6）、跨通道一致性台账 = **Table S6**（旧 S7）、逐鼠文库账本 = **Table S7**、
疾病/病灶层 = **Table S8**（新增，截断交付）。Supplementary Figure S6 仍是 **7 格合成页（A–G，17.00 × 18.86 cm）**，
与 Table S6 是两回事；原先"合成页 G 格 vs Fig6 图注 Fig S6g"的字母撞车，按稿件口径解决：
支持图统一称 `Figure S6A–G`，Fig 6 的支持内容落在 Table S5/S6 与 `Figure S7`。
唯一口径表：`dev/TABLES.md`；推送前检查 (`dev/push_preflight.sh` G9) 会扫旧说法。

⚠ 本目录的 `fig6_v21_legend_final.md` / `fig7_v13_legend_final.md` 是 2026-11-05 的快照，
其中 "Fig S6g / S6h / S6m" 一类字母仍是旧口径 → 需从稿件当前版回改（并行环境工单 §3）。

---

## 2026-11-17 定稿收录（fig6 v36 / figS6 v16 / fig7 v17）

| 文件 | 内容 | 权威源 |
|---|---|---|
| `fig6_v36_legend_final.md` | Figure 6 主图（v36 三格横排 17.00 × 4.89 cm：A score / B UMAP / C cross-channel ladder；Manhattan 删、cross-chem 移 S6G） | 本仓 `scripts/CTM` 真源 FIG6（2026-11-17 重写） |
| `figS6_v16_legend_final.md` | Figure S6（v16 四行七格 17.00 × 18.91 cm：A loci / B permutation / C UMAP bins / D top-1000 landing / E cis-window landing / F per-nucleus boxes / G cross-chemistry Δ1 pair；D/E 面上共享竖版图例） | 同上 FIGS6 |
| `fig7_v17_legend_final.md` | Figure 7（v17 双栏 17.00 × 14.34 cm：A donors / B size-deconfounded ranks / C three axes / D drug-target × cell-type chip，含供体级 S1PR5 复检句） | `fig7/LEGEND_FIG7_final.md`（v17 定稿） |
| `figS6_v11_legend_final.md`（= 包内 `figS6/LEGEND_S6_v11.md` 时代） | 历史：v11 版 S6 图注，被 v16 取代 | `figS6/` |

## 2026-11-18 定稿收录（figS7 v29）

| 文件 | 内容 | 权威源 |
|---|---|---|
| `figS7_v29_legend_final.md` | Figure S7（v29 定稿 17.00 × 14.36 cm：A 双基 rank ladder / B tight pool / C human same-donor / D mouse same-nucleus / E 74×21 ledger / F detection map / G donor-level 22/27；v22 的 C/F/G 撤入 Table S7） | `figS7/LEGEND_S7_final.md` §1（v29） |

- **Fig S7 图注已收录**（v29 定稿，2026-11-18）：v22（A–I）被 v29（A–G）整图取代；旧 S7D→C、S7E→D、S7H→E、S7I→F，
  旧 C/F/G 撤入 Table S7；正文/Fig 7 图注的字母指针、Δ2B 定义（neuron − non-neuron）、Fig 7D 供体句（22/27 +
  8 patients + 3 controls 诊断披露）已同步。
- 已写入投稿管线真源 `scripts/CTM/ctm_legends_text.py`（新增 FIGS7 + LEGENDS/FIGNAME/FIG_FILES 加 S7）并重跑
  `make_ctm_docs.py`：`Manuscript_CTM.docx` 与 `Figure_Legends_Fig4_S4_6_S6_7_S7.docx` 已再生成，
  `check_docx_internal` 未标记汉字 = 0 ✓；正文编辑版 `Manuscript_CTM_v18_revised_20261118.docx` 见
  `<work>/figures/`（回执 `FIGS7_V29_MANUSCRIPT_20261118.md`）。
- `fig6_v21_legend_final.md` / `fig7_v13_legend_final.md` / `figS6_v10_legend_final.md` / `figS7_v29_legend_final.md`
  为历史快照，保留作 provenance，不再随稿；**2026-11-25 起也不入 git**（只留在本工作树，口径见
  `../dev/TABLES.md` §Version policy）。本目录随稿的是 `fig6_v36_legend_final.md`、
  `figS6_v16_legend_final.md`、`fig7_v17_legend_final.md` 与 `figS7/LEGEND_S7_final.md`（v30）。
