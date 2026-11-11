# Results 正文解读稿 — Fig 6 / Fig 7（2026-11-05 整理版）

> 用途：可直接粘进 `Manuscript_CTM*.docx` Results 的英文正文段；**所有结果解读放这里**，
> 图注只保留图面看不到的细节（定义、参数、caveat）。整理依据：
> ① 旧总稿 `online/figures/DRAFT_FIG6S6_7S7_LEGENDS_RESULTS_20261103.md` §5（v0，数字已过期）；
> ② fig6 定稿 `FIG6_v21/fig6_v21_legend_final.md`（log1p 修正后数值，权威）；
> ③ fig7 定稿 `FIG7_v13/fig7_v13_legend_final.md` + 统一表 `seven4_confound_ALL_cohorts.csv`。
> 本轮相对 v0 草稿的**数字变化**（作者核对用）：
> · Fig 6 ladder 判据结论 由 "2 of 3 pass、HT 0.518 FAIL" 更新为修正后 **"3 of 3 pass"**
>   （0.896 / 0.824 / 0.630），MS 例外名次由 7/8、4/55、11/47 更新为 **7/37、4/49、11/41**；
>   桥接由 Δ = −0.273 更新为 **0.748 / 0.734 vs 0.819，Δ = −0.071（仍超 ±0.05 容差）**。
> · Fig 7 donor 级阳性由 "48/49, P = 8.9e-14"（旧 A4 donor 子集）更新为统一表实测
>   **54 of 55 units, binomial P = 1.6e-15**（唯一负向 = GSE277430 donor 2006-041，其 OPC 臂仅 53 核）。
> · Fig 7 MS 名次由 2/4/5/2 → 1/1/2/1（mean 1.2）更新为统一表 **2/4/7/2 → 1/1/3/1（mean 1.5）**；
>   confound ρ 更新为 **−0.84 / −0.69 / −0.78 / −0.87**；匹配窗更新为 **2/40、3/40、5/40、2/40**。
> · v0 的 Fig 7 D/E 格已退出主图（v11/v13 版式决定），其内容全部由本正文段 + Table S8 承载。

---

## 1. Fig 6 对应正文段（替换稿件 §3.9 的 Fig 6 段）

### 1.1 [P1 → Fig 6A, B] The axis in our own data

We first established the mature–precursor MS risk axis in our own dsciCAT-seq data. Applying scDRS
with the multiple-sclerosis (MS) GWAS gene set to the expression channel of adult mouse brain nuclei
(5,962 oligodendrocyte-lineage nuclei) raised the score from the precursor state COP/NFOL to the most
mature states MOL2 and MOL6 (Fig 6A); all five subtype-versus-COP/NFOL contrasts were significant at
uncorrected *P* ≤ 0.039, and the group-level permutation test survived Benjamini–Hochberg adjustment
over the family of five for MOL3/6 (*q* = 0.040) and MOL6 (*q* = 0.047). These are different tests and
are reported separately: at the level of individual nuclei, scDRS's own multiple-testing threshold was
met by **0 of 5,962 cells** (FDR < 0.05/0.1/0.2 all zero), which is why the score is used throughout as
a ranking statistic and never as a per-cell call. The gradient is continuous on the lineage UMAP
(Fig 6B) and is not produced by binning (smoothing controls in Fig S6b).

### 1.2 [P2 → Fig 6C, S6g] The same nuclei, two channels

In the same nuclei the identical contrast can be computed from the chromatin gene-score channel
(Fig 6C). The two channels rank the 74 GWAS traits consistently: Spearman ρ = **0.819**
(permutation *P* = 5.0×10⁻⁴, 10,000 permutations, seed 20261102), **0.823 [0.739–0.874]** after
residualising on gene-set regulatory-domain size, and 0.822–0.831 under alternative covariate
adjustments — **93% of the split-half noise ceiling of 0.878** (Fig S6g b). The chromatin channel is
therefore informative at the level of the trait ladder. Three MS-specific readings, however, cut
against a naïve "same signal in both channels" story and are reported together: (i) MS has the
largest expression-versus-chromatin off-diagonal *z* among the 74 traits (*z* = +2.09, rank 1/74) —
its chromatin placement is much weaker than its expression placement; (ii) the six-subtype profile of
MS is anticorrelated across channels (ρ = −0.143, *P* = 0.787; subtype-mean family — with subtype medians, the family consistent with Δ1, MS = −0.486, *P* = 0.329, Fig S6D), and this reading is source-dependent
(+0.257 on the matrices used in Fig 6D; Fig S6g d), so we do not treat it as established; (iii) at
the single-nucleus level the two channels barely co-vary (ρ = 0.079 for MS; 0.111 SCZ; 0.061 HEIGHT).
Absolute trait levels are unrelated between channels (ρ = 0.029, *P* = 0.81) whereas the contrast is
(ρ = 0.679 after removing trait level): **what the channels share is the contrast, not the level**.
Two method-independent scorers (AUCell, AddModuleScore; Table S6) detect the mature–precursor axis
in neither channel — they are ~10-fold less sensitive than scDRS here and therefore neither support
nor refute the scDRS gene-score reading. Because the own-data chromatin ranks are raw ranks within
the 74 sets, traits such as HEIGHT ("chromatin #1") are not read as trait-specific: under the
domain-adjusted ranks the same trait is 33/42 in our own data.

### 1.3 [P3 → Fig 6D, S6h] Cross-chemistry replication

We then asked whether the concordance is a property of our platform or of the biology, rerunning the
same pre-registered Δ1 cross-channel test on two further library chemistries from two separate
adult mouse-brain studies (10x Multiome; HT-scCAT-seq) through an identical pipeline (Fig 6D). The
ten libraries enter the analysis as four library-preparation groups, and **chemistry is inseparable
from source study** (each chemistry comes from a single preparation experiment); this constrains the
claim to "reproduced for each of three chemistries", not "replicated across independent
experiments", and is stated again in Limitations. Fig 6D is drawn as paired heatmaps whose row order
uses expression-channel ranks only — the smooth appearance of the expression rank columns is by
construction and is not a result — and whose Δ1 colours are within-column *z*-scores, so absolute
magnitudes are never compared across technologies. Under the pre-registered criterion (reproduced
iff ρ_adj ≥ 0.6 and permutation *P* < 0.05), **all three arms reject the null and all three now meet
the effect-size bar**: dsciCAT-seq ρ_adj = **0.896** (95% CI 0.562–0.832 → bootstrap interval
0.806–0.944, *P* = 1.0×10⁻⁴), 10x Multiome **0.824** (0.700–0.902, *P* = 1.0×10⁻⁴), HT-scCAT-seq
**0.630** (0.451–0.768, *P* = 1.0×10⁻⁴) — while noting transparently that the HT bootstrap interval
still spans the 0.6 threshold and rests on 513 oligodendroglial-lineage nuclei (293 Oligo / 220 OPC),
against 5,813/1,488 in ours and 540/298 in 10x (Fig S6h b). The MS exception — chromatin rank worse
than expression rank with a positive off-diagonal *z* — is met in all three arms, with very different
strength: 10x Multiome #4 → #49 (Δ1 +0.661 → −0.473, *z* = +2.62, rank 1/74), HT-scCAT-seq #11 → #41
(+0.339 → +0.064, *z* = +0.89), and our own chemistry +0.397 → −0.153 (#7 → #37, *z* = +2.48; formerly
marginal at +0.20 on the un-logged input before the 2026-11-04 correction). We report this as a
**cross-technology parallel observation**: it is neither a proof of, nor a rebuttal to, the own-data
result, and no platform is benchmarked head-to-head or compared by *P* value.

### 1.4 [P4 → bridge and limits（保留，必须存在）]

The two runs of the own-data contrast are deliberately reported as separate values rather than
merged. Fig 6C is the O1 run (5,962 nuclei, raw ranks, ρ = 0.819) and Fig 6D is the D5 rerun (own
contrast on the 30,559-gene intersection, domain-adjusted ranks, ρ = 0.896). Recomputing the C
contrast on the D matrices gives ρ = **0.748** (0.734 adjusted), Δ = **−0.071** against the delivered
0.819 and outside the pre-registered ±0.05 bridge tolerance; we attribute this to the chromatin
channel's wider control pool (5,962 → 7,523 nuclei) and **disclose it rather than reconcile it**, so
0.819 and 0.896 are not a family of values, are not averaged, and are not presented as a trend. The
expression channel used throughout is the measured Seurat matrix (CP10k + log1p), not the ArchR
GIM/CCA-imputed matrix, and cell types for all three chemistries follow the harmonised 33-class
annotation of Fig 2 rather than the per-channel labels of the source studies (whose 10x ATAC labels
are RNA-transferred and would be circular).

### 1.5 [P5 → cis-window one-liner]

As an orthogonal genome-level proximity check, 38 of 331 platform-identified oligodendrocyte genes
(11.5%) fell within ±500 kb of the 178 genome-wide significant MS sentinel variants, against 10.3%
for the protein-coding background (permutation *P* = 0.26): the overlap is illustrative, not
evidence of enrichment (Fig S6a; Table S6).

---

## 2. Fig 7 对应正文段（替换稿件的 Fig 7 段；主图现为 A/B/C 三格，v0 草稿的 D/E 内容由本段+Table S8 承载）

### 2.1 [P1 → Fig 7A] Donor-level replication in five external studies

The axis itself is reproducible outside our data. In four independent human MS post-mortem cohorts
(GSE277430, GSE118257, GSE301908, GSE231585) and in the mouse EAE spinal-cord multiome study
(GSE250589), per-donor (or per-library) Δ1 — mature-OL minus OPC median MS scDRS score under the
frozen 74-trait protocol — is positive in **54 of 55 human units** (binomial *P* = 1.6×10⁻¹⁵; the
single negative donor, GSE277430 2006-041, has only 53 OPC nuclei) and in **19 of 19 mouse
libraries**. Three of the four human cohorts use the original authors' own cell-type labels; for
GSE231585, where no per-nucleus labels are distributed, we assigned clusters following the authors'
published marker panel and eight-class scheme (Methods), reproducing our previous partition exactly.
**The between-donor MS-versus-control comparison, in contrast, is null in all four human cohorts**
(Δ = +0.226, *P* = 0.500; +0.058, *P* = 1.000; +0.230, *P* = 0.091; +0.316, *P* = 0.393): the axis is
a reproducible feature of oligodendroglial biology, not a per-disease diagnostic marker. In the mouse
the EAE libraries sit *below* controls overall (Δ = −0.419, *P* = 0.0092), but the stage split shows
this is an attenuation during active disease with late reversal — every stage keeps a positive
absolute Δ1 (1.08–1.44), and the Late columns have EAE medians above Ctrl (δ = +0.200) — so the
correct statement is "lesion state moves the level of the axis", not "mice are the opposite of
humans". These disease-state contrasts are confounded with composition (EAE n_MOL/n_OPC 2.78 vs
4.01, *P* = 0.035; the RNA-channel EAE effect after composition adjustment is *P* = 0.387, n.s.; the
ATAC channel shows no effect, +0.011, *P* = 0.775), so state and composition cannot be separated in
these data and we do not interpret the direction mechanistically (Table S8; per-library detail). At
the trait level within the mouse cohort, MS is the third-ranked of the 74 sets (Δ1 = +1.178,
**22/22 libraries positive**, sign *P* = 2.4×10⁻⁷), while SCZ and SWB rank 70 and 74 — the same
specificity pattern as in human cohorts. We additionally probed the GSE250589 paired multiome
chromatin arm directly (scDRS on gene scores, 12 targeted traits, EAE vs Ctrl within-library): the
ATAC channel showed no usable cross-channel trait concordance (all Δ1 negative, a systematic
OPC-greater-than-mature accessibility bias; RNA-vs-ATAC Spearman ρ = −0.12 raw, 0.06 after
domain-size adjustment), which is why the cross-channel claim is confined to the own same-nucleus
data (Fig 6C) and the mouse chromatin arm enters this figure only through the supplement.

### 2.2 [P2 → Fig 7B, S7a] Where MS sits among 74 traits — after confronting the confound

Across the four human cohorts the same frozen 74-gene-set panel can be scored cohort by cohort
(Fig 7B). Read naively, MS is the top-ranked neurological trait in three of four cohorts (raw ranks
2 / 4 / 7 / 2). Before taking that as a result we confronted the obvious confound: a trait's Δ1
scales with the regulatory-domain size of its gene set (Spearman ρ of Δ1 with ln n_peak = −0.84 /
−0.69 / −0.78 / −0.87 across the four cohorts; Fig S7a), so we residualised Δ1 on ln n_peak and
re-ranked. **After adjustment MS is rank 1 / 1 / 3 / 1 — the only trait of 74 that remains at or
below rank 3 in every cohort** (mean adjusted rank 1.50). Two further robustness checks are reported
in Fig S7a rather than on the figure face: size-matched-window sensitivity (MS is in the top two of
40/40, 3/40, 5/40 and 2/40 traits matched to MS's domain size to within ±15%) and the snATAC chromatin-channel rank (MS 36 → 38 after adjustment, i.e. adjustment does not
rescue, and is not needed for, the RNA-channel OL-axis placement). Because ranks are per-cohort permutations of the same 74
traits, the cross-cohort consistency we claim is **on the RNA channel only**; the previously
suggested "74-trait ranking is highly consistent across channels" claim was tested and withdrawn
after residualisation (RNA–RNA rank concordance 0.77–0.96 raw, 0.73–0.93 adjusted, retained;
RNA–chromatin 0.53–0.75 raw but 0.04–0.22 adjusted, collapsed — Fig S7c; Table S7-8).

### 2.3 [P3 → Fig 7C] The placement is specific to the oligodendrocyte maturation axis

Replacing Δ1 by two independent cell-identity axes — Δ2A (astrocyte − immune) and Δ2B
(neuron − glia) — within the same cells and the same 74 traits drops MS outside the top 40 on both
control axes in all four cohorts (#65/#44, #71/#59, #66/#60, #67/#74), while it remains #2/#4/#7/#2
on the OL-maturation axis. Schizophrenia behaves as the mirror image and, deliberately, as the
methodological calibration: SCZ, whose heritability is known to enrich in neuronal lineages (Skene
2018; Finucane 2018), ranks last on the OL axis (#67/#73/#63/#67) yet near the top on astro−immune
(#2/#10/#5/#4). If the ranking procedure were simply noisy, neither pattern would be so consistent;
that two orthogonal identity axes move in opposite directions for the two traits — and cross-axis
correlations ρ(Δ1,Δ2A) = −0.60/−0.66, ρ(Δ1,Δ2B) = −0.70/−0.87 — removes the axis-selection objection.
The placement is therefore OL-lineage maturation specific, not a general glial or neuronal signal.

### 2.4 [P4 → 结论句（可作 Fig 7 段末句）]

Taken together: the genetic-risk gradient from precursor to mature oligodendrocyte is a reproducible
feature across cohorts, platforms and species; the MS-specific top placement within the 74-trait
panel survives domain-size adjustment in every human cohort and is confined to the OL-maturation
axis; and its *level* — but not its *existence* — tracks lesion state rather than diagnosis, with
cell composition an inseparable co-varying explanation we decline to resolve in these data.

---

## 3. 图注瘦身清单（本轮已执行 / 待作者拍板）

**Fig 7（已按此原则更新 `fig7_v13_legend_final.md`）**，移入本正文的文件：
- (A) "The stage-wise mouse pattern is an attenuation...it is not 'mice are the opposite of
  humans'" → 解读，移正文 §2.1。
- (C) "SCZ — the methodological calibration...showing that the MS pattern is axis specific rather
  than ranking noise" 的论证部分 → 移正文 §2.3；图注只留中性一句（SCZ plotted as control trait with
  previously reported neuronal enrichment）。
- (C) "showing that the MS pattern is..." 与 ⟹ 全句（top-7 / out-of-top-40 汇总解读）→ 正文；
  图注保留读数清单本身（#k 名次，图面可见性的补充）。
- 探针（ATAC 无可用信号）一句 → 正文 §2.1 末；不属图注内容。

**Fig 6（`fig6_v21_legend_final.md` 数字无需更新——已是 log1p 修正定稿）**。按"图注只写看不到的
细节"原则，下列句子按现行 final 属"可辩论移入正文"，**因 final 与 draft/verify 61 句门绑定，
未擅动，等作者点头后再删**（删时 draft/final 同步 + 重跑 `A6_fig6_v21_verify.py`）：
- D 段 "Ladder verdict for (D)... 3 of 3 chemistries reproduce the ladder; the ladder is a
  pre-registered effect-size screen..."（判据定义留图注、结论句移正文）。
- D 段 MS 并行观察段末 "The chromatin channel ranks the MS exception lower... so the counterexample
  is reproduced in every technology"（解读句；数值句保留）。
- C 段 "Why HEIGHT #1 on chromatin is not read as trait-specific..."（整段移正文/补充口径）。
- S6g/S6h 面板标题为图面已印的论证式标题（属图面，不能只删图注；若瘦身须连图重做，不建议动）。
