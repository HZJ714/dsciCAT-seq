# Figure S5 — SCENIC+ support: eRSS ranking, state networks, the Mef2c → *Nrxn3* locus, downsampling stability

`Figure_S5.ipynb` reproduces the content of the four Fig S5 blocks from the same analysis
objects as `Figure_5.ipynb`. Compositing (18.00 × 46.75 cm canvas, Arial substitution, panel
letters, the `retext_scenic2.py` text-layer rebuild) is done by hand and is not part of this
notebook.

| Block | What it shows | Notebook section |
|---|---|---|
| A | eRegulon specificity (eRSS) per state, top 10 labelled, 6 states × 96 filtered eRegulons | §1 |
| B | TF–target networks of MOL2, MOL3_6, MOL6 (Fig 5B is MOL3) | §2 |
| C | *Nrxn3* locus: per-state ATAC coverage, four Mef2c-eRegulon enhancers + MEF2 motif hits, *Nrxn3* expression | §3 |
| D | Downsampling stability of the filtered eRegulon set (1k / 3k / full) | §4 |

The provenance of `scplus_obj_final.pkl` — RDS → R export of the binarised peak matrix →
cisTopic LDA at **60 topics** → three binarisations + DARs → pyCisTarget (`v10nr_clust_public`)
→ `create_SCENICPLUS_object` → `run_scenicplus` (10 eRegulons only) → `build_grn` →
`score_eRegulons` → pseudobulk correlation → std filter → eRSS — is documented once, in
**`GitHub/fig5/Figure_5.ipynb` §0** (14 dated steps). It is deliberately not repeated here;
this notebook reads the same object.

## Inputs (`FIGS5_DATA`, default `./data`)

| File | Produced by | Used by |
|---|---|---|
| `scplus_obj_final.pkl` | the chain above (8.0 GB) | A, B, C (expression) |
| `NIHMS70326-supplement-Table_S1.xlsx` | Marques et al. 2016, 24 states × 50 up-regulated genes | B |
| `nrxn3_expr.npz` | `Nrxn3` column of `X_EXP` grouped by state — §3.1 recomputes it and compares | C |
| `nrxn3_tracks.npz` | ArchR **GroupCoverages** of the lineage subset (`FY/oligo/3kdsci/archr-subset/GroupCoverages/cell_type`), FY13+FY14 merged per state, **per-base counts normalised to each state's library size**, 200-bp bins over chr12:88.55–89.60 Mb (mm10), decoded once by `grn_mef2c_nrxn3_track.py` | C |
| `nrxn3_mef2c_enhancers.fa` | the four Mef2c-eRegulon regions in that window, mm10 sequence | C |
| `eregulon_stability_table.csv` | `stability_downsample.py`, seed 123, one full SCENIC+ re-run per depth | D |
| `FIGS5_MOTIF` (default `<work>/aertslab_motif_colleciton/v10nr_clust_public`) | Aerts lab cisTarget release, `singletons/metacluster_179.4.cb` = JASPAR MA0052.4 | C |

## Run

```bash
export FIGS5_DATA=/path/to/objects FIGS5_OUT=figS5_out
jupyter nbconvert --to notebook --execute --inplace --allow-errors Figure_S5.ipynb
```

Needs only the python SCENIC+ environment (`scenicplus`, `pyscenic`, `ctxcore`, `plotnine`,
`dill`) — **no R and no Signac**, unlike Fig 5 §3a: the coverage signal is read from the decoded
`.npz`. `--allow-errors` matters because nbconvert 7 discards *all* cell outputs if one cell
fails. ~10 min, most of it loading the 8 GB object; §2 (three networks) is the slow part.

## Outputs (`FIGS5_OUT`)

`S5a_rss_scores_96.csv` · `S5a_rss_top10.csv` · `figS5a_rss_grid.pdf` ·
`figS5b_grn_{mol2,mol3_6,mol6}.pdf` + `S5b_grn_*_{nodes,edges}.csv` ·
`S5c_nrxn3_expression_by_state.csv` · `S5c_mef2_motif_hits.csv` · `figS5c_nrxn3_locus.pdf` ·
`S5d_stability_table.csv` · `figS5d_stability.pdf` · `requirements.txt`

## Validation against the delivered figure

| Check | Delivered / archived | Notebook | Result |
|---|---|---|---|
| eRSS table (6 × 96) | `pi/plots-scenic+/rss_scores_filtered_96.csv` | recomputed | exact: max abs Δ = `0.000e+00` |
| top-1 eRSS per state | Glis1 0.621 (COP_NFOL), Pax5 0.343 (MFOL), Usf1 0.378 (MOL2), Nr1h3 0.356 (MOL3), Creb5 0.331 (MOL3_6), Creb5 0.513 (MOL6) | same six values | ✔ |
| network sizes / hubs | MOL2 58/97 (Klf7 34), MOL3_6 **56/73** (Mllt10 36), MOL6 **49/65** (Mbnl2 37) | 2026-10-22 run: 58/97, 56/73, 49/65 — edge **sets identical** to the archived Cytoscape exports (97/73/65, undirected set comparison) | ✔ (row corrected 2026-10-29: the archived table had MOL3_6 and MOL6 swapped — the delivered figure face and `figures/scripts/README_FigS5.md` §3 already had it right) |
| *Nrxn3* medians | MOL3 ≈ 7.0 vs 3.1–3.8 elsewhere, one-sided MWU p < 1e-22 | — | to be filled |
| MEF2 motif hits | `pi/plots-scenic+/nrxn3_mef2c_motif_hits.tsv`, 19 hits, best E1 12.94 | — | to be filled |
| stability | 196 / 294 / 292 eRegulons, Jaccard 0.110 / 0.280 / 1.000, top-10 recovery 1.00 | table read as-is | ✔ (D is a re-plot of the archived table) |

## One registered discrepancy — not fixed

The delivered **legend** for the stability block says "1k: 208 → 6k: 292" and "0.11 → 3k: 0.26
→ 1.00". The delivered **figure face** and `eregulon_stability_table.csv` say **196** and **0.28**.
The panel was drawn from that table (`figures/scripts/figs5_stability.py`), so the legend numbers
are stale. `figures.docx` / `figure_captions.md` were signed off on 2026-10-21, so nothing was
edited here; the notebook states both numbers in §4 and waits for the author's call.

## Not reproduced here

Layout and compositing (`compose_figS5_noc3.py`, `figs5_repanels.py`, panel letters, shared
legend, Arial), the full downsampling re-run (one complete SCENIC+ run per depth, ~30–60 min
each — §4 reads its output table), and the ArchR coverage extraction itself (R + ArchR; §3.3
reads the decoded `.npz` and cites the log-verified origin). Figures are saved to disk and only
their paths are printed: a screen render from this notebook is not the manuscript file, and
comparing fonts/sizes would be comparing things this notebook does not control.

## 2026-10-29 之后的两处改动（Codex #90/#91 = 图侧 G3/G4）

1. **§2 的机器可读补表换了导出对象。** 旧版手写两列 CSV，两个坑：
   `G.nodes[n].get('TF', False)` 永远是 `False`（图上的节点属性不叫 `TF`），所以导出的
   `S5b_grn_*_nodes.csv` 里 58 行 TF 列全 False；`G.edges` 只有顶点对，丢了
   `TF2G_importance` 与方向 —— 而图里的**边宽正是这个值**，且 `create_nx_graph` 用
   `from_pandas_edgelist` 建的是**无向图**，`source/target` 顺序是任意的。
   现在直接落 `create_nx_graph` 返回的 `edge_tables` / `node_tables`（就是画图用的那两张，逐列自解释）。
   本轮**没有**为验证这一改动重跑 8 GB 对象；交付用的表仍是
   `pi/plots-scenic+/grn_*_cytoscape_{edges,nodes}.csv`（含 `weight,importance,regulation,eRegulon,P_val,TF_rho,TF_peak`，
   `source` 全为 TF，实测每 TF 边数 = hub degree）。
2. **§1 的图与交付版面分家了。** 交付件 `FINAL_SET_20261021/Figure_S5.pdf` 的 A 面板自 2026-10-29 起
   由 `figures/scripts/figs5_rss_unified.py` 重出（六格统一 y 量程 + eRSS=0.5 基准线 + top-4 标签），
   本笔记本 §1 画的仍是"每格独立量程、top-10 标注"的原始 `plot_rss` 复刻版 ——
   **数值同源**（同一份 eRSS 矩阵，最大绝对差 0.000e+00），只是版式按审稿意见改了。
   B 面板另加了一个边方向/权重图例块（`figures/scripts/network_key.py`，内容实测自成品矢量）。

## What changed on the published figure (2026-10-29, reviewer round 2)

The notebook still reproduces the *analysis*; the four items below are figure-side and live in
`assembly/` (they edit the rendered PDFs / compose the page, they do not re-run SCENIC+).

| Item | Why | Script |
|---|---|---|
| Block A redrawn: all six facets share one y-scale (0.15–0.65), each facet shows its own ticks, an `eRSS = 0.5` reference line, labels cut to top-4 | reviewers read the unlabelled right-hand facets against the neighbouring facet's scale and reported "MOL6 only reaches 0.33"; the underlying values were never wrong (max abs delta vs the notebook output = 0.000e+00) | `assembly/figs5_rss_unified.py` |
| Full top-10 eRSS lists moved out of the figure into **Table S5** (61 rows × 5 columns; this was numbered S6 before the manuscript v16 renumbering) | same reason: definitions belong in the legend/table, not in six crowded panels | `assembly/make_table_s6_rss.py` |
| Network key added to block B (regulator diamond → target diamond, edge-width samples, colour = source eRegulon) | the networks encode direction/weight/colour with no key on the page | `assembly/compose_figS5_noc4.py` + `../fig5/assembly/network_key.py` |
| Block C: two wrong statements on the panel removed | ① the title claimed the tracks were **RPM-normalized**; they are not (the plotted value is per-base, library-size-normalised insertions per 200-bp bin, each row on its own scale — maxima 2.9e5–6.2e5). ② the panel printed `MWU p<1e-22` while the text says `6.2e-80`: both are correct but they are *different comparisons* (weakest pairwise vs COP/NFOL, n = 203 — versus MOL3 against the other five states pooled). The panel now reads `p<1e-79` (pooled, conservative) and the legend gives both | `assembly/patch_figS5_c_units.py` |

⚠ Block C is patched **in the delivered PDF**, not re-rendered: `assembly/figs5_panel_track_box.py`
no longer reproduces the published panel (re-running it unchanged moves 3.62 % of the pixels), so a
re-render would silently change the picture while we were only trying to fix three lines of text.
The re-render is kept aside as `replot/S5c_nrxn3_track_box_rerender_probe.pdf` upstream.

`assembly/figs5_repanels.py` is the older A-panel re-panel script (kept for provenance; it rewrites
notebook cell sources, so it is deliberately excluded from the path-porting tool — see
`../00_Figure_Tools/port_paths.py`).
