# Figure S6 — cross-channel concordance: sensitivity, units, replication (v16/v18 FINAL)

> **2026-11-18**：交付页升级为四行版 `Figure_S6_v16.pdf`（17.00 × 18.91 cm，
> `routeA_out/S6_v16_20261117/`，`assembly/figS6v16_run_all.sh` + `figS6v18_row3.sh` 重渲）。
> 当日审计（`FIGS_OVERALL_AUDIT_20261118.md` §2-①）发现 v16 的 G 面板曾误用 log1p 修正前的
> tables_v16 CSV，已换修正表重渲（仓内 `tables_v16/` 即修正版，`tables_v16_pre_log1p_backup_20261118/`
> 是修正前的本地留档，不入 git）。
> **2026-11-19 → 24**：**Table S6 台账升级为 v13**（`tables_v12/Table_S6_v13.csv/.pdf`，
> 209 数据行 × 12 列，文字列全部改写为面向读者的措辞、数值零改动，残留内部措辞单元格 = 0）。
> 交付树 11-23 那份 `Table_S6.pdf`（5 pp，`f76bef3ea9ec…`）仍是 297 行的 v12 旧文字版；
> 两边已按 82694f94 的仓内 v12 对齐（见 `dev/TABLES.md` 与本目录 §Table S6）。

`Figure_S6.ipynb` reproduces the seven delivered panels of **Supplementary Figure S6** in their
current **v16/v18** form from the corrected (log1p, 2026-11-04) multiome analysis objects. The
notebook reads only tables vendored in this package (`tables_v16/`, `figS6_out/caches/`,
`upstream/`), prints each panel's measured size and asserts the face numbers against the ledger.
Journal assembly (the four-row 17.00 × 18.91 cm sheet, panel letters, shared legend placement) is
done outside the notebook by `assembly/figS6v16_compose.py`; `assembly/` holds the delivered
pipeline — panel scripts, compose, the Table S6/S7 generators, the ledger statistic generators and
the **reverse-trace gate** the delivered page passes — while `upstream/` holds the scripts that
produced the input tables.

This directory supersedes the 2026-10-23 notebook (the old SUPPORT pipeline, whose sheet carried
S6a–S6j). That version is archived in the authors' tree at
`routeA_out/SUPPORT_v2/pdf/superseded/github_figS6_v7_20261023/`.

## Panels (current v16/v18 delivery — sizes measured on the shipped PDFs)

| Panel | File | Size | Content |
|---|---|---|---|
| A | `figS6_out/figS6v16_A_loci.pdf` | 10.78 × 2.79 cm | four example loci (PHLDB1, PHGDH, PLLP, BCAS1) ± 500 kb around the lead SNP; hit gene red; lead-SNP *P*; full hit list in Table S6 |
| B | `figS6_out/figS6v16_B_perm.pdf` | 5.38 × 3.05 cm | window permutation of the hit count: 38/331 = 11.5 %; v13 recomputes the null from source (10,000 draws, seed 20261104, add-one) → *P* = 0.2767723 |
| C | `figS6_out/figS6v16_C_umapbins.pdf` | 17.00 × 3.00 cm | bin-size robustness of the Fig 6B RNA UMAP (single cells / 15 / 25 / 40, shared colour scale) |
| D | `figS6_out/fig6v32_C_top1000_landing_s6v18D.pdf` | 7.45 × 5.75 cm | gene-level landing of the MS MAGMA top-1000, six subtypes × both channels; deliberately unlabelled (null panel) |
| E | `figS6_out/fig6v31_C_landing_s6v18E.pdf` | 7.45 × 5.75 cm | the 38 cis-window genes land on MOL6: mean Δ = +0.145 (RNA) / +0.097 (gene score), *P* = 0.0017 / 0.0007 |
| D+E legend | `figS6_out/figS6v18_DE_legend_v.pdf` | 1.62 × 5.75 cm | the one vertical swatch legend shared by D and E (clip triangles, colour ramp) |
| F | `figS6_out/figS6v16_F_cellboxes.pdf` | 8.08 × 5.75 cm | per-nucleus OPC-vs-mature-OL contrast, 2 channels × 3 chemistries, MWU *P* per panel (asserted against the ledger) |
| G | `figS6_out/fig6D_heatmap_pair_compact_s6v16G.pdf` | 8.64 × 5.39 cm | 74-trait Δ1 pair (expression × chromatin) across the three chemistries; size-adjusted cross-channel ρ = 0.726 / 0.763 / 0.518, 2-of-3-chemistries verdict |

The composed sheet itself is **not** versioned here (see the root README, *What is in each
package*): readers rebuild it with `bash assembly/figS6v16_run_all.sh` and
`assembly/figS6v16_compose.py`. The v11 sheet (`fig6_out/figS6_v11_FINAL.pdf`, 17.00 × 18.86 cm — **local working tree, not
in git**) and the v10/v11 panel PDFs it was built from are superseded and stay in this working tree only.

Numbers that used to sit on the face (Δ1, the gate rho + noise ceilings, mouse 18/18) live in
**Table S6**. Current version = **v13**: `tables_v12/Table_S6_v13.csv` (209 data rows × 12 columns,
md5 `aae38a32133d645b344cd821c5cb6ac0`) and `tables_v12/Table_S6_v13.pdf` (4 pp, A4 landscape
29.7 × 21.0 cm, md5 `175edf1425ebb2df6757d77d5ff1ece6`), produced by `dev/make_table_S6_v13.py`,
which reads `tables_v12/Table_S6_v12.csv` (repo copy, md5 `82694f942067…`) and rewrites **text
cells only** — every number is carried through unchanged; re-running reproduces the published v13
byte-for-byte. Per-mouse rows are in **Table S7** (`fig6_out/Table_S7.pdf`, 1 pp, md5
`65a156b9293711e49f8b44ab5eb032df`, text-sha256 `1f9be8d44be8d26a`, 18 rows × 14 cols).
The typesetter stamps the PDF, so **md5 is per-run**: re-running a generator on the same inputs
yields the same extracted text and geometry but a different md5 — ledger signatures should quote
the text-sha (or pin one file). The v11/v12 table faces (`Table_S6_v11.pdf` md5 `b4cdfba1bcad…`,
`Table_S6_v12.pdf` md5 `35f219d03452…`) are superseded and are not versioned; their machine-readable
twins stay in git because they are generator inputs: `tables_v11/Table_S6_v11.csv` (md5
`171c290a18ce0a870c4445922f59f3a3`, byte-identical to what `assembly/make_table_S6_S7.py` writes)
and `tables_v12/Table_S6_v12.csv`. The legend that ships with the current page is
`../legends/figS6_v16_legend_final.md`; the v11-era snapshot `LEGEND_S6_v11.md` stays in this working
tree only, together with the v10/v11 faces it describes.

## Pipeline

> **Script generations inside `assembly/`.** The v16/v18 line (`figS6v16_panels.py`,
> `figS6v16_compose.py`, `figS6v18_legend_right.py`, `figS6v16_run_all.sh`, `figS6v18_row3.sh`,
> `make_figS6v16_G_heatmap.py`) is the **current delivery pipeline**. The v10/v11 line
> (`make_S6_v10_panels.py`, `make_S6_v11_f_npeak.py`, `S6v10_*.py`, `verify_S6_v10.py`,
> `compose_support_S6S7.py`) is kept as **number provenance** for the 2026-11-06 sheet: it reads
> and writes the analysis tree (`routeA_out/SUPPORT_v2/…`), and the faces it produced are not
> versioned — running it is not part of the reproduction contract for the current version.
```
own dsciCAT-seq cohort (CNP0010263) + GEO cohorts
  └── upstream/: scDRS scoring, O1 same-nucleus tables, M2b per-mouse contrasts,
        locus windows, log1p-corrected ATAC rerun (delivered 2026-11-04, audited by A6)
              └── assembly/: ledger statistics (window perm, permnull, bootstrap, MWU, sixpoint)
                    ├── Figure_S6.ipynb  (this notebook: panels A–G PDFs)  ──> fig6_out/
                    ├── make_S6_v10_panels.py + compose  (v10 sheet + number index)
                    ├── S6v11_bootstrap_npeak.py + make_S6_v11_f_npeak.py  (v11 panel F + sheet)
                    ├── make_table_S6_S7.py  (Table_S6_v11: n_hit relabel + n_peak row)
                    └── verify_S6_v10.py  (reverse-trace gate: 8/8, strict index mode)
```

## Inputs

All inputs are read through environment variables; nothing is downloaded by the notebook.

| Variable | Object | Produced by |
|---|---|---|
| `S6_ROOT` | `routeA_out/SUPPORT_v2/tables_v8/*.csv` | `assembly/` ledger scripts |
| `S6_ROOT` | `SUPPORT_v2/tables_v11/figS6v11_{unit_summary_ext,bootstrap_nhit,bootstrap_npeak,F_number_index}.csv` | `assembly/S6v11_bootstrap_npeak.py` (panel F v11 rungs) |
| `S6_ROOT` | `SUPPORT_v1/tables/deconfound_policy_table.csv` | `upstream/D10_deconfound_policy_table.py` (B = 10,000 n_peak CI in Table S6 v11) |
| `S6_ROOT` | `O1_own_multiome_paired/tables/{T4_delta1_74traits_both_channels,cell_level_table}.csv` | `upstream/O1_own_same_nucleus.py` |
| `S6_ROOT` | `M2b_authoranno_v1/tables/M2b_per_mouse.csv` | `upstream/M2b_formal_authoranno.py` (GSE250589, author annotation) |
| `S6_ROOT` | `D_crosschemistry/{scdrs,scdrs_log1p_20261104}/**/PASS_Multiple_sclerosis.full_score.gz`, `tables_log1p_20261104/{cells_*.csv,fig6D_delta1_74traits_per_tech.csv}` | cross-chemistry scDRS runs; log1p rerun delivered 2026-11-04 (audited by `upstream/A6_fig6_input_transform_audit.py`) |
| `S6_FY` | `oligo/scdrs_ms_out/ms_locus_overlap_table.csv`, `gencode.v45.basic.annotation.gtf.gz`, `oligo/3k_dsci_oligo_subset.h5ad` | own cohort objects + `upstream/ms_locus_overlap.py` |
| `S6_OUT` | panel PDFs | this notebook |
| `ARIAL_TTF` | `ARIAL.TTF` | Arial face (registered when absent from the OS; cosmetic only) |

The large objects (score `.gz`, h5ad) are not in the repository; the notebook reads them in
seconds-to-minutes (G ≈ 30 s, F ≈ 10 s, full run ≈ 2 min).

## Running

```bash
cd figS6
S6_ROOT=/path/to/routeA_out S6_FY=/path/to/FY \
  jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.timeout=-1 Figure_S6.ipynb
```

Executed end-to-end on 2026-11-05 with **0 errors** (python 3.9.23, numpy 1.26, pandas 2.x,
scipy 1.13, matplotlib 3.8, anndata 0.10, PyMuPDF; scdrs not needed by the notebook — the
score tables were produced upstream under `$DSCI_PY` (scdrs 1.0.2).
Regenerate the notebook itself with `python3 ../../dev/build_figS6_notebook.py` (byte-stable cells).

## Validation checklist

- [x] panel sizes match the delivered cells (printed by `savecell`)
- [x] Panel F recomputes all 9 bootstrap rungs (v11) and **asserts** them against the delivered
      `figS6v11_unit_summary_ext.csv` (same seed ⇒ identical draws)
- [x] self-check cell: every required number string is present on its regenerated panel
      (forward half of the gate); no banned pre-correction value reappears
- [x] delivered sheet passes the strict reverse-trace gate `assembly/verify_S6_v10.py`:
      8/8, every face number registered in `figS6v10_number_index.csv`, 0 unsourced,
      0 unregistered, 0 missing-key, 0 old-lineage values
- [x] fabricated-token probes blocked (`assembly/S6v10_gate_tolerance_probe.py`:
      4/4 blocked, 3/3 real values pass)
- [x] width 17.00 cm ≤ 17 cm two-column limit; Arial only; visible text ≥ 5.5 pt
- [x] font gate run on the final sheet 2026-11-06 (`fig6/assembly/font_gate.py`):
      `figS6_v11_FINAL.pdf` 17.00 × 18.86 cm, 176 spans, 5.50–8.00 pt, faces Arial-BoldMT/ArialMT
      only — PASS (Tables S6/S7 likewise Arial-only, 6.20–8.40 pt; their 29.7 cm page width is the
      landscape table format, not a figure limit)
- [x] v11 number audit 2026-11-06: all nine F-rung face numbers match
      `figS6v11_F_number_index.csv` / `figS6v11_unit_summary_ext.csv` digit-for-digit; B/D/E/G legend
      numbers re-verified against `tables_v8` sources; caption-patch CI for the n_peak rung corrected
      0.547–0.815 → **0.547–0.813** (B = 1,000 delivered draws; 0.550–0.815 is the D10 B = 10,000 run)
