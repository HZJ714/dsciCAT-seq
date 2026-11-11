# Table inventory — Table 1, Table 2 and Tables S1–S8 (matches manuscript v16/v17)

The manuscript now carries **10 tables**: 2 main (Table 1, Table 2) and 8 supplementary
(Tables S1–S8). The supplementary numbering was re-ordered when Table S8 was added
(old S6 → S5, old S7 → S6, old S8 → S7), so this file pins **number ↔ content ↔
generator ↔ delivered artefact** in one place.

| No. | content | carried by | generator in this repo | delivered artefact (analysis tree, not in this repo) |
|---|---|---|---|---|
| Table 1 | datasets and libraries generated in this study | embedded in manuscript | statistics cells in `fig1/`, `00_Preprocessing/` (same numbers as the Fig 1/2 legends) | — |
| Table 2 | overlap of MS risk loci with platform-discovered genes (38/331) | embedded in manuscript | locus-overlap step in `fig6/assembly/` | — |
| Table S1 | barcoded primer configurations (99 rows) | embedded in manuscript | `00_Preprocessing/` (primer design side) | — |
| Table S2 | scDRS score statistics per oligodendrocyte subtype + per-library composition (11 columns) | embedded in manuscript | scoring in `fig6/upstream/`, composition in `figS2/` | — |
| Table S3 | fraction of cells with scDRS FDR < 0.2 per modality and subtype (13 rows) | embedded in manuscript | `fig6/upstream/` | — |
| Table S4 | mean MS scDRS norm-score per cell type, four datasets (33 rows; header `ct,10x Multiome,HT-scCAT-seq,dsciCAT-seq-100k,dsciCAT-seq-50k`) | embedded in manuscript (caption: "Mean MS scDRS norm-score by cell type…") | norm-score aggregation in `fig6/upstream/` + `figS2/` | `figures/TABLE_S5_source_allcells_type_means.csv` md5 `adf92abdff9b`; `figures/TABLE_S5_full_rows.json` md5 `150753cb20bb` (both are **norm-score table** sources; the `TABLE_S5_*` name predates the S6→S5 renumbering and now refers to this table) |
| **Table S5** | top-10 eRSS per lineage state (61 rows) | embedded in manuscript | `figS5/assembly/make_table_s5_rss.py` (renamed 2026-11-11; it was make_table_s6_rss.py) | eRSS source tables: `rss_scores_filtered_96.csv` and `figS5/figS5_out/S5a_rss_scores_96.csv` (in this repo); the generator's header self-describes them and asserts 0.000e+00 consistency with the delivered rows |
| **Table S6** | cross-channel / cross-technology consistency ledger — **current version: v13, 209 data rows × 12 columns** (v12 213 rows minus the 4 rows dropped in the 2026-11-19/24 consolidation; the 297-row figure quoted earlier counts the pre-trim v11 sheet). Repo artefacts: `figS6/tables_v12/Table_S6_v13.csv` md5 `aae38a32133d645b344cd821c5cb6ac0`, `Table_S6_v13.pdf` md5 `175edf1425ebb2df6757d77d5ff1ece6` (4 pp, A4 landscape, 0 internal-wording cells, verified by `receipts` V13–V15). Generator: `dev/make_table_S6_v13.py`, which reads `figS6/tables_v12/Table_S6_v12.csv` (repo copy md5 `82694f942067…`, the version both sides reconciled to) and rewrites only text cells — numbers untouched; re-running reproduces the published v13 byte-for-byte. | delivered CSV/PDF | `figS6/assembly/make_table_S6_S7.py` (v11/v12 base) + `dev/make_table_S6_v13.py` (v13 text-column pass) + `figS6/assembly/render_table_S6.py` | superseded faces `Table_S6_v11.pdf`, `Table_S6_v12.pdf` and the analysis-tree delivery copy (`SI_TABLE_DELIVERY_20261123/Table_S6.pdf`, md5 `f76bef3ea9ec…`, 5 pp, text-sha `81d19bde6bc5`, still the 297-row v12 wording) stay in this working tree only, not in git |
| **Table S7** | mouse EAE per-library ledger (18 rows × 14 columns) | delivered CSV/PDF | per-library statistics in `figS7/upstream/` + `figS6/assembly/make_table_S6_S7.py` | `SUPPORT_v2/tables/Table_S7.csv` md5 `5316d6240402`; `SUPPORT_v2/pdf/Table_S7.pdf` md5 `65a156b92937` (1 page, A4 landscape) |
| **Table S8** | disease / lesion layer: 4 human cohorts + mouse stages (28 rows × 12 columns) | delivered CSV/PDF | `figures/S8_DELIVERY_20261111/make_table_S8.py` (trim rules and value-by-value provenance in `RECEIPT_V16_S8_20261111.md`, same directory) | `S8_DELIVERY_20261111/Table_S8.csv` md5 `58b64351f9c4`; `Table_S8.pdf` md5 `42c659ded255` (1 page, A4 landscape, Arial) |

Conventions
-----------
1. The manuscript cites **only Table 1/2 and S1–S8**. Older phrasings such as "S6a–S6g,
   ten panels + Tables S6/S7" are obsolete: Supplementary Figure S6 is the current **v16/v18
   four-row sheet (A–G, 17.00 × 18.91 cm; the v11 sheet at 18.86 cm is superseded)** and is a
   different object from Table S6. Figure S7 is v30 (A–G + K/RB/L2, 17.00 × 14.28 cm, ten panel
   PDFs); "S7 A–I" is the obsolete v22 wording.
2. The caption-less 23-row abbreviation/cell-type-label table that previously sat between
   Table S3 and Table S4 in the manuscript has **no number and is unreferenced** (audit
   2026-11-11, D1): cell-type labels and abbreviations are defined inline in the main
   text. Decision: delete the table from the manuscript (parallel-env revision v17) and
   keep S1–S5 + S8 as the printed set, matching the SI statement "Tables S1–S5 and S8
   are printed below".
2. Table S8 is a **trimmed** delivery: every value is traceable to a row of the internal
   `FIG7_v12/tables/fig7D_Table_S8.csv`; values that were not delivered are written
   `n.d.` and never recomputed here. The full 78-row internal table stays out of the
   public tree.
3. Changing a table number means changing it in five places at once: the manuscript,
   this table, the root README, `legends/*_legend_final.md`, and the comments inside the
   generating script. `dev/push_preflight.sh` scans for stale numbers.

## Version policy (2026-11-25)

Only the **currently final** version of every figure panel, table face and legend is
versioned. Superseded sheets (`_archive_*`, `_superseded_*`, `tables_v16_pre_log1p_backup_*`,
the v10/v11 Figure S6 faces, the v21/v13/v10/v29 legend snapshots, `dev/backups/`) remain in
the authors' working tree as a fallback and are gitignored; each affected directory keeps a
`README.md` in git saying what was archived and what replaced it. Composed journal pages are
assembled outside the repository (see the root README, *What is in each package*).
