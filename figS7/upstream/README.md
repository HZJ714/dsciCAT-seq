# upstream/ — provenance scripts for `../data`

These scripts **produced the input tables**. They are shipped for provenance and are *not*
called by the figure: `Figure_S7.ipynb` runs entirely from `../data`, so a reviewer can redraw
and re-audit the page without any of these.

Two consequences of that role:

1. They read the authors' analysis tree — the EAE probe directory, the ArchR objects, the
   GWAS gene sets — through one variable, `DSCI_WORK_ROOT` (defaults in each script are
   relative to the current directory, so an unset variable means "run from the tree root").
   Re-running them needs those objects; re-running the *figure* does not.
2. Some delivered tables carry the authors' internal run labels in their `source` / `pair`
   columns (e.g. `O1 own dsciCAT-seq same nucleus`, `M3_probe_delta1_3libs.csv`). The tables ship
   byte-identical to the delivery so their md5 stays checkable; the labels are defined in the
   glossary in `../README.md`, and none of them reaches the figure face.

| script | writes into `../data` |
|---|---|
| `confound_all_cohorts.py` | `FIG7_v1/tables/seven4_confound_ALL_cohorts.csv` |
| `matched_window_table.py` | `SUPPORT_v1/tables/fig7Sb_matched_window.csv` |
| `crosschannel_pairs_table.py` | `SUPPORT_v1/tables/figS7c_crosschannel_pairs.csv` |
| `axis_rho_and_rank_matrix.py` | `SUPPORT_v1/tables/fig7Sc_axis_pair_rho.csv`, `fig7Sc_rank_matrix_4cohorts.csv` |
| `delta1_vs_peakdomain.py` | `A6_GSE231584_atac_v1/tables/delta1_vs_peakdomain_allchannels.csv` |
| `human_atac_threeaxis.py` | `SUPPORT_v1/tables/S7_atac_crossaxis_human_full.csv` |
| `mouse_atac_threeaxis.py`, `mouse_atac_threeaxis_composition_fix.py` | `SUPPORT_v1/tables/S7_atac_crossaxis_mouse_full.csv` |
| `mouse_rna_threeaxis.py`, `mouse_rna_gene_stats.py` | `SUPPORT_v1/tables/S7_rna_crossaxis_mouse_full.csv`, `S7_rna_fullstats_GENE_STATS.csv` |
| `export_mouse_same_nucleus_delta1.py` | `SUPPORT_v1/tables/S7_mouse_same_nucleus_delta1.csv` (+ `.json`) |
| `target_stats_GSE118257_v7.py` | `A3_GSE118257_v1/tables/A3_GSE118257_target_stats.csv` |
| `build_manifest.py` | `../data/MANIFEST.csv` (md5 / rows / producer / panel map) |

---

The panel map is regenerated, not maintained by hand: `python3 build_manifest.py` rewrites
`../data/MANIFEST.csv` (md5, rows, producer, cohort, panels, content). Since v29 the letters are
A-G plus the legend pieces K/RB/L2; the withdrawn v22 panels C/F/G are recorded in
`../LEGEND_S7_final.md` section 2 and their numbers now live in Table S7.
