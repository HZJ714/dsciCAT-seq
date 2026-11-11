# assembly/ — the delivered Figure S6 v10 pipeline

These scripts produced the delivered sheet, tables and ledgers that `Figure_S6.ipynb`
reproduces. Run order (the shared helpers live under `$DSCI_WORK_ROOT/scripts/routeA` in the
authors' tree; any python ≥ 3.8 with numpy/pandas/scipy/matplotlib/PyMuPDF — see the root
README *Requirements* — runs them from this directory):

## 1. Ledger statistics (inputs of the notebook and of the panels)

| script | writes |
|---|---|
| `S6v8_window_permutation.py` | `tables_v8/figS6v8_window_permutation{,_null}.csv` (Panel B) |
| `S6v10_perm_null.py` | `tables_v8/figS6v10_permnull_{o1,summary}.csv` (Panel E) |
| `S6v10_bootstrap_units.py` | `tables_v8/figS6v10_bootstrap_units.csv`, `figS6v10_unit_summary.csv` (Panel F core 7 rungs) |
| `S6v11_bootstrap_npeak.py` | `tables_v11/figS6v11_bootstrap_{npeak,nhit}.csv`, `figS6v11_unit_summary_ext.csv` (Panel F v11: the two adjustment rungs, independent re-seeded streams) |
| `S6v10_G_mwu_p.py` | `tables_v8/figS6v10_G_opc_oligo_mwu.csv` (Panel G) |
| `S6v8_derived_stats.py` | `tables_v8/figS6v8_sixpoint_all74.csv`, `figS6v8_state_median_74x6x2.csv` (Panel D), ladder variants |
| `S6v10_window_gene_hits.py` | window hit counts feeding the Table S6 locus block |

## 2. Figure + tables

| script | writes |
|---|---|
| `make_S6_v10_panels.py` | the seven cell PDFs, `figS6_v10_draft.pdf` (17.00 × 18.86 cm, via `compose_support_S6S7.py`), `tables_v8/figS6v10_number_index.csv` (every printed number registered with its source table) |
| `make_S6_v11_f_npeak.py` | `S6f_unit_violins_v11.pdf` + `figS6_v11_candidate.pdf` (A–E, G = v10 cells read-only; F = v11) + `tables_v11/figS6v11_F_number_index.csv`; asserts 9 rungs and the F face numbers |
| `make_table_S6_S7.py` | `tables/Table_S6.csv/.pdf` (v10), **`Table_S6_v11.csv/.pdf`** (n_hit/n_peak relabel + 0.708 row) and **`Table_S7.csv/.pdf`**; run with `SKIP_S7=1` to leave the delivered Table_S7 untouched — the 2026-11-11 rerun deliberately did *not* set it, so both table PDFs were rebuilt together after the group-name fix (see `figS6/README.md`) |
| `S6v10_update_table_face_tags.py` | rewrites the `in_figure` column of Table S6 from the actual PDF faces (strict retag, 2026-11-05) + audit CSV |

Shared helpers: `support_common.py` (paths, palettes, logging), `arial_style.py` (Arial
registration, compose layer).

## 3. QA gates (must pass before delivery)

| script | role |
|---|---|
| `verify_S6_v10.py` | **reverse-trace gate**, strict index mode: every number token on every panel/draft must be registered in the number index (verbatim printed string, or equal to a registered value rounded to the token's significant digits) or be a whitelisted structural value (ticks, coordinates, GSE accession). Categories: NO-SOURCE / UNREGISTERED / MISSING-key / BANNED-old-value. `S6_GATE_LENIENT=1` downgrades UNREGISTERED to a warning (legacy mode). |
| `S6v10_gate_tolerance_probe.py` | adversarial probe: fabricated tokens (`9.9e-77`, `4.0e-10`, near-miss `0.817`, random `0.5123`) must be blocked; three real values must pass. Guards the gate itself (the old 1e-9 absolute tolerance was vacuous for tiny P-values). |
| `S6v10_gate_allow_probe.py` | saturation audit of the ledger-value pool (documents why the global-allowlist gate was replaced by the index gate). |

Delivered run (2026-11-05): gate 8/8 OK, probes 4/4 blocked + 3/3 pass.
