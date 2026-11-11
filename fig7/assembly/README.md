# assembly/ — panel code and the sheet-level steps that stay outside the notebook

**Current (2026-11-17, v17 FINAL).**

```
make_fig7_v17_panels.py    # the four panels, one PDF each ("# %% PANEL A/B/C/D" blocks)
arial_style.py             # the manuscript-wide font policy (Arial / Liberation Sans, two tiers)
```

`make_fig7_v17_panels.py` is the single source of the panel code: `build_fig7_notebook.py` (one level up)
splits it on the `# %%` markers and emits `Figure_7.ipynb`, so the notebook and the script cannot drift.
Edit code here, prose in the builder, then re-run `python3 build_fig7_notebook.py`.

```bash
cd fig7
DSCI_WORK_ROOT=/path/to/work W_A=6.36 H_A=6.30 W_B=4.00 H_B=5.98 W_C=5.90 H_C=5.80 FIG7D_FW=17.00 \
  python3 assembly/make_fig7_v17_panels.py
# fig7_out/fig7v17_{A_donors,B_adjranks,C_threeaxes,D_targets}.pdf
```

What this directory deliberately does **not** contain: the composed sheet, panel letters, and the QA
gates. Those run on the composed page and live with the working tree
(`figures/SUPPL_ANALYSES_20261112/s38_compose_fig7_2col.py` → `fig7_composed_v17_2col.pdf`,
17.00 × 14.34 cm, A | B | C on row 1 and D on row 2; gates `font_gate.py --strict`,
`collide_gate.py`, `verify_text_on_data.py`, all zero findings on the delivered page).

Fonts: the compose scripts resolve Arial in the order `ARIAL_TTF` (env) → `$DSCI_WORK_ROOT/fonts/ARIAL.TTF`
→ `./fonts/ARIAL.TTF` → system Arial; Arial is not redistributed by this repository.

## Decision history (kept because the numbers changed, not only the layout)

| date | decision | consequence |
|---|---|---|
| 2026-11-05 | Panel B = a single raw rank band (author decision, v16) | B colour was the unadjusted rank; MS 2 / 4 / 7 / 2 |
| 2026-11-17 | Panel B re-cut as a **portrait** heatmap on the ranks obtained after residualising Δ1 on ln(regulatory-domain size) | B colour now MS 1 / 1 / 3 / 1; the raw ranking, the ±15 % matched-window control and the ρ strip move to Fig S7; the word for the adjustment never appears on the face (terminology gate), it is stated in the caption |
| 2026-11-17 | the GSE231585 `*` moved off the tile title onto the neuron−glia violin it qualifies; wording corrected from "thinned to 717" to "the cohort's 717 annotated neurons" | see `LEGEND_FIG7_final.md` (C) and `upstream/README.md` Step 2 |
| 2026-11-17 | Panel D (18 audited targets) joined Figure 7 and its inputs were vendored to `tables_v17/` | the notebook gained section 4. The donor-level re-test of the S1PR5 detection rule was computed the same day and, on the author's ruling, is **not** part of Figure 7: it was handed to the supplement (`figS7/handoff/donor_s1pr5_20261117/`) |
| 2026-11-18 | caption-only correction in `LEGEND_FIG7_final.md` (D): the sentence "both over the three cohorts that carry the required annotations" was wrong for the left block | the delivered panel draws **mature OL vs OPC over four cohorts** (`DS4`) and only **lineage vs non-lineage over three** (`DS3`), because GSE277430 is an OL-subset release with no non-lineage cells; the hatched (`not in matrix`) and pale (q >= 0.05) keys are now also stated in words. No figure file changed |

## Superseded (2026-11-06 and earlier, archived, not pushed)

`_archive_v13_v16/assembly_old/` (this working tree only — not versioned since 2026-11-25) holds
`make_fig7_v13_column.py`, `make_fig7_B_v16_rawband.py` and
`compose_fig7_v13_column.py`, which produced the single-column sheet (7.93 × 13.98 cm, limit 8.0 cm).
Their outputs are in `_archive_v13_v16/fig7_out_old/`. Copies of the same files remain in the working
tree (`scripts/routeA/`, `figures/GITHUB_PUSH_READINESS_20261111/`), so nothing is lost; the archive is
excluded from the push by the repository `.gitignore` (`fig*/_archive*/`).
