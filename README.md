# dsciCAT-seq analysis pipeline

Reproducible analysis pipeline for dsciCAT-seq / ISSAAC-seq dual-omics data:
the upstream FASTQ → count-matrix WDL workflows and RNA preprocessing
(`00_Preprocessing/`), plus fifteen self-contained Jupyter notebooks
(14 figure + 1 preprocessing) that reproduce every main- and
supplementary-figure panel of the paper **from the delivered data** — each one
runs end-to-end on a fresh machine (R ≥ 4.1 + Python ≥ 3.8) and writes the
figure files to its own `figX_out/` directory.

| directory | notebook | panels | executed end-to-end |
|---|---|---|---|
| `fig1/` | `Figure_1.ipynb` | 1a–1f | ✅ 2026-10-23, 0 errors |
| `fig2/` | `Figure_2.ipynb` | 2a–2e | ✅ 2026-10-23, 0 errors |
| `fig3/` | `Figure_3.ipynb` | 3a–3f | ✅ 2026-10-23, 0 errors |
| `fig4/` | `Figure_4.ipynb` | 4a–4d | ✅ 2026-11-05 (synced to the 2026-10-28 final panels), 0 errors |
| `fig5/` | `Figure_5.ipynb` | 5a–5f | ✅ 2026-10-22, 0 errors |
| `fig6/` | `Figure_6.ipynb` | 6a–6c (one row, 17.00 × 4.89 cm) | ✅ 2026-11-17 (**v36 final**: A score / B UMAP / C cross-channel ladder; Manhattan withdrawn to Table 2, cross-chemistry panel moved to Fig S6G; row composed by `assembly/fig6v36_abd_final.py`) |
| `fig7/` | `Figure_7.ipynb` | 4 single panels (A–D); the 17.00 × 14.42 cm two-column sheet is composed outside this repo | ✅ 2026-11-19 (**v17d final**: row 1 = A donors | B size-deconfounded rank heatmap | C three axes, row 2 = D drug-target × cell-type chip; the rotated-label re-anchoring that makes 17d differs from 17c; v13/v16 sheets live in this working tree under `fig7/_archive_v13_v16/`, README only in git) |
| `figS1/` | `Figure_S1.ipynb` | S1a–S1e | ✅ 2026-10-22, 0 errors |
| `figS2/` | `Figure_S2.ipynb` | S2a–S2d | ✅ 2026-10-23, 0 errors |
| `figS3/` | `Figure_S3.ipynb` | S3a–S3g | ✅ 2026-10-23, 0 errors |
| `figS4/` | `Figure_S4.ipynb` | S4a–S4d | ✅ 2026-11-05 (synced to the 2026-10-28 final panels), 0 errors |
| `figS5/` | `Figure_S5.ipynb` | S5a–S5e | ✅ 2026-10-22, 0 errors |
| `figS6/` | `Figure_S6.ipynb` | S6a–S6g (seven panels) + the Table S6 ledger | ✅ 2026-11-17 (**v16 final**, four-row sheet 17.00 × 18.91 cm composed outside this repo: A loci / B permutation / C UMAP bins / D top-1000 landing / E cis-window landing / F per-nucleus boxes / G cross-chemistry Δ1 pair with shared vertical legend; v11 2026-11-06 sheet superseded) |
| `figS7/` | `Figure_S7.ipynb` | S7 A–G + K/RB/L2 (one 17.00 × 14.28 cm sheet) = 10 single-panel PDFs | ✅ 2026-11-18 (**v30 final**, ten asserted panel PDFs in `figS7_out/panels_v30/` + the composed v30 sheet, **not versioned**; the v22 A–I sheet and the v29 intermediate are superseded and kept locally in `figS7/_archive_v22/` / `_archive_v29/`; `handoff/donor_s1pr5_20261117/` holds the generator of the per-donor S1PR5 tables quoted in the Fig 7D legend) |
| `00_Preprocessing/` | `RNA_soupx_singlet.ipynb` + WDL workflows | upstream: FASTQ → matrices → RNA QC/doublets | executed upstream (KP12 demo; FY13/14 via env vars) |

Each notebook is regenerated from exactly **one** live generator (regenerate → notebook,
byte-for-byte cell sources): most live in `dev/`, several stay in their package root and
`figS7` has none — the authoritative table is
[`dev/LIVE_GENERATORS.md`](dev/LIVE_GENERATORS.md). Table numbering (Tables 1, 2, S1–S8)
and where each delivered table comes from is pinned in [`dev/TABLES.md`](dev/TABLES.md).
Each `figX/README.md`
documents the panel → input → output map, the run recipe (env vars, versions) and
a validation checklist.

## Data availability

All input data are public at **CNGBdb, project CNP0010263**. The notebooks
only read what the delivered analysis produced and never bundle raw data into
this repository; large intermediate objects (`*.RDS`, `*.mtx*`, …) are
gitignored and rebuilt when a notebook runs. See each figure's README for the
exact input list and sizes.

## Data layout and environment variables

Nothing in this repository is configured with a machine-specific absolute path.
Every script and notebook resolves its inputs through environment variables whose
defaults are **relative to the current working directory**, so a plain `cd figX` run
reads and writes inside the checkout, and a reviewer who downloaded the delivered
tables elsewhere exports one variable instead (no variable may hold a hard-coded
author path — see the policy at the end of this section):

| variable | points at | used by |
|---|---|---|
| `DSCI_WORK_ROOT` | the root of the analysis tree (`FY/`, `figures/`, `routeA_probe/`, `fonts/`, `scripts/`) | upstream/provenance scripts, notebooks' external-data defaults |
| `DSCI_ONLINE_ROOT` | the same analysis root under its original name, used by the **fig6 upstream** steps that read the ArchR project (`<root>/FY/archr-subset-copnfol`); default `.` | `fig6/upstream/D1_export_genescore_archr.R`, `D2d_read_rna_from_arrows.py`, `D2f_pack_rna_seurat.py` |
| `DSCI_FY_ROOT` | `<work>/FY` — the oligo/dsciCAT-seq processing tree | fig6, figS6, upstream scoring scripts |
| `DSCI_FIG_ROOT` | `<work>/figures` — delivered figure/analysis output (`routeA_out/…`) | fig6, figS6, cross-chemistry scripts |
| `DSCI_SEURAT_ROOT` | `<work>/pi/plots-seurat` (Seurat replot outputs) | fig5, figS2 |
| `DSCI_SCENIC_ROOT`, `DSCI_SCENIC_SCRIPTS` | `<work>/pi/plots-scenic+`, `<work>/pi/scripts-scenic+` | figS5 and its assembly steps |
| `DSCI_SCRIPTS_ROOT` | `<work>/scripts` (analysis scripts, incl. the MS GWAS sweeps) | fig6 / fig7 upstream steps |
| `DSCI_PROBE_ROOT`, `DSCI_REF_ROOT`, `DSCI_MOTIF_ROOT` | `routeA_probe/`, `referrence/`, `aertslab_motif_colleciton/` | fig6 probes, liftOver and homology inputs |
| `DSCI_FONTS_NAMED` | `<work>/fonts_named` (one file per weight) | `00_Figure_Tools/arialize.py` |
| `FIG1_DATA` | raw count matrices (`ResultData/…`), **not distributed** | fig1, `00_Preprocessing/` |
| `DSCI_REFDATA_ROOT` | STAR/chromap indexes, GTF, TSS BED, bead whitelist | `00_Preprocessing/` WDLs |
| `DSCI_REPO_ROOT` | this checkout (default `.`) | scripts that address another package |
| `DSCI_LEGACY_ROOT` | the author tree whose scripts `figS5/assembly/figs5_repanels.py` patches by text | that one script only |
| `DSCI_HOMOLOGENE` | `homologene.data` (human↔mouse orthology, panel D of fig6) | `fig6/Figure_6.ipynb` |
| `DSCI_PYLIBS` | optional extra site-packages directory (nothing is bound to a private venv) | `fig6/assembly` gates, `figS6/assembly/support_common.py` |
| `DSCI_PY` | the interpreter of the analysis/scoring environment (python 3.9 + scdrs 1.0.2) — the authors' own env is not published and is never named in the tree | `figS6/assembly/figS6v18_row3.sh`, upstream scDRS/probe steps, every `scdrs` call in the READMEs |
| `ARIAL_TTF`, `FIG7_FONT`, `FONT_DIR`, `DSCI_FONTS` | an Arial TrueType file / a directory of them | every figure; order is env → `./fonts/ARIAL.TTF` → system Arial |
| `FIG*_OUT`, `S6_OUT`, `FIG7_MASTER`, `S*_DATA` | per-figure output / delivered-table overrides | the figure's own notebook |

External data are never hard-coded: `FY/oligo`, the ArchR `ArrowFiles`, `homologene.data`
and the GEO probe downloads have to be supplied by the user (CNGBdb **CNP0010263** plus the
GEO accessions named in each README), and the path is given through the variables above.
The authors' own runs set `DSCI_WORK_ROOT` to their analysis tree; the delivered artifacts
in `figX_out/` do not depend on it.

Check what is still missing on your machine:

```bash
python3 00_Figure_Tools/env_check.py     # required/optional variables, whether set, what they resolve to
```

Large input objects are not in git either. `fig5/data/` and `figS5/data/` used to be
absolute symlinks into the authors' disk; they are now created on demand:

```bash
export DSCI_WORK_ROOT=/path/to/your/analysis/tree
bash fig5/data/link_data.sh
bash figS5/data/link_data.sh
```

## Reproducible environment (`dsci-repro`)

The delivered figures were produced under **python 3.9.23 / scanpy 1.10.3 /
anndata 0.10.9 / pandas 2.3.1 / numpy 1.26.4 (import-time)** with conda R 4.4.3 —
see `environment.yml` and `requirements-dsci-repro.txt` at the repository root for
the full pinned set (231 conda + 49 pip packages; the `.explicit` conda export used
to build it is kept in the authors' tree).

```bash
conda env create -f environment.yml        # creates an env named `dsci-repro`
conda activate dsci-repro
python -m ipykernel install --user --name dsci-repro   # register the kernel
```

Re-running any notebook must use the **`dsci-repro`** kernel
(`nbconvert --execute --ExecutePreprocessor.kernel_name=dsci-repro`): the drawing
output depends on the versions above, and pandas 1.5 → 2.3 API differences would
silently change numbers under the wrong kernel. The one R notebook
(`00_Preprocessing/RNA_soupx_singlet.ipynb`, `ir` kernel) uses the same conda R
4.4.3.

Provenance, sizes and md5 prefixes of every such object are tabulated in
`fig5/data/README.md` and `figS5/data/README.md`.

**Policy (2026-11-11).** No author-machine absolute path may appear anywhere in this
repository, not even as a default value: `00_Figure_Tools/port_paths.py` rewrites them
into the variables above, and `dev/push_preflight.sh` fails the push if a new one is
introduced (the only exempt file is `port_paths.py` itself — see
`dev/PATH_LEAK_ALLOWLIST.txt`).

## Running a notebook

```bash
# example: figS2 (see each figX/README.md for that figure's env vars)
cd figS2
S2_DATA=/path/to/delivered/data S2_OUT=figS2_out \
  jupyter nbconvert --to notebook --execute --inplace --allow-errors \
  --ExecutePreprocessor.timeout=-1 Figure_S2.ipynb
```

## Requirements

- **R ≥ 4.1**: Seurat 4.x, Signac, harmony, ArchR (fig1/fig2/fig3/S1/S5),
  aricode ≥ 1.0.1 (figS2), ggtext (figS2), ggplot2/dplyr/tidyr/patchwork/pals
- **Python ≥ 3.8**: numpy, pandas, scipy, matplotlib, seaborn, nbformat,
  jupyter/nbconvert
- Versions actually used are printed by each notebook's §0.x versions cell.

## Layout

- `figX/Figure_X.ipynb` — the executed, self-checking notebook (0 errors)
- `dev/LIVE_GENERATORS.md` — which single generator produces which notebook (incl. the figS7 exception)
- `dev/TABLES.md` — Tables 1, 2, S1–S8: content, generator, delivered file, md5
- `dev/push_preflight.sh` — the pre-push gate (identity / paths / numbering / syntax / sizes)
- `figX/README.md` — panel/input/output map + run recipe + validation
- `figX/data/` — symlinks to large non-distributed inputs, created by
  `bash figX/data/link_data.sh` (`fig5/data/README.md`, `figS5/data/README.md` list them)
- `figX/figX_out/` — the **currently final** figure artifacts (panel PDFs + the CSVs the
  notebook writes and runs) and the R scripts it calls. Superseded sheets are not
  versioned: they stay in this working tree under `figX/_archive_*/` (whose `README.md`
  says what was archived and why) and under `_superseded_20261023/`. Composed journal
  pages are assembled outside the repository, so the repository ships panels plus the
  compose scripts that rebuild them.
- `figX/assembly/` — the scripts that assemble, patch and check the *published*
  PDF for that figure (page splicing, in-place text surgery, pixel self-checks).
  They do not re-run any analysis.
- `figX/upstream/` — the scripts that produced the *input tables* the notebook
  reads (scDRS score sweeps, cross-chemistry statistics, locus windows, per-mouse
  contrasts). Each one's docstring states the criterion, formula, seed and reason
  for its CSV. They read the analysis tree through `DSCI_WORK_ROOT` / `DSCI_FY_ROOT`
  / `DSCI_FIG_ROOT` (see *Data layout*), so they are auditable and portable; they are
  provenance evidence and are **not** called by the figure — the runnable entry
  point is always the notebook. The data themselves stay out of
  git (CNGBdb CNP0010263 + public GEO).
- `00_Figure_Tools/` — the shared figure-QA tooling (`final_audit.py`,
  `inkcheck.py`, `arialize.py`, `pdfix.py`, `figkit.py`, `diaudit.py`) plus the
  gates every delivered page passes and the four traps that are now guarded in
  code: [README](00_Figure_Tools/README.md)
- `.gitignore` — large binaries (raw data / `*.RDS` / `*.mtx*` / `*.pkl`), superseded
  archives (`*/_superseded_20261023/`, README kept) and run logs stay out of git

## Reviewer-round-2 status (2026-10-29)

Figures **1, 2, 3, 5, S1, S2, S3, S5** are authored here; **4, S4, 6, S6 and the
graphical abstract** come from a parallel environment. What is folded in now,
and what is deliberately not:

- **folded in**: Figure S3 panels H–J (a pre-registered test of the "100k
  COP/NFOL splits into two clusters" claim, `figS3/s3c_cop2/`, with the frozen
  criteria and the unedited result side by side), Figure 5/S5 network keys and
  the unified eRSS panel (Table S5 in the current numbering, formerly S6), the Figure S5C unit and p-value corrections.
- **Figure 6 is final** as of 2026-11-03 and its notebook is rewritten for the delivered design
  (A scDRS boxplot, B UMAP hexbin, C cross-channel Δ1 over 74 traits, D the same criterion in three
  library chemistries). The 2026-10-23 notebook (A/B/Manhattan/74-trait strip) is kept in
  `fig6/_superseded_20261023/`; the Manhattan and the strip now belong to the supplement. The
  reproduction is bit-identical to the delivered tables for every statistic except the split-half
  noise ceiling of panel C (±0.002, a 20-draw median whose delivered value came from an already
  consumed RNG stream) — see `fig6/README.md` for the full checklist.
- **Figures 4 / S4**: the 2026-10-28 sheets are accepted as final (the 2026-11-02
  re-lettered re-composition in `figures/FINAL_SET_v2_20261028/` is the delivered
  pair; md5 `0b98…`/`703a…`). On 2026-11-05 both notebooks were folded in to that
  final state: Fig 4 C/D panels (signed 2ΔlogL axis + negative-count line, `93.6%
  of 1798`, `x/114` switch denominators, example genes Anln/Cemip2/Tcf7l2, p99
  display clipping, post-switch shading removed) and S4 (model-2 heatmap end label
  MOL6 → MOL3/6, `latent time (0 → 1)`, FDR underflow printed as `<1e-300`, widened
  non-cubic 3D boxes with the tc50/ts50 definition line). Output directory of the
  S4 notebook renamed `fig4_out/` → `figS4_out/`. The house text fixes of the
  delivered PDF (`patch_fig4_house.py` upstream) are noted in `fig4/README.md`.
  Ten `LRT` labels in Figure S4 still need a source re-render rename to `2ΔlogL`
  (in-place PDF patching is unsafe — overlapping bboxes); tracked in
  `figS4/README.md`, exact strings in `REPLY_FIG4_S4_AUDIT_20261029.md` upstream.
- **Figure 1 is now final** as `Fig1_v4.pdf` (items G0/G1 installed; the two remaining reviewer
  asks — recolouring the A-panel scheme and slimming the figure — were declined by the authors,
  see `fig1/README.md`). It is the one delivered figure that knowingly does **not** pass the font
  gates: 1.84 pt minimum inside the AI-master schematic, and the waiver is recorded item by item
  rather than left as a silent FAIL.
- **Figure S6 is final** as of 2026-11-06 (**v11**): the 10-panel 2026-10-23 notebook was
  replaced by the v10 seven-panel design (17.00 × 18.86 cm sheet + Tables S6/S7,
  strict index-level reverse-trace gate 8/8, fabricated-token probes blocked —
  see `figS6/README.md`). It now passes the font gates too (Arial only, minimum
  visible span 5.5 pt); the old S6's failures (5.4 pt minimum, one un-subsetted
  face) went with it. v11 adds the two trait-level adjustment rungs to panel F
  (ln n_hit 0.826, ln n_peak 0.708, B = 1,000 seed 20261104, independent re-seeded streams —
  the seven core rungs bit-for-bit unchanged) and the Table S6 v11 relabel (the "domain-size
  adjusted 0.823" row corrected to gene-set-size / ln n_hit; new ln n_peak row 0.708
  [0.547, 0.813]); notebook re-executed end-to-end 2026-11-06 with 0 errors, all 9 F rungs
  asserted against the delivered ledger.
- **Figure S7 is final** as of 2026-11-18 (**v30**): one 17.00 × 14.28 cm two-column sheet
  (assembled from the ten versioned panels; the sheet PDF itself is **not** versioned) = panels A–G plus the channel key band K and the legend bands
  RB/L2, drawn one asserted PDF at a time by `assembly/make_figS7_v30_panels.py`
  (`figS7_out/panels_v30/`); the notebook never composes, and the page splice plus the house
  gates (`render_v30.py`, `gates_v30.py`, `verify_v30_ink.py`) run on the analysis tree, as
  recorded in `figS7/README.md`. The relayout withdrew the old C/F/G panels into Table S7
  (`assembly/make_table_S7.py`) and re-lettered the survivors; the v22 sheet (A–I,
  17.00 × 16.81 cm) and the v29 intermediate are kept in this working tree under
  `_archive_v22/` / `_archive_v29/` and are **not** versioned — only `_archive_v22/README.md`
  is, as the pointer to what they were. Statistics are read from `data/` (row-by-row md5 in
  `data/MANIFEST.csv`) and from `tables_v29/S41_*.csv` (generated by
  `handoff/donor_s1pr5_20261117/s41_donor_all_human_cohorts.py`), never recomputed in the
  drawing code; the numbers deliberately kept off the face are recorded in
  `figS7/LEGEND_S7_final.md`, the on-face register is
  `assembly/register_figS7_faces.py` → `figS7_out/figS7v30_face_numbers.csv`, the geometry
  frames are `figS7_out/{panel_frames_v30,sheet_v30_frames}.json`, and the text-on-data gate
  is `figS7/verify_text_on_data.py`.
- **open**: the graphical abstract (2.21 pt minimum, un-subsetted face inside the
  schematic). Figure 6 passes (minimum 5.5 pt, Arial only — `fig6/assembly/font_gate.py`).

## Reproducibility hygiene (2026-11-11)

What changed while preparing the first public push, and how to keep it that way:

- **Paths**: 1245 author-machine paths rewritten (175 files) into the `DSCI_*` contract
  above; everything inside the repo is relative to `DSCI_REPO_ROOT` or the working directory.
- **Identity**: nothing in this tree spells out a person, a machine account, a cluster
  project number or the site's absolute mount prefix. Private Docker registry tags and the
  cluster project/account strings in the four WDLs are masked (`<registry-namespace>`,
  `<image-tag-masked-for-upload>`, `<project-id>`, `<user>`; restore procedure in
  `00_Preprocessing/README.md`), the institutional HPC mount root is written as
  `<hpc-mount>/`, the conda environment name that encoded a machine account is printed as
  `<conda-env>` (call it through `$DSCI_PY`), and notebook outputs and READMEs use
  `<work>`/`<repo>`/`<raw>` placeholders. Gate G1 enforces this by **shape** (see
  `dev/identity_tokens.example.txt`); the verbatim list used for extra strictness inside the
  analysis environment lives in `dev/identity_tokens.txt`, which is gitignored on purpose.
- **Sizes**: the 70 MB liftOver chain and the 64 MB Figure S5C vector PDF left git
  (md5-pinned regeneration notes in `figS6/README.md`); the 7.5 GB SCENIC+ object and the
  two Harmony `.RDS` files are reached through `link_data.sh`.
- **Self-consistency**: every local module imported by repo code is now inside the repo
  (`figaudit.py`, `font_census.py`, `compose_figS5_noc3.py` vendored 2026-11-11), and
  `dev/push_preflight.sh` re-checks imports, Python/R/shell syntax, notebook JSON, stale
  table numbers, symlinks and file sizes.
- **Not resolved here** (needs the authors): push credentials and the commit-identity /
  LICENSE wording (G10). The institutional mount root question, open in
  `dev/TIDY_RECEIPT_20261111.md`, was closed on 2026-11-25 in favour of masking: the
  repository went public and is being archived on Zenodo, so the mount prefix is now the
  `<hpc-mount>/` placeholder everywhere.

## License

MIT (see `LICENSE`). The input data remain under their own terms (CNGBdb,
project CNP0010263).
