# Repo tidy-up receipt — 2026-11-11 (this environment, branch `main`)

Commits: `502ac28` (path port), `c2b636d` (identity / parallel-env sync / docs / gate),
`b4f14cc` (fig6 config fix), `21e3783`/`9e2f82a` (this file and its follow-up).
State after: **486 tracked files, 19 commits, clean working tree,
`bash dev/push_preflight.sh` → PASS 9 / WARN 2 / FAIL 0.**

## 1. Author-machine paths: gone (1245 rewrites across 175 files)

`00_Figure_Tools/port_paths.py` was rewritten (v3) and run over the whole repository.
Policy now: **no absolute author-machine path may appear anywhere, not even as a default value.**

| form in code | after |
|---|---|
| `'<work>/FY/oligo/x.h5ad'` (literal) | `os.path.join(os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT',''), 'FY')), 'oligo/x.h5ad')` |
| `os.environ.get('DSCI_PYLIBS', <literal default>)` | `os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT',''), 'pylibs'))` |
| `${DSCI_FIG_ROOT:-<literal>}/routeA_out` | `${DSCI_FIG_ROOT:-${DSCI_WORK_ROOT:-}/figures}/routeA_out` |
| `Sys.getenv("X", <literal>)` | `Sys.getenv("X", file.path(Sys.getenv("DSCI_WORK_ROOT",""), "…"))` |
| `/home/<account>/.local/bin/scdrs` | `$HOME/.local/bin/scdrs` / `os.path.expanduser('~/.local/bin/scdrs')` |
| docs, notebook outputs, markdown cells | `<work>` / `<repo>` / `<raw>` / `<refdata>` placeholders |

The variable contract was **unified on `DSCI_WORK_ROOT`**, the name already documented in the
root README; a second root name introduced by an earlier pass was renamed away. Notebooks were
rewritten through `json`, so code cells got code-form rewrites while markdown cells and stored
outputs got placeholders. Python / R / shell syntax and notebook JSON were re-parsed afterwards
(0 failures; 179 notebook code cells checked).

Two files are on the port tool's `SKIP` list and were hand-ported instead:
`figS5/assembly/figs5_repanels.py` (it rewrites *other* scripts by exact text match, so its
search strings must stay literal — the legacy prefixes now live in one `DSCI_LEGACY_ROOT`
variable) and `00_Preprocessing/RNA_soupx_singlet.ipynb` (R `glue()` strings).

The only file still containing machine paths is `port_paths.py` itself, because its MAP keys
*are* the search strings; it is the first entry of `dev/PATH_LEAK_ALLOWLIST.txt`, together with
the gate and its token list (same reason: inspection only, no I/O on those paths).

## 2. Identity: masked

- WDLs: 34 private image tags → `<registry-namespace>/<image-tag-masked-for-upload>:latest`;
  cluster project id and the two personal accounts → `<project-id>` / `<user>`.
  Restore procedure is documented in `00_Preprocessing/README.md`. The institutional mount root
  is kept (it identifies the facility, not a person) — flagged for the authors in §4.
- The conda environment name that encoded a machine account → `<conda-env>` (12 files; the
  python/scdrs version strings next to it are kept, they are the reproducible part).
- 33 workspace paths in READMEs/scripts → `<work>`; notebook stderr lines that printed a
  per-user `site-packages` path → `[site-packages]`.
- Gate G1 reads its forbidden list from `dev/identity_tokens.txt` — currently 0 hits in the tree.
  **(2026-11-25 update:** the repository became public that day, so the verbatim list moved out of
  git — it is now gitignored and local-only; the versioned gate input is the *structural* pattern
  list `dev/identity_tokens.example.txt`, which catches accounts/project numbers/mount prefixes by
  shape instead of naming them.)

## 3. Sync from the parallel environment (its 2026-11-10 snapshot)

`merge_from_parallel.sh --apply` was executed here: **96 same-name files taken from that side,
103 new files added (~18 MB)**, including the previously missing `figS7/` package and `dev/`.
The incoming `.gitignore` was *not* taken (it had lost 15 binary-exclusion rules); the local one
was kept and extended. Layout migration honoured: root generators for `fig4`/`figS4`/`figS6`
removed in favour of `dev/`, `fig6`/`fig7` keep theirs. Everything merged was then re-ported
(§1) and re-masked (§2), so the merged tree is clean rather than pending.

## 4. Structure and self-consistency

- 6 tracked absolute symlinks replaced by `fig5/data/link_data.sh` + `figS5/data/link_data.sh`
  plus provenance tables with sizes and md5 prefixes (`fig5/data/README.md`, `figS5/data/README.md`).
- Untracked, files kept on disk: `*/_superseded_20261023/` binaries (including the 70 MB liftOver
  chain), 9 `regen.log`, the 64 MB `figS5c_nrxn3_locus.pdf`. Only the two archive READMEs stay tracked.
  `.gitignore` gained `_superseded*/`, `regen.log`, `*.pkl` and the two `data/` directory rules.
- Vendored 3 helper modules that repo code imported but that existed only outside the repo:
  `00_Figure_Tools/figaudit.py`, `00_Figure_Tools/font_census.py`,
  `figS5/assembly/compose_figS5_noc3.py`, with `sys.path` bootstrap so `figkit`/`network_key` resolve.
- `fig6/upstream/A5_fig6_v15_stats.py`: deleted the `from config import CONFIG` line that shadowed
  the inline `CONFIG` with a module that ships nowhere (`ModuleNotFoundError` on a fresh clone).
- New docs: `dev/LIVE_GENERATORS.md` (one live generator per notebook + the figS7 exception),
  `dev/TABLES.md` (Tables 1, 2 and S1–S8: content, generator, delivered file, md5),
  `dev/PATH_LEAK_ALLOWLIST.txt`, `dev/identity_tokens.txt` (local-only since 2026-11-25; the shipped
  counterpart is `dev/identity_tokens.example.txt`), `00_Figure_Tools/env_check.py`.
- `figS5/assembly/make_table_s6_rss.py` → `make_table_s5_rss.py` and its 4 references updated, so
  the file name no longer contradicts the manuscript numbering; stale numbering phrasing fixed in
  the root README and in `figS5`/`figS6` docs.

## 5. Gate

`bash dev/push_preflight.sh [--verbose]` runs G1 identity tokens, G2 machine paths,
G3 outside-repo symlinks, G4 junk/caches, G5 + G5b file sizes, G6 syntax (Python, notebook JSON
and code cells, R, shell), G7 local-import completeness, G8a/G8b the environment-variable
contract, G9 stale supplementary-table numbering, G10 commit identity. Any FAIL exits non-zero.

Accepted WARNs (not blockers): two tracked files in the 10–50 MB band
(`fig3/fig3_out/joint_meta.csv.gz` 24.5 MB, `figS2/figS2_out/S2a_batch_shadow_facets.pdf` 12.0 MB)
and G10 (commit identity + commit messages mentioning machine paths).

## 6. What cannot be done from this environment

1. **Push**: no credentials here (no `GH_TOKEN`/`GITHUB_TOKEN`, no SSH key, no upstream on `main`).
   The first push is `git push -u origin main`.
2. **History decision**: 18 commits, all authored by the same real name/email, `.git` ≈ 223 MB,
   never pushed. This is the cheapest moment to squash into one clean initial commit, or to run
   `git filter-repo` to rewrite author/committer and drop machine paths from commit messages.
   The choice — signed with real identity, or anonymous via the GitHub noreply address — is the
   authors'. Once decided, this environment can execute it before the first push.
3. **Re-run verification** belongs to whoever owns the analysis tree: see
   `RECEIPT_PARALLEL_ENV_SCRIPTS_20261111.md` (§1) for the notebook list and acceptance criteria.

## 6. History squash (2026-11-11, decision ② of work order §4-3)

All 22 commits (20 pre-push + this round's fixes) were squashed into **one clean
initial commit** after this receipt was written; commit author/committer is
`HZJ714 <HZJ714@users.noreply.github.com>` (GitHub noreply; no private email).
The pre-squash history is preserved byte-for-byte in a git bundle held outside
the repository (authors' tree: `github_prep/history_full_20261111.bundle` (analysis-tree root, outside the repo)),
which also proves the audit trail below. The gate's G10 shows exactly two items
by design: the noreply identity and the LICENSE line — both are author sign-off
items, not leaks.

Pre-squash `git log` (22 commits, oldest first):
  01d2f94 dsciCAT-seq: reproducible figure notebooks (12 figures + 1 preprocessing) and upstream dsciCAT/ISSAAC WDL workflows
  634fb83 add heavy intermediates (1/3): figS5c_nrxn3_locus.pdf (63.8 MB genome-browser locus PDF)
  7eae15b add heavy intermediates (2/3): joint_meta, batch-shadow/S2a, HT_harmony_meta, figS3ab umaps, joint_metrics_sub, fig2 umaps, tenx_meta, truth_annotCOP, fig3abc umaps
  7dd08fc README: title -> dsciCAT-seq analysis pipeline; drop stale gitignore entries for re-added intermediates
  224de33 Add MIT license
  fc74e07 figS6 Panel I: interval-overlap fix (10,322 -> 10,323, matching the delivered panel and manuscript); pin chain/peak md5s in README
  10d1b7a Fold in reviewer round 2: Figure S3 H-J (pre-registered COP/NFOL split test), Fig5/S5 network keys + unified eRSS + Table S6, S5C unit/p-value corrections; add 00_Figure_Tools (figure-QA gates + traps) and 01_Legend_Edits (anchored docx surgery); de-hardcode machine paths via DSCI_* env vars
  05273ed Figure 1 finalised as v4 (G0/G1 installed; A-panel recolour and figure slimming declined by authors); add fig1/assembly chain, sync the stale figures.docx Figure 1 legend, extend port_paths to the FIG1_DATA root and fix its idempotency guard (was matching path prefixes, so nothing ported)
  b668f89 Text-side finalisation round (2026-10-30): manuscript v7->v8 legend/body fixes ([0076]x4, [1311], [1319], [1320], [1348]/[1352], [1355], S4C columns), sync figures.docx Figure 2/3/4/5/S1-S5 blocks to manuscript at run level (keeps docx styles, carries italics); Fig 6/S6 blocks deferred to finalisation
  17aba0a Trim repo: drop duplicate PNGs (keep vector PDFs; readers regenerate rasters from code) and remove 01_Legend_Edits (docx-only legend surgery, not in the figure reproduction chain); sync README references
  4ba1cd9 Sync final Fig4/S4/6/7 code, panels, legends from parallel env (2026-11-05)
  70d46cc Fix Fig7 hard flaws #1-#3 in the unified copies (2026-11-05)
  a0e0a26 Sync finalized Fig S6 code + legend from parallel env (2026-11-05)
  15e0e82 Fix S6-family disclosures + write parallel-env action receipt
  502ac28 Port all author-machine paths out of the repo (port_paths v2 + 1236 rewrites)
  c2b636d Unify repo for first public push: single path contract, identity masking, generator/table docs, preflight gate
  b4f14cc fig6: drop the shadowing 'from config import CONFIG' (no such module ships; A5 could not run on a fresh clone)
  21e3783 dev: add this-side tidy receipt (paths/identity/sync/structure/gate, 2026-11-11)
  9e2f82a docs: resolve stale numbering decisions, finish DSCI_WORK_ROOT rename, gate self-exemption
  011fc98 env_check: RESULTDATA optional; gate: exempt numbering-mapping docs
  d5f40c8 post-audit unification: pin table inventory (S4 norm-score/S5 eRSS sources, D1/D2), mark M2b + v15 sheet title internal (D7/D9), drop 2 heavy regenerable intermediates (joint_meta, S2a facets), add dsci-repro environment files, merge S6 numbering-mapping legend header, clean ignored checkpoints
  003e6bf legend: avoid G9 false positive in numbering-mapping header
