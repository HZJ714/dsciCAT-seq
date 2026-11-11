# dev/ — notebook generators and other build-time provenance

Not part of the reproduction path. A reviewer never needs anything in here: each
`figX/Figure_X.ipynb` runs end-to-end on its own, and the delivered artifacts in
`figX/figX_out/` were produced by those notebooks plus the scripts in `figX/assembly/`
and `figX/upstream/`.

What is here:

| file | builds | notes |
|---|---|---|
| `build_fig4_notebook.py` | `../fig4/Figure_4.ipynb` | |
| `build_fig6_notebook.py` | — (2026-10-23 layout, superseded) | the live generator is `../fig6/build_fig6_notebook.py` |
| `build_fig7_notebook.py` | — (2026-11-05 draft) | the live generator is `../fig7/build_fig7_notebook.py` |
| `build_figS4_notebook.py` | `../figS4/Figure_S4.ipynb` | |
| `build_figS6_notebook.py` | `../figS6/Figure_S6.ipynb` | |
| `identity_tokens.example.txt` | — (gate input) | the **structural** patterns G1 matches by shape (home/users paths, project-number shape, mount-prefix shape, non-noreply e-mails). The verbatim list lives in `identity_tokens.txt`, which is deliberately gitignored: it names real accounts, so it must not end up in a public tree or in the Zenodo archive. G1 uses both when the local one exists.
| `preflight_checks.py` | — (gate) | the PASS/WARN/FAIL report run before every push; `push_preflight.sh` wraps it |
| `backups/` | — | pre-relabel copy of the superseded v13 compose script; **not versioned** (2026-11-25): it stays in this working tree as a fallback, `dev/TABLES.md` §Version policy |

Why they were moved out of the repository root (upload fixlist, 2026-11-10 §3): they are
generators, not figure code, and keeping them at the top level made the public tree look
like it needed them. They are kept rather than deleted because every notebook in this
repository is reproducible from its generator — cell sources, markdown and all — and the
byte-for-byte claim in the root README is checked by regenerating.

Paths: the generators write the notebook into its own figure directory, resolved from
`__file__`, and the cell text they embed reads the analysis tree through `DSCI_WORK_ROOT`
exactly like the notebooks do (root README, *Data layout and environment variables*). Run
them from this directory:

```bash
cd dev && python3 build_figS6_notebook.py      # -> ../figS6/Figure_S6.ipynb
```
