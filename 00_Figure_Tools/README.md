# 00_Figure_Tools — the figure-QA tooling the paper's PDFs were built and checked with

Nothing here draws a result. These are the tools that make a figure *publishable and checkable*:
they set the fonts, splice pages, and — the part worth publishing — **prove that a change did what
it says and nothing else**.

| File | What it is for |
|---|---|
| `figkit.py` | the house bootstrap: locate the persistent PyMuPDF (`DSCI_PYLIBS`) and the Arial TTFs (`DSCI_FONTS`), expose `FONT`/`CM`/`write_spans`/`use_matplotlib_arial()`, and the font audit used below. Falls back to metric-compatible Liberation Sans when no Arial is present, and says so. |
| `final_audit.py` | per-figure gate table: page size, **minimum font size**, panel letters, non-Arial fonts, un-embedded fonts, **hidden characters** (ink without text / text without ink), subsetting, Type 3. |
| `inkcheck.py` | second-pass overlap check: renders at 600 dpi, takes each span's *tight ink* box, and lists pairs overlapping > 0.3 pt. A candidate list, not a verdict — see its docstring for why the connectivity test alone is too strict for 6 pt type. |
| `pdfix.py` | locate a run of text in a page and rewrite it in place (used where the text lives in the page stream). |
| `arialize.py` | replace non-Arial text with Arial **inside each Form XObject**, in place, preserving z-order and occlusion; `strip_everywhere()` removes the text layer from a page *and* all nested forms. |
| `diaudit.py` | drawing/vector audit of a page (paths, clips, ExtGState) — how we established facts like "no alpha is applied to network edges". |
| `port_paths.py` | the mechanical pass that removed hard-coded machine paths from every script in this repo (see below). |

## The gates every delivered figure passes

- 18.00 cm wide unless the sheet was delivered otherwise; minimum font **≥ 5.5 pt**; Arial only,
  embedded, **subsets**, no Type 3; panel letters 8 pt regular; no hidden characters.
- Every edit is accompanied by a pixel self-check: render before and after at 150 and 300 dpi, mask
  the intended footprint(s), and assert that **nothing outside the mask changed** (0 px). A
  grey-level threshold of `|Δ| > 24` is used rather than ink/non-ink, because 249-vs-250
  anti-aliasing flips are not visible differences.

## Four traps that cost real time, and are now guarded in code

1. **Never write Arial through PyMuPDF's `TextWriter`/`insert_text` with the full TTF.** The
   generated `/ToUnicode` maps spaces to U+00A0 and hyphen-minus to U+00AD, so the text layer reads
   back wrong (`pdfix.rewrite` produces exactly that). Either subset the font first with fontTools
   and `insert_text(fontfile=subset)`, or use `figkit.write_spans`, which appends word by word.
2. **`apply_redactions` decides by bbox intersection, not containment.** Two lines whose bboxes
   overlap by 0.05 cm cannot be redacted separately — deleting the middle one also deletes its
   neighbours (measured: 16 extra lines in Figure S4D; 0.02 cm of pad on a redact box was enough to
   eat the line above in the same figure). Measure the gap before choosing redaction over re-render.
3. **Compare against a freshly opened file.** Rendering the in-memory page object you just edited is
   not the same as rendering it before: after a redaction on a pasted page, the source page's scatter
   panel re-rendered 931 px differently in-process. A self-check that compares a product with itself
   passes vacuously.
4. **Paste offsets must sit on the rasterisation lattice.** Off-grid pastes shift glyph phase and
   make untouched panels "change" under a pixel diff. Snap to a multiple of 0.72 pt (1/100 inch),
   which is common to 200 and 300 dpi.

Also: a matplotlib script that forgets `figkit.use_matplotlib_arial()` silently falls back to DejaVu
Sans, whose taller boxes make unrelated labels collide in `inkcheck`.

## `port_paths.py`

Every script in this repo resolves its roots from environment variables, with the author's machine
as the default:

```
DSCI_PYLIBS   DSCI_FONTS   DSCI_FIG_ROOT   DSCI_SEURAT_ROOT
DSCI_SCENIC_ROOT   DSCI_SCENIC_SCRIPTS   DSCI_FY_ROOT   DSCI_REPO_ROOT
```

Run it as `python3 00_Figure_Tools/port_paths.py <dir>...`; it is idempotent (a file that already
mentions `DSCI_` is skipped) and refuses to port itself — its own source contains the literals it
replaces, and the first version of this tool rewrote its own replacement templates into
`getenv('X', getenv('X', …))`. Files that rewrite *other* people's source (e.g. notebook generators
that `str.replace` a hard-coded path out of a generated cell) are listed in `SKIP`, because porting
the search string would break the match.
