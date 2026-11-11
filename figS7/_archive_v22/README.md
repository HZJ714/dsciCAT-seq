# Superseded deliveries, kept for provenance, not pushed

Everything here is the **v22** Fig S7 (17.00 x 16.81 cm, panels A-I, signed 2026-11-10) and the
scripts that made it. It is kept because the delivered manuscript supplement and several audit
receipts reference those PDFs and those face numbers, and because the v29 legend's section 2 quotes
numbers that only exist in the v22 panels.

| Item | What it was |
|---|---|
| `figS7_out_v22/` | the v22 sheet, its nine standalone panel PDFs, the face-number CSV and the frames JSON |
| `compose_figS7.py` | the v22 sheet composer (ten panels + LEG band, letters baked into the panels) |
| `verify_figS7.py`, `check_figS7_repro.py` | the v22 gate battery and the point-by-point reproduction check |
| `LEGEND_S7_v22.md` | the v22 legend and face-number contract, verbatim |

The current version is **v30**: the ten panels in
`../figS7_out/panels_v30/`, drawn by `../assembly/make_figS7_v30_panels.py`; the composed v30 sheet
is assembled from those panels and is not versioned. The v29 intermediate
is itself archived next to this directory (`../_archive_v29/`). The GitHub root `.gitignore` excludes `fig*/_archive_*/` contents (only these archive `README.md`
files are versioned), so nothing in this directory ships to the public repository.
