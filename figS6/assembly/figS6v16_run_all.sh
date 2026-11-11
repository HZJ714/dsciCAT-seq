#!/usr/bin/env bash
# figS6v16_run_all.sh - one-shot re-render of the delivered Figure S6 (v16/v18): four rows, seven cells + vertical legend + compose + three gates
# ---------------------------------------------------------------------------------------------
# Why this script exists: row-3 geometry of v16 is SOLVED (not eyeballed) and row-4 G's column widths are a string of env
# tunings; both previously lived only in one-off command lines that the receipt could not reproduce. The commands now live here.
#
# Default output is IN-REPO (figS6_out/); inputs are read-only tables_v16/ + figS6_out/caches/,
# so a clone can run straight to the figure without the 6 GB scDRS tree. To point back at the analysis tree:
#   unset S6_OUT S6_V8 S6_CACHE S6_SCMS S6_SHEET and set DSCI_WORK_ROOT correctly (missing caches rebuild automatically).
#
# Width conservation (house-spec §17.9: to land exactly on 17.00 spread the gaps, never scale panels):
#   PW = 1.78 + 12 x CW (1.78 = gene-name column 1.54 + right pad 0.24)
#   2 x PW + 2 x 0.24 + LW = 17.00, LW = the legend script's own "width needed" estimate + 0.11
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export DSCI_WORK_ROOT=${DSCI_WORK_ROOT:-}
PY=${PY:-python3}
# Official escape hatch when PyMuPDF is missing from the conda base env: DSCI_PYLIBS points at a site-packages containing fitz.
# (This used to default to PYTHONPATH=${DSCI_WORK_ROOT:-}, a private path that cannot go into the repo.)
if [ -n "${DSCI_PYLIBS:-}" ]; then export PYTHONPATH=$DSCI_PYLIBS${PYTHONPATH:+:$PYTHONPATH}; fi
$PY -c "import fitz, matplotlib, pandas, numpy" 2>/dev/null || {
  echo "missing dependencies: $PY needs fitz(PyMuPDF)/matplotlib/pandas/numpy."
  echo "  either pip install pymupdf, or export DSCI_PYLIBS=<a directory containing these packages>"; exit 1; }
A=$REPO/assembly
export S6_OUT=${S6_OUT:-$REPO/figS6_out}
export S6_V8=${S6_V8:-$REPO/tables_v16}
export S6_CACHE=${S6_CACHE:-$S6_OUT/caches}
export S6_SCMS=${S6_SCMS:-$REPO/tables_v16}
export S6_SHEET=${S6_SHEET:-$S6_OUT}
mkdir -p "$S6_OUT" "$S6_CACHE" "$S6_SHEET"

echo "== 1/5 rows 1+2+4F: A locus tracks / B window permutation / C bin UMAP / F per-nucleus boxes (merged cell)"
$PY $A/figS6v16_panels.py

echo "== 2/5 row 3: measure the legend ink width -> solve CW -> D, E, vertical legend"
COMMON=(PANEL_H=5.75 PANEL_BODY_B=1.35 PANEL_ZMODE=1 PANEL_LEGEND=0 PANEL_TITLE=0
        PANEL_DCOLUMN=0 PANEL_CONTRAST=mol6)
run_common() { for kv in "${COMMON[@]}"; do export "${kv%%=*}=${kv#*=}"; done; }
MEAS=$($PY $A/figS6v18_legend_right.py | grep "legend needs width")
LW=$(echo "$MEAS" | grep -oE "[0-9]+\.[0-9]+" | head -1)
GEO=$($PY -c "LW=$LW+0.11; CW=(17.00-2*0.24-LW-3.56)/24.0; PW=1.78+12*CW; tot=2*PW+2*0.24+LW; assert abs(tot-17.00)<0.02, tot; print('%.4f %.2f %.2f'%(CW,PW,LW))")
read -r CW PW LW2 <<<"$GEO"
echo "   legend ink $LW -> column width LW=$LW2; solved CW=$CW => PW=$PW; row width 2 x PW + 2 x 0.24 + LW = 17.00 OK"
run_common
PANEL_W=$PW PANEL_CW=$CW PANEL_NAMES=0 PANEL_TAG=_s6v18D $PY $A/p0_8c_fig6c_top1000.py
PANEL_W=$PW PANEL_CW=$CW PANEL_PRESPEC=PHLDB1,PHGDH,BCAS1,PLLP PANEL_TAG=_s6v18E $PY $A/p0_8c_fig6c_landing.py
LEG_W=$LW2 LEG_H=5.75 $PY $A/figS6v18_legend_right.py

echo "== 3/5 row 4G: 74-trait dual heatmap (moved from Fig 6 D; narrowed only by env, zero vector scaling)"
COMPACT=1 SUF=_s6v16G WD=0.40 WK=0.40 WG=0.12 LM=1.55 GAPC=0.42 RMC=1.10 \
  $PY $A/make_figS6v16_G_heatmap.py

echo "== 4/5 compose (four rows; letters stamped by the compose layer; the self-calibrated left margin EDGE measured from the source pages' text boxes)"
$PY $A/figS6v16_compose.py

echo "== 5/5 gates: collision / two-size type + font whitelist / text on data"
S=$S6_SHEET/Figure_S6_v16.pdf
$PY $A/collide_gate.py $S
$PY $A/font_gate.py $S --strict | head -3
$PY $A/verify_text_on_data.py $S
echo "[done] $S"
