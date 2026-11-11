#!/usr/bin/env bash
# figS5/data/link_data.sh — 把分析盘上的输入对象挂到本目录
# ==========================================================
# 仓库里不放这些大对象（见 ../../.gitignore）；Figure S5 的 notebook 只按相对路径
# data/<name> 读取，clone 之后先跑一次本脚本即可。真身位置全部由环境变量解析。
#
#   DSCI_WORK_ROOT   分析盘 online 根（必填）
#   DSCI_FY_ROOT       默认 $DSCI_WORK_ROOT/FY
#   DSCI_SCENIC_ROOT   默认 $DSCI_WORK_ROOT/pi/plots-scenic+
#   DSCI_FIG_ROOT      默认 $DSCI_WORK_ROOT/figures
#
# 用法: bash figS5/data/link_data.sh [--force]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ONLINE="${DSCI_WORK_ROOT:-}"
FY="${DSCI_FY_ROOT:-${ONLINE%/}/FY}"
SCENIC="${DSCI_SCENIC_ROOT:-${ONLINE%/}/pi/plots-scenic+}"
FIGS="${DSCI_FIG_ROOT:-${ONLINE%/}/figures}"
REPO="$(cd "$HERE/../.." && pwd)"

if [[ -z "$ONLINE" ]]; then
  echo "先 export DSCI_WORK_ROOT=<你的分析盘 online 根>（README《Environment variables》）" >&2
  exit 1
fi

NAMES=()
SRCS=()
add() { NAMES+=("$1"); SRCS+=("$2"); }

add "scplus_obj_final.pkl"                 "${FY}/oligo/3kdsci/GRN/scplus/scplus_obj_final.pkl"
add "eregulon_stability_table.csv"         "${SCENIC}/stability/eregulon_stability_table.csv"
add "nrxn3_mef2c_enhancers.fa"             "${SCENIC}/nrxn3_mef2c_enhancers.fa"
add "nrxn3_tracks.npz"                     "${FIGS}/scripts/cache/nrxn3_tracks.npz"
add "nrxn3_expr.npz"                       "${FIGS}/scripts/cache/nrxn3_expr.npz"
add "NIHMS70326-supplement-Table_S1.xlsx"  "$REPO/fig5/data/NIHMS70326-supplement-Table_S1.xlsx"

rc=0
for i in "${!NAMES[@]}"; do
  name="${NAMES[$i]}"; src="${SRCS[$i]}"
  if [[ ! -e "$src" ]]; then
    echo "MISS  $name  <- 找不到 $src"; rc=1; continue
  fi
  if [[ -L "$HERE/$name" && "${1:-}" != "--force" ]]; then
    echo "OK    $name  (已链接)"; continue
  fi
  ln -sfn "$src" "$HERE/$name"
  echo "LINK  $name -> $src"
done
exit $rc
