#!/usr/bin/env bash
# fig5/data/link_data.sh — 把分析盘上的输入对象挂到本目录
# ==========================================================
# 仓库里不放这些大二进制（见 ../../.gitignore）；Figure 5 的 notebook 只按相对路径
# data/<name> 读取，所以 clone 之后先跑一次本脚本。真身位置全部由环境变量解析，
# 脚本本身不含任何机器绝对路径（对照表见 fig5/data/README.md）。
#
#   DSCI_WORK_ROOT   分析盘 online 根（必填）
#   DSCI_FY_ROOT       默认 $DSCI_WORK_ROOT/FY
#
# 用法: bash fig5/data/link_data.sh [--force]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ONLINE="${DSCI_WORK_ROOT:-}"
FY="${DSCI_FY_ROOT:-${ONLINE%/}/FY}"

if [[ -z "$ONLINE" ]]; then
  echo "先 export DSCI_WORK_ROOT=<你的分析盘 online 根>（README《Environment variables》）" >&2
  exit 1
fi

NAMES=()
SRCS=()
add() { NAMES+=("$1"); SRCS+=("$2"); }

add "NIHMS70326-supplement-Table_S1.xlsx" "${ONLINE%/}/NIHMS70326-supplement-Table_S1.xlsx"
add "scplus_obj_final.pkl"                "${FY}/oligo/3kdsci/GRN/scplus/scplus_obj_final.pkl"
add "10x_subset_harmony.RDS"              "${FY}/tenx/10x_subset_harmony.RDS"
add "HT_harmony.RDS"                      "${FY}/HT/HT_harmony.RDS"

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
