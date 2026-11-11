#!/usr/bin/env bash
# figS6v18_row3.sh —— FigS6 第 3 行（D | E | 竖版图例）一键重渲
# ---------------------------------------------------------------------------
# 作者 2026-11-17：20.11 cm 太高，"DE 整体左移，图例塞最右边一列"。
# 省下的 1.20 cm = 横版图例带行高 1.05 + 贴附间距 0.15。
#
# 宽度不是拍的，是**解出来的**（与 Fig6 摊栏间距同一思路，SKILL §17.9）：
#   面板页宽 PW = 1.78 + 12·CW          （1.78 = 基因名列 1.54 + 右垫 0.24，CW = 亚型列宽）
#   行总宽 2·PW + 2·GAP + LW = 17.00    （GAP=0.24，LW = 竖版图例实测墨迹宽 + 0.11）
#   ⇒ CW = (17.00 − 2·GAP − LW − 3.56) / 24
# 先跑一遍图例脚本拿 LW，再反解 CW 渲 D/E，最后用 LW 重渲图例一次（画布=墨迹宽）。
# 字号一律不动（8.0/6.0 两档），行高/体表高度一律不动（RH 0.0182 保持）。
set -euo pipefail
export DSCI_WORK_ROOT=${DSCI_WORK_ROOT:-}
export PYTHONPATH=${PYTHONPATH:-${DSCI_WORK_ROOT:-}/pylibs}
PY=${PY:-${DSCI_PY:?先 export DSCI_PY=<分析环境 python（3.9 + scdrs 1.0.2）>}}
A=$DSCI_WORK_ROOT/GitHub/figS6/assembly

# 与上一轮完全相同的显示/统计口径，只改几何
COMMON=(PANEL_H=5.75 PANEL_BODY_B=1.35 PANEL_ZMODE=1 PANEL_LEGEND=0 PANEL_TITLE=0
        PANEL_DCOLUMN=0 PANEL_CONTRAST=mol6)
run_common() { for kv in "${COMMON[@]}"; do export "${kv%%=*}=${kv#*=}"; done; }

# ① 量图例宽度（画布给足 2.00，脚本自己报"实际需要宽"）
MEAS=$($PY $A/figS6v18_legend_right.py | grep "实际需要宽")
LW=$(echo "$MEAS" | grep -oE "[0-9]+\.[0-9]+" | head -1)
echo "图例墨迹宽 $LW cm"
GEO=$($PY -c "LW=$LW+0.11; CW=(17.00-2*0.24-LW-3.56)/24.0; PW=1.78+12*CW; tot=2*PW+2*0.24+LW; assert abs(tot-17.00)<0.02, tot; print('%.4f %.2f %.2f'%(CW,PW,LW))")
read -r CW PW LW2 <<<"$GEO"
echo "解得 CW=$CW → PW=$PW，图例列 $LW2 ⇒ 行宽校验通过（=17.00）"

# ② D：top-1000 可测 223 行，**不点名**（阴性面板，名字会引导读者找规律）
run_common
PANEL_W=$PW PANEL_CW=$CW PANEL_NAMES=0 PANEL_TAG=_s6v18D \
  $PY $A/p0_8c_fig6c_top1000.py
# ③ E：cis38，只标 S6A 已按顺式窗口独立解析的 4 个基因（命名依据在图上看不到 ⇒ 不循环）
PANEL_W=$PW PANEL_CW=$CW PANEL_PRESPEC=PHLDB1,PHGDH,BCAS1,PLLP PANEL_TAG=_s6v18E \
  $PY $A/p0_8c_fig6c_landing.py
# ④ 图例按最终列宽重渲（画布 = 墨迹宽，才能与 D/E 顶对齐共行）
LEG_W=$LW2 LEG_H=5.75 $PY $A/figS6v18_legend_right.py
echo "[完成] 第 3 行三格已重渲；接着跑 figS6v16_compose.py 拼页"
