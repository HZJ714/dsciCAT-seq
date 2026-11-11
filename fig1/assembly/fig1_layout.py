#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig1_layout.py — Figure 1 版式几何的**唯一真源**（replot_v3.py 与 compose_fig1_v3.py 都从这里取数）
================================================================================
2026-10-20 按作者意见重排（改前 17.00 x 26.60 → 改后 17.00 x 25.49，省 1.11 cm）：
  1) A 与 B/C 的行距：A 流程图墨迹底 6.956 → B 字母墨迹顶，视觉空档 0.74 cm **压到 0.35 cm**；
  2) B 与 C **等高**（ROWH = 6.53 = C 在可用宽 11.85 cm 下按比例的最大高度；C 不能再高，
     再高就要超出可用宽 / 非等比拉伸，所以是 B 缩到 C，而不是 C 涨到 B）；
  3) D 的**上行 = 下行**（原 4.04 / 3.00 不等高 → 都取 H_ROW = 3.50，从 add_axes 真改箱体，
     不是只改外框）；
  4) E 的总高 = D 上行高（H_ROW），且 E 改为**自绘**（原来直接嵌 22x12 的 plots-cellline 原图，
     缩放 0.437，高度被原图白边钉死，无法与 D 对齐）；
  5) F 的方块边长 = H_ROW、标签带 = LBL，于是 F 的绘图框与 D 下行**逐边对齐**；
  6) D 的 HEK293T / NIH3T3 标记从上行左上角**移到下行左下角**（与 E 同一约定：
     transAxes (0.02, 0.02) va=bottom ha=left，加粗 8pt）。

A 的来源变更：作者批准的 A 带是**在 Illustrator 里重排过**的版本（`Workflow_combined.pdf` 是重排前的
旧布局：'Integrated matrix' 框内三行互相压线），唯一载体是 `ai/Fig1_v3.ai`（实测与交付件逐像素一致）。
所以 A 现在按**矢量裁贴**从 AI 母版取（CLIP_A，1:1），再经 `p1_pkg_fix.py` 摘掉带外幽灵文字，
成品 `replot/F1A_workflow.pdf`；面板字母 A 已在带内，compose 不再另画（与 C 同理）。

改尺寸只改这里，然后**必须**依次重跑：
  python3 scripts/extract_fig1_A.py            # 若交付件 A 带有变
  python3 scripts/replot_v3.py                 # 出 replot/{B,D,E,F}_v3_*.pdf
  python3 scripts/compose_fig1_v3.py           # 出 Fig1_v3.pdf
"""
CM = 72 / 2.54

PW = 17.00                                    # 页宽（本室 SI 用 18.00，Fig1 历史上是 17.00，不动）

# ---- A：流程图（矢量裁贴自 AI 母版 ai/Fig1_v3.ai，**含面板字母 A**）-------------------------
# 坐标在**交付件页坐标**下量得：字母 A 墨迹顶 0.344 / 流程图墨迹底 6.956（旧版槽位 0.42-7.42）
CLIP_A = (0.00, 0.20, 17.00, 7.05)            # 带 17.00 x 6.85 -> replot/F1A_workflow.pdf
A_H = round(CLIP_A[3] - CLIP_A[1], 3)         # 6.85
A_INK_TOP_OFF = round(0.344 - CLIP_A[1], 3)   # 0.144 带内：字母 A 墨迹顶
A_INK_BOT_OFF = round(6.956 - CLIP_A[1], 3)   # 6.756 带内：流程图墨迹底
A_TOP = 0.10                                  # 带上沿到页顶

# ---- B/C 行 -------------------------------------------------------------------------------
GAP_A_BC_VIS = 0.38                           # A 墨迹底 → B 字母墨迹顶 的视觉空档（原 0.74）
LTR_OFF = 0.29                                # 面板字母参数 y = 槽位顶 - 0.29（墨迹顶≈槽位顶-0.25）
# 字母**基线 = 面板顶 + 0.44**（墨迹顶 = 面板顶 + 0.04）——这是参考图里 C 字母的真实做法（交付版
# C 基线 8.53 vs 面板顶 8.05）；旧版 B/D/E/F 用 LTR_OFF 把字母顶到面板上方 0.32 cm，
# 于是 B 和 C 差 0.37 cm（作者指出）。改用 LTR_DY 统一跟 C 对齐。
LTR_DY = 0.10
# ---- 面板字母偏移（作者 2026-10-20 22:0x："和图靠太近了，模仿 A/C 离图有点距离；
#      B 和 F 左上一点，D E 就上一点"）----------------------------------------------------
# 元组 = (dx, dy)，作用于旧位置之上；y 是**字顶**（compose 里再加 ascent 得基线），
# 负 dy = 上抬。D 的 x 已在页边 0.00，不能再左，故只上抬；E/F 同列，一起左移保持对齐。
# F 的 y 偏移含一段补偿：v3.5 方块 2.90→2.78 使 Y_F 下移 0.16，字母若只跟槽位
# 就离 D 下行顶脊只剩 0.068 cm（实测），故 -0.22 再减 0.10 => 间隙回到 0.15 cm。
LTR = {'B': (-0.20, -0.16), 'D': (0.00, -0.16), 'E': (0.00, -0.16), 'F': (-0.20, -0.32)}
#   ⚠ E 不能左移：D 上行右边框在 x=8.42，字母会压上去（实测重叠）；F 可以左移，因为它抬到
#     y 18.79-19.18 正好在 D 上下两行之间的空档（D 下行从 19.30 才开始），左边无墨。
LTR_BOLD = False                              # 作者要求：BDEF 换常规体，与 A/C（ArialMT）一致
ROWH = 6.53                                   # B/C 行高 = C 按宽 11.85 等比得到的实际高
BW, BX = 3.82, 0.50                           # B 槽位宽 / x
CX, CW = 5.05, 11.85                          # C 槽位 x / 宽
# v4 (2026-10-28, Codex #237): **C 不再从参考图矢量裁贴**，改由 scripts/fig1c_v4.py 从原始数据
#   整块重画成 replot/C_v4_tracks.pdf（尺寸就是 CW x ROWH = 11.85 x 6.53，compose 时 1:1 放）。
#   CLIP_C / C_ASPECT 只保留给 v3 回退用；v4 的 C 字母也不再需要白覆盖（墨迹不在裁贴里了）。
CLIP_C = (4.38, 3.18, 16.10, 9.64)            # 参考图 `Figure 1.pdf` 里的 tracks 区：11.72 x 6.46
C_ASPECT = (CLIP_C[2] - CLIP_C[0]) / (CLIP_C[3] - CLIP_C[1])
Y_BC = round(A_TOP + A_INK_BOT_OFF + GAP_A_BC_VIS + LTR_OFF, 3)     # 7.53

# ---- B 图内部：两个**正方形**散点子图（行高变小时自动缩边，不再手调分数）-----------------
# 子图间距必须容纳**上行子图的两行 x 轴标签**（~0.51 cm）+ 下行子图右上角图例的余量，
# 0.55 会压字（实测），故取 0.72；上下边距相应收紧，方形子图 2.51 cm。
# 子图间距 0.85：72 会压字（实测），95 会把下行的两行 x 标签顶到图例行上；0.85 + 底边距 0.80
# 是同时满足「上下拉开」与「不撞图例行」的组合（方形子图边长仍 2.39）。
B_L, B_R, B_T, B_B, B_GAP = 1.00, 0.18, 0.10, 0.80, 0.85
B_SIDE = min(BW - B_L - B_R, (ROWH - B_T - B_B - B_GAP) / 2)


def b_adj(W=None, H=None):
    """返回 fig.subplots_adjust(**...) 的图形分数参数（B 的两个方形子图）。"""
    W = W or BW; H = H or ROWH
    side = min(W - B_L - B_R, (H - B_T - B_B - B_GAP) / 2)
    return dict(left=B_L / W, right=(B_L + side) / W, top=1 - B_T / H,
                bottom=B_B / H, hspace=B_GAP / side)


# ---- 技术图例行（D/E/F 的配色钥匙） ---------------------------------------------------------
LEG_OFF = 0.36                                # 图例基线在 B/C 行之下的距离
LEG_Y = round(Y_BC + ROWH + LEG_OFF, 3)       # 14.37

# ---- 下排：D 两行 + E/F 右列 ----------------------------------------------------------------
H_ROW = 3.50                                  # D 上行 = D 下行 = E 总高 = F 方块边长
GAP_D = 0.80                                  # D 两行之间 = E/F 之间
LBL = 2.20                                    # 旋转全称标签带（D 下行 / F 方块之下）——实测墨迹只到 2.10
DW = 8.45                                     # D 宽
EW = 8.35                                     # E/F 宽
EX = 8.55                                     # E/F x
Y_D = round(LEG_Y + 0.58, 3)                  # 15.00  D/E 槽位顶
# ---- F 内部纵向版式（v3.5：标题带 → 横向色条带 → 方块 → 旋转全称标签带）----
# 方块之上这一整带 = TITLE_F = 1.00 cm，内含：标题行 0.30 + (间隙0.12+渐变条0.18+间隙0.22) = 0.70
# ⚠ 纵向预算是量出来的：8pt 标题大写高 0.20、6.5pt 数字行高 0.23、旋转 6pt 全称最长 2.01+pad。
#   上一版把标题放在 H-0.25 处，实测 bbox 顶 -0.06 => 被页顶裁掉，且与 0.2/1.0 数字压在一起。
F_TITH, F_CBAND, F_SIDE, F_LBLH = 0.26, 0.70, 2.78, 2.32
F_SW_D, F_SW_W = 0.10, 0.16                   # 行色块：距方块右缘 / 宽
F_SWH, F_SWD2 = 0.14, 0.06                    # 列色块条：高 / 距方块下缘（作者：只有右边有读者对不上，下面也要）
#   F_LBLH 里含：0.06 间隙 + 0.14 色条 + 0.04 标签 pad + 2.02 旋转全称墨迹 + 0.06 底余量 = 2.32
#   => 标签墨底 = 页 25.06，与 D 的标签墨底 25.00 齐平（再挤就得把 F 方块缩到 2.6 以下）
F_BAR_H = 0.18                                # 渐变条自身高
F_BAR_GAP = 0.12                              # 方块顶 → 渐变条底
TITLE_F = round(F_TITH + F_CBAND, 3)          # 0.96  方块之上的带（标题 + 横向色条）
D_TOP = 0.16                                  # D 图顶留白（上行 y 轴「5」刻度字的上半截 + 外 whisker）
D_H = round(D_TOP + 2 * H_ROW + GAP_D + LBL, 3)   # 10.16
E_TOP, E_BOT = 0.00, 0.16                     # E 上/下边距：下行要留出 y 轴「0」刻度字的下半截
E_H = round(H_ROW + E_TOP + E_BOT, 3)         # 3.66  E 总高 = 绘图框(H_ROW) + 上下边距
F_H = round(F_TITH + F_CBAND + F_SIDE + F_LBLH, 3)   # 6.06
Y_E = Y_D                                     # 槽位顶 = D 槽位顶；E_TOP=0 => **绘图框顶**也与 D 上行对齐
Y_F = round((Y_D + 2 * H_ROW + GAP_D) - (TITLE_F + F_SIDE), 3)   # 18.90  使 F **方块底边** = D 下行底边 22.80
PH = round(Y_D + D_H - D_TOP + 0.19, 3)       # 25.19  页高（D 绘图框底 22.80 + 标签带 2.20 + 下边距 0.19）

# D 图内部四轴（cm，自图底）：**行间距 GAP_D 必须落在两行之间**，不是图顶（这是上一版贴一起的原因）
D_LX, D_RX, D_AW = 0.70, 4.735, 3.685
D_AXPOS = [(D_LX, LBL + H_ROW + GAP_D, D_AW, H_ROW), (D_RX, LBL + H_ROW + GAP_D, D_AW, H_ROW),   # 上行 UMIs
           (D_LX, LBL,                 D_AW, H_ROW), (D_RX, LBL,                 D_AW, H_ROW)]  # 下行 genes
# E 图内部两轴：绘图框高 = H_ROW（与 D 上行逐边对齐），上下边距另给 E_TOP/E_BOT
D_YLIM_TOP = 5.5                              # RNA lg(UMIs) 上行 ylim 顶
D_YLIM_BOT = 4.8                              # RNA lg(detected genes) 下行 ylim 顶
E_YLIM_TOP = 5.9                              # ATAC lg(fragments) 数据最大 5.586，留头
E_YLIM_BOT = -0.30                            # ⚠ 不能从 0 起：seaborn 默认 cut=2 会把 KDE 下尾外插到
                                              #   数据最小值以下，ylim 底设 0 会把最后一根（sciCAR-seq，
                                              #   数据最低 ≈0.4）的下尾**削平在底框上**（作者指出"超出边框"）。
                                              #   原 notebook 没写 set_ylim ⇒ 自动缩放本就留了 ≈-0.3 的底，
                                              #   实测旧交付版 0 刻度悬在框内 0.30 个单位，这里照抄。

# E/F 图内部几何（cm）：E 绘图框 = H_ROW（上边距 0、下边距 E_BOT），F 方块边长 = H_ROW
# v4 (2026-10-28, Codex #222/#237 同族问题): E 右半（NIH3T3）原来 labelleft=False + left=False
#   => 整个右半**一个 y 刻度都没有**，只能靠眼睛横向比左半。两半本来就共用 ylim，所以正确修法是
#   右半把刻度搬到自己的**右脊柱**上（左半不动），代价是右边距 0.05 -> 0.32，两个子图各窄 0.13 cm。
E_LBL_W, E_RGT, E_GAP = 0.96, 0.32, 0.00      # E 左边距(ylabel+刻度) / 右边距 / 两轴间隙（共用一条竖线）
F_X0 = [0.85, 4.55]                           # 两方块左缘（v3.5 整体右移 0.57；原 [0.28, 4.02]）
F_CB_X, F_CB_W = 7.58, 0.20                   # （旧：右侧共享竖色条，v3.5 已改到方块上方，保留仅供回退）
#   （F_TITH/F_CBAND/F_SIDE/F_LBLH/F_SW_*/F_BAR_H 与 TITLE_F/F_H/Y_F 一起定义在上面 D/E/F 段）

SLOT = {'B': (BW, ROWH), 'D': (DW, D_H), 'E': (EW, E_H), 'F': (EW, F_H)}
SLOTS = {'B': (BX, Y_BC, BW, ROWH),
         'D': (0.00, round(Y_D - D_TOP, 3), DW, D_H),     # 槽位顶上移 D_TOP，使**绘图框**顶 = Y_D
         'E': (EX, Y_E, EW, E_H),
         'F': (EX, Y_F, EW, F_H)}
C_SLOT = (CX, Y_BC, CW, ROWH)
# C 字母原本在参考图裁贴内部（实测字顶 7.572、x 5.076）：白覆盖后按 B 同一约定重画，
# 这样 B/C 既同行同基线、又一起离图。白覆盖区实测除字母外零墨迹（600 dpi 校验）。
LTR_C = (CX - 0.20, Y_BC + LTR_DY - 0.16)
LTR_C_COVER = (CX - 0.03, Y_BC - 0.03, CX + 0.37, Y_BC + 0.51)
LETTERS = {'B': (BX + LTR['B'][0], Y_BC + LTR_DY + LTR['B'][1]),
           'C': LTR_C,
           'D': (0.00 + LTR['D'][0], Y_D + LTR_DY + LTR['D'][1]),
           'E': (EX + LTR['E'][0], Y_E + LTR_DY + LTR['E'][1]),
           'F': (EX + LTR['F'][0], Y_F + LTR_DY + LTR['F'][1])}

if __name__ == '__main__':
    print('页 %.2f x %.2f cm' % (PW, PH))
    print('F 内部: 方块 %s 宽 %.2f -> 右缘 %.2f；色条 %.2f-%.2f (+标签~0.50) <= %.2f' % (
        F_X0, H_ROW, F_X0[1] + H_ROW, F_CB_X, F_CB_X + F_CB_W, EW))
    print('A  y %.2f-%.2f (墨迹 %.2f-%.2f)' % (A_TOP, A_TOP + A_H, A_TOP + A_INK_TOP_OFF,
          A_TOP + A_INK_BOT_OFF))
    print('B/C 行 y %.2f-%.2f  高 %.2f  视觉空档 %.2f' % (Y_BC, Y_BC + ROWH, ROWH,
          Y_BC - LTR_OFF - (A_TOP + A_INK_BOT_OFF)))
    print('图例基线 %.2f' % LEG_Y)
    print('D  槽位 y %.2f-%.2f  绘图框: 上行 %.2f-%.2f  下行 %.2f-%.2f  行间隙 %.2f' % (
        Y_D - D_TOP, Y_D - D_TOP + D_H, Y_D, Y_D + H_ROW,
        Y_D + H_ROW + GAP_D, Y_D + 2 * H_ROW + GAP_D, GAP_D))
    print('E  y %.2f-%.2f (绘图框 %.2f-%.2f = D 上行)   F y %.2f-%.2f  方块 %.2f-%.2f' % (
        Y_E, Y_E + E_H, Y_E + E_TOP, Y_E + E_TOP + H_ROW, Y_F, Y_F + F_H,
        Y_F + TITLE_F, Y_F + TITLE_F + H_ROW))
    print('C 等比校验: 宽 %.2f -> 高 %.4f (ROWH=%.2f)  缩放 %.4f' % (
        CW, CW / C_ASPECT, ROWH, CW / (CLIP_C[2] - CLIP_C[0])))
