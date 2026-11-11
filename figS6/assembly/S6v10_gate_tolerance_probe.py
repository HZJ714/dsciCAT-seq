# -*- coding: utf-8 -*-
"""反向溯源门的**容差体检**（只读探针，不改任何交付物）。

背景：verify_S6_v10.py 的 sourced() 判定是
        abs(round(a, dec) - f) < 1e-9   （a=台账里的真实值，f=图面 token）
当图面写的是极小数（P = 1.7e-33 这类），round(a, dec) 与 f 都≈0，
于是**任何** 绝对值 < 1e-9 的 token 都会被判"有出处" —— 反向门对小数位 P 失效。

本探针把四种编造数字分别写成单页 PDF，交给真正的门（verify_S6_v10.py <pdf>）跑，
打印每种写法是"放过"还是"拦下"。严格版门上线后再跑一遍，应当：
  * 8 个真图面仍然 8/8 OK（不误伤）
  * 4 个编造 token 全部被拦（不再靠 1e-9 蒙混）
"""
import os
import subprocess
import sys

import fitz

OUT = "/tmp/s6v10_gate_probe"
os.makedirs(OUT, exist_ok=True)

# 作者 2026-11-05：探针字体改用内置 helv。原先用 ARIAL.TTF 子集，fitz 把指数里的 ASCII 连字符
# 提取成软连字符 U+00AD（"9.9e\xad77"），token 被拆成 "9.9"/"77" 两个真台账值 → 假阴性“auto-pass”。
# 严格门下拦截 token 可能报 NO-SOURCE 或 UNREGISTERED，两者都算拦下。

# (标签, 图面文字, 期望：strict 应拦/放)
CASES = [
    ("fake_tiny_P",      "P = 9.9e-77",      "block"),   # 完全编造的极小 P
    ("fake_tiny_P2",     "P = 4.0e-10",      "block"),   # 真值 4e-09 附近的编造
    ("near_miss_rho",    "rho = 0.817",      "block"),   # 真值 0.819 的"手滑版"
    ("random_rho",       "rho = 0.5123",     "block"),   # 无中生有
    ("real_P",           "P = 1.7e-33",      "pass"),    # 台账里真有（G rna dsciCAT-seq）
    ("real_P_small",     "P = 5e-26",        "pass"),    # 台账里真有（A PHLDB1）
    ("real_rho",         "rho = 0.819",      "pass"),    # 台账里真有
]

gate = os.path.join(os.path.dirname(os.path.abspath(__file__)), "verify_S6_v10.py")
print("%-16s %-16s %-8s %s" % ("case", "face text", "want", "gate verdict"))
for tag, line, want in CASES:
    pdf = "%s/%s.pdf" % (OUT, tag)
    d = fitz.open()
    pg = d.new_page(width=120, height=60)
    pg.insert_text((6, 30), line, fontname="helv", fontsize=6)
    d.save(pdf)
    d.close()
    r = subprocess.run([sys.executable, gate, pdf], capture_output=True, text=True)
    caught = ("NO-SOURCE" in r.stdout) or ("UNREGISTERED" in r.stdout)
    verdict = "BLOCKED" if caught else "auto-pass"
    ok = (verdict == "BLOCKED") == (want == "block")
    print("%-16s %-16s %-8s %-10s %s" % (tag, line, want, verdict, "" if ok else "  <== 与期望不符"))
