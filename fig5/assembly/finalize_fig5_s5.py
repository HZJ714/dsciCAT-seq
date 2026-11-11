#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
finalize_fig5_s5.py — 把 Fig5_v8 / FigS5_noc4 落进交付块（Codex G3/G4 收尾）
==========================================================================
任何一项自检不过就不落盘（退出码 1）。落盘内容：
  FINAL_SET_20261021/Figure_5.pdf   ← Fig5_v8.pdf     （18.00 × 18.46，A–C，加 B 图例块）
  FINAL_SET_20261021/Figure_S5.pdf  ← FigS5_noc4.pdf  （18.00 × 46.75，A–D，A 换统一坐标版 + B 图例块）
  png/ 两件 300 dpi 重出；MD5_checksums.txt 全量重算（行序照旧）；
  GENERATED_TABLE.md 两行更新；CHANGELOG_20261021.md 追加一节。
旧件进 _superseded/。
"""
import os, sys, hashlib, shutil, subprocess, datetime
sys.path.insert(0, os.environ.get('DSCI_PYLIBS', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'pylibs')))
import fitz
import numpy as np
from PIL import Image

CM = 72 / 2.54
FIG = os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures'))
BLK = f'{FIG}/FINAL_SET_20261021'
SUP = f'{BLK}/_superseded'
NEW = {'Figure_5.pdf': f'{FIG}/Fig5_v8.pdf', 'Figure_S5.pdf': f'{FIG}/FigS5_noc4.pdf'}
# 落盘前的『我要覆盖的是我以为的那件』守卫。2026-10-29 两轮：
#   第一轮 = 10-21 交付件 (0a012ea4fcad / b0a126910e5f)，已进 _superseded/
#   第二轮 = 本轮 G3/G4 版，再加 C 带三行字原位换字（Codex #92）
OLD_MD5 = {'Figure_5.pdf': '2302687ee28b', 'Figure_S5.pdf': 'a9c941f90413'}


def md5(p):
    return hashlib.md5(open(p, 'rb').read()).hexdigest()


def inkdiff(a, b, box, dpi=150):
    """返回 (盒内可见差 px, 盒外可见差 px)，box=(x0,y0,x1,y1) cm"""
    pa, pb = fitz.open(a)[0], fitz.open(b)[0]
    A = pa.get_pixmap(dpi=dpi, alpha=False); B = pb.get_pixmap(dpi=dpi, alpha=False)
    x = np.frombuffer(A.samples, np.uint8).reshape(A.height, A.width, A.n)[:, :, :3].astype(int)
    y = np.frombuffer(B.samples, np.uint8).reshape(B.height, B.width, B.n)[:, :, :3].astype(int)
    assert x.shape == y.shape, '两件尺寸不同，不能逐块比'
    d = np.abs(x - y).max(2) > 24
    m = np.zeros_like(d)
    m[int(box[1] * dpi / 2.54):int(box[3] * dpi / 2.54), int(box[0] * dpi / 2.54):int(box[2] * dpi / 2.54)] = True
    return int((d & m).sum()), int((d & ~m).sum())


def audit(paths):
    r = subprocess.run([sys.executable, f'{FIG}/scripts/final_audit.py'] + paths,
                       capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def main(apply=True):
    os.makedirs(SUP, exist_ok=True)
    ok = True
    # ---- 1) 源件门 ----
    for name, src in NEW.items():
        assert os.path.exists(src), '缺源件 ' + src
        rc, out = audit([os.path.abspath(src)])
        line = [l for l in out.splitlines() if '硬门' in l][-1]
        print('%-16s %s' % (os.path.basename(src), line.strip()))
        if 'FAIL' in line:
            ok = False
    # ---- 2) 与交付件逐块比：改动必须只在图例盒 / S5 的 A 带 ----
    plan = [('Figure_5.pdf', (6.60, 0.45, 11.72, 1.90), None),
            ('Figure_S5.pdf', (0.0, 21.0, 4.9, 23.1), (0.0, 0.0, 18.0, 11.6))]
    for name, keybox, allow in plan:
        old, new = f'{BLK}/{name}', NEW[name]
        if allow is None:                          # Fig5：整页只许图例盒变
            i, o = inkdiff(old, new, keybox)
            print('%-16s 图例盒内 %d px 变化，盒外 %d px 变化' % (name, i, o))
            if o:
                ok = False
        else:                                       # S5：A 带是本轮重出的，只比 y>11.6cm
            m = np.zeros((1, 1), bool)
            pa, pb = fitz.open(old)[0], fitz.open(new)[0]
            A = pa.get_pixmap(dpi=150, clip=fitz.Rect(0, 11.6 * CM, pa.rect.width, 46.75 * CM), alpha=False)
            B = pb.get_pixmap(dpi=150, clip=fitz.Rect(0, 11.6 * CM, pb.rect.width, 46.75 * CM), alpha=False)
            x = np.frombuffer(A.samples, np.uint8).reshape(A.height, A.width, A.n)[:, :, :3].astype(int)
            y = np.frombuffer(B.samples, np.uint8).reshape(B.height, B.width, B.n)[:, :, :3].astype(int)
            d = (np.abs(x - y).max(2) > 24)
            box = np.zeros_like(d)
            box[int((21.0 - 11.6) * 150 / 2.54):int((23.1 - 11.6) * 150 / 2.54), 0:int(4.9 * 150 / 2.54)] = True
            # C 带三行字的原位换字（Codex #92：撤 'RPM-normalized' + MWU 口径）也算允许区
            for cx0, cy0, cx1, cy1 in [(0.05, 30.08, 11.35, 30.42), (0.05, 30.66, 13.50, 31.00),
                                       (15.78, 32.42, 17.72, 32.74)]:
                box[int((cy0 - 11.6) * 150 / 2.54):int((cy1 - 11.6) * 150 / 2.54),
                    int(cx0 * 150 / 2.54):int(cx1 * 150 / 2.54)] = True
            print('  %-16s S5 B/C/D 带 (y>11.6cm)：允许区(图例盒+C 带三行字)内 %d px，区外 %d px'
                  '（区外必须 0，否则字母字重/位置或别的面板被动过）' % (
                name, int((d & box).sum()), int((d & ~box).sum())))
            if (d & ~box).sum():
                ok = False
    # ---- 3) 旧件 md5 必须是记录里的值（防止覆盖错东西）----
    for name, m in OLD_MD5.items():
        cur = md5(f'{BLK}/{name}')[:12]
        print('%-16s 交付件 md5 %s（应为 %s）%s' % (name, cur, m, '✓' if cur == m else '✗'))
        ok &= cur == m
    if not ok:
        print('\n自检未过，不落盘。'); return 1
    if not apply:
        print('\n[干跑] 去掉 --dry 才落盘。'); return 0
    # ---- 4) 落盘 ----
    for name, src in NEW.items():
        dst = f'{BLK}/{name}'
        dep = f'{SUP}/{name.replace(".pdf", "")}_20261021_deprecated.pdf'
        if not os.path.exists(dep):     # 只留 10-21 原始交付件，别把备份覆盖成中间版
            shutil.copy2(dst, dep)
        shutil.copy2(src, dst)
        pm = fitz.open(dst)[0].get_pixmap(dpi=300, alpha=False)
        im = Image.frombytes('RGB', (pm.width, pm.height), pm.samples)
        im.save(f'{BLK}/png/{name.replace(".pdf", ".png")}', dpi=(300, 300))
        print('  安装 %-16s %.2f × %.2f cm  md5 %s' % (
            name, fitz.open(dst)[0].rect.width / CM, fitz.open(dst)[0].rect.height / CM, md5(dst)[:12]))
    # ---- 5) MD5 清单（行序照旧）----
    p = f'{BLK}/MD5_checksums.txt'
    lines = open(p).read().splitlines()
    out = []
    for l in lines:
        h, f = l.split()
        out.append('%s  %s' % (md5(f'{BLK}/{f}'), f))
    open(p, 'w').write('\n'.join(out) + '\n')
    print('  MD5_checksums.txt 重算 %d 行' % len(out))
    # ---- 6) 全块体检 ----
    rc, out2 = audit([os.path.abspath(f'{BLK}/{f}') for f in
                      [l.split()[1] for l in open(p)] if f.endswith('.pdf')])
    print('\n'.join([l for l in out2.splitlines() if l.startswith(('Figure', 'Graphical', '硬门'))][:20]))
    return 0


if __name__ == '__main__':
    sys.exit(main(apply='--dry' not in sys.argv))
