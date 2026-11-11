# -*- coding: utf-8 -*-
"""Fig 6 · D rank colour map: quantifying "can you tell the levels apart" (author 11-03: with the digits gone, the rank colours need more contrast)

Method: sample each candidate ramp into 74 levels (= the 74 traits, matching the drawing script: rank 1 -> 1.0, rank 74 -> 0.0,
so both ends are taken full, no 0.02/0.98 trim), then compute the Lab ΔE*76 of ADJACENT levels. JND ~ 2.3 is only an order-of-magnitude reference for resolvability;
the real criterion is the median adjacent ΔE and the fraction of adjacent pairs above JND. The old greyblue is 1.09 / 0%,
so "they all look the same" is decided by the colour map, not by taste.

Usage: python3 A6_fig6_D_rank_contrast.py [hex1,hex2,...]     # append a custom ramp
Output: a stdout table sorted by median ΔE; `--pick NAME` prints that ramp's 74 level colours (for colour-bar verification)
"""
import sys
import numpy as np
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

N = 74
JND = 2.3
M = np.array([[0.4124564, 0.3575761, 0.1804375],
              [0.2126729, 0.7151522, 0.0721750],
              [0.0193339, 0.1191920, 0.9503041]])


def lab(rgb):
    c = np.atleast_2d(np.asarray(rgb, float))
    c = np.where(c > 0.04045, ((c + 0.055) / 1.055) ** 2.4, c / 12.72)
    X, Y, Z = (c @ M.T).T
    f = lambda t: np.where(t > 0.00885675, np.cbrt(t), 7.787 * t + 16 / 116)
    fx, fy, fz = f(X / 0.95047), f(Y), f(Z / 1.08883)
    return np.stack([116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz)], 1)


def dE(a, b):
    return np.linalg.norm(lab(a) - lab(b), axis=-1)


def stat(cm):
    q = cm(np.linspace(0, 1, N))[:, :3]          # same sampling as the drawing script: ends taken full
    de = np.array([dE(q[i], q[i + 1]) for i in range(N - 1)])
    return dict(med=np.median(de), lo=de.min(), hi=de.max(),
                above=100.0 * (de > JND).mean(), total=float(dE(q[0], q[-1])), q=q)


CAND = {"greyblue": LinearSegmentedColormap.from_list("x", ["#F4F6F8", "#8FA0B8", "#1B2A45"]),  # the old delivery
        "ink": LinearSegmentedColormap.from_list("x", ["#FFFFFF", "#000000"])}
for k in ["cividis_r", "viridis_r", "magma_r", "inferno_r", "plasma_r", "YlGnBu_r", "PuBu_r", "Greys_r",
          "Blues_r", "rocket_r", "crest_r"]:
    if k in mpl.colormaps:
        CAND[k] = mpl.colormaps[k]
for a in sys.argv[1:]:
    if a.startswith("--"):
        continue
    CAND["custom:" + a] = LinearSegmentedColormap.from_list("x", a.split(","))

if "--pick" in sys.argv:
    name = sys.argv[sys.argv.index("--pick") + 1]
    q = stat(CAND[name])["q"]
    print("rank 1  -> %s   ...   rank 74 -> %s" % (np.round(q[0] * 255).astype(int), np.round(q[-1] * 255).astype(int)))
    print("74-level RGB:")
    for i in range(0, N, 6):
        print("  #%02d " % (i + 1) + "  ".join("%3d %3d %3d" % tuple(v) for v in q[i:i + 6]))
    sys.exit(0)

print("%-12s %8s %8s %8s %10s %8s" % ("ramp", "medianΔE", "minΔE", "maxΔE", ">JND%", "totalΔE"))
for k, cm in sorted(CAND.items(), key=lambda x: -stat(x[1])["med"]):
    s = stat(cm)
    print("%-12s %8.2f %8.2f %8.2f %9.0f%% %8.1f" % (k, s["med"], s["lo"], s["hi"], s["above"], s["total"]))
print("\nJND~%.1f; 74 levels = 74 traits; sampling takes both ends full (rank1=1.0 darkest, rank74=0.0 lightest)" % JND)
