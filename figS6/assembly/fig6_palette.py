# -*- coding: utf-8 -*-
"""The single colour source for the three named Fig 6 traits (MS / SCZ / HEIGHT).

Why a dedicated file: C and D label THE SAME three traits, so the colours must match byte for byte.
Previously each side re-declared the constants; editing one and forgetting the other produced "looks aligned, actually not".
Now both scripts source from here, and the acceptance script A6_fig6_v21_verify.py reads the colours back from the RENDERED PDF
and compares C's named-trait colours against D's (comparing the printed ink, not the source).
"""
import matplotlib.colors as mpl

# identical to the constants in make_fig6_v17_A_B_C_panels.py (these three values since v14)
RED = "#D51F26"      # MS      - farthest from the diagonal (the chromatin channel ranks MS lower)
BLU = "#272E6A"      # SCZ     - negative in both channels, but more so in chromatin
MUT = "0.42"         # HEIGHT  - consistent across channels (on the identity line)

TRAIT_COLOR = {"MS": RED, "SCZ": BLU, "HEIGHT": MUT}


def rgb(name):
    """'MS' -> (r,g,b) 0..255 ints, for comparing PDF span colours."""
    r, g, b = mpl.colors.to_rgb(TRAIT_COLOR[name])
    return tuple(int(round(v * 255)) for v in (r, g, b))
