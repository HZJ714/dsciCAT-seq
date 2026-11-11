#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The manuscript-wide font policy (hard convention since 2026-10-13)

Convention: EVERY re-run/new figure uses the Arial family. This machine has no Arial; the user provides
$DSCI_WORK_ROOT/fonts/LiberationSans-*.ttf (Ascender, OFL 1.1) - it is METRICALLY IDENTICAL to Arial
(unitsPerEm 2048 / ascender 1854 / descender -434, byte-identical to Arial in all three),
so widths, wrapping and spacing all compute as Arial; at submission, Illustrator's Find Font
"LiberationSans" -> "Arial" is a ZERO-DISPLACEMENT swap.

Why in code rather than installed into the system: /opt/conda pip installs do not persist (fitz has been lost twice),
and $HOME's matplotlib font cache can be cleared. Importing this module in any script takes effect at once, with no environment dependency.

Known pitfall: neither Liberation Sans nor Arial HAS U+207B superscript minus ⁻ ("5 × 10⁻⁸" would lose the glyph),
the other 89 glyphs used in the figures (≥ ± · – … → − ₀ ₁ ⁸ etc.) all exist.
**mathtext is not the answer**: measured - with `mathtext.fontset` = "custom" or "dejavusans",
`$5\times10^{-8}$` hands the whole math expression to DejaVu Sans, so the PDF embeds two fonts.
The right way is `sup_text()`: base and superscript become two text objects; the superscript uses the body font's own U+2212 "−", raised and set small,
measured: the PDF then embeds ONLY the Arial family (Arial itself when Liberation is absent). Never type the Unicode "⁻" directly in code again.
"""
import os
import matplotlib

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

FONT_DIR = WORK + "/fonts"
CM = 1 / 2.54

# Prefer the true Arial: Monotype v7.00; the user uploaded ARIAL/ARIALBD/ARIALI/ARIALBI.TTF from this machine on 2026-10-13.
# Liberation Sans is only a fallback - it is not Arial, but upm/ascender/descender are identical
# (2048 / 1854 / -434); with identical metrics a font swap cannot move a laid-out page.
_CAND = {
    "normal":     ["ARIAL.TTF",   "LiberationSans-Regular.ttf"],
    "bold":       ["ARIALBD.TTF", "LiberationSans-Bold.ttf"],
    "italic":     ["ARIALI.TTF",  "LiberationSans-Italic.ttf"],
    "bolditalic": ["ARIALBI.TTF", "LiberationSans-BoldItalic.ttf"],
}
FONTS = {}
for _k, _names in _CAND.items():
    FONTS[_k] = next((os.path.join(FONT_DIR, n) for n in _names
                      if os.path.exists(os.path.join(FONT_DIR, n))), None)
    if FONTS[_k] is None:
        raise FileNotFoundError("missing font %s (looked at %s)" % (_k, _names))
FAMILY = "Arial" if os.path.basename(FONTS["normal"]).upper() == "ARIAL.TTF" else "Liberation Sans"

# ---- FONTS_PDF: renamed copies for fitz insert_text (2026-10-15; their final review hard-asserts §6.3) ----
# When fitz embeds a font it writes the INTERNAL name (family(+subfamily)) into the PDF's /BaseFont, which journal preflight reads.
# True Arial's family + " Bold" would be written as "Arial Bold", outside their whitelist {ArialMT, Arial-BoldMT, ...}.
# The copy approach: nameID1<-PostScript name, nameID2<-empty, nameID4<-PostScript name, so fitz emits BaseFont=Arial-BoldMT.
# Glyphs/metrics untouched (name table only); matplotlib's use_arial still uses the original FONTS, so family matching is unaffected.
FONTS_PDF = {}
_PSNAME = {"normal": "ArialMT", "bold": "Arial-BoldMT",
           "italic": "Arial-ItalicMT", "bolditalic": "Arial-BoldItalicMT"}
_NAMED_DIR = WORK + "/fonts_named"


def _make_named_copy(src, ps, dst):
    from fontTools.ttLib import TTFont
    f = TTFont(src)
    for rec in f["name"].names:
        if rec.nameID in (1, 4):
            rec.string = ps.encode("utf-16-be") if rec.platformID == 3 else ps.encode("latin-1")
        elif rec.nameID == 2:
            rec.string = ("" if rec.platformID == 3 else b"").encode("utf-16-be") if rec.platformID == 3 else b""
    f.save(dst)


if FAMILY == "Arial":
    os.makedirs(_NAMED_DIR, exist_ok=True)
    for _k, _ps in _PSNAME.items():
        _dst = os.path.join(_NAMED_DIR, "%s.ttf" % _ps)
        if not os.path.exists(_dst):
            _make_named_copy(FONTS[_k], _ps, _dst)
        FONTS_PDF[_k] = _dst
else:  # Liberation fallback: not forced (only used by non-delivered figures); keep the original names
    FONTS_PDF = dict(FONTS)

# Journal norm: body annotations 6-8 pt; the floor is 5 pt (dense spots like heatmap row names)
BASE = {
    "font.family": FAMILY,
    "mathtext.fontset": "custom",   # see the file header: keep mathtext on the body font, never DejaVu
    "pdf.fonttype": 42, "ps.fonttype": 42,        # live type, editable in AI
    "svg.fonttype": "none",                        # SVG keeps live type too
    "font.size": 7.0,
    "axes.titlesize": 7.0, "axes.labelsize": 7.0,
    "xtick.labelsize": 6.0, "ytick.labelsize": 6.0,
    "legend.fontsize": 6.0, "axes.titleweight": "normal",
    "axes.linewidth": 0.6, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
    "xtick.major.size": 2.0, "ytick.major.size": 2.0,
    "xtick.minor.size": 1.0, "ytick.minor.size": 1.0,
    "xtick.major.pad": 1.6, "ytick.major.pad": 1.6,
    "lines.linewidth": 0.8, "patch.linewidth": 0.5,
    "savefig.bbox": None, "savefig.dpi": 600, "figure.dpi": 120,
    "axes.spines.top": True, "axes.spines.right": True,
}


def use_arial(extra=None):
    """Registers the fonts and applies the policy. Must be called before any Figure is created."""
    for f in sorted(os.listdir(FONT_DIR)):
        if f.lower().endswith((".ttf", ".otf")):
            try:
                fm.fontManager.addfont(os.path.join(FONT_DIR, f))
            except Exception as e:
                print("arial_style: registration failed %s: %s" % (f, e))
    plt.rcParams.update(BASE)
    if extra:
        plt.rcParams.update(extra)
    for p in FONTS.values():                       # explicit registration, not directory-scan order
        fm.fontManager.addfont(p)
    got = {f.name for f in fm.fontManager.ttflist if f.fname.startswith(FONT_DIR)}
    assert FAMILY in got, "%s not registered: %s" % (FAMILY, got)
    return plt


def cm2fig(w_cm, h_cm):
    return plt.figure(figsize=(w_cm * CM, h_cm * CM))


def align_axis_labels(ax, pad_x=None, pad_y=None):
    """For axis-title alignment: keeps a fixed gap between x/y titles and tick text, reproducible across panels."""
    if pad_x is not None:
        ax.xaxis.label.set_label_coords(0.5, -pad_x)
    if pad_y is not None:
        ax.yaxis.label.set_label_coords(-pad_y, 0.5)


def sup_text(ax, x, y, base, sup, fontsize=None, transform=None, **kw):
    """Single-font true superscript (Liberation Sans / Arial lack U+207B ⁻, and mathtext falls back to DejaVu).

    All measurement happens in display (= point) space, then transforms back to the target coordinate system AS A WHOLE.
    Two pits we hit: feeding data-space y into the inverse transform, or adding fontsize directly to data-space y -
    both push the superscript off-canvas and it vanishes SILENTLY - the figure looks fine but is missing a glyph, more dangerous than an error.
    """
    fs = fontsize or plt.rcParams["font.size"]
    tr = transform or ax.transData
    t1 = ax.text(x, y, base, fontsize=fs, va="baseline", transform=tr, **kw)
    fig = ax.figure
    fig.canvas.draw()
    bb = t1.get_window_extent(renderer=fig.canvas.get_renderer())
    xp, yp = tr.inverted().transform((bb.x1, bb.y0 + 0.74 * (bb.y1 - bb.y0)))
    t2 = ax.text(xp, yp, sup, fontsize=fs * 0.70, va="baseline", transform=tr, **kw)
    return t1, t2


def ten_pow_sup(ax, x, y, coef="5", exp="8", prefix="p = ", fontsize=7, **kw):
    """Draws 'p = 5 × 10⁻⁸' (single font). kw may carry transform / color / ha etc."""
    return sup_text(ax, x, y, "%s%s \u00d7 10" % (prefix, coef), "\u2212%s" % exp,
                    fontsize=fontsize, **kw)


def open_spines(ax, top=True, right=True):
    """Statistical panels share an open frame (top/right spines off); spatial figures (UMAP/manhattan colour frames) keep all four spines."""
    if top:
        ax.spines["top"].set_visible(False)
    if right:
        ax.spines["right"].set_visible(False)


def frame(ax, lw=0.6):
    """Panels needing all four spines (UMAP) declare it explicitly, so "some framed, some not" can never be an accident."""
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(lw)



def sub_text(ax, x, y, base, sub, tail="", fontsize=None, rotation=0, transform=None,
             center=False, **kw):
    """True subscript, horizontal and 90-degree vertical (manhattan's y-axis title -log10(P) is vertical).

    Why not the ready-made "10" subscript characters: Arial LACKS U+2080/U+2081 (Liberation Sans has them),
    and they would be DROPPED SILENTLY - no error, no DejaVu, invisible in the font table, and the axis title simply loses its base digit.
    mathtext falls back as a whole too. So it is typeset in three pieces: base + a small sunk sub + tail, all Arial.

    All geometry is computed in display(point) space, then inverse-transformed to the target coordinate system; vertical text uses rotation_mode="anchor"
    so ha/va act on the unrotated glyph and the baseline position is predictable (reading direction = display +y, below baseline = +x).
    The translation must be written `tr + Affine2D().translate(...)`: reversed, the shift is treated as axes-fraction units,
    scaled hundreds of times, and the whole string flies off the canvas (hit).
    """
    assert rotation in (0, 90), "this function only handles 0/90 degrees"
    fs = fontsize or plt.rcParams["font.size"]
    tr = transform or ax.transData
    fig = ax.figure
    ex = dict(kw); ex.pop("ha", None); ex.pop("va", None); ex.pop("rotation_mode", None)
    gap, drop = 0.06 * fs, 0.22 * fs
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    A = tuple(tr.transform((x, y)))

    def place(pt, txt, size):
        t = ax.text(*tr.inverted().transform(pt), txt, fontsize=size, rotation=rotation,
                    rotation_mode="anchor" if rotation else None,
                    ha="left", va="baseline", transform=tr, **ex)
        fig.canvas.draw()
        return t, t.get_window_extent(renderer=r)

    t1, b1 = place(A, base, fs)
    if rotation == 0:
        p_sub = (b1.x1 + gap, A[1] - drop)
    else:                                   # vertical: advance along +y, the subscript sinks toward +x
        p_sub = (A[0] + drop, A[1] + b1.height + gap)
    t2, b2 = place(p_sub, sub, fs * 0.70)
    t3 = None
    if tail:
        p3 = (b2.x1 + 0.5 * gap, A[1]) if rotation == 0 else (A[0], p_sub[1] + b2.height + gap)
        t3, _ = place(p3, tail, fs)
    if center:                              # centre on the anchor along the reading direction (vertical axis titles need centring)
        from matplotlib.transforms import Bbox, Affine2D
        ts = [t for t in (t1, t2, t3) if t is not None]
        u = Bbox.union([t.get_window_extent(renderer=r) for t in ts])
        off = (A[0] - (u.x0 + u.x1) / 2.0, A[1] - (u.y0 + u.y1) / 2.0)
        for t in ts:
            t.set_transform(tr + Affine2D().translate(*off))
            fig.canvas.draw()
    return t1, t2, t3



def vsub_label(ax, y, base, sub, tail="", fontsize=7, pad_pt=6.0, drop_pt=None,
               dx_pt=None, sub_scale=0.70, **kw):
    """True subscript for vertical axis titles: base + a small sunk sub + tail (for -log10(P)).

    Why not the ready-made subscript characters: Arial LACKS U+2080/U+2081 (Liberation Sans has them),
    writing "-log\u2081\u2080(P)" directly gets DROPPED SILENTLY - no error, no DejaVu, nothing odd in the font table,
    and the axis title loses its base digit (hit, measured). mathtext falls back as a whole too.

    Placement conventions:
      - horizontal: derive the axes fraction from pad_pt (how many pt from the PAGE left), never hard-code a negative fraction -
        manhattan's left=0.058; a hard-coded -0.105 would fly the whole string off-canvas without an error (hit).
      - all three pieces ha=va="center"; lengths measured in pixels and divided by the axes pixel size to get fractions;
        horizontal offsets divide by width, vertical by height - mixing them once skews the whole layout (hit).
      - no whole-bbox translation (the union measured by the center=True variant was unreliable).
    """
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bb = ax.get_window_extent(renderer=r)
    fs = fontsize
    gap = r.points_to_pixels(0.10 * fs)
    drop = r.points_to_pixels(0.30 * fs) if drop_pt is None else r.points_to_pixels(drop_pt)
    if dx_pt is not None:          # offset relative to the AXIS left edge (when a multi-axis ylabel is not at the page edge)
        x = r.points_to_pixels(dx_pt) / bb.width
    else:
        x = (r.points_to_pixels(pad_pt) - bb.x0) / bb.width

    def adv(txt, size):
        t = ax.text(0, 0, txt, fontsize=size, rotation=90, ha="center", va="center",
                    transform=ax.transAxes, **kw)
        fig.canvas.draw()
        h = t.get_window_extent(renderer=r).height      # vertical: text length = bbox height (px)
        t.remove()
        return h

    a1, a2 = adv(base, fs), adv(sub, fs * sub_scale)
    a3 = adv(tail, fs) if tail else 0.0
    total = a1 + gap + a2 + (gap + a3 if tail else 0.0)
    yc = y - total / (2.0 * bb.height) + a1 / (2.0 * bb.height)
    out = [ax.text(x, yc, base, fontsize=fs, rotation=90, ha="center", va="center",
                   transform=ax.transAxes, **kw)]
    yc += (a1 + gap) / bb.height
    out.append(ax.text(x + drop / bb.width, yc, sub, fontsize=fs * sub_scale, rotation=90,
                       ha="center", va="center", transform=ax.transAxes, **kw))
    if tail:
        yc += (a2 + gap) / bb.height
        out.append(ax.text(x, yc, tail, fontsize=fs, rotation=90, ha="center", va="center",
                           transform=ax.transAxes, **kw))
    fig.canvas.draw()
    return out


def align_xlabel(ax, fig, y_fig=0.040):
    """Pins x-axis titles at THE SAME HEIGHT across the whole page (cross-panel alignment).

    The default labelpad is relative to the axis, so panels WITH tick text push their titles down:
    Fig6C's UMAP has no tick digits and Fig6D has them, differing by ~0.3 cm. Switched to
    a blended transform (x in axes fraction, y in figure fraction), aligning same-row panels on one baseline.
    set_label_coords also disables _autolabelpos (else every draw recomputes the position).
    """
    from matplotlib.transforms import blended_transform_factory as btf
    # note: the transform must be passed as an argument to set_label_coords.
    # calling ax.xaxis.label.set_transform(...) first and then set_label_coords(x, y) does nothing -
    # in the source, transform=None makes it RESET to transAxes, so y reverts to axes fraction,
    # and the cross-panel "one baseline" is gone (hit; D/E measured 4.2 pt apart).
    ax.xaxis.set_label_coords(0.5, y_fig,
                              transform=btf(ax.transAxes, fig.transFigure))


def patch_fitz_arial(fitz_mod):
    """Makes fitz's insert_text / get_text_length substitute the base-14 "helv"/"hebo" with EMBEDDED Arial.

    Why: PyMuPDF's fontname="helv"/"hebo" is a PDF standard-14 font and is NOT EMBEDDED in the PDF.
    Without true Helvetica on this machine, viewers substitute another font (measured: FigS4 has 466 glyphs on this path,
    ~17% of that figure), and journals flag "font not embedded" at submission.
    A monkeypatch rather than editing call sites: S4's scripts have 10+ insert_text calls; fixing them one by one invites misses,
    and get_text_length must change in sync, or measured widths disagree with the actual typesetting and text overlaps.
    """
    _ins = fitz_mod.Page.insert_text

    def _resolve(fontname, fontfile):
        if fontfile is None and fontname in ("helv", "hebo", "Helvetica"):
            bold = fontname == "hebo"
            return ("Arial-Bold" if bold else "Arial-Regular"), (FONTS_PDF["bold"] if bold else FONTS_PDF["normal"])
        return fontname, fontfile

    def _insert(self, point, text, fontname="helv", fontfile=None, **kw):
        fontname, fontfile = _resolve(fontname, fontfile)
        return _ins(self, point, text, fontname=fontname, fontfile=fontfile, **kw)

    fitz_mod.Page.insert_text = _insert

    # fitz.get_text_length() only knows base-14 and has NO fontfile argument, so measuring Arial widths
    # requires fitz.Font(fontfile=...).text_length(). helv measurements cannot stay:
    # Helvetica and Arial differ in width; under-measuring makes text collide.
    _cache = {}

    def _font(bold):
        key = "b" if bold else "r"
        if key not in _cache:
            _cache[key] = fitz_mod.Font(fontfile=FONTS_PDF["bold"] if bold else FONTS_PDF["normal"])
        return _cache[key]

    def _gtl_patched(text, fontname="helv", fontsize=11, **kw):
        if fontname in ("helv", "hebo", "Helvetica"):
            return _font(fontname == "hebo").text_length(text, fontsize=fontsize)
        ff = kw.pop("fontfile", None)
        if ff:
            return fitz_mod.Font(fontfile=ff).text_length(text, fontsize=fontsize)
        return fitz_mod.get_text_length(text, fontname=fontname, fontsize=fontsize, **kw)

    fitz_mod.get_text_length = _gtl_patched
    return fitz_mod


def hsub_label(ax, x, y, base, sub, tail="", fontsize=7, align="center",
               dx_pt=None, dy_pt=None, sub_scale=0.70, **kw):
    """Horizontal true subscript: base + a small sunk sub + tail (for subscripts in axis titles/annotations like t_sw1).

    Same approach as vsub_label: all three pieces ha=va="center", lengths measured in pixels then divided by the axes pixel size,
    NO mathtext and NO Unicode subscript digits (Arial lacks U+2080/U+2081; they would be silently swallowed).
    with align="left", x is the whole string's left end; dx_pt/dy_pt are in points relative to the axis left/bottom edges,
    where in multi-axis figures the title sits off the page edge, points are more stable than axes fractions.
    sub_scale: the subscript's size relative to the base (default 0.70; pass 1.0 in small annotations that must keep the floor size).
    """
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bb = ax.get_window_extent(renderer=r)
    W, H = bb.width, bb.height
    fs = fontsize
    gap = r.points_to_pixels(0.08 * fs)
    drop = r.points_to_pixels(0.24 * fs)
    if dx_pt is not None:
        x = r.points_to_pixels(dx_pt) / W
    if dy_pt is not None:
        y = r.points_to_pixels(dy_pt) / H

    def adv(txt, size):
        t = ax.text(0, 0, txt, fontsize=size, ha="center", va="center",
                    transform=ax.transAxes, **kw)
        fig.canvas.draw()
        w = t.get_window_extent(renderer=r).width      # horizontal: length = bbox width
        t.remove()
        return w

    a1, a2 = adv(base, fs), adv(sub, fs * sub_scale)
    a3 = adv(tail, fs) if tail else 0.0
    total = a1 + gap + a2 + (gap + a3 if tail else 0.0)
    x0 = x * W if align == "left" else x * W - total / 2.0
    out = [ax.text((x0 + a1 / 2.0) / W, y, base, fontsize=fs, ha="center", va="center",
                   transform=ax.transAxes, **kw)]
    x0 += a1 + gap
    out.append(ax.text((x0 + a2 / 2.0) / W, y - drop / H, sub, fontsize=fs * sub_scale,
                       ha="center", va="center", transform=ax.transAxes, **kw))
    if tail:
        x0 += a2 + gap
        out.append(ax.text((x0 + a3 / 2.0) / W, y, tail, fontsize=fs, ha="center",
                           va="center", transform=ax.transAxes, **kw))
    fig.canvas.draw()
    return out


FINAL_DIR = WORK + "/FY/oligo/plots-merged"


def out_dir(env, default, guard=FINAL_DIR):
    """Output directory for figures. Defaults to the draft directory; writing to the FINAL directory requires an explicit unlock.

    Why: today's two incidents were both "the script's default pointed at plots-merged and one run forgot the env var",
    so the Arial panels overwrote the final ones (looked fine, no error; only an md5 comparison caught it later).
    Flipping the default and adding an assert makes such overwrites deliberate-only.
    """
    d = os.environ.get(env, default)
    if os.path.abspath(d) == os.path.abspath(guard):
        assert os.environ.get("ALLOW_FINAL_OVERWRITE") == "1", (
            "refusing to write the final directory %s. To overwrite, back up first and explicitly set ALLOW_FINAL_OVERWRITE=1." % guard)
        print("arial_style: writing to the final directory %s" % guard)
    return d


def fix_subset_prefix_consistency(doc):
    """S1 (ACK_BATCH3_SUPPLEMENT) + B6 (ACK_BATCH4): font-naming repair after subset_fonts().

    S1: a Type0/CIDFont's /BaseFont must carry the same XXXXXX+ tag as its descriptor's /FontName
        (ISO 32000-1 section 9.6.3).
    B6: fitz generates tags by (document, source font) - when each panel container subsets its own copy, the same
        "same name, different face" arises (one GZROMK+ArialMT shared by 4-6 programs with different glyph counts),
        veraPDF's §6.3.2 same-name consistency flags it, and a RIP deduplicating by name renders with the wrong program.
        Fix: regenerate the tags from the DECOMPRESSED STREAM's md5 (base32, first 6 chars, uppercase) -
        different faces must differ in name and identical faces share one. Object-tree-only rename; the rendering byte path is untouched.
    """
    import re as _re
    import base64 as _b64
    import hashlib as _hashlib

    def _get(obj, key):
        m = _re.search(r"/" + key + r"\s+(\d+)\s+0\s+R", obj)
        return int(m.group(1)) if m else None

    def _stream_md5(xref):
        try:
            return _hashlib.md5(doc.xref_stream(xref)).hexdigest()
        except Exception:
            return None

    # ---- pass 1: each FontDescriptor gets its new tag from the embedded stream's md5
    fd_new = {}                                     # fd xref -> new FontName
    seen = {}                                       # md5 -> tag (identical faces share a tag)
    for xref in range(1, doc.xref_length()):
        obj = doc.xref_object(xref, compressed=True) or ""
        if "/FontDescriptor" not in obj or "/FontName" not in obj:
            continue
        sx = None
        for k in ("FontFile2", "FontFile", "FontFile3"):
            sx = _get(obj, k)
            if sx:
                break
        if sx is None:
            continue                                # no embedded stream: leave it (a separate gate covers non-embedded)
        h = _stream_md5(sx)
        if h is None:
            continue
        mn = _re.search(r"/FontName\s*/([^\s/]+)", obj)
        base = mn.group(1).split("+")[-1]
        tag = seen.get(h)
        if tag is None:
            tag = _b64.b32encode(_hashlib.md5(bytes.fromhex(h)).digest()).decode()[:6].upper()
            seen[h] = tag
        want = f"{tag}+{base}"
        fd_new[xref] = want
        if mn.group(1) != want:
            doc.xref_set_key(xref, "FontName", "/" + want)

    # ---- pass 2: align BaseFont along the Type0 -> DescendantFonts -> FontDescriptor and SimpleFont chains
    fixed = []
    for xref in range(1, doc.xref_length()):
        obj = doc.xref_object(xref, compressed=True) or ""
        if "/BaseFont" not in obj:
            continue
        mb = _re.search(r"/BaseFont\s*/([^\s/]+)", obj)
        if not mb:
            continue
        fdx = _get(obj, "FontDescriptor")
        if fdx is None:
            md = _re.search(r"/DescendantFonts\s*\[\s*(\d+)\s+0\s+R", obj)
            if md:
                fdx = _get(doc.xref_object(int(md.group(1)), compressed=True) or "",
                           "FontDescriptor")
        if fdx is None or fdx not in fd_new:
            continue
        want = fd_new[fdx]
        if mb.group(1) != want:
            doc.xref_set_key(xref, "BaseFont", "/" + want)
            fixed.append(want)
        if "/Type0" in obj:                          # sync the CIDFont sub-layer
            md = _re.search(r"/DescendantFonts\s*\[\s*(\d+)\s+0\s+R", obj)
            if md:
                cdx = int(md.group(1))
                cd = doc.xref_object(cdx, compressed=True) or ""
                mcb = _re.search(r"/BaseFont\s*/([^\s/]+)", cd)
                if mcb and mcb.group(1) != want:
                    doc.xref_set_key(cdx, "BaseFont", "/" + want)
    return fixed


# ---- the two-tier overlay (house spec §9; added 2026-11-12) ----
# Why an "overlay function" instead of changing the 7.0 in BASE to 6.0:
#   BASE is part of HISTORICAL REPRODUCIBILITY (§14.13). Old scripts must be able to re-render the delivered layouts; move the default
#   and every old figure (signed ones included) becomes unrereaderable. Old scripts stay untouched; only NEW/THIS-ROUND figures explicitly call
#   use_two_tier() to tighten to two tiers. font_gate.py re-audits each PDF to confirm only 8.0/6.0 remain.
TWO_TIER_SIZES = (8.0, 6.0)          # 8.0 = panel letters + row-head titles; 6.0 = everything else


def use_two_tier(letter_cm=(1.0, 9.5)):
    """Tightens rcParams to the two-tier policy; returns the previous rcParams so a panel can temporarily revert.

    Rules (§9 / §11.3 / §14.7):
      - every readable text (axis titles, ticks, corner notes, legend, table cells, footnotes) = 6.0
      - only panel letters and row-head titles = 8.0 (added by the compose layer; the panel scripts never print letters)
      - 5.5/5.6/5.7/5.8 all raise to 6.0; 7.0/6.8/6.5/6.4/6.2 all lower to 6.0
      - when space runs out NEVER shrink type: follow §9's order: breathe -> two lines -> trim -> move below the panel
    """
    import matplotlib as mpl
    saved = {k: mpl.rcParams[k] for k in ("font.size", "axes.titlesize", "axes.labelsize",
                                          "xtick.labelsize", "ytick.labelsize", "legend.fontsize",
                                          "axes.titleweight")}
    mpl.rcParams.update({
        "font.size": 6.0, "axes.titlesize": 6.0, "axes.labelsize": 6.0,
        "xtick.labelsize": 6.0, "ytick.labelsize": 6.0, "legend.fontsize": 6.0,
        "axes.titleweight": "normal",
    })
    # hard self-check: prove the two tiers are live before drawing (§14.13 requires an rcParams print assert)
    assert mpl.rcParams["font.size"] == 6.0 and mpl.rcParams["axes.labelsize"] == 6.0, "two-tier policy not in effect"
    return saved


def snap_fontsize(fs, role="body"):
    """Buckets font sizes scattered through old code into the two tiers (role="letter"/"rowhead" -> 8.0, all else -> 6.0).
    For rewriting old scripts line by line while keeping semantics: after the edit, re-render and pass font_gate."""
    return TWO_TIER_SIZES[0] if role in ("letter", "rowhead") else TWO_TIER_SIZES[1]
