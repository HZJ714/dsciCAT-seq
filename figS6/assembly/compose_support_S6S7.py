# -*- coding: utf-8 -*-
"""把 S6 / S7 的格子拼成整页附图（矢量 show_pdf_page，不栅格化）。

输出落在本会话目录（草稿，字母号待与主图会话统一）：
  SUPPORT_v1/pdf/figS6_composed_v10_draft.pdf/.png
  SUPPORT_v1/pdf/figS7_composed_v1_draft.pdf/.png
布局规则：按行给入面板；同一行的面板按宽度并排，行高取最高者，行间距 0.22 cm。
"""
import os, sys

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
sys.path.insert(0, WORK + "/scripts/routeA")
sys.path.insert(0, WORK + "/scripts/MS_GWAS/fig6_repro/deps")
import fitz
from arial_style import patch_fitz_arial, FONTS_PDF
patch_fitz_arial(fitz)

CM = 72 / 2.54
GAP = 0.22
OLD = WORK + "/FY/oligo/fig6_rework_v2"
SUP = WORK + "/figures/routeA_out/SUPPORT_v1/pdf"
D7 = WORK + "/figures/routeA_out/FIG7_v1/pdf"
D6 = WORK + "/figures/routeA_out/FIG6_v20/tables"


def size(p):
    d = fitz.open(p); r = d[0].rect; d.close()
    return r.width, r.height


A4_H = 29.7          # 交付件实测：Figure_S6_FINAL_v7.pdf = 17.00 x 29.55 cm 单页（A4 竖版内）


def compose_pages(pages, out_pdf, out_png=None, title=None):
    """pages = [ [row, ...], ... ]，row = [(path, letter), ...]。

    11-04 统一尺寸（用户「图片大小和面板字母也要统一」）：照交付件定的规矩 ——
    页宽一律取最大行宽（S6/S7 都是 17.00 cm），页高取各页内容高的最大值（页页一致），
    内容高超过 A4 可用高度就**必须分页**，不许再出 65 cm 的长条页。
    字母由本层统一落笔：ArialMT 8.0，origin = (格左 + 1.0 pt, 格顶 + 11.5 pt)。
    """
    sized_pages = []
    for rows in pages:
        sized, ph, rwmax = [], 0.0, 0.0
        for row in rows:
            items = []
            for p, k in row:
                w, h = size(p)
                items.append((p, k, w / CM, h / CM))
            rw = sum(i[2] for i in items) + GAP * (len(items) - 1)
            rh = max(i[3] for i in items)
            sized.append((items, rw, rh)); ph += rh + GAP
            rwmax = max(rwmax, rw)
        sized_pages.append((sized, rwmax, ph - GAP))
    W = max(rw for _, rw, _ in sized_pages)
    H = max(ph for _, _, ph in sized_pages)
    if H > A4_H:
        raise SystemExit("[尺寸门] 页高 %.2f cm 超过 A4 可用 %.1f cm —— 必须分页" % (H, A4_H))
    doc = fitz.open()
    for sized, _, _ in sized_pages:
        page = doc.new_page(width=W * CM, height=H * CM)
        if title:
            page.insert_text((4, 11), title, fontname="Arial-Regular", fontfile=FONTS_PDF["normal"],
                             fontsize=8, color=(0, 0, 0))
        y = 0.0
        for items, rw, rh in sized:
            x = 0.0
            for p, k, w, h in items:
                src = fitz.open(p)
                page.show_pdf_page(fitz.Rect(x * CM, y * CM, (x + w) * CM, (y + h) * CM), src, 0)
                src.close()
                if k:
                    page.insert_text((x * CM + 1, y * CM + 11.5), k, fontname="Arial-Regular",
                                     fontfile=FONTS_PDF["normal"], fontsize=8, color=(0, 0, 0))
                x += w + GAP
            y += rh + GAP
    doc.save(out_pdf)
    d = fitz.open(out_pdf)
    print("%-46s %.2f x %.2f cm  (%d 页 / %d 行)" % (os.path.basename(out_pdf), W, H,
                                                    d.page_count, sum(len(p) for p in pages)))
    if out_png and os.environ.get("PNGOUT"):
        for i, pg in enumerate(d):
            suf = "" if d.page_count == 1 else "_p%d" % (i + 1)
            pg.get_pixmap(matrix=fitz.Matrix(150 / 72, 150 / 72)).save(
                os.path.join(os.environ["PNGOUT"],
                             os.path.splitext(os.path.basename(out_png))[0] + suf + ".png"))
    d.close()
    return out_pdf


def compose(rows, out_pdf, out_png=None, title=None):
    """rows = [[(path, letter), ...], ...] —— 单页包装（S6 走这条）。"""
    return compose_pages([rows], out_pdf, out_png, title)

S6 = [
    [(OLD + "/figS6a_strips.pdf", "a")],
    [(OLD + "/figS6b_perm.pdf", "b"), (SUP + "/S6c_mol3_refit.pdf", "c")],
    [(OLD + "/figS6d_umap_cells.pdf", "d"), (OLD + "/figS6e_distance.pdf", "e")],
    [(OLD + "/figS6f_gridsize.pdf", "f")],
    [(SUP + "/S6g_support_C.pdf", "")],      # g/h 自带 "Fig S6g/S6h" 标题行，拼版不再叠字母
    [(SUP + "/S6h_support_D.pdf", "")],
]
S7 = [
    [(SUP + "/S7a_provenance.pdf", "A")],     # 11-04：字母统一 Arial 常规 8 pt 大写（与交付 Figure_6/S6 同一样式），
    [(SUP + "/S7b_rank_matrix.pdf", "B")],    # 格内 "Fig S7x …" 标题行已撤下图面、移交图注
    [(SUP + "/S7c_three_axes.pdf", "C")],
    [(SUP + "/S7d_disease_support.pdf", "D")],
]

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    if which in ("s6", "both"):
        compose(S6, SUP + "/figS6_composed_v11_draft.pdf", SUP + "/figS6_composed_v11_draft.png")
    if which in ("s7", "both"):
        # 每格一页：内容高 12.60 / 17.40 / 17.40 / 17.40 cm，合成 4 页 17.00 x 17.40
        compose_pages([[r] for r in S7], SUP + "/figS7_composed_v1_draft.pdf",
                      SUP + "/figS7_composed_v1_draft.png")
