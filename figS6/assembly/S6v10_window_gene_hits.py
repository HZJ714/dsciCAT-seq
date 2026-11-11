# -*- coding: utf-8 -*-
"""S6 v10：A/B 格逐基因命中清单落附表（作者指令：A 主要展示哪些基因在 TSS 窗口内，
大部分未富集 -> 明细进 Table S6，图面只留 4 个例子位点）。

复用 S6v8_window_permutation.py 的同一判据（IMSGC lead SNP ±500 kb, hg38, GENCODE v45 basic
protein_coding），对 331 个平台基因逐个输出：是否命中 + 最近 lead SNP 的 chr:pos 与距离。
产出 tables_v8/figS6v10_window_gene_hits.csv（331 行）。
"""
import os, sys, gzip, re
import numpy as np, pandas as pd

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

sys.path.insert(0, WORK + "/scripts")
import ms_locus_overlap as mlo

FY = WORK + "/FY/oligo"
OUT = WORK + "/figures/routeA_out/SUPPORT_v2/tables_v8"
WINDOW = mlo.WINDOW_KB * 1000

snps = mlo.load_ms_snps(f"{FY}/scdrs_ms_out/gwas-catalog-associations-full.zip")
gtf = f"{FY}/scdrs_ms_out/gencode.v45.basic.annotation.gtf.gz"
recs = []
with gzip.open(gtf, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] != "gene":
            continue
        nm = re.search(r'gene_name "([^"]+)"', c[8])
        tp = re.search(r'gene_type "([^"]+)"', c[8])
        if not nm or tp.group(1) != "protein_coding":
            continue
        chrom = c[0].replace("chr", "")
        tss = int(c[3]) if c[6] == "+" else int(c[4])
        recs.append((chrom, tss, nm.group(1).upper()))
G = pd.DataFrame(recs, columns=["chr", "tss", "name"]).drop_duplicates("name")

plat_file = f"{FY}/scdrs_ms_out/genes_of_interest.txt"
P331 = [l.strip().upper() for l in open(plat_file) if l.strip()]
gmap = {n: (c, t) for c, t, n in G.itertuples(index=False) for n in [n]}
snp_by_chr = {c: np.sort(g["pos"].values) for c, g in snps.groupby(snps["chr"].astype(str))}

rows = []
for gname in P331:
    key = gname.upper()
    if key not in gmap:
        rows.append(dict(gene=gname, in_window="NA(gene not in universe)", chr="", tss="",
                         nearest_snp="", dist_kb="")); continue
    ch, tss = gmap[key]
    if ch not in snp_by_chr:
        rows.append(dict(gene=gname, in_window=0, chr=ch, tss=tss, nearest_snp="none", dist_kb="")); continue
    a = snp_by_chr[ch]
    i = int(np.searchsorted(a, tss))
    cand = [j for j in (i - 1, i, min(i + 1, len(a) - 1)) if 0 <= j < len(a)]
    j = min(cand, key=lambda j: abs(a[j] - tss))
    d = abs(int(a[j]) - tss)
    rows.append(dict(gene=gname, in_window=int(d <= WINDOW), chr=ch, tss=tss,
                     nearest_snp="%s:%d" % (ch, int(a[j])), dist_kb=round(d / 1000.0, 1)))
T = pd.DataFrame(rows)
T.to_csv(OUT + "/figS6v10_window_gene_hits.csv", index=False)
print("[hits] 331 平台基因：%d 命中窗口（<=500kb），%d 未命中，%d 不在基因全集"
      % ((T.in_window == 1).sum(), (T.in_window == 0).sum(),
         (T.in_window.astype(str) == "NA(gene not in universe)").sum()))
print("written:", OUT + "/figS6v10_window_gene_hits.csv")

# ---- 追加进 Table_S6（先备份）----
TAB = WORK + "/figures/routeA_out/SUPPORT_v2/tables/Table_S6.csv"
BK = TAB.replace(".csv", "_v8_backup.csv")
if not os.path.exists(BK):
    import shutil; shutil.copy(TAB, BK)
old = pd.read_csv(TAB)
old = old[~old["check"].astype(str).str.contains("B gene-level window hits", na=False)]
hits = T[T.in_window == 1].sort_values("dist_kb")
add = pd.DataFrame([dict(check="B gene-level window hits (331 platform genes)",
                         statistic=g.gene, value="%.1f kb" % g.dist_kb,
                         ci95_or_IQR="in window (<=500 kb)",
                         reference="%s (chr%s TSS %d)" % (g.nearest_snp, g.chr, g.tss),
                         rejects_zero="", meets_effect="", n="331 platform genes",
                         source=OUT + "/figS6v10_window_gene_hits.csv",
                         source_column="gene, in_window, dist_kb, nearest_snp",
                         in_figure="Fig S6v10 A/B (list on figure? no - Table only)",
                         in_main_text="") for g in hits.itertuples()])
miss = T[~(T.in_window == 1)]
add = pd.concat([add, pd.DataFrame([dict(
    check="B gene-level window hits (331 platform genes)",
    statistic="NON-hit genes (summary row)",
    value="%d of %d not in any window" % (len(miss), len(T)),
    ci95_or_IQR="full per-gene list in source csv",
    reference="background rate 10.34%; perm P = 0.277",
    rejects_zero="no", meets_effect="", n="331 platform genes",
    source=OUT + "/figS6v10_window_gene_hits.csv", source_column="in_window",
    in_figure="Fig S6v10 B (counts only)", in_main_text="")])], ignore_index=True)
pd.concat([old, add], ignore_index=True).to_csv(TAB, index=False)
print("[Table_S6] %d -> %d 行（+38 命中明细 +1 汇总）" % (len(old), len(old) + len(add)))
