#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
MS GWAS gene-locus overlap analysis (the optional Fig 5d panel)
=============================================================================
Purpose: check whether the genes your own platform discovered (oligodendrocyte-lineage master regulators, MOL subtype markers,
mouse data -> human orthologue symbols) fall near MS GWAS genome-significant loci
(index SNP ± 500 kb cis-window).

Significance: MS-associated genetic variation physically near the gene gives genetic corroboration that the gene - and its cell state -
is disease-relevant. Complementary to scDRS (transcriptome-wide enrichment).
Notes:
  - this is "locus-overlap" evidence, not causation (one Methods sentence suffices);
  - the overlap set must be your own platform genes, never the scDRS MS gene set itself
    (that set was built from MS genetics - circular).

Pipeline:
  1) download the full GWAS Catalog associations table (includes the MS 2019 Science study,
     PMID 31604244, ~333 association records)
  2) keep that study's genome-significant SNPs (default p <= 5e-8; the code includes the equality)
  3) download the GENCODE hg38 gene annotation and take each gene's TSS
  4) for each significant SNP, find genes whose TSS lies within SNP ±500 kb -> intersect with your gene list
  5) output: an overlap table (CSV) + two figures (a genome-wide manhattan annotation + a per-gene
     locus strip plot)

-----------------------------------------------------------------------------
Input files (all auto-downloadable, or point at local copies)
-----------------------------------------------------------------------------
  [automatic] GWAS Catalog associations:
         https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/
         gwas-catalog-associations-full.zip   (~68 MB; cache it locally)
  [automatic] GENCODE hg38 annotation:
         https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/
         release_45/gencode.v45.basic.annotation.gtf.gz (~50 MB)
  [required-config] your gene list: GENES_FILE points at a file, one human symbol per line
         (recommended: human orthologues of the lineage TFs + MOL subtype markers);
         empty = the built-in default list (oligodendrocyte-lineage genes).

-----------------------------------------------------------------------------
Environment
-----------------------------------------------------------------------------
- same as the scDRS scripts: Linux server, Python 3.9+
- dependencies: numpy, pandas, matplotlib (pure CPU, finishes within tens of minutes)
- the first run downloads two files; afterwards it runs offline

Usage:
  python ms_locus_overlap.py            # the full pipeline
  python ms_locus_overlap.py --help
Output (under OUT_DIR):
  ms_locus_overlap_table.csv   per (significant SNP, your gene) window overlaps
  ms_manhattan.png/.pdf        genome-wide significant loci, hit genes in red + gene names
  ms_locus_strips.png/.pdf     a locus strip (SNP ± 500 kb) per hit gene
=============================================================================
"""

import argparse
import bisect
import gzip
import io
import os
import re
import sys
import urllib.request
import zipfile

import numpy as np
import pandas as pd

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

# ============================== configuration ==================================
MS_PUBMEDID = "31604244"      # IMSGC 2019 Science (47,429 cases / 68,375 controls; GWAS Catalog record)
#   convention check: python3 $DSCI_WORK_ROOT/scripts/MS_GWAS/check_fig6_facts.py
P_THRESH = 5e-8               # genome-significance threshold (EQUALITY INCLUDED: the code drops p > P_THRESH, so of the 178
#                                SNPs, 4 have p exactly 5e-8; strict p<5e-8 leaves 174)
WINDOW_KB = 500               # cis-window: TSS within SNP ± 500 kb
BUILD = "hg38"                # verified against Ensembl that the catalogue is GRCh38
GENES_FILE = WORK + "/FY/oligo/scdrs_ms_out/genes_of_interest.txt"  # platform-discovered genes (mouse->human orthologues)
OUT_DIR = WORK + "/FY/oligo/scdrs_ms_out"

# download URLs (existing local files skip the download)
ASSOC_URL = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations-full.zip"
GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.basic.annotation.gtf.gz"

# the built-in default gene list (lineage master regulators + MOL markers, human symbols; replace with your own)
DEFAULT_GENES = [
    "SOX10", "OLIG2", "NKX2-2", "TCF7L2", "TNR", "MAL", "MBP", "PLP1",
    "MOG", "PTGDS", "CA2", "PDGFRA", "CSPG4", "ANLN", "CNP", "MOBP",
    "NKX6-2", "GPR17", "QKI", "KLK6",
]
# =========================================================================


def log(msg):
    print(f"[ms-locus] {msg}", flush=True)


def fetch(url, dst):
    if os.path.exists(dst) and os.path.getsize(dst) > 1_000_000:
        log(f"using local cache: {dst}")
        return dst
    log(f"downloading {url} ...")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=1800) as resp, open(dst, "wb") as f:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
    log(f"saved: {dst} ({os.path.getsize(dst) / 1e6:.1f} MB)")
    return dst


# -------------------------- 1) MS significant SNPs -----------------------------
def load_ms_snps(assoc_zip):
    """Extract the MS study's significant SNPs from the GWAS Catalog associations."""
    log("parsing the full associations table, filtering for the MS study...")
    rows = []
    with zipfile.ZipFile(assoc_zip) as z:
        name = z.namelist()[0]
        with z.open(name) as raw:
            f = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
            header = next(f).rstrip("\n").split("\t")
            idx = {h: i for i, h in enumerate(header)}
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) <= idx["CHR_POS"]:
                    continue
                if cols[idx["PUBMEDID"]] != MS_PUBMEDID:
                    continue
                try:
                    p = float(cols[idx["P-VALUE"]])
                except ValueError:
                    continue
                if p > P_THRESH:
                    continue
                if not cols[idx["CHR_POS"]].strip() or not cols[idx["CHR_ID"]].strip():
                    continue
                rows.append({
                    "snp": cols[idx["SNP_ID_CURRENT"]],
                    "chr": str(cols[idx["CHR_ID"]]).replace("chr", ""),
                    "pos": int(cols[idx["CHR_POS"]]),
                    "p": p,
                    "mapped_gene": cols[idx["MAPPED_GENE"]],
                })
    df = pd.DataFrame(rows).drop_duplicates(subset=["snp", "chr", "pos"])
    log(f"MS (PMID {MS_PUBMEDID}) genome-significant SNPs: {len(df)}")
    if len(df) == 0:
        sys.exit("no significant SNPs found: check MS_PUBMEDID or P_THRESH")
    return df


# -------------------------- 2) GENCODE gene annotation ------------------------
def load_genes(gtf_gz):
    """Parse the GENCODE GTF's gene features, returning (chr, tss, name)."""
    log("parsing the GENCODE gene annotation...")
    genes = []
    with gzip.open(gtf_gz, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            attr = cols[8]
            m = re.search(r'gene_name "([^"]+)"', attr)
            if not m:
                continue
            name = m.group(1)
            chrom = cols[0].replace("chr", "")
            start, end = int(cols[3]), int(cols[4])
            tss = start if cols[6] == "+" else end
            genes.append((chrom, tss, name))
    log(f"genes: {len(genes)}")
    return genes


# -------------------------- 3) window overlaps --------------------------------
def window_overlap(ms_snps, genes, window_bp):
    """For each significant SNP find the genes whose TSS is inside the window; returns a long DataFrame."""
    log(f"computing SNP±{window_bp/1000:.0f}kb window overlaps...")
    by_chr = {}
    for chrom, tss, name in genes:
        by_chr.setdefault(chrom, []).append((tss, name))
    for c in by_chr:
        by_chr[c].sort()

    recs = []
    for _, s in ms_snps.iterrows():
        arr = by_chr.get(s["chr"], [])
        if not arr:
            continue
        tss_list = [t for t, _ in arr]
        lo = bisect.bisect_left(tss_list, s["pos"] - window_bp)
        hi = bisect.bisect_right(tss_list, s["pos"] + window_bp)
        for tss, name in arr[lo:hi]:
            recs.append({
                "gene": name, "chr": s["chr"], "snp": s["snp"],
                "snp_pos": s["pos"], "gene_tss": tss,
                "distance_bp": tss - s["pos"], "pvalue": s["p"],
            })
    out = pd.DataFrame(recs)
    if len(out):
        out = out.sort_values("pvalue").drop_duplicates(
            subset=["gene", "snp"]).reset_index(drop=True)
    return out


def main():
    ap = argparse.ArgumentParser(description="MS GWAS gene-locus overlap analysis")
    ap.add_argument("--genes-file", default=GENES_FILE,
                    help="your gene-list file (one human symbol per line)")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    genes_file = args.genes_file

    # gene list
    if genes_file and os.path.exists(genes_file):
        genes_of_interest = [l.strip().upper() for l in open(genes_file)
                             if l.strip() and not l.startswith("#")]
    else:
        genes_of_interest = DEFAULT_GENES
        log("no gene-list file given; using the built-in oligodendrocyte-lineage default list")
    log(f"target genes: {len(genes_of_interest)}")

    # data
    assoc_zip = fetch(ASSOC_URL, os.path.join(OUT_DIR, "gwas-catalog-associations-full.zip"))
    gtf_gz = fetch(GTF_URL, os.path.join(OUT_DIR, "gencode.v45.basic.annotation.gtf.gz"))

    ms_snps = load_ms_snps(assoc_zip)
    genes = load_genes(gtf_gz)
    gene_set = {n for _, _, n in genes}

    # window overlaps (all genes -> the hit list)
    ov_all = window_overlap(ms_snps, genes, WINDOW_KB * 1000)
    n_snp_with_gene = ov_all["snp"].nunique() if len(ov_all) else 0
    log(f"sanity check: {n_snp_with_gene}/{len(ms_snps)} significant SNPs have a gene in their window "
        f"({n_snp_with_gene/len(ms_snps)*100:.0f}%); >80%% means the annotation is right")

    # intersect with the target genes
    hit = ov_all[ov_all["gene"].isin(genes_of_interest)].copy()
    missing = [g for g in genes_of_interest if g not in gene_set]
    if missing:
        log(f"note: these genes are absent from the GENCODE annotation (check spelling / human symbols): {missing}")

    # cross-check against the catalogue's own MAPPED_GENE column
    mapped = set()
    for _, s in ms_snps.iterrows():
        for g in str(s["mapped_gene"]).split(","):
            if g.strip().upper() in genes_of_interest:
                mapped.add(g.strip().upper())
    log(f"target genes hit by GWAS Catalog's own MAPPED_GENE column: {sorted(mapped)}")

    if len(hit):
        # keep only the most significant SNP per gene
        hit = (hit.sort_values("pvalue")
                  .groupby("gene", as_index=False).first())
        hit.to_csv(os.path.join(OUT_DIR, "ms_locus_overlap_table.csv"), index=False)
        log("overlap table saved: ms_locus_overlap_table.csv")
        log("\nhit genes (most significant SNP):")
        for _, r in hit.iterrows():
            d = r["distance_bp"] / 1000
            print(f"  {r['gene']:8s} chr{r['chr']}  SNP {r['snp']}  "
                  f"locus distance {d:+.0f} kb   p={r['pvalue']:.1e}")
    else:
        log("no target gene falls inside an MS significant-locus window -> report 0 hits honestly,"
            "or relax WINDOW_KB / change the gene list; never force a hit")

    # ---- figure 1: genome-wide manhattan (hit SNPs in red + gene names) ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ms_snps = ms_snps.sort_values(["chr", "pos"])
    chrs = [c for c in sorted(ms_snps["chr"].unique(), key=lambda x: (len(x), x))]
    offsets, cum = {}, 0
    for c in chrs:
        offsets[c] = cum
        cum += ms_snps.loc[ms_snps["chr"] == c, "pos"].max() * 1.08
    ms_snps["x"] = ms_snps.apply(lambda r: offsets[r["chr"]] + r["pos"], axis=1)
    hit_pos = set(zip(hit["chr"], hit["snp_pos"])) if len(hit) else set()

    fig, ax = plt.subplots(figsize=(12, 4.2), constrained_layout=True)
    colors = plt.cm.tab20(np.linspace(0, 1, len(chrs)))
    for i, c in enumerate(chrs):
        sub = ms_snps[ms_snps["chr"] == c]
        ax.scatter(sub["x"], -np.log10(sub["p"]), s=7, c=[colors[i % 20]])
    # hit SNPs in red
    for i, (_, r) in enumerate(ms_snps.iterrows()):
        if (r["chr"], r["pos"]) in hit_pos:
            gene = hit.loc[hit["chr"] == r["chr"], "gene"].tolist()
            ax.scatter(r["x"], -np.log10(r["p"]), s=40, c="red", zorder=5)
            ax.annotate("/".join(gene), (r["x"], -np.log10(r["p"])),
                        xytext=(4, 6), textcoords="offset points", fontsize=8, color="red")
    ax.axhline(-np.log10(5e-8), color="grey", lw=0.8, ls="--")
    ax.text(0.995, -np.log10(5e-8), "genome-wide", va="bottom", ha="right",
            fontsize=7, color="grey")
    ax.set_xticks([offsets[c] + (ms_snps.loc[ms_snps["chr"] == c, "x"].max()
                                 - offsets[c]) / 2 for c in chrs])
    ax.set_xticklabels(chrs, fontsize=7)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("-log10(P)")
    ax.set_title(f"MS GWAS significant loci (PMID {MS_PUBMEDID}); "
                 f"red = within ±{WINDOW_KB}kb of platform-identified genes")
    fig.savefig(os.path.join(OUT_DIR, "ms_manhattan.png"), dpi=300)
    fig.savefig(os.path.join(OUT_DIR, "ms_manhattan.pdf"))
    log("figure saved: ms_manhattan.png/.pdf")

    # ---- figure 2: a locus strip plot per hit gene ----
    if len(hit):
        n = len(hit)
        fig2, axes = plt.subplots(n, 1, figsize=(9, 1.35 * n), constrained_layout=True)
        if n == 1:
            axes = [axes]
        # per gene, take all SNP-gene records inside its window (non-target genes included, as background)
        ov_gene = ov_all[ov_all["gene"].isin(hit["gene"])]
        for ax, (_, h) in zip(axes, hit.iterrows()):
            g = h["gene"]
            lo, hi = h["snp_pos"] - WINDOW_KB * 1000, h["snp_pos"] + WINDOW_KB * 1000
            bg = ov_all[(ov_all["chr"] == h["chr"]) &
                        (ov_all["snp_pos"].between(lo, hi))]
            # background genes: all genes in the window (approximated by the SNP records), TSS drawn deduplicated
            bg_genes = bg.drop_duplicates("gene")
            for _, r in bg_genes.iterrows():
                if r["gene"] == g:
                    continue
                ax.plot([r["gene_tss"], r["gene_tss"]], [0.4, 0.6],
                        color="grey", lw=1.2)
                ax.text(r["gene_tss"], 0.28, r["gene"], fontsize=6,
                        color="grey", rotation=90, ha="center")
            # the hit gene
            ax.plot([h["gene_tss"], h["gene_tss"]], [0.4, 0.6], color="red", lw=2)
            ax.text(h["gene_tss"], 0.28, h["gene"], fontsize=8, color="red",
                    fontweight="bold", rotation=90, ha="center")
            # SNP
            ax.plot([h["snp_pos"], h["snp_pos"]], [0.65, 0.85], color="black", lw=1.5)
            ax.text(h["snp_pos"], 0.9, f"SNP {h['snp']}\np={h['pvalue']:.1e}",
                    fontsize=7, ha="center")
            ax.set_xlim(lo, hi)
            ax.set_ylim(0, 1)
            ax.set_yticks([])
            ax.set_title(f"{g} — chr{h['chr']}:{lo:,}-{hi:,} "
                         f"(±{WINDOW_KB}kb around sentinel SNP)", fontsize=9)
        fig2.savefig(os.path.join(OUT_DIR, "ms_locus_strips.png"), dpi=300)
        fig2.savefig(os.path.join(OUT_DIR, "ms_locus_strips.pdf"))
        log("figure saved: ms_locus_strips.png/.pdf")

    log("done. If the manhattan's SNP positions clearly disagree with chromosome density, switch BUILD/GTF to "
        "grch37 and rerun (currently verified as hg38).")


if __name__ == "__main__":
    main()
