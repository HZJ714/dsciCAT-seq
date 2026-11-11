#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S6v8_window_permutation.py - recomputes S6b's "locus-overlap permutation null" from the SOURCE DATA, for real.

Background (receipt ROUTEA_METHODS_RECEIPT / last round's A2-class issues): v6/v7's null was DIGITISED BACKWARDS from the delivered PDF's
vector rectangles (bar height -> counts); the raw permutation vectors are not in the repo, hence:
  - the draw count could only be written as "n = 1,000" (binning loses the true count);
  - the face P = 0.26 could not align with the source-convention P = 0.27.
This script recomputes with ms_locus_overlap.py's original criteria and writes the TRUE null vectors to a table; drawing only reads tables:

  statistic = among the 331 platform genes, the fraction whose TSS falls inside the ±500 kb cis-window of any MS genome-significant SNP (IMSGC 2019, PMID 31604244,
            p <= 5e-8, hg38 CHR_POS)
  null      = random gene sets of the same size (331) drawn from the GENCODE v45 basic PROTEIN-CODING genes,
            counted by the same criterion; N_PERM draws, seed fixed
  reference = Fisher 2x2 (hit vs not-hit x platform vs genome-wide) - the same criterion as the historically delivered 0.27

Outputs under SUPPORT_v2/tables_v8/:
  figS6v8_window_permutation.csv        summary (one quantity per row + source)
  figS6v8_window_permutation_null.csv   the true null vector (n_perm rows)
Usage: python3 S6v8_window_permutation.py
"""
import os, sys, gzip, re
import numpy as np, pandas as pd
from scipy.stats import fisher_exact

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

sys.path.insert(0, WORK + "/scripts")
import ms_locus_overlap as mlo

FY = WORK + "/FY/oligo"
OUT = WORK + "/figures/routeA_out/SUPPORT_v2/tables_v8"
os.makedirs(OUT, exist_ok=True)
WINDOW = mlo.WINDOW_KB * 1000
N_PERM = int(os.environ.get("NPERM", 10000))
SEED = 20261104

# ---------- 1) MS significant SNPs (cached zip, no re-download) ----------
zip_path = f"{FY}/scdrs_ms_out/gwas-catalog-associations-full.zip"
snps = mlo.load_ms_snps(zip_path)
print("[1] MS lead SNPs (p<=5e-8, hg38): %d" % len(snps), flush=True)

# ---------- 2) GENCODE v45 basic protein-coding TSSs (built-in parser, gene_type kept) ----------
gtf = f"{FY}/scdrs_ms_out/gencode.v45.basic.annotation.gtf.gz"
recs = []
with gzip.open(gtf, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] != "gene":
            continue
        at = c[8]
        nm = re.search(r'gene_name "([^"]+)"', at)
        tp = re.search(r'gene_type "([^"]+)"', at)
        if not nm or tp.group(1) != "protein_coding":
            continue
        chrom = c[0].replace("chr", "")
        tss = int(c[3]) if c[6] == "+" else int(c[4])
        recs.append((chrom, tss, nm.group(1).upper()))
G = pd.DataFrame(recs, columns=["chr", "tss", "name"]).drop_duplicates("name")
print("[2] protein-coding gene universe: %d (historical convention 20,049)" % len(G), flush=True)

# ---------- 3) the window-hit gene set ----------
hit = set()
by_chr = {c: g.sort_values("tss") for c, g in G.groupby("chr")}
tss_arr = {c: g["tss"].values for c, g in by_chr.items()}
name_arr = {c: g["name"].values for c, g in by_chr.items()}
for _, s in snps.iterrows():
    ch = str(s["chr"])
    if ch not in tss_arr:
        continue
    arr = tss_arr[ch]
    lo, hi = np.searchsorted(arr, s["pos"] - WINDOW), np.searchsorted(arr, s["pos"] + WINDOW, "right")
    hit.update(name_arr[ch][lo:hi])
n_uni, n_hit = len(G), len(hit)
print("[3] genes in the windows: %d / %d = %.2f%% (delivered 2,071 / 20,049 = 10.3%%)"
      % (n_hit, n_uni, 100.0 * n_hit / n_uni), flush=True)

# ---------- 4) the 331 platform genes ----------
# True source = ms_locus_overlap.GENES_FILE (human orthologue symbols, one per line).
# The previous version wrongly used scmega_p2g_S1markers.csv: a MOUSE symbol table (Apod/Timp3...),
# intersecting the historical 331 in only 306, with 21 symbols absent from the GENCODE v45 protein-coding universe
# (RIK clone IDs / old symbols), undercounting hits by one (37 vs 38) and misnaming the denominator.
plat_file = mlo.GENES_FILE
P331 = set(l.strip().upper() for l in open(plat_file) if l.strip())
missing = [g for g in P331 if g not in set(G["name"])]
obs_set = P331 & hit
obs = len(obs_set)
hist = set(pd.read_csv(f"{FY}/scdrs_ms_out/ms_locus_overlap_table.csv")["gene"].astype(str).str.upper())
print("[4] platform genes %d | absent from the universe: %d | window hits %d = %.1f%% (historical %d / %.1f%%)"
      % (len(P331), len(missing), obs, 100.0 * obs / len(P331), len(hist), 100.0 * len(hist) / len(P331)), flush=True)
print("    vs the historical hit table: reproduced %d/%d | extra this run %s | missing %s"
      % (len(obs_set & hist), len(hist), sorted(obs_set - hist), sorted(hist - obs_set)), flush=True)

# ---------- 5) Fisher + permutation ----------
OR, pf = fisher_exact([[obs, len(P331) - obs], [n_hit - obs, n_uni - n_hit - (len(P331) - obs)]])
rng = np.random.default_rng(SEED)
names = G["name"].values
null = np.empty(N_PERM, dtype=int)
for i in range(N_PERM):
    sel = rng.choice(names, size=len(P331), replace=False)
    null[i] = np.fromiter((g in hit for g in sel), dtype=int, count=len(sel)).sum()
P_perm = (1 + np.sum(null >= obs)) / (1 + N_PERM)
pd.DataFrame([dict(stat="MS lead SNPs (p<=5e-8, hg38)", value=len(snps), source="GWAS Catalog PMID 31604244 via ms_locus_overlap.load_ms_snps"),
              dict(stat="window half-width (kb)", value=mlo.WINDOW_KB, source="ms_locus_overlap.WINDOW_KB"),
              dict(stat="universe: protein-coding genes (GENCODE v45 basic)", value=n_uni, source=gtf),
              dict(stat="genes inside any window", value=n_hit, source="computed here"),
              dict(stat="background window rate (%)", value=round(100.0 * n_hit / n_uni, 3), source="computed here"),
              dict(stat="platform genes (denominator)", value=len(P331), source=plat_file),
              dict(stat="historical hits reproduced", value=len(obs_set & hist), source="FY/oligo/scdrs_ms_out/ms_locus_overlap_table.csv"),
              dict(stat="hits only in recompute", value=";".join(sorted(obs_set - hist)) or "none", source="computed here"),
              dict(stat="hits only in historical table", value=";".join(sorted(hist - obs_set)) or "none", source="computed here"),
              dict(stat="platform genes in window (observed count)", value=obs, source="computed here"),
              dict(stat="observed (%)", value=round(100.0 * obs / len(P331), 3), source="computed here"),
              dict(stat="Fisher OR", value=round(float(OR), 4), source="computed here (2x2 platform vs universe)"),
              dict(stat="Fisher P", value=float(pf), source="computed here"),
              dict(stat="permutation draws", value=N_PERM, source="computed here, seed=%d" % SEED),
              dict(stat="permutation mean count", value=float(null.mean()), source="computed here"),
              dict(stat="permutation 2.5%% count", value=float(np.quantile(null, .025)), source="computed here"),
              dict(stat="permutation 97.5%% count", value=float(np.quantile(null, .975)), source="computed here"),
              dict(stat="permutation P (one-sided, observed vs null)", value=float(P_perm), source="computed here"),
              dict(stat="observed (%) in null units", value=round(100.0 * obs / len(P331), 3), source="for the marker on the panel"),
              ]).to_csv(f"{OUT}/figS6v8_window_permutation.csv", index=False)
pd.DataFrame(dict(null_hit_count=null, null_pct=np.round(100.0 * null / len(P331), 4))).to_csv(
    f"{OUT}/figS6v8_window_permutation_null.csv", index=False)
print("[5] obs %d/331 = %.1f%% | background %.2f%% | Fisher P=%.3f (OR %.2f) | perm P=%.4f (n=%d, seed=%d) | null mean %.1f [%.0f, %.0f]"
      % (obs, 100.0 * obs / len(P331), 100.0 * n_hit / n_uni, pf, OR, P_perm, N_PERM, SEED,
         null.mean(), np.quantile(null, .025), np.quantile(null, .975)), flush=True)
print("    null vector written to figS6v8_window_permutation_null.csv (%d rows)" % N_PERM, flush=True)
