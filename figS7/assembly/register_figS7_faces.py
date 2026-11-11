#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""register_figS7_faces.py -- face-number register for Fig S7 v30 (house spec SS3, SS14.13).

Every number printed on the figure must trace to a delivered table row, and every table row the figure
claims must be on the page.  This script does both directions:

  forward  - each register entry is computed HERE from the vendored table (not copied from the figure),
             formatted the way the panel prints it, and then looked for on the composed page;
  reverse  - every numeric token on the page must be either a registered value or an explicitly listed
             axis tick / colour-bar endpoint / dataset accession.  An unexplained number is a failure,
             not a warning: that is how a stale annotation survives a relayout.

Outputs: figS7_out/figS7v30_face_numbers.csv  (panel, on_face, source_table, source_row)

Run from the figS7/ directory with the figure kernel:
  DSCI_WORK_ROOT=<work> python3 assembly/register_figS7_faces.py
"""
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import fitz

CM = 72 / 2.54
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DATA = os.environ.get("FIGS7_DATA", os.path.join(REPO, "data"))
TABD = os.environ.get("FIGS7_TABLES", os.path.join(REPO, "tables_v29"))
PDF = os.environ.get("FIGS7_PDF", os.path.join(REPO, "figS7_out", "figS7_v30.pdf"))
OUT = os.path.join(REPO, "figS7_out", "figS7v30_face_numbers.csv")

MS = "PASS_Multiple_sclerosis"
COH = ["A2", "A3", "A4", "A7"]
RHO = "\u03c1"          # the face uses the Greek letter, not "rho"
ARROW = "\u2192"          # the rank arrow printed on the face
PAIRARROW = "\u2194"      # the channel-pair arrow used by the pairs table
D1 = "\u03941"

rows = []


def reg(panel, on_face, table, source_row):
    rows.append(dict(panel=panel, on_face=on_face, source_table=table, source_row=source_row))


def read(rel, **kw):
    p = os.path.join(DATA, rel)
    assert os.path.exists(p), "missing vendored input %s" % p
    return pd.read_csv(p, **kw)


# ------------------------------------------------------------------ A  the eight MS ranks
S = read("FIG7_v1/tables/seven4_confound_ALL_cohorts.csv").set_index("trait")
iMS = int(S.index.get_loc(MS))
raw = [int(S["rank_raw_%s" % c].iloc[iMS]) for c in COH]
adj = [int(S["rank_adj_%s" % c].iloc[iMS]) for c in COH]
assert raw == [2, 4, 7, 2] and adj == [1, 1, 3, 1], (raw, adj)
for c, v in zip(COH, raw):
    reg("A", str(v), "FIG7_v1/tables/seven4_confound_ALL_cohorts.csv",
        "%s rank_raw_%s = %d (MS row, as computed)" % (MS, c, v))
for c, v in zip(COH, adj):
    reg("A", str(v), "FIG7_v1/tables/seven4_confound_ALL_cohorts.csv",
        "%s rank_adj_%s = %d (MS row, size trend removed)" % (MS, c, v))

# ------------------------------------------------------------------ C  human same-donor RNA x ATAC
S4 = read("FIG7_v1/tables/seven4_confound_ALL_cohorts.csv")
PK = read("A6_GSE231584_atac_v1/tables/delta1_vs_peakdomain_allchannels.csv")
HUM = S4.merge(PK[["trait", "median_groupdiff"]], on="trait")      # exactly the plotted columns
CP = read("SUPPORT_v1/tables/figS7c_crosschannel_pairs.csv")
i = int(HUM.index[HUM.trait == MS][0])
x, y = HUM["A7"].values, HUM["median_groupdiff"].values
rank_x, rank_y = int((x > x[i]).sum()) + 1, int((y > y[i]).sum()) + 1
assert (rank_x, rank_y) == (2, 36), (rank_x, rank_y)
u = CP[CP.pair.astype(str).str.contains("GSE231585 " + PAIRARROW + " GSE231584", regex=False)].iloc[0]
assert abs(float(u.rho_raw) - spearmanr(x, y).statistic) < 1e-4, "the plotted cloud is not the rho"
reg("C", str(rank_x), "FIG7_v1/tables/seven4_confound_ALL_cohorts.csv",
    "MS rank on the plotted RNA column (A7) = %d" % rank_x)
reg("C", str(rank_y), "A6_GSE231584_atac_v1/tables/delta1_vs_peakdomain_allchannels.csv",
    "MS rank on the plotted ATAC column (median_groupdiff) = %d; the donor-paired table says 21 and "
    "the size-matched series says 36, and the face must show the plotted series" % rank_y)
reg("C", "%.3f" % u.rho_raw, "SUPPORT_v1/tables/figS7c_crosschannel_pairs.csv",
    "pair GSE231585 %s GSE231584(snATAC): rho_raw (recomputed from the plotted columns above)" % PAIRARROW)
reg("C", "%.3f" % u.rho_adj, "SUPPORT_v1/tables/figS7c_crosschannel_pairs.csv",
    "same pair: rho_adj after residualising both channels on peak-domain size")
reg("C", "%.2f" % u.P_adj, "SUPPORT_v1/tables/figS7c_crosschannel_pairs.csv",
    "same pair: P_adj printed as 'P %.2f'" % u.P_adj)

# ------------------------------------------------------------------ D  mouse same-nucleus pair
SNU = read("SUPPORT_v1/tables/S7_mouse_same_nucleus_delta1.csv")
import json
SNUJ = json.load(open(os.path.join(DATA, "SUPPORT_v1/tables/S7_mouse_same_nucleus_delta1.json")))
rho = spearmanr(SNU.delta1_rna.values, SNU.delta1_atac.values)
assert abs(rho.statistic - SNUJ["spearman_rho"]) < 1e-6 and abs(rho.pvalue - SNUJ["spearman_p"]) < 1e-9
reg("D", "%.3f" % SNUJ["spearman_rho"], "SUPPORT_v1/tables/S7_mouse_same_nucleus_delta1.json",
    "spearman_rho over the 74 same-nucleus traits (recomputed from the CSV)")
reg("D", "%.2g" % SNUJ["spearman_p"], "SUPPORT_v1/tables/S7_mouse_same_nucleus_delta1.json",
    "spearman_p, printed as 'P %.2g'")
for col, lab in [("delta1_rna", "RNA"), ("delta1_atac", "ATAC")]:
    v = int((SNU[col].values > SNU[col].values[int(SNU.index[SNU.trait == MS][0])]).sum()) + 1
    reg("D", str(v), "SUPPORT_v1/tables/S7_mouse_same_nucleus_delta1.csv",
        "MS rank on the plotted %s column = %d (json says rna %d / atac %d)"
        % (lab, v, SNUJ["ms_rank_rna"], SNUJ["ms_rank_atac"]))
    assert v == (SNUJ["ms_rank_rna"] if lab == "RNA" else SNUJ["ms_rank_atac"])

# ------------------------------------------------------------------ F  population denominators
TS = read("A3_GSE118257_v1/tables/A3_GSE118257_target_stats.csv")
for g, n in TS.groupby("group").cells.max().items():
    reg("F", "{:,}".format(int(n)), "A3_GSE118257_v1/tables/A3_GSE118257_target_stats.csv",
        "population %s: nuclei in the column = the denominator every %% in that column is taken over" % g)

# ------------------------------------------------------------------ G  per-donor rule tallies
S41 = pd.read_csv(os.path.join(TABD, "S41_summary_by_cohort.csv"), index_col=0)
for coh, r in S41.iterrows():
    if coh == "ALL":
        continue
    reg("G", "%d/%d" % (int(r.rule_met), int(r.donors_evaluable)), "S41_summary_by_cohort.csv",
        "%s: rule met / donors evaluable" % coh)

# the red dashed arm in G is the margin half of the frozen rule; it is a printed quantity, so it is
# registered here and the rule is re-derived from the per-donor table rather than taken on trust
DN = pd.read_csv(os.path.join(TABD, "S41_S1PR5_per_donor.csv"))
E = DN[DN.evaluable].copy()
E["re"] = (E.margin_pp >= 1.0) & ((E.fold >= 1.5) | (E.binding_pct == 0))
assert len(E) == 27 and bool((E["re"] == E["rule"]).all()), \
    "the printed +1.0 pp arm no longer matches the rule column of S41_S1PR5_per_donor.csv"
reg("G", "+1.0", "S41_S1PR5_per_donor.csv",
    "margin arm of the frozen rule (margin_pp >= 1.0 AND (fold >= 1.5 OR binding population = 0)); "
    "re-derived here from all %d evaluable donors and matched to the rule column" % len(E))

# ------------------------------------------------------------------ colour bar endpoints
reg("RB", "1", "FIG7_v1/tables/seven4_confound_ALL_cohorts.csv", "rank scale lower end")
reg("RB", "74", "FIG7_v1/tables/seven4_confound_ALL_cohorts.csv",
    "rank scale upper end = n traits in the sweep = %d" % len(S))

REG = pd.DataFrame(rows, columns=["panel", "on_face", "source_table", "source_row"])

# ------------------------------------------------------------------ the two gates
doc = fitz.open(PDF)
pg = doc[0]
words = [(w[0] / CM, (pg.rect.height - w[3]) / CM, w[4]) for w in pg.get_text("words")]
page_text = " ".join(w[2] for w in words)

# (1) forward: registered strings must be on the page
missing = [r for r in REG.itertuples() if r.on_face not in page_text]

# (2) reverse: every numeric token must be registered or an allowed tick / label
TICKS = {
    # axis ticks and colour-bar endpoints, each stated in the panel code that draws it
    "\u22125": "C x-axis tick (RNA d1 range)", "\u22120.50": "C y-axis tick (ATAC d1)",
    "\u22120.25": "C y-axis tick", "0.00": "C y-axis tick", "0.25": "C y-axis tick",
    "0.50": "C y-axis tick", "0": "tick / rank-0 label", "\u22121": "D x-axis tick (nucleus d1)",
    "1": "colour-bar / tick", "\u22122.5": "D y-axis tick", "0.0": "D x-axis tick",
    "2.5": "D x-axis tick", "4": "G x-axis tick (comparator nuclei)", "8": "G x-axis tick",
    "12": "G x-axis tick", "1%": "F bubble size legend", "10%": "F bubble size legend",
    "50%": "F bubble size legend", "2.2": "F expression colour-bar endpoint",
    "4.4": "F expression colour-bar endpoint", "74)": "RB '(1 of 74)' closing token",
    "(1": "RB '(1 of 74)' opening token",
}
# A numeric token is a quantity, not a label that happens to contain digits: CD52, ITGA4, NFE2L2 and
# Oligo1-6 are gene / population names and belong to the word ledger, not to the number register.
NUMTOK = re.compile(r"^[+\u2212-]?\d{1,3}(,\d{3})*(\.\d+)?%?\)?$", re.UNICODE)
unreg = set()
regset = set(REG.on_face)
for _, _, t in words:
    if not NUMTOK.match(t):
        continue
    core = t.strip("()")
    if core in regset or t in TICKS or core in TICKS:
        continue
    if re.match(r"^GSE\d+$", core):      # dataset accessions: labels, registered in MANIFEST/README
        continue
    unreg.add(t)

REG.to_csv(OUT, index=False)
print("register rows %d -> %s" % (len(REG), OUT))
print("page %.2f x %.2f cm | numeric tokens %d"
      % (pg.rect.width / CM, pg.rect.height / CM,
         sum(1 for w in words if NUMTOK.match(w[2]))))
print("FORWARD  printed strings not found on page: %d %s" % (len(missing), [m.on_face for m in missing][:8]))
print("REVERSE  numeric tokens with no register entry and no tick role: %d %s" % (len(unreg), sorted(unreg)[:10]))
if missing or unreg:
    sys.exit("register gate FAILED")
print("register gate PASS")
