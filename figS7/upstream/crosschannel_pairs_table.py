# -*- coding: utf-8 -*-
"""crosschannel_pairs_table.py — cross-channel / cross-cohort d1 concordance table: the only numbers panel D and the legend may quote.

make_S7c_pair_table.py - the face-number contract for Fig S7c (dual-omics / cross-channel concordance).

Discipline: no hand-typed numbers on the face. This CSV is the complete set of numbers that the four
S7c sub-panels may show; A7_figS7_v3_verify.py reads every cell back from it, and reads its upstream

Rows = one "pair" (one Spearman between two d1 series):  6 RNA-RNA pairs (between human cohorts, the
positive control: retained after domain adjustment);  4 RNA-chromatin pairs (external human cohorts vs
GSE231584 snATAC, **not paired** - same donors, separately profiled libraries);  1 RNA-chromatin pair
(our own dsciCAT-seq, **same nuclei**, 74 traits);  1 pair mouse EAE multiome **same nuclei**, MS
single trait x 22 libraries (different unit, drawn separately on the face);  1 pair mouse chromatin
arm, 12-trait probe (1 library, stored in routeA_probe, copied into SUPPORT_v1/tables/). Authoritative
rho values come from D_crosschemistry/tables/D3_crossconcordance_domainadjust.csv (SEED 20261102,
10,000 permutations + 10,000 bootstraps; adjustment = ln n_peak, own same-nucleus pair residualised on ln n_hit_atac)."""
import sys, os

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
import os
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
sys.path.insert(0, WORK + "/scripts/routeA")
from support_common import *
import numpy as np, pandas as pd
from scipy.stats import spearmanr

D3 = read("D_crosschemistry/tables/D3_crossconcordance_domainadjust.csv")
MOUSE = read("M2_mouse_multiome_paired/tables/cross_channel_delta1_same_nuclei.csv")
PROBE = "SUPPORT_v1/tables/M3_probe_delta1_3libs.csv"
assert os.path.exists(os.path.join(R, PROBE)), "copy routeA_probe/m3_probe/tables/M3_probe_delta1_3libs.csv into SUPPORT_v1/tables/ first"
PRB = read(PROBE)

rows = []
for _, r in D3[D3.kind.isin(["RNA↔RNA", "RNA↔chromatin(separately profiled)", "RNA↔chromatin(same nucleus)"])].iterrows():
    rows.append(dict(
        pair=r.pair, kind=r.kind, channel_unit="trait (74)", n=int(r.n),
        rho_raw=r.rho_raw, P_raw=r.P_raw, rho_adj=r.rho_adj, P_adj=r.P_adj,
        perm_P_adj=r.perm_P_adj, retained=r.rho_retained,
        # reading rule: "kept" needs a significant permutation P **and** adjusted rho >= 0.5
        verdict=("kept" if (r.perm_P_adj <= 0.05 and r.rho_adj >= 0.5) else
                 ("weak but sig" if r.perm_P_adj <= 0.05 else "collapsed")),
        source="D3_crossconcordance_domainadjust.csv"))

# mouse same nucleus: MS single trait x 22 libraries (the unit is "library", not "trait", so it is
rr, pp = spearmanr(MOUSE.rna_delta1, MOUSE.chr_delta1_MSpk)
rows.append(dict(pair="Mouse GSE250589 same nucleus (MS single trait per library)", kind="RNA↔chromatin(same nucleus)",
                 channel_unit="library (22)", n=22, rho_raw=rr, P_raw=pp, rho_adj=np.nan, P_adj=np.nan,
                 perm_P_adj=np.nan, retained=np.nan,
                 verdict="n.s. (P>=0.05)" if pp > 0.05 else "p<0.05",
                 source="cross_channel_delta1_same_nuclei.csv (Spearman over 22 libraries)"))

# mouse chromatin-arm 12-trait probe (1 library EarlyCtrl4; magnitude and sign diagnostic)
r12, p12 = spearmanr(PRB.rna, PRB.atac)
ra12, pa12 = spearmanr(PRB.rna_adj, PRB.atac_adj)
rows.append(dict(pair="Mouse chromatin-arm probe 12 traits (1 library EarlyCtrl4)", kind="RNA↔chromatin(probe)",
                 channel_unit="trait (12)", n=12, rho_raw=r12, P_raw=p12, rho_adj=ra12, P_adj=pa12,
                 perm_P_adj=np.nan, retained=ra12 / r12,
                 verdict="collapsed (ATAC 12/12 negative)", source="M3_probe_delta1_3libs.csv (copied from routeA_probe)"))

out = pd.DataFrame(rows)
p = os.path.join(TABDIR, "figS7c_crosschannel_pairs.csv")
out.to_csv(p, index=False)
say("wrote %s  %d rows" % (p, len(out)))
print(out[["pair", "kind", "channel_unit", "n", "rho_raw", "rho_adj", "verdict"]].to_string(index=False))
