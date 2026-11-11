#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""D8: the "four-layer negative disclosure" table of Fig 6 C - the four sentences printed on the v17 face, each traceable to a data row.

Background (handover sheet §7 / review receipt §10.2): v16's A cell printed "chromatin channel ranks the same 6 subtypes
ρ=0.943, P=0.005", but that 0.943 came from T7-pooled "one trait misplaced once" configuration (with n=6, Σd²=2),
NOT MS's six-point ordering; MS's own six-point cross-channel ordering is reversed. With that sentence deleted,
C must give all four layers at once, else printing only ρ=0.819 is selective reporting:
  L1 ladder (74-trait Δ1 cross-channel)  ρ = 0.819          <- concordant
  L2 MS off-diagonal residual z / rank  z = +2.09, rank 1  <- MS breaks out only in the expression layer
  L3 six-subtype mean profiles (MS six-point) ρ = -0.143, P=0.79 <- not concordant at the subtype level
  L4 per-nucleus cross-channel rho (same-nucleus 5,962)  mean ρ = 0.079 <- almost nothing at the single-nucleus level
L3 is computed on the fly from O1's cell_level_table.csv (mean per cell_type group -> six-point Spearman),
with SCZ / HEIGHT as the two controls (the two values named in the review receipt), no cherry-picking.
"""
import os, numpy as np, pandas as pd
from scipy.stats import spearmanr

# Root of the author's analysis tree (FY/, figures/, routeA_probe/, fonts/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")

O1 = os.path.join(os.environ.get('DSCI_FIG_ROOT', WORK + "/figures"), 'routeA_out/O1_own_multiome_paired/tables')
OUT = os.path.join(os.environ.get('DSCI_FIG_ROOT', WORK + "/figures"), 'routeA_out/FIG6_v17/tables')
os.makedirs(OUT, exist_ok=True)
ORDER = ["COP/NFOL", "MFOL", "MOL2", "MOL3", "MOL3/6", "MOL6"]
SAY = lambda fmt, *a: print("[%s] %s" % (pd.Timestamp.now().strftime("%H:%M:%S"), fmt % a if a else fmt), flush=True)

t6 = pd.read_csv(O1 + "/T6_cross_channel_summary.csv").set_index("item")["value"]
t1 = pd.read_csv(O1 + "/T1_same_nucleus_crosschannel_rho.csv")
cl = pd.read_csv(O1 + "/cell_level_table.csv")
SAY("cell_level_table: %s | subtypes %s", cl.shape, dict(cl.cell_type.value_counts()))

rows = [dict(line="L1 ladder (74 traits, Delta1 = median(MOL6)-median(COP/NFOL), cross-channel)",
             value=float(t6["cross-channel rho (raw)"]), P=np.nan, n=74,
             source="O1_own_multiome_paired/tables/T6_cross_channel_summary.csv",
             note="noise ceiling %.3f | attenuation-corrected %.3f | independent-halves %.3f"
                  % (t6["noise ceiling sqrt(relRNA*relATAC)"], t6["attenuation-corrected rho"],
                     t6["independent-nucleus-halves rho (median)"])),
        dict(line="L2 MS off-diagonal (robust regression residual of Delta1_rna on Delta1_atac)",
             value=float(t6["MS off-diagonal residual z"]), P=np.nan, n=74,
             source="same as above", note="rank %d of 74 (the largest of the 74 traits)" % int(t6["MS off-diagonal rank of 74"]))]

# L3: six-point cross-channel (MS + the two controls + the all-74-trait median)
trait_cols = [c for c in cl.columns if c.endswith("_rna") and c.startswith(("PASS_", "UKB_"))]
six = cl.groupby("cell_type").mean(numeric_only=True)
six = six.reindex([o for o in ORDER if o in six.index])
per_trait = []
for c in trait_cols:
    r, p = spearmanr(six[c].values, six[c.replace("_rna", "_atac")].values)
    per_trait.append(dict(trait=c[:-4], n_points=len(six), rho=r, P=p,
                          mean_rna=float(six[c].mean()), mean_atac=float(six[c.replace("_rna", "_atac")].mean())))
PT = pd.DataFrame(per_trait).sort_values("rho", ascending=False)
for nm, key in [("MS", "PASS_Multiple_sclerosis"), ("SCZ", "PASS_Schizophrenia_Pardinas2018"),
                ("HEIGHT", "UKB_460K.body_HEIGHTz")]:
    r = PT[PT.trait == key].iloc[0]
    rows.append(dict(line="L3 six-subtype mean profile, cross-channel (trait=%s)" % nm, value=float(r.rho),
                     P=float(r.P), n=int(r.n_points),
                     source="fig6C_sixpoint_all74.csv (computed on the fly from cell_level_table.csv)",
                     note="six points = %s" % "/".join(ORDER)))
rows.append(dict(line="L3 median across the %d traits stored per-nucleus (six-subtype profiles)" % len(PT),
                 value=float(PT.rho.median()), P=np.nan, n=len(PT),
                 source="fig6C_sixpoint_all74.csv", note="frac positive = %.3f" % (PT.rho > 0).mean()))

rows.append(dict(line="L4 per-nucleus cross-channel rho (same 5,962 nuclei, 74-trait vector)",
                 value=float(t1.mean_rho.iloc[0]), P=float(t1.P_emp.iloc[0]), n=int(t1.n_cell.iloc[0]),
                 source="O1_own_multiome_paired/tables/T1_same_nucleus_crosschannel_rho.csv",
                 note="median %.4f | frac positive %.3f | null mean %.4f | Z vs null %.1f"
                      % (t1.median_rho.iloc[0], t1.frac_pos.iloc[0], t1.null_mean.iloc[0], t1.Z_vs_null.iloc[0])))

D = pd.DataFrame(rows)
D["value"] = D.value.astype(float).round(4)
PT.round(4).to_csv(OUT + "/fig6C_sixpoint_all74.csv", index=False)
D.round(4).to_csv(OUT + "/fig6C_fourline_disclosure.csv", index=False)
SAY("\n" + D[["line", "value", "P", "n"]].to_string(index=False))
SAY("\ncontrols (audit §10.2 expectations MS -0.143 / SCZ +0.771 / HEIGHT +0.943): here MS=%.3f SCZ=%.3f HEIGHT=%.3f",
    PT.set_index("trait").loc["PASS_Multiple_sclerosis", "rho"],
    PT.set_index("trait").loc["PASS_Schizophrenia_Pardinas2018", "rho"],
    PT.set_index("trait").loc["UKB_460K.body_HEIGHTz", "rho"])
SAY("written %s", OUT)
