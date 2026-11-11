#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""s32_flags_v2.py -- evidence flags, author-side re-audit (review receipt §9.2/§9.3).

The parallel environment could not open openFDA/DailyMed/DrugBank, so 9 flag squares were drawn
half-filled ("reported elsewhere, not re-read here").  The author-side environment re-read all 18
rows x 3 flags against openFDA label text + Europe PMC and returned a value for every square, so
the half state disappears.  This script stores those values verbatim, together with the quotable
source, and writes the diff against the previous (local-audit) values.

Rule A (label/MOA naming standard): a label names the target if it uses the receptor/complex name
("type I interferon receptor" -> IFNAR1/2, "melanocortin receptors" -> MC2R, "S1P receptors 1 and 5"
-> S1PR1/5); a drug-CLASS name does not count ("Glucocorticoids" -> NR3C1 = 0).
Rule B (MS flag): MS = 1 if an approved OR a phase II/III MS trial acts through this target.
"""
import os
import pandas as pd

# Root of the author's analysis tree (figures/, routeA_out/, routeA_probe/).
# Defaults are cwd-relative; export DSCI_WORK_ROOT to point it at your copy of
# that tree (repository README, "Data layout").
WORK = os.environ.get("DSCI_WORK_ROOT", ".")
R = WORK + "/figures/SUPPL_ANALYSES_20261112"
REV = "author-side openFDA + Europe PMC re-audit, 2026-11-17 (receipt §9.2)"

V2 = [
    # gene, label, MOA, MS, source
    ("S1PR5",   1, 1, 1, 'openFDA MAYZENT, mechanism_of_action: "binds with high affinity to S1P receptors 1 and 5"'),
    ("S1PR1",   1, 1, 1, 'openFDA MAYZENT/PLEGRIDY MOA + ponesimod PMID 42561748'),
    ("S1PR2",   0, 0, 0, 'openFDA GILENYA label: "S1P receptor modulator", no subtype named -> target not named'),
    ("S1PR3",   0, 0, 0, 'openFDA GILENYA label, no subtype named (S1P3 = off-target/safety, not claimed as target)'),
    ("S1PR4",   0, 0, 0, 'openFDA GILENYA label, no subtype named'),
    ("IFNAR1",  1, 1, 1, 'openFDA AVONEX pharmacodynamics: "type I interferon receptor ... binding to specific receptors" (Rule A)'),
    ("IFNAR2",  1, 1, 1, 'openFDA AVONEX pharmacodynamics, same sentence; IFNAR1/IFNAR2 scored symmetrically (Rule A)'),
    ("CD52",    1, 1, 1, 'openFDA LEMTRADA indications: "CD52-directed cytolytic monoclonal antibody"'),
    ("DCK",     1, 1, 1, 'openFDA MAVENCLAD clinical_pharmacology: "phosphorylated to Cd-AMP by deoxycytidine kinase"'),
    ("ITGA4",   1, 1, 1, 'openFDA TYSABRI indications: "integrin receptor antagonist" (MOA names a4b1-integrin/VLA-4)'),
    ("MS4A1",   1, 1, 1, 'openFDA OCREVUS indications: "CD20-directed cytolytic antibody"'),
    ("DHODH",   1, 1, 1, 'openFDA AUBAGIO mechanism_of_action: "inhibits dihydroorotate dehydrogenase"'),
    ("NFE2L2",  1, 1, 1, 'openFDA TECFIDERA mechanism_of_action: "activation of Nrf2 pathway in vitro and in vivo"'),
    ("NR3C1",   0, 0, 1, 'openFDA SOLU-MEDROL: no "glucocorticoid receptor" in label, class name "Glucocorticoids" only (Rule A -> 0/0); MS = 1 from 18 SOLU-MEDROL labels with MS indication'),
    ("TOP2A",   1, 1, 1, 'openFDA NOVANTRONE: "potent inhibitor of topoisomerase II" + MS-worsening indication'),
    ("MC2R",    1, 1, 1, 'openFDA ACTHAR mechanism_of_action: "bind to melanocortin receptors" (Rule A)'),
    ("BTK",     0, 1, 1, 'no approved label (correctly 0); MOA PMID 42075951; MS = 1 by Rule B, tolebrutinib MS review PMID 42894057 (DOI 10.1007/s10072-026-09459-6)'),
    ("CD40LG",  0, 1, 1, 'no approved label (correctly 0); MOA frexalimab "CD40L Antagonist" phase 2 MS trial PMID 42203212 (DOI 10.1111/dom.70785); MS = 1 by Rule B'),
]
F = pd.DataFrame(V2, columns=["gene", "label_v2", "MOA_v2", "MS_v2", "source_v2"])
F["reviewed_by"] = REV
F.to_csv(R + "/tables/S29_flags_v2.csv", index=False)

OLD = pd.read_csv(R + "/tables/S26_rowset_C12.csv")[["gene", "label", "MOA", "MS"]]
rows = []
for r in F.itertuples():
    o = OLD[OLD.gene == r.gene].iloc[0]
    for k, new in [("label", r.label_v2), ("MOA", r.MOA_v2), ("MS", r.MS_v2)]:
        old = int(o[k])
        if old != new:
            rows.append(dict(gene=r.gene, flag=k, old=old, new=new,
                             kind=("half -> resolved" if old == 2 else
                                   ("upgrade 0 -> 1" if new == 1 else "downgrade -> 0")),
                             source=r.source_v2))
D = pd.DataFrame(rows)
D.to_csv(R + "/tables/S29_flags_diff.csv", index=False)

print("v2 flags: %d rows x 3 = %d squares" % (len(F), 3 * len(F)))
print("half squares remaining:", int((F[["label_v2", "MOA_v2", "MS_v2"]] == 2).sum().sum()))
print("changed squares: %d\n" % len(D))
print(D[["gene", "flag", "old", "new", "kind"]].to_string(index=False))
print("\ntally v2:", F[["label_v2", "MOA_v2", "MS_v2"]].sum().to_dict(), "of", 3 * len(F))
