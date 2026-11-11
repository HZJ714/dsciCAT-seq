#!/usr/bin/env bash
# Build the two .h5ad files that Figure_6.ipynb reads as its *own data* inputs.
# Both were verified to reproduce the delivered files bit-for-bit (X / obs / var all identical).
#
#   (1) expression channel   3k_dsci_oligo_subset.RDS      -> 3k_dsci_oligo_subset.h5ad
#       5,977 nuclei x 31,263 genes, X = Seurat RNA `data` slot = log1p(counts/colsum*1e4)
#       (the Fig 2 lineage object; obs carries the six subtype labels used by panels A/B/C)
#   (2) chromatin channel    3k_dsci_ht_10x_gs_subset.RDS  -> genescore_oligo_5977.h5ad
#       5,977 nuclei x 24,333 genes, X = Signac `gene_score` data slot = log1p(ArchR gene score)
#       same nuclei, same order as (1); subtype labels joined from (1) by barcode
#
# Env overrides (defaults match the working tree):
#   DSCI_FY_ROOT        default $DSCI_WORK_ROOT/FY/oligo  (repository README)
#   DSCI_GS_RDS         default <online>/FY/3k_dsci_ht_10x_gs_subset.RDS   (13 GB)
#   DSCI_OWN_RNA_H5AD   default $DSCI_FY_ROOT/3k_dsci_oligo_subset.h5ad
#   DSCI_WORK           scratch directory for the dumps (default ./own_h5ad_build)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
FY="${DSCI_FY_ROOT:-${DSCI_WORK_ROOT:-.}/FY/oligo}"
WORK="${DSCI_WORK:-$PWD/own_h5ad_build}"
mkdir -p "$WORK"

# --- (1) RNA: RDS -> dump -> h5ad -------------------------------------------------------------
Rscript "$HERE/rds_to_dump.R" "$FY/3k_dsci_oligo_subset.RDS" "$WORK/rna_dump"
python3 "$HERE/build_h5ad.py" "$WORK/rna_dump" "$WORK/3k_dsci_oligo_subset.h5ad"

# --- (2) chromatin: RDS -> mtx dump -> h5ad ---------------------------------------------------
export DSCI_GS_PROBE_DIR="$WORK/gs_dump" DSCI_OWN_RNA_H5AD="$WORK/3k_dsci_oligo_subset.h5ad"
mkdir -p "$DSCI_GS_PROBE_DIR"
# the nucleus list is exactly the obs_names of the RNA object, in that order:
python3 - "$DSCI_OWN_RNA_H5AD" "$DSCI_GS_PROBE_DIR/lineage_cells_5977.txt" <<'PY'
import sys, anndata as ad
open(sys.argv[2], "w").write("\n".join(ad.read_h5ad(sys.argv[1], backed="r").obs_names) + "\n")
PY
Rscript "$HERE/export_genescore.R"
python3 "$HERE/build_genescore_h5ad.py" "$DSCI_GS_PROBE_DIR" "$WORK/genescore_oligo_5977.h5ad"

echo "done -> $WORK"
echo "panel C reads these two files; see README.md for the scale of each and for why the"
echo "panel-D cross-chemistry gene-score matrices must be log1p'd before scDRS."
