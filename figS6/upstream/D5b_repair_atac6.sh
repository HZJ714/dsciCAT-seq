#!/usr/bin/env bash
# repair the 6 traits missing from dsciCAT_seq_atac (run one .gs at a time, dodging D5's "skip if the dir has files" bug)
set -u
SCDRS=${SCDRS:-$HOME/.local/bin/scdrs}
GSALL=${DSCI_WORK_ROOT:-.}/FY/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs
H5AD=${DSCI_WORK_ROOT:-.}/figures/routeA_out/D_crosschemistry/h5ad/dsciCAT_seq_atac.h5ad
COV=${DSCI_WORK_ROOT:-.}/figures/routeA_out/D_crosschemistry/h5ad/dsciCAT_seq_atac.cov.tsv
O=${DSCI_WORK_ROOT:-.}/figures/routeA_out/D_crosschemistry/scdrs/dsciCAT_seq_atac/repair
mkdir -p "$O"
MISSING="PASS_DrinksPerWeek_Liu2019 PASS_FastingGlucose_Manning PASS_Primary_biliary_cirrhosis PASS_ReactionTime_Davies2018 PASS_Type_1_Diabetes PASS_Type_2_Diabetes"
for t in $MISSING; do
  d="$O/$t"; mkdir -p "$d"
  [ -f "$d/$t.score.gz" ] && { echo "skip $t"; continue; }
  { echo -e "TRAIT\tGENESET"; grep -P "^\Q$t\E\t" "$GSALL"; } > "$O/${t}.gs"
  "$SCDRS" compute-score --h5ad-file "$H5AD" --h5ad-species mouse --gs-file "$O/${t}.gs" \
      --gs-species human --out-folder "$d" --flag-filter-data True --flag-raw-count False \
      --n-ctrl 1000 --flag-return-ctrl-raw-score False --flag-return-ctrl-norm-score True \
      --cov-file "$COV" > "$d/compute.log" 2>&1 \
    && echo "done $t $(date +%H:%M:%S)" || echo "FAIL $t $(date +%H:%M:%S)"
done
echo REPAIR_DONE $(date +%H:%M:%S)
