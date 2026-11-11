#!/usr/bin/env bash
# D5: 3 chemistries x 2 channels x 74 traits = 444 scDRS runs.
# Parameters are copied verbatim from own-C's run_full74_genescore.sh (same runner, same gene set, same n_ctrl, same cov column count),
# only the h5ad changes for D - so cross-chemistry rho differences cannot come from "different parameters".
#
# Usage: D5_run_scdrs.sh [smoke|full] ["tag_ch list"]
#
# 2026-11-04 additions (input-transform defect, see upstream/README.md "Input transform"):
#   * --flag-raw-count False makes scDRS expm1() the input, so the matrix must be on a log scale.
#     A preflight gate now refuses to launch a job whose matrix violates that promise, and every
#     batch is post-checked for the numpy overflow warnings that this defect produced. Both are
#     hard failures unless explicitly skipped (SKIP_PREFLIGHT=1 / SKIP_POSTCHECK=1).
#   * H5 / SCDRS_OUT are overridable so a corrected rerun can live in a parallel tree without
#     touching the delivered one.
#   * NICE / IONICE keep CPU+IO away from other jobs on the box (this run is not the only tenant).
set -u
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SCDRS=${SCDRS:-$HOME/.local/bin/scdrs}
PY=${PY:-python3}
GSDIR=${GSDIR:-${DSCI_FY_ROOT:-${DSCI_WORK_ROOT:-.}/FY}/oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_1000.74_traits.gs.batch}
OUT=${OUT:-${DSCI_FIG_ROOT:-${DSCI_WORK_ROOT:-.}/figures}/routeA_out/D_crosschemistry}
H5=${H5:-$OUT/h5ad}
SCDRS_OUT=${SCDRS_OUT:-$OUT/scdrs}
LOG=${LOG:-$OUT/logs}
mkdir -p "$LOG"
MODE=${1:-full}
JOBS=${2:-"10x_rna 10x_atac HT_scCAT_seq_HT_rna HT_scCAT_seq_HT_atac dsciCAT_seq_rna dsciCAT_seq_atac"}
NP=${NP:-8}
NICE=${NICE:-10}                    # 0-19; 10 keeps CPU for the other tenant on this box
IONICE=${IONICE:-2}                 # 2 = best-effort
IONICE_N=${IONICE_N:-7}             # 0-7 within best-effort; 0 disables ionice
PRIO="nice -n $NICE"
[ "$IONICE_N" != "0" ] && PRIO="$PRIO ionice -c $IONICE -n $IONICE_N"

run_one() {            # $1=h5ad $2=cov $3=outdir $4=gsfile
  local h5ad="$1" cov="$2" out="$3" gs="$4"
  local d="$out/$(basename "$gs" .gs)"; mkdir -p "$d"
  if ls "$d"/*.score.gz >/dev/null 2>&1; then echo "skip $(basename "$gs")"; return; fi
  $PRIO "$SCDRS" compute-score --h5ad-file "$h5ad" --h5ad-species mouse \
      --gs-file "$gs" --gs-species human --out-folder "$d" \
      --flag-filter-data True --flag-raw-count False --n-ctrl 1000 \
      --flag-return-ctrl-raw-score False --flag-return-ctrl-norm-score True \
      --cov-file "$cov" > "$d/compute.log" 2>&1 \
    && echo "done $(basename "$gs") $(date +%H:%M:%S)" || echo "FAIL $(basename "$gs") $(date +%H:%M:%S)"
}
export -f run_one; export SCDRS GSDIR OUT H5 SCDRS_OUT LOG PRIO

for job in $JOBS; do
  H5AD="$H5/${job}.h5ad"; COV="$H5/${job}.cov.tsv"; O="$SCDRS_OUT/$job"; mkdir -p "$O"
  if [ ! -s "$H5AD" ] || [ ! -s "$COV" ]; then echo "!! missing input: $H5AD / $COV"; continue; fi

  # --- preflight: is ${job}.h5ad the scale that --flag-raw-count False promises? ---
  if [ "${SKIP_PREFLIGHT:-0}" != "1" ]; then
    if ! "$PY" "$HERE/preflight_h5ad.py" "$H5AD"; then
      echo "!! PREFLIGHT FAILED: $job -> not starting. Run D4b_pack_atac_log1p.py first, or switch to --flag-raw-count True."
      continue
    fi
  fi

  if [ "$MODE" = smoke ]; then GS=$(ls "$GSDIR"/*.gs | head -1); echo "== SMOKE $job using $(basename $GS)"; run_one "$H5AD" "$COV" "$O" "$GS"; continue; fi
  echo "== $job start $(date +%H:%M:%S) | $(ls -la $H5AD | awk '{print $5}') B | NP=$NP"
  ls "$GSDIR"/*.gs | xargs -P "$NP" -I{} bash -c 'run_one "$@"' _ "$H5AD" "$COV" "$O" {} \
     > "$LOG/D5_${job}.log" 2>&1
  echo "== $job end $(date +%H:%M:%S) | produced $(find $O -name '*.score.gz' | wc -l) score.gz files"

  # --- postcheck: scDRS warns instead of failing; a warning here means the scores are unusable ---
  NG=$(ls "$GSDIR"/*.gs | wc -l)                       # batch .gs count (25)
  NEXP=$(cat "$GSDIR"/*.gs | awk -F'\t' 'NR>1 && $1!="" {print $1}' | sort -u | wc -l)   # total traits (74)
  NS=$(find "$O" -name '*.score.gz' 2>/dev/null | wc -l)   # one score.gz per trait (full_score.gz does not match *.score.gz)
  NBDONE=0; for bd in "$O"/magma_10kb_1000.batch*; do [ -d "$bd" ] && ls "$bd"/*.score.gz >/dev/null 2>&1 && NBDONE=$((NBDONE+1)); done
  NF=$(grep -c '^FAIL' "$LOG/D5_${job}.log" 2>/dev/null || true)
  echo "   batch .gs $NG | batches done $NBDONE | trait score.gz $NS/$NEXP | FAIL lines ${NF:-0}"
  if [ "$NS" -ne "$NEXP" ] || [ "$NBDONE" -ne "$NG" ] || [ "${NF:-0}" -ne 0 ]; then
      echo "!! POSTCHECK FAILED: $job output incomplete ($NS/$NEXP traits, $NBDONE/$NG batches). See $LOG/D5_${job}.log and the compute.log files"
  fi
  if [ "${SKIP_POSTCHECK:-0}" != "1" ]; then
    N=$(grep -ril "overflow encountered" "$O" --include=compute.log 2>/dev/null | wc -l)
    T=$(find "$O" -name compute.log 2>/dev/null | wc -l)
    if [ "$N" -gt 0 ]; then
      echo "!! POSTCHECK FAILED: $job has $N/$T batch logs with 'overflow encountered' (expm1 or its following multiplication)"
    else
      echo "   postcheck OK: $T batch logs without overflow warnings"
    fi
  fi
done
echo "ALL_FINISHED $(date +%H:%M:%S)"
