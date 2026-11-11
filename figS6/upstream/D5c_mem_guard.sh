#!/usr/bin/env bash
# D5c_mem_guard.sh - a co-running guardian: other sessions are running on this machine (currently 4 shards of 4.9G),
# we have a single cgroup (memory.max=256G, cpu.max=32 cores), and the OOM killer picks the highest-RSS process.
# This script does exactly one thing: when memory is tight, SIGSTOP OUR OWN scDRS process group, and SIGCONT it once pressure eases.
# Suspending is not failing: scDRS writes .score.gz per batch, and D5's resume skips completed batches.
#
# Usage: nohup bash D5c_mem_guard.sh > /dev/null 2>&1 &
#   PAUSE_AT_GB   suspend once cgroup memory.current exceeds this (default 200)
#   RESUME_AT_GB  resume only below this (default 160; hysteresis)
#   AVAIL_MIN_GB  also suspend when /proc/meminfo MemAvailable falls below this (default 60)
set -u
PAUSE_AT_GB=${PAUSE_AT_GB:-200}
RESUME_AT_GB=${RESUME_AT_GB:-160}
AVAIL_MIN_GB=${AVAIL_MIN_GB:-60}
LOG=${LOG:-${DSCI_WORK_ROOT:-.}/figures/routeA_out/D_crosschemistry/logs_log1p/memguard.log}
mkdir -p "$(dirname "$LOG")"
STOPPED=0
echo "$(date +%H:%M:%S) guard start limit=$(awk '{printf "%.0f",$1/1073741824}' /sys/fs/cgroup/memory.max)GB pause>${PAUSE_AT_GB}GB resume<${RESUME_AT_GB}GB" >> "$LOG"
while true; do
  CUR=$(awk '{printf "%.0f", $1/1048576/1024}' /sys/fs/cgroup/memory.current 2>/dev/null)
  AVA=$(awk '/MemAvailable/{printf "%.0f", $2/1048576}' /proc/meminfo)
  RSS=$(ps -eo rss,cmd | awk '/compute-score/ && !/awk/ {s+=$1} END{print int(s/1048576+0.5)}')
  NRUN=$(pgrep -fc "compute-score" 2>/dev/null || echo 0)
  if [ "$STOPPED" = "0" ] && { awk -v c="$CUR" -v l="$PAUSE_AT_GB" 'BEGIN{exit !(c>l)}' || awk -v a="$AVA" -v m="$AVAIL_MIN_GB" 'BEGIN{exit !(a<m)}'; }; then
    pkill -SIGSTOP -f "compute-score" 2>/dev/null
    STOPPED=1; echo "$(date +%H:%M:%S) PAUSE cur=${CUR}GB avail=${AVA}GB scdrs=${NRUN} rss=${RSS}GB" >> "$LOG"
  elif [ "$STOPPED" = "1" ] && awk -v c="$CUR" -v l="$RESUME_AT_GB" 'BEGIN{exit !(c<l)}'; then
    pkill -SIGCONT -f "compute-score" 2>/dev/null
    STOPPED=0; echo "$(date +%H:%M:%S) RESUME cur=${CUR}GB avail=${AVA}GB" >> "$LOG"
  fi
  echo "$(date +%H:%M:%S) cur=${CUR}GB avail=${AVA}GB scdrs=${NRUN} rss=${RSS}GB stopped=$STOPPED" >> "$LOG"
  # once we have finished, stand down - do not leave the guardian on the machine
  if [ "$NRUN" = "0" ] && ! pgrep -f "D5_run_scdrs.sh" >/dev/null 2>&1; then
    [ "$STOPPED" = "1" ] && pkill -SIGCONT -f "compute-score" 2>/dev/null
    echo "$(date +%H:%M:%S) guard exit (no scdrs left)" >> "$LOG"; break
  fi
  sleep 20
done
