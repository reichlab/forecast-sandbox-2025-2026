#!/bin/bash
# Shared per-task body for the WISAR6_fourthroot season array jobs (adaptive_t variant)
# (submit-unity-season-*.sh). Not submitted directly: each season script
# defines `dates=(...)` and then sources this file from src/WISAR6.
#
# Copied from ../WISAR6_fourthroot/wisar6-task.sh (same env vars, progress log and
# heartbeat); only the output path differs. chains == cores, single-threaded BLAS, a
# progress log with MCMC/stage markers and a CPU/RSS heartbeat.
#
# Extras for the full run:
#   - Skips a date whose output CSV already exists (set WISAR6_OVERWRITE=1
#     to refit it), so resubmitting a season only redoes missing dates.
#   - Gives each task its own PyTensor compile cache on node-local disk.
#     The test tasks ran one at a time; here up to 10 run at once, and a
#     shared ~/.pytensor cache on /home would mean lock contention and
#     metadata traffic. Compiling costs well under a minute per task.

set -euo pipefail

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export WISAR6_MCMC_CHAINS="${SLURM_CPUS_PER_TASK:-1}"
export WISAR6_MCMC_CORES="${SLURM_CPUS_PER_TASK:-1}"
export WISAR6_PROGRESS_EVERY=250
HEARTBEAT_SECS=300

date="${dates[$SLURM_ARRAY_TASK_ID]}"
out_csv="../../model-output/UMass-WISAR6_fourthroot_adaptive_t/${date}-UMass-WISAR6_fourthroot_adaptive_t.csv"
progress_log="logs/${SLURM_JOB_NAME}-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.progress"
export WISAR6_PROGRESS_LOG="$progress_log"

echo "=========================================="
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Running forecast for date: $date"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}  Memory: ${SLURM_MEM_PER_NODE:-?}M"
echo "Start time: $(date)"
echo "Progress log: $progress_log"
echo "=========================================="

if [ -s "$out_csv" ] && [ "${WISAR6_OVERWRITE:-0}" != "1" ]; then
    echo "Output already exists, skipping: $out_csv"
    exit 0
fi

module list 2>&1 || true

if [ ! -d ".venv" ]; then
    echo "Error: .venv not found. Run 'uv sync' once from a login node first."
    exit 1
fi
# Activate the synced env directly rather than `uv run`, which may try to
# re-sync .venv and would race across array tasks.
source .venv/bin/activate
echo "Python: $(which python) ($(python --version))"

compiledir="${TMPDIR:-/tmp}/${USER}-pytensor-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$compiledir"
trap 'rm -rf "$compiledir"' EXIT
export PYTENSOR_FLAGS="compiledir=${compiledir}"
# ArviZ's once-a-day FutureWarning writes a stamp file under the user cache
# dir (~/.cache/arviz) via a fixed temp name, so tasks starting in the same
# second on the first run after midnight race on it and crash at import
# (job 65229735 tasks 6, 8, 9). Give each task its own cache dir instead.
export XDG_CACHE_HOME="${compiledir}/xdg-cache"

echo "[$(date '+%F %T')] job started on $(hostname) for $date" > "$progress_log"

python main.py --today_date="$date" &
py_pid=$!

# Heartbeat: CPU% and RSS for the python process plus its chain subprocesses.
(
  while kill -0 "$py_pid" 2>/dev/null; do
    sleep "$HEARTBEAT_SECS"
    procs=$(ps --no-headers -o pid=,pcpu=,rss= -p "$py_pid" --ppid "$py_pid" 2>/dev/null || true)
    [ -z "$procs" ] && break
    summary=$(echo "$procs" | awk '{n++; cpu+=$2; rss+=$3} END {printf "%d procs, %.0f%% CPU total (lifetime avg per ps), %.2f GB RSS total", n, cpu, rss/1048576}')
    echo "[$(date '+%F %T')] heartbeat: $summary" >> "$progress_log"
  done
) &
hb_pid=$!

status=0
wait "$py_pid" || status=$?
kill "$hb_pid" 2>/dev/null || true

if [ "$status" -eq 0 ]; then
    echo "Successfully completed forecast for $date"
    echo "[$(date '+%F %T')] job finished OK" >> "$progress_log"
else
    echo "Error: Forecast failed for $date (exit $status)"
    echo "[$(date '+%F %T')] job FAILED (exit $status)" >> "$progress_log"
fi

echo "End time: $(date)"
echo "=========================================="
exit "$status"
