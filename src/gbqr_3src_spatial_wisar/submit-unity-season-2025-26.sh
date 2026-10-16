#!/bin/bash
#SBATCH -J gbqr3ssw_2025_26
#SBATCH -N 1
#SBATCH -c 12                           # gbqr/gbqr_us/nssp run concurrently (main.py caps each to c/3
                                         # threads) alongside wisar's own MCMC chains; see main.py
#SBATCH --mem=16G                       # 4 model processes resident at once instead of 1
#SBATCH -p cpu
#SBATCH -t 6:00:00                      # gbqr_3src_spatial alone took up to 3h50m serially; the 4
                                         # components now run concurrently, so wall time should track
                                         # the slowest single component rather than their sum
#SBATCH --array=0-27%10                 # 28 reference dates in the 2025-26 season, at most 10 running at once
#SBATCH -o logs/%x-%A_%a.out
#SBATCH -e logs/%x-%A_%a.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=trobacker@umass.edu
#
# gbqr_3src_spatial_wisar(_nssp) forecasts for the 2025-26 season. Submit
# from src/gbqr_3src_spatial_wisar after `mkdir -p logs` and the one-time
# environment setup in README.md:
#   sbatch submit-unity-season-2025-26.sh            # all 28 dates
#   sbatch --array=0 submit-unity-season-2025-26.sh  # smoke test, one date
#
# Dates already in model-output/UMass-gbqr_3src_spatial_wisar_nssp/ (the
# combined submission file) are skipped (set ENSEMBLE_OVERWRITE=1 to rerun
# all four components and reensemble), so resubmitting only redoes missing
# dates.
#
# Dates are this season's valid reference dates in ../../hub-config/tasks.json
# (Saturdays; main.py rounds --today_date forward to Saturday, so each maps
# to itself). gbqr_4src_nssp (the ED-visits component) additionally skips
# 2025-11-29, 2025-12-27, and 2026-01-03 (holiday reporting lag); main.py
# publishes the hosp ensemble alone under the combined name for those dates.

set -euo pipefail

dates=(
  "2025-11-22"
  "2025-11-29"
  "2025-12-06"
  "2025-12-13"
  "2025-12-20"
  "2025-12-27"
  "2026-01-03"
  "2026-01-10"
  "2026-01-17"
  "2026-01-24"
  "2026-01-31"
  "2026-02-07"
  "2026-02-14"
  "2026-02-21"
  "2026-02-28"
  "2026-03-07"
  "2026-03-14"
  "2026-03-21"
  "2026-03-28"
  "2026-04-04"
  "2026-04-11"
  "2026-04-18"
  "2026-04-25"
  "2026-05-02"
  "2026-05-09"
  "2026-05-16"
  "2026-05-23"
  "2026-05-30"
)

cd "$SLURM_SUBMIT_DIR"

# main.py sets OMP/MKL/OPENBLAS_NUM_THREADS itself for the gbqr/gbqr_us/nssp
# subprocesses (splitting $SLURM_CPUS_PER_TASK three ways; see its
# _gbqr_step_env()), so don't export them here too. wisar's chain count does
# need to come from the environment, sized to leave room for the other three
# components running at the same time.
export OMP_WAIT_POLICY=PASSIVE
export WISAR6_MCMC_CHAINS=2
export WISAR6_MCMC_CORES=2

date="${dates[$SLURM_ARRAY_TASK_ID]}"
out_csv="../../model-output/UMass-gbqr_3src_spatial_wisar_nssp/${date}-UMass-gbqr_3src_spatial_wisar_nssp.csv"

echo "=========================================="
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Running gbqr_3src_spatial_wisar(_nssp) forecast for date: $date"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-?}  Memory: ${SLURM_MEM_PER_NODE:-?}M"
echo "Start time: $(date)"
echo "=========================================="

if [ -s "$out_csv" ] && [ "${ENSEMBLE_OVERWRITE:-0}" != "1" ]; then
    echo "Output already exists, skipping: $out_csv"
    exit 0
fi

module list 2>&1 || true

R_MODULE="r/4.4.1-4oqw6bx"
if module load $R_MODULE 2>/dev/null; then
    echo "Loaded R module: $R_MODULE"
else
    echo "ERROR: Could not load R module: $R_MODULE"
    module spider r 2>&1 | head -20
    exit 1
fi
export PATH
echo "R: $(which Rscript) ($(Rscript --version 2>&1 | head -1))"

if [ ! -d ".venv-gbqr" ]; then
    echo "Error: .venv-gbqr not found. See README.md for one-time setup."
    exit 1
fi
if [ ! -d "wisar/.venv" ]; then
    echo "Error: wisar/.venv not found. Run 'uv sync' from wisar/ first (see README.md)."
    exit 1
fi
if [ ! -d ".venv-nssp" ]; then
    echo "Error: .venv-nssp not found. See README.md for one-time setup."
    exit 1
fi

status=0
.venv-gbqr/bin/python main.py --today_date="$date" || status=$?

if [ "$status" -eq 0 ]; then
    echo "Successfully completed forecast for $date"
else
    echo "Error: Forecast failed for $date (exit $status)"
fi

echo "End time: $(date)"
echo "=========================================="
exit "$status"
