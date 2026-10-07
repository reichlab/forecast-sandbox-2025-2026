#!/bin/bash
#SBATCH -J gbqr3ss_2025_26
#SBATCH -N 1
#SBATCH -c 4                            # LightGBM threads; 8 threads stalled on busy nodes (job 65325493)
#SBATCH --mem=8G                         # v2.1 smoke test (job 65285718) MaxRSS 1.3 GB
#SBATCH -p cpu
#SBATCH -t 12:00:00                     # job 65325493: mean 3h22m, max 3h50m; 3 of 12 hit 6h with stalled bagging
#SBATCH --array=0-27%10                 # 28 reference dates in the 2025-26 season, at most 10 running at once
#SBATCH -o logs/%x-%A_%a.out
#SBATCH -e logs/%x-%A_%a.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=trobacker@umass.edu
#
# gbqr_3src_spatial forecasts for the 2025-26 season. Submit from
# src/gbqr_3src_spatial after `mkdir -p logs`:
#   sbatch submit-unity-season-2025-26.sh            # all 28 dates
#   sbatch --array=0 submit-unity-season-2025-26.sh  # smoke test, one date
# Dates already in model-output/UMass-gbqr_3src_spatial/ are skipped (set
# GBQR_OVERWRITE=1 to refit them), so resubmitting only redoes missing dates.
#
# Uses the idmodels v2.1.0 venv in ../gbqr_3src (same pin as requirements.txt
# here; ./.venv is the old pre-v2.1 install without wave features). Training
# data is downloaded from S3 at run time, so the compute node needs internet
# access.
#
# Dates are this season's valid reference dates in ../../hub-config/tasks.json
# (Saturdays; main.py rounds --today_date forward to Saturday, so each maps to
# itself).

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

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
# don't let idle OpenMP threads spin and compete for the allocated cores
export OMP_WAIT_POLICY=PASSIVE

date="${dates[$SLURM_ARRAY_TASK_ID]}"
out_csv="../../model-output/UMass-gbqr_3src_spatial/${date}-UMass-gbqr_3src_spatial.csv"

echo "=========================================="
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Running forecast for date: $date"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-?}  Memory: ${SLURM_MEM_PER_NODE:-?}M"
echo "Start time: $(date)"
echo "=========================================="

if [ -s "$out_csv" ] && [ "${GBQR_OVERWRITE:-0}" != "1" ]; then
    echo "Output already exists, skipping: $out_csv"
    exit 0
fi

module list 2>&1 || true

venv="../gbqr_3src/.venv"  # idmodels v2.1.0, same pin as requirements.txt here
if [ ! -d "$venv" ]; then
    echo "Error: Virtual environment not found at $venv"
    exit 1
fi
source "$venv/bin/activate"
echo "Python: $(which python) ($(python --version))"

status=0
python main.py --today_date="$date" || status=$?

if [ "$status" -eq 0 ]; then
    echo "Successfully completed forecast for $date"
else
    echo "Error: Forecast failed for $date (exit $status)"
fi

echo "End time: $(date)"
echo "=========================================="
exit "$status"
