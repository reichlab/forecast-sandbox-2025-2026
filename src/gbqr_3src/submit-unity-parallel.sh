#!/bin/bash
#SBATCH -J gbqr3src
#SBATCH -N 1
#SBATCH -c 4                            # LightGBM threads; 8 threads stalled on busy nodes (job 65325494)
#SBATCH --mem=8G                        # job 65325494 MaxRSS ~0.7-1.0 GB
#SBATCH -p cpu
#SBATCH -t 08:00:00                     # job 65325494: mean 52 min, max 2h10m; 6 of 31 hit 3h with stalled bagging
#SBATCH --array=0-79%10                 # all 80 reference dates in ../../hub-config/tasks.json, at most 10 at once
#SBATCH -o logs/%x-%A_%a.out
#SBATCH -e logs/%x-%A_%a.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=trobacker@umass.edu
#
# gbqr_3src forecasts (all 53 locations, fit jointly) for every hub
# reference date, 2023-24 through 2025-26. Used for the national (US) GBQR
# component of the gbqr_3src_spatial x AR ensembles, as in
# flusion_spatial2_prod. Submit from src/gbqr_3src:
#   sbatch --array=0 submit-unity-parallel.sh   # smoke test, one date
#   sbatch submit-unity-parallel.sh             # all 80 dates
# Dates already in model-output/UMass-gbqr_3src/ are skipped (set
# GBQR_OVERWRITE=1 to refit them).
#
# Uses .venv built from requirements.txt (idmodels v2.1.0). Training data
# is downloaded from S3 at run time, so the compute node needs internet
# access.

set -euo pipefail

dates=(
  "2023-10-21"
  "2023-10-28"
  "2023-11-04"
  "2023-11-11"
  "2023-11-18"
  "2023-11-25"
  "2023-12-02"
  "2023-12-09"
  "2023-12-16"
  "2023-12-23"
  "2023-12-30"
  "2024-01-06"
  "2024-01-13"
  "2024-01-20"
  "2024-01-27"
  "2024-02-03"
  "2024-02-10"
  "2024-02-17"
  "2024-02-24"
  "2024-03-02"
  "2024-03-09"
  "2024-03-16"
  "2024-03-23"
  "2024-03-30"
  "2024-04-06"
  "2024-04-13"
  "2024-04-20"
  "2024-04-27"
  "2024-05-04"
  "2024-11-30"
  "2024-12-07"
  "2024-12-14"
  "2024-12-21"
  "2024-12-28"
  "2025-01-04"
  "2025-01-11"
  "2025-01-18"
  "2025-02-01"
  "2025-02-08"
  "2025-02-15"
  "2025-02-22"
  "2025-03-01"
  "2025-03-08"
  "2025-03-15"
  "2025-03-22"
  "2025-03-29"
  "2025-04-05"
  "2025-04-12"
  "2025-04-19"
  "2025-05-03"
  "2025-05-10"
  "2025-05-17"
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
out_csv="../../model-output/UMass-gbqr_3src/${date}-UMass-gbqr_3src.csv"

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

venv=".venv"
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
