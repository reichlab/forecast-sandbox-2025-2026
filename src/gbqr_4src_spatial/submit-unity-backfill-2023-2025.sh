#!/bin/bash
#SBATCH -J gbqr_4src_spatial_backfill   # Job name
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 8                            # Number of cores per task
#SBATCH --mem=32G                       # Memory per node (same budget as the current-season run)
#SBATCH -p cpu                          # Partition name
#SBATCH -t 02:30:00                     # Time limit (2.5 hours; per-bag time varies a lot under node contention)
#SBATCH --array=0-51%4                  # Array indices (52 dates total; throttled to 4 concurrent to avoid node contention)
#SBATCH -o logs/slurm-%A_%a.out         # Output file (%A=job ID, %a=array index)
#SBATCH -e logs/slurm-%A_%a.err         # Error file
#SBATCH --mail-type=FAIL,TIME_LIMIT_80  # Email on failure or 80% time reached
#SBATCH --mail-user=lshandross@umass.edu

# Array of dates to process (2023-24 and 2024-25 seasons; current 2025-26 season is
# handled by submit-unity-parallel.sh)
dates=(
  "2023-10-18"
  "2023-10-25"
  "2023-11-01"
  "2023-11-08"
  "2023-11-15"
  "2023-11-22"
  "2023-11-29"
  "2023-12-06"
  "2023-12-13"
  "2023-12-20"
  "2023-12-27"
  "2024-01-03"
  "2024-01-10"
  "2024-01-17"
  "2024-01-24"
  "2024-01-31"
  "2024-02-07"
  "2024-02-14"
  "2024-02-21"
  "2024-02-28"
  "2024-03-06"
  "2024-03-13"
  "2024-03-20"
  "2024-03-27"
  "2024-04-03"
  "2024-04-10"
  "2024-04-17"
  "2024-04-24"
  "2024-05-01"
  "2024-11-27"
  "2024-12-04"
  "2024-12-11"
  "2024-12-18"
  "2024-12-25"
  "2025-01-01"
  "2025-01-08"
  "2025-01-15"
  "2025-01-29"
  "2025-02-05"
  "2025-02-12"
  "2025-02-19"
  "2025-02-26"
  "2025-03-05"
  "2025-03-12"
  "2025-03-19"
  "2025-03-26"
  "2025-04-02"
  "2025-04-09"
  "2025-04-16"
  "2025-04-30"
  "2025-05-07"
  "2025-05-14"
)

# Get the date for this array task
date="${dates[$SLURM_ARRAY_TASK_ID]}"

echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Running forecast for date: $date"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "=========================================="

# Set up Python virtual environment
if [ -d ".venv" ]; then
    source .venv/bin/activate
else
    echo "Error: Virtual environment not found at .venv"
    exit 1
fi

echo "Python: $(which python)"
echo "Python version: $(python --version)"

python main.py --today_date="$date"

if [ $? -eq 0 ]; then
    echo "Successfully completed forecast for $date"
else
    echo "Error: Forecast failed for $date"
    exit 1
fi

echo "End time: $(date)"
echo "=========================================="
