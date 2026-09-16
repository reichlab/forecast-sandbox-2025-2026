#!/bin/bash
#SBATCH -J gbqr_5src_mobsgleam_otid10             # Job name
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 8                            # Number of cores per task
#SBATCH --mem=32G                       # Memory per node
#SBATCH -p cpu                          # Partition name
#SBATCH -t 01:00:00                     # Time limit (1 hours)
#SBATCH --array=0-27                    # Array indices (28 dates total; NSSP requires as_of >= 2025-09-17)
#SBATCH -o logs/slurm-%A_%a.out         # Output file (%A=job ID, %a=array index)
#SBATCH -e logs/slurm-%A_%a.err         # Error file
#SBATCH --mail-type=FAIL,TIME_LIMIT_80  # Email on failure or 80% time reached
#SBATCH --mail-user=lshandross@umass.edu

# Array of dates to process
dates=(
  "2025-11-19"
  "2025-11-26"
  "2025-12-03"
  "2025-12-10"
  "2025-12-17"
  "2025-12-24"
  "2025-12-31"
  "2026-01-07"
  "2026-01-14"
  "2026-01-21"
  "2026-01-28"
  "2026-02-04"
  "2026-02-11"
  "2026-02-18"
  "2026-02-25"
  "2026-03-04"
  "2026-03-11"
  "2026-03-18"
  "2026-03-25"
  "2026-04-01"
  "2026-04-08"
  "2026-04-15"
  "2026-04-22"
  "2026-04-29"
  "2026-05-06"
  "2026-05-13"
  "2026-05-20"
  "2026-05-27"
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
