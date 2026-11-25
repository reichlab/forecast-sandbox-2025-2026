#!/bin/bash
#SBATCH -J AR6_pooled_holiday           # Job name
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 8                            # Number of cores per task
#SBATCH --mem=32G                       # Memory per node
#SBATCH -p cpu                          # Partition name
#SBATCH -t 01:00:00                     # Time limit (1 hour)
#SBATCH --array=0-144                   # 5 variants x 29 dates
#SBATCH -o logs/slurm-%A_%a.out         # Output file (%A=job ID, %a=array index)
#SBATCH -e logs/slurm-%A_%a.err         # Error file
#SBATCH --mail-type=FAIL,TIME_LIMIT_80  # Email on failure or 80% time reached

# Submit from this folder, which needs .venv and logs/ (Slurm doesn't create logs/):
#   mkdir -p logs && sbatch submit-unity-parallel.sh
# 2023/24 reference dates only; later seasons are held out for now.
variants=(none shared partial xm1 xm1_ny1)
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
)

variant="${variants[$((SLURM_ARRAY_TASK_ID / ${#dates[@]}))]}"
date="${dates[$((SLURM_ARRAY_TASK_ID % ${#dates[@]}))]}"
echo "Job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID: variant=$variant date=$date node=$SLURM_NODELIST start=$(date)"

if [ -d ".venv" ]; then
    source .venv/bin/activate
else
    echo "Error: Virtual environment not found at .venv"
    exit 1
fi

python main.py --variant="$variant" --today_date="$date" || { echo "Error: forecast failed for $variant $date"; exit 1; }
echo "End time: $(date)"
