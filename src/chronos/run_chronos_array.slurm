#!/bin/bash
#SBATCH --job-name=chronos_array
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --array=0-9                  # Array indices (adjust based on number of dates)
#SBATCH --output=logs/chronos_array_%A_%a.out
#SBATCH --error=logs/chronos_array_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=trobacker@umass.edu

################################################################################
# Chronos Array Job - Run Multiple Forecast Dates in Parallel
#
# Usage:
#   1. Edit the DATES array below with your forecast dates
#   2. Update --array parameter to match number of dates (e.g., 0-9 for 10 dates)
#   3. Submit: sbatch run_chronos_array.slurm
#
# Each array task runs one forecast date independently.
################################################################################

set -e

# Define forecast dates (one per array task)
# Edit this array with your desired forecast dates
DATES=(
    "2024-11-30"
    "2024-12-07"
    "2024-12-14"
    "2024-12-21"
    "2024-12-28"
    "2025-01-04"
    "2025-01-11"
    "2025-01-18"
    "2025-01-25"
    "2025-02-01"
)

# Configuration
MODEL_SIZE=${MODEL_SIZE:-small}
CONTEXT_LEN=${CONTEXT_LEN:-104}
NUM_SAMPLES=${NUM_SAMPLES:-200}

# Get the date for this array task
TODAY_DATE=${DATES[$SLURM_ARRAY_TASK_ID]}

echo "=================================================="
echo "Chronos Array Forecast Job"
echo "=================================================="
echo "Array Job ID:  $SLURM_ARRAY_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Job ID:        $SLURM_JOB_ID"
echo "Node:          $SLURM_NODELIST"
echo "Date:          $(date)"
echo "Forecast Date: $TODAY_DATE"
echo "Model Size:    $MODEL_SIZE"
echo "=================================================="

# Create logs directory
mkdir -p logs

# Activate virtual environment
source .venv/bin/activate

# Set threading
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Starting forecast for $TODAY_DATE..."
start_time=$(date +%s)

# Run forecast
python main.py \
    --today_date="$TODAY_DATE" \
    --model_size="$MODEL_SIZE" \
    --context_length="$CONTEXT_LEN" \
    --num_samples="$NUM_SAMPLES"

end_time=$(date +%s)
elapsed=$((end_time - start_time))

echo ""
echo "=================================================="
echo "Array task $SLURM_ARRAY_TASK_ID completed!"
echo "Date: $TODAY_DATE"
echo "Runtime: ${elapsed} seconds ($((elapsed / 60)) minutes)"
echo "=================================================="
