#!/bin/bash
#SBATCH --job-name=chronos_forecast
#SBATCH --partition=cpu              # Use 'cpu' partition (adjust based on your Unity allocation)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8            # Chronos benefits from multi-threading
#SBATCH --mem=16G                    # 16GB should be sufficient for small/base models
#SBATCH --time=00:30:00              # 30 minutes (adjust: base=15min, small=10min, tiny=5min)
#SBATCH --output=logs/chronos_%j.out
#SBATCH --error=logs/chronos_%j.err
#SBATCH --mail-type=END,FAIL         # Email on completion or failure
#SBATCH --mail-user=trobacker@umass.edu

# For GPU acceleration (optional - uncomment and adjust partition):
# #SBATCH --partition=gpu
# #SBATCH --gres=gpu:1
# #SBATCH --mem=32G

################################################################################
# Chronos Time Series Foundation Model - SLURM Job Script
#
# Usage:
#   sbatch run_chronos.slurm
#   sbatch --export=TODAY_DATE=2025-01-04,MODEL_SIZE=base run_chronos.slurm
#
# Environment variables:
#   TODAY_DATE   - Forecast date (YYYY-MM-DD), defaults to today
#   MODEL_SIZE   - Model size (tiny/mini/small/base/large), defaults to small
#   CONTEXT_LEN  - Historical weeks to use, defaults to 104 (2 years)
#   NUM_SAMPLES  - Samples for quantile estimation, defaults to 200
################################################################################

# Exit on error
set -e

# Default parameters (can be overridden via --export)
TODAY_DATE=${TODAY_DATE:-$(date +%Y-%m-%d)}
MODEL_SIZE=${MODEL_SIZE:-small}
CONTEXT_LEN=${CONTEXT_LEN:-104}
NUM_SAMPLES=${NUM_SAMPLES:-200}

echo "=================================================="
echo "Chronos Forecast Job"
echo "=================================================="
echo "Job ID:        $SLURM_JOB_ID"
echo "Node:          $SLURM_NODELIST"
echo "Date:          $(date)"
echo "Forecast Date: $TODAY_DATE"
echo "Model Size:    $MODEL_SIZE"
echo "Context:       $CONTEXT_LEN weeks"
echo "Samples:       $NUM_SAMPLES"
echo "CPUs:          $SLURM_CPUS_PER_TASK"
echo "Memory:        $SLURM_MEM_PER_NODE MB"
echo "=================================================="

# Create logs directory if it doesn't exist
mkdir -p logs

# Load any necessary modules (adjust based on Unity's available modules)
# Example: module load python/3.9 or module load anaconda
# module load python/3.11

# Activate virtual environment
# If .venv doesn't exist, create it first (see setup instructions below)
echo "Activating Python virtual environment..."
if [ ! -d ".venv" ]; then
    echo "ERROR: Virtual environment not found at .venv"
    echo "Please run setup first: python -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt"
    exit 1
fi

source .venv/bin/activate

# Verify packages are installed
if ! python -c "import chronos" 2>/dev/null; then
    echo "ERROR: chronos package not found. Installing dependencies..."
    pip install -r requirements.txt
fi

# Set number of threads for optimal CPU performance
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo ""
echo "Starting forecast generation..."
start_time=$(date +%s)

# Run the forecast
python main.py \
    --today_date="$TODAY_DATE" \
    --model_size="$MODEL_SIZE" \
    --context_length="$CONTEXT_LEN" \
    --num_samples="$NUM_SAMPLES"

end_time=$(date +%s)
elapsed=$((end_time - start_time))

echo ""
echo "=================================================="
echo "Job completed successfully!"
echo "Runtime: ${elapsed} seconds ($((elapsed / 60)) minutes)"
echo "Output: model-output/UMass-chronos_${MODEL_SIZE}/"
echo "=================================================="
