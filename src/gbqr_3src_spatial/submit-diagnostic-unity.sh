#!/bin/bash
#SBATCH -J gbqr_spatial_diagnostic     # Job name
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 8                            # Number of cores per task
#SBATCH --mem=32G                       # Memory per node
#SBATCH -p cpu                          # Partition name
#SBATCH -t 06:00:00                     # Time limit (6 hours - spatial features take longer)
#SBATCH --array=0-1                     # Array indices (2 dates)
#SBATCH -o logs/diagnostic-slurm-%A_%a.out  # Output file
#SBATCH -e logs/diagnostic-slurm-%A_%a.err  # Error file
#SBATCH --mail-type=FAIL,END            # Email on failure or completion
#SBATCH --mail-user=nick@umass.edu

# Target dates for diagnostic analysis
dates=(
  "2024-12-18"  # Early season date with high divergence
  "2025-02-12"  # Peak season date with maximum divergence
)

# Get the date for this array task
date="${dates[$SLURM_ARRAY_TASK_ID]}"

echo "=========================================="
echo "DIAGNOSTIC RUN - gbqr_3src_spatial"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Running diagnostic forecast for date: $date"
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

# Verify Python environment
echo "Python: $(which python)"
echo "Python version: $(python --version)"

# Create logs directory if it doesn't exist
mkdir -p logs

# Run the diagnostic forecast with artifact saving
python main_diagnostic.py \
    --today_date="$date" \
    --artifact_dir="artifacts"

# Check exit status
if [ $? -eq 0 ]; then
    echo "Successfully completed diagnostic forecast for $date"
    echo "Artifacts saved to: artifacts/$(date -d "$date" +%Y-%m-%d -I 2>/dev/null || date -j -f "%Y-%m-%d" "$date" +%Y-%m-%d)"
else
    echo "Error: Diagnostic forecast failed for $date"
    exit 1
fi

echo "End time: $(date)"
echo "=========================================="
