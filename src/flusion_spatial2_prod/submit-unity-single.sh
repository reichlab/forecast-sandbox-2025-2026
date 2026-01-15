#!/bin/bash
#SBATCH -J flusion_spatial2_prod       # Job name
#SBATCH -N 1                           # Number of nodes
#SBATCH -c 12                          # Number of cores per task
#SBATCH --mem=64G                      # Memory per node
#SBATCH -p cpu                         # Partition name
#SBATCH -t 10:00:00                    # Time limit
#SBATCH -o logs/slurm-%j.out           # Output file (%j=job ID)
#SBATCH -e logs/slurm-%j.err           # Error file
#SBATCH --mail-type=FAIL               # Email on failure
#SBATCH --mail-user=nick@umass.edu

# Usage: sbatch submit-unity-single.sh YYYY-MM-DD
# Example: sbatch submit-unity-single.sh 2025-01-15

# Check if date argument is provided
if [ -z "$1" ]; then
    echo "Error: No date provided"
    echo "Usage: sbatch submit-unity-single.sh YYYY-MM-DD"
    exit 1
fi

date="$1"

# Validate date format (basic check)
if ! [[ "$date" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]]; then
    echo "Error: Invalid date format. Expected YYYY-MM-DD, got: $date"
    exit 1
fi

echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Running flusion_spatial2_prod ensemble forecast for date: $date"
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

# Load required modules on Unity (after venv activation)
R_MODULE="r/4.4.1-4oqw6bx"
if module load $R_MODULE 2>/dev/null; then
    echo "Loaded R module: $R_MODULE"
else
    echo "ERROR: Could not load R module: $R_MODULE"
    echo "Available R modules on this system:"
    module spider r 2>&1 | head -20
    exit 1
fi

# Explicitly export PATH to ensure R is available to subprocesses
export PATH

# Verify Python environment
echo "Python: $(which python)"
echo "Python version: $(python --version)"

# Verify R is available
echo "R: $(which Rscript)"
echo "R version: $(Rscript --version 2>&1 | head -1)"

# Run the ensemble forecast
python main.py --today_date="$date"

# Check exit status
if [ $? -eq 0 ]; then
    echo "Successfully completed ensemble forecast for $date"
else
    echo "Error: Ensemble forecast failed for $date"
    exit 1
fi

echo "End time: $(date)"
echo "=========================================="
