#!/bin/bash
#SBATCH -J flusion_spatial2_prod       # Job name
#SBATCH -N 1                           # Number of nodes
#SBATCH -c 8                           # Number of cores per task
#SBATCH --mem=32G                      # Memory per node
#SBATCH -p cpu                         # Partition name
#SBATCH -t 08:00:00                    # Time limit (2 hours)
#SBATCH --array=0-51                   # Array indices (55 dates total)
#SBATCH -o logs/slurm-%A_%a.out        # Output file (%A=job ID, %a=array index)
#SBATCH -e logs/slurm-%A_%a.err        # Error file
#SBATCH --mail-type=FAIL               # Email on failure or 80% time reached
#SBATCH --mail-user=nick@umass.edu

# Array of dates to process
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
  # "2024-05-08" # removing due to no new data
  "2024-11-27"
  "2024-12-04"
  "2024-12-11"
  "2024-12-18"
  "2024-12-25"
  "2025-01-01"
  "2025-01-08"
  "2025-01-15"
  # "2025-01-22" # removing due to data errors for this date
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
  # "2025-04-23" # missing data for this date
  "2025-04-30"
  "2025-05-07"
  "2025-05-14"
)

# Get the date for this array task
date="${dates[$SLURM_ARRAY_TASK_ID]}"

echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Running flusion_spatial ensemble forecast for date: $date"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "=========================================="

# Set up Python virtual environment
# Option 1: If using venv (uncomment and adjust path)
# source /path/to/your/venv/bin/activate

# Option 2: If using conda (uncomment and adjust environment name)
# module load miniconda3
# conda activate your_env_name

# Option 3: If using a local venv in the src directory
if [ -d ".venv" ]; then
    source .venv/bin/activate
else
    echo "Error: Virtual environment not found at .venv"
    exit 1
fi

# Load required modules on Unity (after venv activation)
# Try multiple R module versions in order of preference
# Unity has both EasyBuild (R/x.x.x-xxx) and Spack (r/x.x.x) modules
r_loaded=false
for r_module in r/4.4.0 R/4.3.2-gfbf-2023a R/4.2.1-foss-2022a R/4.2.0-foss-2021b R/4.0.0-foss-2020a R/4.4 R/4.3 R/4.2 R/4 R r; do
    if module load $r_module 2>/dev/null; then
        echo "Loaded R module: $r_module"
        r_loaded=true
        break
    fi
done

if [ "$r_loaded" = false ]; then
    echo "ERROR: Could not load any R module"
    echo "Available R modules on this system:"
    module spider R 2>&1 | grep "R/" | head -10
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
