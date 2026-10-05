#!/bin/bash
#SBATCH -J wisar6_fourthroot_adaptive_t_forecast   # Job name
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 8                            # Number of cores per task
#SBATCH --mem=32G                       # Memory per node
#SBATCH -p cpu                          # Partition name
#SBATCH -t 05:00:00                     # Time limit. The sibling log1p
                                         #   WISAR6 model measured ~22 min
                                         #   wall-clock per date on a laptop
                                         #   (8 cores, cores=1/sequential
                                         #   chains, 4th-root transform --
                                         #   that was in fact this model's
                                         #   transform; see model.py). This
                                         #   fourth-root variant's own Unity
                                         #   timing is not yet measured; 5
                                         #   hours is deliberately generous
                                         #   until a real measurement exists
                                         #   -- tighten this once you have one.
#SBATCH --array=0-83                    # Array indices (84 dates total)
#SBATCH -o logs/slurm-%A_%a.out         # Output file (%A=job ID, %a=array index)
#SBATCH -e logs/slurm-%A_%a.err         # Error file
#SBATCH --mail-type=FAIL,TIME_LIMIT_80  # Email on failure or 80% time reached
#SBATCH --mail-user=trobacker@umass.edu # UPDATE to your own address if needed

# Adapted from ../WISAR6/submit-unity-parallel.sh (same model, log1p
# transform instead of fourth-root), which was itself adapted from
# ../gbqr/submit-unity-parallel.sh, updated to use `uv` for environment
# management instead of a manually-created .venv + pip install.
#
# One-time setup before the first `sbatch submit-unity-parallel.sh` (run on
# a Unity login node, from this directory):
#   module load uv        # or: curl -LsSf https://astral.sh/uv/install.sh | sh
#   uv sync                # builds .venv from pyproject.toml/uv.lock
#   mkdir -p logs

# Run PyMC's 2 MCMC chains in parallel OS processes instead of main.py's
# default sequential cores=1. The sequential default exists because of a
# multiprocessing EOFError seen on the macOS development machine (see
# model.fit_model's docstring) -- that failure mode is specific to macOS's
# process-start behavior and has not been observed on Linux, so it's worth
# trying here for the speedup. That said, this has NOT yet been validated
# at full 53-location production scale on Unity itself: watch the first
# array job's logs for a crash (grep for "EOFError" in logs/slurm-*.err)
# before trusting this unattended. If it does crash, remove this line (or
# set it to 1) and resubmit the failed array indices.
export WISAR6_MCMC_CORES=2

# Array of dates to process -- all 84 actual FluSight round (reference)
# dates across the 3 complete historical seasons present in the hub's
# target data as of 2026-09-30 (2023-24, 2024-25, 2025-26), verified
# directly against model-output/FluSight-baseline/'s own submitted
# reference dates (the authoritative source for "which Saturdays were
# real round dates") rather than assumed from a fixed weekly cadence --
# the season boundaries have real gaps (e.g. no submissions between
# 2024-05-04 and 2024-11-23) that a naive date-range generator would get
# wrong. Each date is passed directly as the actual Saturday reference
# date (main.py's --today_date rounds forward to the nearest Saturday via
# relativedelta(weekday=5), so a Saturday maps to itself). Also verified
# that every date here has an at-or-before snapshot in the hub's
# auxiliary-data/target-data-archive/, so main.py's point-in-time-correct
# historical lookup (see README.md's "Historical vs. live data" section)
# works for all of them without silently falling back to the live
# (potentially leaky) target-data file.
#
# Keep this list in sync with ../WISAR6/submit-unity-parallel.sh's (same
# date list, different transform).
dates=(
  "2023-10-14"
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
  "2024-11-23"
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
  "2025-04-26"
  "2025-05-03"
  "2025-05-10"
  "2025-05-17"
  "2025-05-24"
  "2025-05-31"
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

# Verify uv is available and the environment is already synced (uv sync
# should be run once, interactively, before submitting the array job --
# doing it here too would race across array tasks trying to write the same
# .venv simultaneously).
if ! command -v uv &> /dev/null; then
    echo "Error: 'uv' not found on PATH. Run 'module load uv' in your Unity"
    echo "environment or install it (see README.md), then resubmit."
    exit 1
fi

if [ ! -d ".venv" ]; then
    echo "Error: .venv not found. Run 'uv sync' once from a login node before"
    echo "submitting this array job."
    exit 1
fi

echo "uv: $(which uv)"
echo "Python: $(uv run python --version)"

# Run the forecast
uv run python main.py --today_date="$date"

# Check exit status
if [ $? -eq 0 ]; then
    echo "Successfully completed forecast for $date"
else
    echo "Error: Forecast failed for $date"
    exit 1
fi

echo "End time: $(date)"
echo "=========================================="
