#!/bin/bash
#SBATCH -J wisar6frat_2024_25
#SBATCH -N 1
#SBATCH -c 4                            # One core per MCMC chain (test: ~100% CPU per chain process)
#SBATCH --mem=4G                        # Test MaxRSS was 0.48-0.78 GB; 4G leaves headroom
#SBATCH -p cpu
#SBATCH -t 03:00:00                     # Student-t sampling cost unknown; WISAR6_fourthroot took up to 47 min
#SBATCH --array=0-22%10                # 23 reference dates in the 2024-25 season, at most 10 running at once
#SBATCH -o logs/%x-%A_%a.out
#SBATCH -e logs/%x-%A_%a.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=trobacker@umass.edu
#
# WISAR6_fourthroot_adaptive_t forecasts for the 2024-25 season. Submit from src/WISAR6_fourthroot_adaptive_t, either on
# its own (sbatch submit-unity-season-2024-25.sh) or via submit-unity-all-seasons.sh.
# Dates already in model-output/UMass-WISAR6_fourthroot_adaptive_t/ are skipped.
#
# Sizing (4 CPUs, 4G) copied from ../WISAR6_fourthroot. Time limits are ~3-4x the
# slowest WISAR6_fourthroot task (jobs 65186564-6: 31 / 47 / 69 min) since the
# Student-t WIS loss (betaincinv quantiles) is likely slower to sample.
#
# Dates are this season's valid reference dates in ../../hub-config/tasks.json.

dates=(
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
)

cd "$SLURM_SUBMIT_DIR"
source ./wisar6-task.sh
