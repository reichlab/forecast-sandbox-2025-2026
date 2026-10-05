#!/bin/bash
#SBATCH -J wisar6frat_2023_24
#SBATCH -N 1
#SBATCH -c 4                            # One core per MCMC chain (test: ~100% CPU per chain process)
#SBATCH --mem=4G                        # Test MaxRSS was 0.48-0.78 GB; 4G leaves headroom
#SBATCH -p cpu
#SBATCH -t 02:00:00                     # Student-t sampling cost unknown; WISAR6_fourthroot took up to 31 min
#SBATCH --array=0-28%10                # 29 reference dates in the 2023-24 season, at most 10 running at once
#SBATCH -o logs/%x-%A_%a.out
#SBATCH -e logs/%x-%A_%a.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=trobacker@umass.edu
#
# WISAR6_fourthroot_adaptive_t forecasts for the 2023-24 season. Submit from src/WISAR6_fourthroot_adaptive_t, either on
# its own (sbatch submit-unity-season-2023-24.sh) or via submit-unity-all-seasons.sh.
# Dates already in model-output/UMass-WISAR6_fourthroot_adaptive_t/ are skipped.
#
# Sizing (4 CPUs, 4G) copied from ../WISAR6_fourthroot. Time limits are ~3-4x the
# slowest WISAR6_fourthroot task (jobs 65186564-6: 31 / 47 / 69 min) since the
# Student-t WIS loss (betaincinv quantiles) is likely slower to sample.
#
# Dates are this season's valid reference dates in ../../hub-config/tasks.json.

dates=(
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
)

cd "$SLURM_SUBMIT_DIR"
source ./wisar6-task.sh
