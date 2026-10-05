#!/bin/bash
#SBATCH -J wisar6frat_2025_26
#SBATCH -N 1
#SBATCH -c 4                            # One core per MCMC chain (test: ~100% CPU per chain process)
#SBATCH --mem=4G                        # Test MaxRSS was 0.48-0.78 GB; 4G leaves headroom
#SBATCH -p cpu
#SBATCH -t 04:00:00                     # Student-t sampling cost unknown; WISAR6_fourthroot took up to 69 min
#SBATCH --array=0-26%10                # 27 reference dates in the 2025-26 season, at most 10 running at once
#SBATCH -o logs/%x-%A_%a.out
#SBATCH -e logs/%x-%A_%a.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=trobacker@umass.edu
#
# WISAR6_fourthroot_adaptive_t forecasts for the 2025-26 season. Submit from src/WISAR6_fourthroot_adaptive_t, either on
# its own (sbatch submit-unity-season-2025-26.sh) or via submit-unity-all-seasons.sh.
# Dates already in model-output/UMass-WISAR6_fourthroot_adaptive_t/ are skipped.
#
# Sizing (4 CPUs, 4G) copied from ../WISAR6_fourthroot. Time limits are ~3-4x the
# slowest WISAR6_fourthroot task (jobs 65186564-6: 31 / 47 / 69 min) since the
# Student-t WIS loss (betaincinv quantiles) is likely slower to sample.
#
# Dates are identical to ../WISAR6_fourthroot/submit-unity-season-2025-26.sh.

dates=(
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

cd "$SLURM_SUBMIT_DIR"
source ./wisar6-task.sh
