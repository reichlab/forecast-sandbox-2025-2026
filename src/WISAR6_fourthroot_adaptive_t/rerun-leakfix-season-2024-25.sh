#!/bin/bash
#SBATCH -J wisar6frat_rerun_2024_25
#SBATCH -N 1
#SBATCH -c 4                            # One core per MCMC chain (same as submit-unity-season-2024-25.sh)
#SBATCH --mem=4G
#SBATCH -p cpu
#SBATCH -t 02:00:00                     # Previous run of this season: slowest date took 00:55:53
#SBATCH --array=0-19%10                # 20 affected reference dates in 2024-25, at most 10 running at once
#SBATCH -o logs/%x-%A_%a.out
#SBATCH -e logs/%x-%A_%a.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=trobacker@umass.edu
#
# Rerun of the 2024-25 WISAR6_fourthroot_adaptive_t dates whose existing output was
# produced with the reference-date snapshot look-ahead leak (see
# ../../WISAR6_REFERENCE_DATE_FIX.md). Only dates where the old snapshot
# selection (file named <= R) differs from the fixed one (file named < R) are
# listed (valid reference dates in ../../hub-config/tasks.json only);
# unaffected dates are left as-is.
# WISAR6_OVERWRITE=1 makes wisar6-task.sh refit dates that already have output.
# main.py writes the CSV only at the end, so a failed task leaves the old file.
#
# Submit from src/WISAR6_fourthroot_adaptive_t, on its own
# (sbatch rerun-leakfix-season-2024-25.sh) or via submit-rerun-leakfix.sh.

export WISAR6_OVERWRITE=1

dates=(
  "2024-11-30"
  "2024-12-14"
  "2024-12-21"
  "2024-12-28"
  "2025-01-11"
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
