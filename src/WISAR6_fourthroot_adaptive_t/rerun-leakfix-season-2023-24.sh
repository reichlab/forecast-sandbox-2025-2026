#!/bin/bash
#SBATCH -J wisar6frat_rerun_2023_24
#SBATCH -N 1
#SBATCH -c 4                            # One core per MCMC chain (same as submit-unity-season-2023-24.sh)
#SBATCH --mem=4G
#SBATCH -p cpu
#SBATCH -t 02:00:00                     # Previous run of this season: slowest date took 00:41:56
#SBATCH --array=0-27%10                # 28 affected reference dates in 2023-24, at most 10 running at once
#SBATCH -o logs/%x-%A_%a.out
#SBATCH -e logs/%x-%A_%a.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=trobacker@umass.edu
#
# Rerun of the 2023-24 WISAR6_fourthroot_adaptive_t dates whose existing output was
# produced with the reference-date snapshot look-ahead leak (see
# ../../WISAR6_REFERENCE_DATE_FIX.md). Only dates where the old snapshot
# selection (file named <= R) differs from the fixed one (file named < R) are
# listed (valid reference dates in ../../hub-config/tasks.json only);
# unaffected dates are left as-is.
# WISAR6_OVERWRITE=1 makes wisar6-task.sh refit dates that already have output.
# main.py writes the CSV only at the end, so a failed task leaves the old file.
#
# Submit from src/WISAR6_fourthroot_adaptive_t, on its own
# (sbatch rerun-leakfix-season-2023-24.sh) or via submit-rerun-leakfix.sh.

export WISAR6_OVERWRITE=1

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
)

cd "$SLURM_SUBMIT_DIR"
source ./wisar6-task.sh
