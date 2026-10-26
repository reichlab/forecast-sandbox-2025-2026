#!/bin/bash
#SBATCH -J wisar6frat_rerun_2025_26
#SBATCH -N 1
#SBATCH -c 4                            # One core per MCMC chain (same as submit-unity-season-2025-26.sh)
#SBATCH --mem=4G
#SBATCH -p cpu
#SBATCH -t 04:00:00                     # Previous run of this season: slowest date took 02:12:50
#SBATCH --array=0-21%10                # 22 affected reference dates in 2025-26, at most 10 running at once
#SBATCH -o logs/%x-%A_%a.out
#SBATCH -e logs/%x-%A_%a.err
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=trobacker@umass.edu
#
# Rerun of the 2025-26 WISAR6_fourthroot_adaptive_t dates whose existing output was
# produced with the reference-date snapshot look-ahead leak (see
# ../../WISAR6_REFERENCE_DATE_FIX.md). Only dates where the old snapshot
# selection (file named <= R) differs from the fixed one (file named < R) are
# listed (valid reference dates in ../../hub-config/tasks.json only);
# unaffected dates are left as-is.
# WISAR6_OVERWRITE=1 makes wisar6-task.sh refit dates that already have output.
# main.py writes the CSV only at the end, so a failed task leaves the old file.
#
# Submit from src/WISAR6_fourthroot_adaptive_t, on its own
# (sbatch rerun-leakfix-season-2025-26.sh) or via submit-rerun-leakfix.sh.

export WISAR6_OVERWRITE=1

dates=(
  "2025-11-29"
  "2025-12-06"
  "2025-12-13"
  "2025-12-27"
  "2026-01-03"
  "2026-01-10"
  "2026-01-17"
  "2026-01-24"
  "2026-02-07"
  "2026-02-14"
  "2026-02-21"
  "2026-02-28"
  "2026-03-14"
  "2026-03-21"
  "2026-03-28"
  "2026-04-04"
  "2026-04-11"
  "2026-04-18"
  "2026-05-02"
  "2026-05-09"
  "2026-05-16"
  "2026-05-30"
)

cd "$SLURM_SUBMIT_DIR"
source ./wisar6-task.sh
