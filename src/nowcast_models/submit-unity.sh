#!/bin/bash
#SBATCH -J nowcast_retro
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -p cpu
#SBATCH -t 04:00:00
#SBATCH --array=0-335%80                # 12 models x 28 dates, at most 80 at a time
#SBATCH -o logs/slurm-%A_%a.out
#SBATCH -e logs/slurm-%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nick@umass.edu

# Retrospective 2025/26 forecasts of the nowcast comparison models (see main.py). Run from this directory after
# creating .venv (pip install -r requirements.txt) and saving the Delphi archives:
#   python fetch_archives.py archives
#   sbatch submit-unity.sh
# Forecasts are written to ../../model-output/UMass-<model>/. A task whose forecast file already exists is skipped,
# so the script can be resubmitted to fill in failed tasks.

models=(
  gbqr_1src_nc_none gbqr_1src_nc_scale gbqr_1src_nc_state
  gbqr_4src_nc_none gbqr_4src_nc_scale gbqr_4src_nc_state
  AR6_pooled_nc_none AR6_pooled_nc_scale AR6_pooled_nc_state
  AR2_nc_none AR2_nc_scale AR2_nc_state
)
# Wednesday run dates for the 2025/26 reference dates (2025-11-22 to 2026-05-30), as for gbqr_1src
dates=(
  2025-11-19 2025-11-26 2025-12-03 2025-12-10 2025-12-17 2025-12-24 2025-12-31
  2026-01-07 2026-01-14 2026-01-21 2026-01-28 2026-02-04 2026-02-11 2026-02-18
  2026-02-25 2026-03-04 2026-03-11 2026-03-18 2026-03-25 2026-04-01 2026-04-08
  2026-04-15 2026-04-22 2026-04-29 2026-05-06 2026-05-13 2026-05-20 2026-05-27
)

# one task per (model, date); dates vary fastest, so the early tasks cover every date of the first models
model="${models[$((SLURM_ARRAY_TASK_ID / ${#dates[@]}))]}"
date="${dates[$((SLURM_ARRAY_TASK_ID % ${#dates[@]}))]}"
ref_date=$(date -d "$date + 3 days" +%F)
out="../../model-output/UMass-${model}/${ref_date}-UMass-${model}.csv"
echo "task ${SLURM_ARRAY_TASK_ID}: ${model} ${date} (reference date ${ref_date}) on $(hostname), $(date)"
if [ -f "$out" ]; then
  echo "already done: $out"
  exit 0
fi

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
.venv/bin/python main.py --model "$model" --today_date "$date" --archive_dir archives
status=$?
echo "exit ${status}, $(date)"
exit $status
