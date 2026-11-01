#!/bin/bash
# Retrospective ablation: every variant for every other reference date in the
# 2023/24 and 2024/25 seasons (the dates gbqr_3src has outputs for).
#
# Usage: [VARIANTS="none obs"] ./run_eval.sh OUTPUT_ROOT [NUM_BAGS] [N_PARALLEL]
# Requires the data/gefs_*.csv tables built by gefs_data.py.

set -euo pipefail
cd "$(dirname "$0")"
OUT=${1:?output root}
BAGS=${2:-20}
NPAR=${3:-3}
PY=${PY:-../gbqr_3src_weather/.venv/bin/python}
export PYTHONWARNINGS=ignore
export OMP_NUM_THREADS=$(( $(sysctl -n hw.ncpu 2>/dev/null || nproc) / NPAR ))

REF_DATES=$(ls ../../model-output/UMass-gbqr_3src/ | sed 's/-UMass.*//' | awk 'NR % 2 == 1')
mkdir -p "$OUT/logs"
for d in $REF_DATES; do
  today=$(date -j -v-3d -f %Y-%m-%d "$d" +%Y-%m-%d 2>/dev/null || date -d "$d -3 days" +%F)
  for v in ${VARIANTS:-none obs obs_fc obs_oracle}; do
    echo "$today $v"
  done
done | xargs -P "$NPAR" -L 1 bash -c '
  f="'"$OUT"'/model-output/UMass-gbqr_3src_gefs_$1/$(date -j -v+3d -f %Y-%m-%d "$0" +%Y-%m-%d 2>/dev/null || date -d "$0 +3 days" +%F)-UMass-gbqr_3src_gefs_$1.csv"
  [ -f "$f" ] && exit 0
  '"$PY"' main.py --today_date "$0" --variant "$1" --num_bags '"$BAGS"' --short_q \
    --output_root "'"$OUT"'/model-output" > "'"$OUT"'/logs/$0-$1.log" 2>&1 \
    && echo "ok $0 $1" || echo "FAIL $0 $1"
'
