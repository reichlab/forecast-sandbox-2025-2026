#!/bin/bash
# Generate 2025/26 forecasts (the hub's 28 reference dates, 2025-11-22 .. 2026-05-30) for the variants carried into
# 2025/26, sequentially. 2026/27 is still held out.
# Usage: bash run-all-forecasts-2025-2026.sh [variant ...]   (default: none shared xm1)

variants=("$@")
[ ${#variants[@]} -eq 0 ] && variants=(none shared xm1)

dates=(
  "2025-11-19"
  "2025-11-26"
  "2025-12-03"
  "2025-12-10"
  "2025-12-17"
  "2025-12-24"
  "2025-12-31"
  "2026-01-07"
  "2026-01-14"
  "2026-01-21"
  "2026-01-28"
  "2026-02-04"
  "2026-02-11"
  "2026-02-18"
  "2026-02-25"
  "2026-03-04"
  "2026-03-11"
  "2026-03-18"
  "2026-03-25"
  "2026-04-01"
  "2026-04-08"
  "2026-04-15"
  "2026-04-22"
  "2026-04-29"
  "2026-05-06"
  "2026-05-13"
  "2026-05-20"
  "2026-05-27"
)

for variant in "${variants[@]}"
do
  for date in "${dates[@]}"
  do
    echo "Running variant $variant for date: $date"
    python main.py --variant="$variant" --today_date="$date"
  done
done
