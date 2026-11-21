#!/bin/bash
# Generate 2024/25 forecasts (the hub's 23 reference dates, 2024-11-30 .. 2025-05-17; 2025-01-25 and 2025-04-26 are
# not hub rounds) for the variants carried into 2024/25, sequentially. 2025/26 and later are still held out.
# Usage: bash run-all-forecasts-2024-2025.sh [variant ...]   (default: none shared xm1 xm1_ny1)

variants=("$@")
[ ${#variants[@]} -eq 0 ] && variants=(none shared xm1 xm1_ny1)

dates=(
  "2024-11-27"
  "2024-12-04"
  "2024-12-11"
  "2024-12-18"
  "2024-12-25"
  "2025-01-01"
  "2025-01-08"
  "2025-01-15"
  "2025-01-29"
  "2025-02-05"
  "2025-02-12"
  "2025-02-19"
  "2025-02-26"
  "2025-03-05"
  "2025-03-12"
  "2025-03-19"
  "2025-03-26"
  "2025-04-02"
  "2025-04-09"
  "2025-04-16"
  "2025-04-30"
  "2025-05-07"
  "2025-05-14"
)

for variant in "${variants[@]}"
do
  for date in "${dates[@]}"
  do
    echo "Running variant $variant for date: $date"
    python main.py --variant="$variant" --today_date="$date"
  done
done
