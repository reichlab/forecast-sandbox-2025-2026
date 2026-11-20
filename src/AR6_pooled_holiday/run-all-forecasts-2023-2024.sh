#!/bin/bash
# Generate 2023/24 forecasts (reference dates 2023-10-21 .. 2024-05-04) for all variants, sequentially.
# Later seasons are held out for now: do not add 2024/25 or later dates here yet.
# Usage: bash run-all-forecasts-2023-2024.sh [variant ...]   (default: none shared partial xm1 xm1_ny1)

variants=("$@")
[ ${#variants[@]} -eq 0 ] && variants=(none shared partial xm1 xm1_ny1)

dates=(
  "2023-10-18"
  "2023-10-25"
  "2023-11-01"
  "2023-11-08"
  "2023-11-15"
  "2023-11-22"
  "2023-11-29"
  "2023-12-06"
  "2023-12-13"
  "2023-12-20"
  "2023-12-27"
  "2024-01-03"
  "2024-01-10"
  "2024-01-17"
  "2024-01-24"
  "2024-01-31"
  "2024-02-07"
  "2024-02-14"
  "2024-02-21"
  "2024-02-28"
  "2024-03-06"
  "2024-03-13"
  "2024-03-20"
  "2024-03-27"
  "2024-04-03"
  "2024-04-10"
  "2024-04-17"
  "2024-04-24"
  "2024-05-01"
)

for variant in "${variants[@]}"
do
  for date in "${dates[@]}"
  do
    echo "Running variant $variant for date: $date"
    python main.py --variant="$variant" --today_date="$date"
  done
done
