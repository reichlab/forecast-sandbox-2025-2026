#!/bin/bash
# Generate forecasts for this model across every reference_date in hub-config/tasks.json.
# Excludes 2025-11-26, 2025-12-24, 2025-12-31 (Thanksgiving/Christmas/New Year's run dates):
# NSSP's holiday-week reporting lag produces invalid output for these reference dates.

dates=(
  "2025-11-19"
  "2025-12-03"
  "2025-12-10"
  "2025-12-17"
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

for date in "${dates[@]}"
do
  echo "Running for date: $date"
  python main.py --today_date="$date"
done
