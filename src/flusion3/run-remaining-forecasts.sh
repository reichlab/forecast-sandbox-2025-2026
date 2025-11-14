#!/bin/bash

# Generate remaining forecasts with optimized settings (20 GBQR bags instead of 100)
# Skipping already completed forecasts

dates=(
  "2023-12-09"
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
  "2024-05-08"
)

echo "Remaining forecasts to generate: ${#dates[@]}"
echo ""

for date in "${dates[@]}"
do
  echo "============================================================"
  echo "Running flusion3 for date: $date"
  echo "============================================================"
  python main.py --today_date="$date"

  if [ $? -eq 0 ]; then
    echo "✓ Successfully completed forecast for $date"
  else
    echo "✗ Error generating forecast for $date"
    exit 1
  fi

  echo ""
done

echo "============================================================"
echo "All remaining forecasts completed successfully!"
echo "============================================================"
