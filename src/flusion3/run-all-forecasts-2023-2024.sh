#!/bin/bash

# Generate full forecasts for 2023-24 season (no --short_run flag)
# This will use full MCMC samples (2000 warmup + 2000 samples) and 23 quantile levels

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
  "2024-05-08"
)

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
echo "All forecasts completed successfully!"
echo "============================================================"
