#!/bin/bash
# Generate forecasts for every (AR model, GBQR model) ensemble pair, across every
# reference_date in hub-config/tasks.json.
# Excludes 2025-11-26, 2025-12-24, 2025-12-31 (Thanksgiving/Christmas/New Year's run dates):
# NSSP's holiday-week reporting lag produces invalid output for these reference dates.

ar_models=(
  "AR4_fourierP_thetaP_nssp"
  "AR2_fourierP_thetaP_nssp"
  "AR6_fourierP_thetaP_nssp"
)

gbqr_models=(
  "gbqr_2src_flusurvnet_nssp"
  "gbqr_4src_nssp"
  "gbqr_4src_spatial_nssp"
  "gbqr_2src_ilinet_nssp"
  "gbqr_2src_nhsn_nssp"
)

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

for ar_model in "${ar_models[@]}"
do
  for gbqr_model in "${gbqr_models[@]}"
  do
    echo "=== Ensemble: $ar_model + $gbqr_model ==="
    for date in "${dates[@]}"
    do
      echo "Running for date: $date"
      Rscript main.R --ar_model="$ar_model" --gbqr_model="$gbqr_model" --today_date="$date"
    done
  done
done
