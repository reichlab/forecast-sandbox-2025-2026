#!/usr/bin/env Rscript
#
# Compare performance between UMass-flusion3 (equal weights) and
# UMass-flusion3_weighted (learned weights)
#

library(dplyr)
library(readr)
library(tidyr)

# Load target data
target_ts <- readr::read_csv(
  "https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/target-data/target-hospital-admissions.csv",
  show_col_types = FALSE
)

# Test dates to evaluate (where weights were learned)
test_dates <- c("2023-12-16", "2023-12-23", "2023-12-30")

# Function to load forecasts for a model
load_forecasts <- function(model_name, dates) {
  forecasts <- list()
  for (date in dates) {
    file_path <- sprintf("../../model-output/%s/%s-%s.csv", model_name, date, model_name)
    if (file.exists(file_path)) {
      df <- readr::read_csv(file_path, show_col_types = FALSE) %>%
        dplyr::mutate(reference_date = as.Date(date))
      forecasts[[date]] <- df
    } else {
      cat("Warning: File not found:", file_path, "\n")
    }
  }
  dplyr::bind_rows(forecasts)
}

# Function to compute WIS approximation (median absolute error)
compute_wis <- function(forecasts, targets) {
  # Join forecasts with targets
  eval_data <- forecasts %>%
    dplyr::filter(output_type == "quantile") %>%
    dplyr::inner_join(
      targets %>% dplyr::select(location, date, value),
      by = c("location", "target_end_date" = "date")
    ) %>%
    dplyr::rename(forecast = value.x, observed = value.y)

  # Compute WIS approximation: mean absolute error at median
  wis_by_model <- eval_data %>%
    dplyr::filter(abs(as.numeric(output_type_id) - 0.5) < 0.01) %>%
    dplyr::group_by(model_id = reference_date) %>%
    dplyr::summarize(
      n_forecasts = n(),
      wis = mean(abs(forecast - observed), na.rm = TRUE),
      .groups = "drop"
    )

  return(wis_by_model)
}

# Load forecasts
cat("Loading forecasts...\n")
flusion3_forecasts <- load_forecasts("UMass-flusion3", test_dates) %>%
  dplyr::mutate(model = "flusion3 (equal weights)")

flusion3_weighted_forecasts <- load_forecasts("UMass-flusion3_weighted", test_dates) %>%
  dplyr::mutate(model = "flusion3_weighted (learned weights)")

# Combine forecasts
all_forecasts <- dplyr::bind_rows(flusion3_forecasts, flusion3_weighted_forecasts)

cat("\nForecast counts by model:\n")
print(all_forecasts %>%
  dplyr::group_by(model, reference_date) %>%
  dplyr::summarize(n = n(), .groups = "drop"))

# Compute WIS for each model
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("Computing WIS (Weighted Interval Score approximation)...\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

wis_results <- list()

for (test_date in test_dates) {
  cat("Date:", test_date, "\n")

  date_forecasts <- all_forecasts %>%
    dplyr::filter(reference_date == as.Date(test_date))

  for (model_name in c("flusion3 (equal weights)", "flusion3_weighted (learned weights)")) {
    model_forecasts <- date_forecasts %>%
      dplyr::filter(model == model_name)

    if (nrow(model_forecasts) > 0) {
      wis <- compute_wis(model_forecasts, target_ts)
      wis_results[[paste(test_date, model_name)]] <- data.frame(
        date = test_date,
        model = model_name,
        wis = wis$wis,
        n_forecasts = wis$n_forecasts
      )
      cat(sprintf("  %s: WIS = %.3f (n=%d)\n", model_name, wis$wis, wis$n_forecasts))
    } else {
      cat(sprintf("  %s: No forecasts found\n", model_name))
    }
  }
  cat("\n")
}

# Combine results
wis_summary <- dplyr::bind_rows(wis_results)

# Overall summary
cat(paste(rep("=", 60), collapse=""), "\n")
cat("OVERALL SUMMARY\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

overall <- wis_summary %>%
  dplyr::group_by(model) %>%
  dplyr::summarize(
    mean_wis = mean(wis, na.rm = TRUE),
    n_dates = n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(mean_wis)

print(overall)

# Compute improvement
if (nrow(overall) == 2) {
  baseline_wis <- overall$mean_wis[overall$model == "flusion3 (equal weights)"]
  weighted_wis <- overall$mean_wis[overall$model == "flusion3_weighted (learned weights)"]
  improvement_pct <- (baseline_wis - weighted_wis) / baseline_wis * 100

  cat("\nPerformance comparison:\n")
  cat(sprintf("  Baseline (equal weights): WIS = %.3f\n", baseline_wis))
  cat(sprintf("  Weighted (learned weights): WIS = %.3f\n", weighted_wis))
  cat(sprintf("  Improvement: %.2f%%\n", improvement_pct))

  if (improvement_pct > 0) {
    cat("\n✓ Learned weights IMPROVED performance!\n")
  } else {
    cat("\n✗ Learned weights did not improve performance.\n")
  }
}

cat("\n")
