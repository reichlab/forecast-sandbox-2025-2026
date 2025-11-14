#!/usr/bin/env Rscript

# Weighted ensemble for flusion3 using learned weights
# Learns weights based on historical WIS performance

library(dplyr)
library(hubData)
library(hubEnsembles)

args <- commandArgs(trailingOnly = TRUE)
ref_date <- as.Date(args[1])

# Load component forecasts for this reference date
hub_con <- hubData::connect_model_output("intermediate-output/model-output")
component_forecasts <- dplyr::collect(hub_con) |>
  dplyr::filter(reference_date == ref_date, horizon >= 0)

if (nrow(component_forecasts) == 0) {
  stop("No component forecasts found for reference date ", ref_date)
}

# Get list of component models
component_models <- unique(component_forecasts$model_id)
cat("Component models:", paste(component_models, collapse=", "), "\n")

# Load historical performance to compute weights
# We'll use inverse WIS from all previous forecasts
all_past_forecasts <- dplyr::collect(hub_con) |>
  dplyr::filter(reference_date < ref_date, horizon >= 0)

if (nrow(all_past_forecasts) > 0) {
  # Load target data to compute WIS
  target_ts <- readr::read_csv(
    "https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/target-data/target-hospital-admissions.csv",
    show_col_types = FALSE
  )

  # Compute WIS for each model on historical data
  wis_by_model <- list()

  for (model in component_models) {
    model_forecasts <- all_past_forecasts |>
      dplyr::filter(model_id == model, output_type == "quantile")

    if (nrow(model_forecasts) > 0) {
      # Join with targets
      eval_data <- model_forecasts |>
        dplyr::inner_join(
          target_ts |> dplyr::select(location, date, value),
          by = c("location", "target_end_date" = "date")
        )

      if (nrow(eval_data) > 0) {
        # Simple WIS approximation: mean absolute error at median
        median_forecasts <- eval_data |>
          dplyr::filter(abs(as.numeric(output_type_id) - 0.5) < 0.01) |>
          dplyr::mutate(ae = abs(value.x - value.y))

        if (nrow(median_forecasts) > 0) {
          wis_by_model[[model]] <- mean(median_forecasts$ae, na.rm = TRUE)
        }
      }
    }
  }

  # Convert to weights (inverse WIS, normalized)
  if (length(wis_by_model) > 0) {
    # Inverse WIS (lower WIS = higher weight)
    inverse_wis <- sapply(component_models, function(m) {
      if (m %in% names(wis_by_model)) {
        1 / wis_by_model[[m]]
      } else {
        1  # Default if no history
      }
    })

    # Normalize to sum to 1
    weights_vec <- inverse_wis / sum(inverse_wis)
    names(weights_vec) <- component_models

    cat("Learned weights (inverse WIS):\n")
    for (i in seq_along(weights_vec)) {
      cat(sprintf("  %s: %.3f (WIS: %.2f)\n",
                  names(weights_vec)[i],
                  weights_vec[i],
                  if (names(weights_vec)[i] %in% names(wis_by_model)) wis_by_model[[names(weights_vec)[i]]] else NA))
    }

    # Convert to data frame format required by hubEnsembles
    weights <- data.frame(
      model_id = component_models,
      weight = weights_vec[component_models],
      row.names = NULL
    )
  } else {
    # No historical data, use equal weights
    weights_vec <- rep(1/length(component_models), length(component_models))
    names(weights_vec) <- component_models
    cat("No historical data available, using equal weights\n")

    # Convert to data frame format required by hubEnsembles
    weights <- data.frame(
      model_id = component_models,
      weight = weights_vec,
      row.names = NULL
    )
  }
} else {
  # No historical data, use equal weights
  weights_vec <- rep(1/length(component_models), length(component_models))
  names(weights_vec) <- component_models
  cat("No historical data available, using equal weights\n")

  # Convert to data frame format required by hubEnsembles
  weights <- data.frame(
    model_id = component_models,
    weight = weights_vec,
    row.names = NULL
  )
}

# Create weighted ensemble using linear pool
model_out_tbl <- hubEnsembles::linear_pool(
  model_out_tbl = component_forecasts,
  weights = weights,
  model_id = "UMass-flusion3_weighted"
)

# Add categorical target predictions
target_ts <- readr::read_csv(
  "https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/target-data/target-hospital-admissions.csv",
  show_col_types = FALSE
)
location_meta <- readr::read_csv(
  "https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/refs/heads/main/auxiliary-data/locations.csv",
  show_col_types = FALSE
)
bin_endpoints <- idforecastutils::get_flusight_bin_endpoints(
  target_ts = target_ts,
  location_meta = location_meta,
  season = "2024/25"
)
categorical_outputs <- idforecastutils::transform_quantile_to_pmf(
  model_out_tbl = model_out_tbl,
  bin_endpoints = bin_endpoints
) |>
  dplyr::mutate(target = "wk flu hosp rate change") |>
  dplyr::ungroup()

# Save output
reference_date <- model_out_tbl$reference_date[1]
output_dir <- "../../model-output/UMass-flusion3_weighted"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_file <- file.path(
  output_dir,
  paste0(reference_date, "-UMass-flusion3_weighted.csv")
)

utils::write.csv(
  model_out_tbl |> dplyr::select(-model_id),
  file = output_file,
  row.names = FALSE
)

cat("\nWeighted ensemble forecast saved to:", output_file, "\n")
