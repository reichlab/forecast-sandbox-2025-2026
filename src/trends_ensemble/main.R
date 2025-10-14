library(trendsEnsemble)
library(idforecastutils)
library(hubData)
library(hubVis)
library(fs)
library(readr)
library(lubridate)

args <- commandArgs(trailingOnly = TRUE)

run_date <- as.Date(args[1])
reference_date <- lubridate::ymd(run_date) + 3

bucket_name <- "cdcepi-flusight-forecast-hub"
hub_bucket <- s3_bucket(bucket_name)
hub_con <- hubData::connect_hub(hub_bucket, file_format = "parquet", skip_checks = TRUE)
target_data_dates <- bucket_name |>
  aws.s3::get_bucket_df(prefix = "auxiliary-data/target-data-archive/target-hospital-admissions_") |>
  dplyr::mutate(date = ymd(stringr::str_sub(Key, start = -14, end = -5))) |>
  dplyr::pull(date)

locations <- read.csv("https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/refs/heads/main/auxiliary-data/locations.csv")
required_quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)

# load target data
# deals with weeks where target data wasn't archived
file_date <- min(target_data_dates[target_data_dates > reference_date - 7])
target_data <- readr::read_csv(paste0("https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/auxiliary-data/target-data-archive/target-hospital-admissions_", file_date, ".csv")) |>
  dplyr::filter(date <= reference_date - 7, !is.na(value))
target_ts <- target_data |>
  dplyr::select("date", "location", "value") |>
  dplyr::arrange(dplyr::desc(date)) |>
  dplyr::rename(time_index = date, observation = value)

# set up variations of baseline to fit
component_variations <- tidyr::expand_grid(
  transformation = c("none", "sqrt"),
  symmetrize = c(TRUE, FALSE),
  window_size = c(3, 4),
  temporal_resolution = "weekly"
)

# Generate ensemble
outputs_list <- create_trends_ensemble(
  component_variations,
  target_ts,
  reference_date,
  horizons = 0:3,
  target = "wk inc flu hosp",
  quantile_levels = required_quantiles,
  n_samples = 100,
  round_predictions = TRUE,
  return_baseline_predictions = TRUE
)
component_outputs <- outputs_list[["baselines"]] |>
  dplyr::mutate(
    output_type_id = ifelse(
      .data[["output_type"]] == "sample",
      paste0(.data[["location"]], sprintf("%02g", .data[["output_type_id"]])),
      as.character(.data[["output_type_id"]])
    )
  )
model_names <- unique(component_outputs$model_id)

# pmf forecasts
trends_ensemble_raw <- outputs_list[["ensemble"]] |>
  dplyr::mutate(
    output_type_id = ifelse(
      .data[["output_type"]] == "sample",
      paste0(.data[["location"]], sprintf("%02g", .data[["output_type_id"]])),
      as.character(.data[["output_type_id"]])
    )
  )

bin_endpoints <- get_flusight_bin_endpoints(
  target_data,
  locations,
  season = ifelse(reference_date < ymd("2024-09-01"), "2023/24", "2024/25")
)
trends_ensemble_pmf <- trends_ensemble_raw |>
  dplyr::filter(output_type == "quantile") |>
  idforecastutils::transform_quantile_to_pmf(bin_endpoints = bin_endpoints) |>
  dplyr::mutate(target = "wk flu hosp rate change") |>
  dplyr::ungroup()
trends_ensemble_outputs <- trends_ensemble_raw |>
  dplyr::bind_rows(trends_ensemble_pmf) |>
  dplyr::mutate(horizon = as.integer(.data[["horizon"]]))
trendsEnsemble::save_model_out_tbl(trends_ensemble_outputs, path = "../../model-output", extension = "csv")

# open PDF
# model_id <- "UMass-trends_ensemble"
# model_folder <- file.path("../../plots", model_id)
# if (!file.exists(model_folder)) dir.create(model_folder, recursive = TRUE)
# grDevices::pdf(file = paste0(model_folder, "/", reference_date, "-", model_id, ".pdf"), paper = "a4r")
# cat_names <- c("large_decrease", "decrease", "stable", "increase", "large_increase")
# idforecastutils::plot_quantile_pmf_outputs_pdf(
#   trends_ensemble_outputs,
#   target_data, locations,
#   reference_date, cats_ordered = cat_names[5:1],
#   quantile_title = "Inc Flu Hosp", pmf_title = "Flu Hosp Rate Change"
# )
