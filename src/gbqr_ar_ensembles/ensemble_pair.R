# Build an equally weighted quantile-average ensemble of a GBQR model and an AR
# model for every reference date the components share in model-output/.
#
# Usage (from this directory):
#   Rscript ensemble_pair.R <output_model_id> <gbqr_model_id> <ar_model_id> [us_gbqr_model_id]
#
# If us_gbqr_model_id is given, it replaces gbqr_model_id for the national (US)
# forecast, as in flusion_spatial2_prod (state-only gbqr_3src_spatial at state
# level, gbqr_3src for US).

library(dplyr)
library(hubUtils)
library(hubEnsembles)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(3, 4))
ens_model_id <- args[1]
gbqr_id <- args[2]
ar_id <- args[3]
us_gbqr_id <- if (length(args) == 4) args[4] else gbqr_id

state_models <- c(gbqr_id, ar_id)
us_models <- c(us_gbqr_id, ar_id)
component_ids <- unique(c(state_models, us_models))

model_output_root <- "../../model-output"

ref_dates_for <- function(model_id) {
  files <- list.files(file.path(model_output_root, model_id), pattern = "\\.csv$")
  substr(files, 1, 10)
}
ref_dates <- sort(Reduce(intersect, lapply(component_ids, ref_dates_for)))
if (length(ref_dates) == 0) stop("No shared reference dates for ", paste(component_ids, collapse = ", "))

output_dir <- file.path(model_output_root, ens_model_id)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (ref_date in ref_dates) {
  component_dat <- lapply(component_ids, function(model_id) {
    readr::read_csv(
      file.path(model_output_root, model_id, paste0(ref_date, "-", model_id, ".csv")),
      col_types = readr::cols(location = "c", output_type_id = "c", .default = "?")
    ) |>
      mutate(model_id = model_id)
  }) |>
    bind_rows() |>
    filter(horizon >= 0)

  state_dat <- component_dat |> filter(location != "US", model_id %in% state_models)
  us_dat <- component_dat |> filter(location == "US", model_id %in% us_models)

  # every location must have a forecast from both of its components
  n_models <- bind_rows(state_dat, us_dat) |>
    distinct(location, model_id) |>
    count(location)
  if (any(n_models$n != 2)) {
    stop(ref_date, ": locations without both components: ",
         paste(n_models$location[n_models$n != 2], collapse = ", "))
  }

  ens_model <- bind_rows(state_dat, us_dat) |>
    as_model_out_tbl() |>
    simple_ensemble(model_id = ens_model_id) |>
    select(location, reference_date, horizon, target_end_date, target,
           output_type, output_type_id, value)

  utils::write.csv(
    ens_model,
    file = file.path(output_dir, paste0(ref_date, "-", ens_model_id, ".csv")),
    row.names = FALSE
  )
  message("Wrote ", ens_model_id, " ", ref_date, " (", nrow(ens_model), " rows)")
}
