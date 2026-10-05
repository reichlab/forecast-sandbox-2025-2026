# Build a quantile-average ensemble of two existing hub models for every
# reference date both models have in model-output/.
#
# Usage (from this directory):
#   Rscript ensemble_pair.R <output_model_id> <component_model_id> <component_model_id>

library(dplyr)
library(hubUtils)
library(hubEnsembles)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
ens_model_id <- args[1]
component_ids <- args[2:3]

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
    filter(horizon >= 0) |>
    as_model_out_tbl()

  ens_model <- simple_ensemble(component_dat, model_id = ens_model_id) |>
    select(location, reference_date, horizon, target_end_date, target,
           output_type, output_type_id, value)

  utils::write.csv(
    ens_model,
    file = file.path(output_dir, paste0(ref_date, "-", ens_model_id, ".csv")),
    row.names = FALSE
  )
  message("Wrote ", ens_model_id, " ", ref_date, " (", nrow(ens_model), " rows)")
}
