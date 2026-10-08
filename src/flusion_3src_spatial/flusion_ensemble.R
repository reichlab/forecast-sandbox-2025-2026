library(dplyr)
library(readr)
library(hubEnsembles)


args <- commandArgs(trailingOnly = TRUE)
ref_date <- as.Date(args[1])

# Read each component's CSV directly rather than via hubData::connect_model_output,
# which scans every model directory in ../../model-output and can fail with a
# schema-merge error ("incompatible types: double vs int64") caused by unrelated
# models elsewhere in the hub.
read_model_csv <- function(model_abbr) {
  path <- file.path(
    "../../model-output",
    paste0("UMass-", model_abbr),
    paste0(ref_date, "-UMass-", model_abbr, ".csv")
  )
  readr::read_csv(
    path,
    col_types = readr::cols(
      location = readr::col_character(),
      value = readr::col_double(),
      .default = readr::col_guess()
    ),
    show_col_types = FALSE
  ) |>
    dplyr::mutate(
      model_id = paste0("UMass-", model_abbr),
      reference_date = as.Date(reference_date)
    )
}

state_models_to_blend <- c("gbqr_3src_spatial", "AR6_pooled")
us_models_to_blend <- c("gbqr_3src", "AR6_pooled")

state_dat <- bind_rows(lapply(state_models_to_blend, read_model_csv)) |>
  filter(reference_date == ref_date, location != "US", horizon >= 0)

us_dat <- bind_rows(lapply(us_models_to_blend, read_model_csv)) |>
  filter(reference_date == ref_date, location == "US", horizon >= 0)

stopifnot(
  "state-level forecasts missing for one or more component models" =
    length(unique(state_dat$model_id)) == length(state_models_to_blend),
  "US-level forecasts missing for one or more component models" =
    length(unique(us_dat$model_id)) == length(us_models_to_blend)
)

ens_model <- bind_rows(state_dat, us_dat) |>
  simple_ensemble(model_id = "UMass-flusion_3src_spatial")


# save
output_dir <- "../../model-output/UMass-flusion_3src_spatial"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

utils::write.csv(
  ens_model |> dplyr::select(-model_id),
  file = file.path(
    output_dir,
    paste0(ref_date, "-UMass-flusion_3src_spatial.csv")
  ),
  row.names = FALSE
)
