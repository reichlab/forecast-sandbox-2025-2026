library(dplyr)
library(hubData)
library(hubEnsembles)


args <- commandArgs(trailingOnly = TRUE)
ref_date <- as.Date(args[1])

hub_con <- hubData::connect_model_output(
  "intermediate-output/model-output",
  schema = arrow::schema(location = arrow::utf8())
)

# gbqr_3src_spatial has no US-level forecasts (its directional wave features
# are state-level only), so the US blend uses the plain (non-spatial)
# gbqr_3src instead, same as flusion_spatial2_prod's state/US split.
state_models_to_blend <- c("UMass-gbqr_3src_spatial", "UMass-WISAR6_fourthroot_adaptive_t")
us_models_to_blend <- c("UMass-gbqr_3src", "UMass-WISAR6_fourthroot_adaptive_t")

state_dat <- hub_con |>
  filter(
    reference_date == ref_date,
    model_id %in% state_models_to_blend,
    location != "US",
    horizon >= 0
  ) |>
  collect_hub()

us_dat <- hub_con |>
  filter(
    reference_date == ref_date,
    model_id %in% us_models_to_blend,
    location == "US",
    horizon >= 0
  ) |>
  collect_hub()

ens_model <- bind_rows(state_dat, us_dat) |>
  simple_ensemble(model_id = "UMass-gbqr_3src_spatial_wisar") |>
  dplyr::select(-model_id)

# publish the hosp-target ensemble on its own, unchanged
hosp_output_dir <- "../../model-output/UMass-gbqr_3src_spatial_wisar"
if (!dir.exists(hosp_output_dir)) {
  dir.create(hosp_output_dir, recursive = TRUE)
}
utils::write.csv(
  ens_model,
  file = file.path(
    hosp_output_dir,
    paste0(ref_date, "-UMass-gbqr_3src_spatial_wisar.csv")
  ),
  row.names = FALSE
)

# gbqr_4src_nssp is an independent model (not blended), predicting a
# different target (wk inc flu prop ed visits) than the hosp ensemble above.
# When it has run for this reference date, row-bind its forecast onto the
# hosp ensemble into one combined submission file; when it hasn't (it only
# covers 2025/26-forward, unlike the hosp components' full 2023/24-2025/26
# history), publish the hosp ensemble alone under the combined name.
nssp_dat <- hub_con |>
  filter(
    reference_date == ref_date,
    model_id == "UMass-gbqr_4src_nssp",
    horizon >= 0
  ) |>
  collect_hub()

# gbqr_4src_nssp can produce NA for a location it has no current data for
# (e.g. a state with a live NSSP reporting gap); a submission file should
# simply omit that location rather than include null values, so drop NA rows
# here and report what got dropped.
na_rows <- nssp_dat |> dplyr::filter(is.na(value))
if (nrow(na_rows) > 0) {
  na_locations <- sort(unique(na_rows$location))
  message(ref_date, ": dropping ", nrow(na_rows), " NA row(s) from UMass-gbqr_4src_nssp ",
          "for location(s): ", paste(na_locations, collapse = ", "))
  nssp_dat <- nssp_dat |> dplyr::filter(!is.na(value))
}

if (nrow(nssp_dat) > 0) {
  combined <- bind_rows(ens_model, nssp_dat |> dplyr::select(-model_id))
} else {
  message(ref_date, ": no UMass-gbqr_4src_nssp output found, publishing hosp-only")
  combined <- ens_model
}

combined_output_dir <- "../../model-output/UMass-gbqr_3src_spatial_wisar_nssp"
if (!dir.exists(combined_output_dir)) {
  dir.create(combined_output_dir, recursive = TRUE)
}
utils::write.csv(
  combined,
  file = file.path(
    combined_output_dir,
    paste0(ref_date, "-UMass-gbqr_3src_spatial_wisar_nssp.csv")
  ),
  row.names = FALSE
)
