library(dplyr)
library(hubData)
library(hubEnsembles)


args <- commandArgs(trailingOnly = TRUE)
ref_date <- as.Date(args[1])

hub_con <- hubData::connect_model_output("intermediate-output/model-output")

state_models_to_blend <- c("gbqr_3src_spatial2", "AR6_pooled")
us_models_to_blend <- c("gbqr_3src", "AR6_pooled")

state_dat <- hub_con |>
  filter(
    reference_date == ref_date,
    model_id %in% state_models,
    location != "US",
    horizon >= 0
  ) |>
  collect_hub()

us_dat <- hub_con |>
  filter(
    reference_date == ref_date,
    model_id %in% us_models,
    location == "US",
    horizon >= 0
  ) |>
  collect_hub()

ens_model <- bind_rows(state_dat, us_dat) |>
  simple_ensemble(model_id = "UMass-flusion_spatial2_prod")


# save
output_dir <- "../../model-output/UMass-flusion_spatial2_prod"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

utils::write.csv(
  model_out_tbl |> dplyr::select(-model_id),
  file = file.path(
    output_dir,
    paste0(ref_date, "-UMass-flusion_spatial2_prod.csv")
  ),
  row.names = FALSE
)
