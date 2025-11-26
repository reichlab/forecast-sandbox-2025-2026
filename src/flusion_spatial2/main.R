## ad hoc script to build ensemble model from existing model files

library(hubEnsembles)
library(hubData)
library(dplyr)
library(readr)

## build new model directory
new_model <- "UMass-flusion_spatial2"
if (!dir.exists(file.path("model-output", new_model))) {
  dir.create(file.path("model-output", new_model))
}

## obtain lists of existing files
state_models_to_blend <- c("gbqr_3src_spatial2", "AR6_pooled")
us_models_to_blend <- c("gbqr_3src", "AR6_pooled")
if (length(state_models_to_blend) != 2) {
  stop("only testing exists for two models.")
}

file_list <- lapply(state_models_to_blend, FUN = function(x) {
  list.files(path = file.path("model-output", paste0("UMass-", x)))
})

## check that dates match
date_list <- lapply(file_list, FUN = function(x) {
  substr(x, 0, 10)
})
if (!all.equal(date_list[[1]], date_list[[2]])) {
  stop("Dates for the models are not equivalent.")
}

## loop through and ensemble
hub_con <- connect_hub(hub_path = ".", skip_checks = TRUE)
us_models <- paste0("UMass-", us_models_to_blend)
state_models <- paste0("UMass-", state_models_to_blend)

for (i in 1:length(file_list[[1]])) {
  state_dat <- hub_con |>
    filter(
      reference_date == date_list[[1]][i],
      model_id %in% state_models,
      location != "US"
    ) |>
    collect_hub()
  us_dat <- hub_con |>
    filter(
      reference_date == date_list[[1]][i],
      model_id %in% us_models,
      location == "US"
    ) |>
    collect_hub()

  all_dat <- bind_rows(state_dat, us_dat)

  ens_model <- simple_ensemble(
    all_dat,
    model_id = new_model,
    agg_fun = "median"
  ) |>
    select(-model_id)

  write_csv(
    ens_model,
    file = file.path(
      "model-output",
      new_model,
      paste0(date_list[[1]][i], "-", new_model, ".csv")
    )
  )
}
