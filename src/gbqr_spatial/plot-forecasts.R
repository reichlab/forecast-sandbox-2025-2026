library(hubData)
library(hubVis)
library(dplyr)

## load two gbqr model files
hub_con <- connect_hub(".", skip_checks = TRUE)
dat <- hub_con |>
  filter(model_id %in% c("UMass-gbqr", "UMass-gbqr_spatial")) |>
  collect_hub()

## use hubvis to plot forecasts to compare
dat |>
  filter(location == "25") |>
  plot_step_ahead_model_output(
    target_data = dat,
    plot_target = FALSE,
    use_median_as_point = TRUE,
    interactive = FALSE,
    x_col_name = "target_end_date"
  )
