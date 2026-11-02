suppressMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(hubEnsembles)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  val <- sub(paste0("^--", flag, "="), "", grep(paste0("^--", flag, "="), args, value = TRUE))
  if (length(val) != 1) stop(paste0("Missing or duplicated --", flag, " argument"))
  val
}

ar_model <- get_arg("ar_model")
gbqr_model <- get_arg("gbqr_model")
today_date_arg <- sub("^--today_date=", "", grep("^--today_date=", args, value = TRUE))
today_date <- if (length(today_date_arg) == 1) ymd(today_date_arg) else Sys.Date()

# Roll the run date forward to the following Saturday reference date, matching the
# convention used by the Python models (main.py's `relativedelta(weekday=5)`).
iso_wday <- lubridate::wday(today_date, week_start = 1)
reference_date <- today_date + ((6 - iso_wday) %% 7)

# Ensemble model name: AR and GBQR component names (each ending in "_nssp") joined, with
# "_ens_nssp" replacing the trailing "_nssp" so the combined name isn't doubly suffixed.
ar_short <- sub("_nssp$", "", ar_model)
gbqr_short <- sub("_nssp$", "", gbqr_model)
model_name <- paste0(ar_short, "_", gbqr_short, "_ens_nssp")

component_models <- c(ar_model, gbqr_model)
component_dfs <- lapply(component_models, function(m) {
  path <- file.path("..", "..", "model-output", paste0("UMass-", m),
                     paste0(reference_date, "-UMass-", m, ".csv"))
  read_csv(path, show_col_types = FALSE) |>
    mutate(model_id = paste0("UMass-", m))
})
combined <- bind_rows(component_dfs)

ensemble <- simple_ensemble(
  combined,
  agg_fun = median,
  model_id = paste0("UMass-", model_name),
  task_id_cols = c("location", "reference_date", "horizon", "target_end_date", "target")
) |>
  select(location, reference_date, horizon, target_end_date, target, output_type, output_type_id, value)

output_dir <- file.path("..", "..", "model-output", paste0("UMass-", model_name))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_path <- file.path(output_dir, paste0(reference_date, "-UMass-", model_name, ".csv"))
write_csv(ensemble, output_path)
cat("Wrote", nrow(ensemble), "rows to", output_path, "\n")
