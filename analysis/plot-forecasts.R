#!/usr/bin/env Rscript
# Plot forecasts from a model output file against truth data
# Usage: Rscript plot-forecasts.R <forecast_file> [truth_file]
# Example: Rscript plot-forecasts.R ../model-output/UMass-flusion_spatial2_prod/2026-01-10-UMass-flusion_spatial2_prod.csv

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(lubridate)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript plot-forecasts.R <forecast_file> [truth_file]")
}

forecast_file <- args[1]

# Default truth file location
if (length(args) >= 2) {
  truth_file <- args[2]
} else {
  # Try relative path from forecast-sandbox to Flusight-forecast-hub
  script_dir <- tryCatch({
    dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
  }, error = function(e) getwd())
  truth_file <- normalizePath(file.path(script_dir, "../../Flusight-forecast-hub/target-data/time-series.csv"))
}

cat(sprintf("Forecast file: %s\n", forecast_file))
cat(sprintf("Truth file: %s\n", truth_file))

# Load data
cat("\nLoading data...\n")

forecasts <- read_csv(forecast_file, show_col_types = FALSE) |>
  mutate(
    reference_date = as.Date(reference_date),
    target_end_date = as.Date(target_end_date),
    output_type_id = as.numeric(output_type_id)
  )

truth_data_raw <- read_csv(truth_file, show_col_types = FALSE) |>
  mutate(
    target_end_date = as.Date(target_end_date),
    as_of = as.Date(as_of)
  )

# Filter to most recent as_of date
most_recent_as_of <- max(truth_data_raw$as_of)
cat(sprintf("Using most recent as_of date: %s\n", most_recent_as_of))

truth_data <- truth_data_raw |>
  filter(as_of == most_recent_as_of)

# Get reference date from forecast
ref_date <- unique(forecasts$reference_date)[1]
cat(sprintf("Reference date: %s\n", ref_date))
cat(sprintf("Loaded %d forecast rows for %d locations\n",
            nrow(forecasts), n_distinct(forecasts$location)))

# Extract model name from forecast file path
model_name <- tryCatch({
  basename(dirname(forecast_file)) |>
    gsub("UMass-", "", x = _)
}, error = function(e) "model")
cat(sprintf("Model: %s\n", model_name))

# Determine current season year (season runs roughly Aug-Jul)
current_season_year <- if (month(ref_date) >= 8) year(ref_date) else year(ref_date) - 1
cat(sprintf("Current season: %d/%d\n", current_season_year, current_season_year + 1))

# Create historical season data aligned to current year
create_historical_seasons <- function(target_df, ref_date) {
  current_season_year <- if (month(ref_date) >= 8) year(ref_date) else year(ref_date) - 1

  # Add season info to target data
  target_with_season <- target_df |>
    mutate(
      season_year = if_else(month(target_end_date) >= 8,
                            year(target_end_date),
                            year(target_end_date) - 1),
      # Calculate days since Aug 1 of season year for alignment
      season_start = as.Date(paste0(season_year, "-08-01")),
      day_of_season = as.numeric(target_end_date - season_start)
    )

  # Get past seasons (exclude current season)
  past_seasons <- target_with_season |>
    filter(season_year < current_season_year) |>
    mutate(
      # Shift dates to align with current season
      current_season_start = as.Date(paste0(current_season_year, "-08-01")),
      aligned_date = current_season_start + day_of_season,
      season_label = paste0(season_year, "/", (season_year + 1) %% 100)
    )

  return(past_seasons)
}

historical_data <- create_historical_seasons(truth_data, ref_date)
n_historical_seasons <- n_distinct(historical_data$season_year)
cat(sprintf("Historical seasons available: %d\n", n_historical_seasons))

# Pivot forecasts to wide format
forecasts_wide <- forecasts |>
  filter(output_type == "quantile") |>
  select(location, reference_date, target_end_date, target, output_type_id, value) |>
  pivot_wider(
    names_from = output_type_id,
    values_from = value,
    names_prefix = "q"
  )

# Output directory
script_dir <- tryCatch({
  dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
}, error = function(e) getwd())

output_dir <- file.path(script_dir, "plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Get unique locations from forecasts
forecast_locations <- unique(forecasts_wide$location)

# Get location names from truth data
location_names <- truth_data |>
  select(location, location_name) |>
  distinct()

cat(sprintf("\nFound %d forecast locations\n", length(forecast_locations)))

# Function to plot a single location
plot_location <- function(loc, forecasts_df, target_df, historical_df, location_names_df, ref_date) {

  loc_name <- location_names_df |>
    filter(location == loc) |>
    pull(location_name) |>
    first()
  if (is.na(loc_name)) loc_name <- loc

  # Get forecast data for this location
  loc_forecasts <- forecasts_df |>
    filter(location == loc)

  if (nrow(loc_forecasts) == 0) return(NULL)

  # Get target type
  target_type <- unique(loc_forecasts$target)[1]

  # Get target data for this location
  loc_targets <- target_df |>
    filter(location == loc, target == target_type)

  # Date range
  min_date <- min(loc_forecasts$target_end_date) - 84  # ~12 weeks before

max_date <- as.Date(paste0(if (month(ref_date) >= 8) year(ref_date) + 1 else year(ref_date), "-04-01"))

  loc_targets_filtered <- loc_targets |>
    filter(target_end_date >= min_date, target_end_date <= max_date)

  # Prepare historical data
  loc_historical <- historical_df |>
    filter(location == loc, target == target_type) |>
    filter(aligned_date >= min_date, aligned_date <= max_date)

  n_hist_seasons <- n_distinct(loc_historical$season_label)

  # Build plot
  p <- ggplot()

  # Add historical season curves (grey, in background)
  if (nrow(loc_historical) > 0) {
    p <- p + geom_line(
      data = loc_historical,
      aes(x = aligned_date, y = observation, group = season_label),
      color = "grey70", linewidth = 0.4, alpha = 0.6
    )
  }

  # Current season observations (black)
  if (nrow(loc_targets_filtered) > 0) {
    p <- p +
      geom_line(
        data = loc_targets_filtered,
        aes(x = target_end_date, y = observation),
        color = "black", linewidth = 0.6
      ) +
      geom_point(
        data = loc_targets_filtered,
        aes(x = target_end_date, y = observation),
        color = "black", size = 1.5
      )
  }

  # 95% prediction interval
  if ("q0.025" %in% names(loc_forecasts) && "q0.975" %in% names(loc_forecasts)) {
    p <- p + geom_ribbon(
      data = loc_forecasts,
      aes(x = target_end_date, ymin = q0.025, ymax = q0.975),
      fill = "steelblue", alpha = 0.2
    )
  }

  # 50% prediction interval
  if ("q0.25" %in% names(loc_forecasts) && "q0.75" %in% names(loc_forecasts)) {
    p <- p + geom_ribbon(
      data = loc_forecasts,
      aes(x = target_end_date, ymin = q0.25, ymax = q0.75),
      fill = "steelblue", alpha = 0.3
    )
  }

  # Median forecast line
  if ("q0.5" %in% names(loc_forecasts)) {
    p <- p + geom_line(
      data = loc_forecasts,
      aes(x = target_end_date, y = q0.5),
      color = "steelblue", linewidth = 0.8
    ) +
    geom_point(
      data = loc_forecasts,
      aes(x = target_end_date, y = q0.5),
      color = "steelblue", size = 2
    )
  }

  # Reference date line
  p <- p + geom_vline(
    xintercept = as.numeric(ref_date),
    linetype = "dashed", color = "gray50", alpha = 0.5
  )

  p <- p +
    labs(
      title = sprintf("%s (%s)", loc_name, loc),
      x = "Date", y = "Hospitalizations",
      caption = sprintf("Black: current season | Grey: %d historical season%s | Blue: forecast",
                       n_hist_seasons, if(n_hist_seasons != 1) "s" else "")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      plot.caption = element_text(size = 7, color = "gray50")
    ) +
    scale_x_date(date_labels = "%b %d", date_breaks = "2 weeks") +
    scale_y_continuous(limits = c(0, NA))

  return(p)
}

# Generate multi-page PDF with grid of locations
cat("\nGenerating PDF...\n")
pdf_file <- file.path(output_dir, sprintf("%s-%s-forecasts.pdf", ref_date, model_name))

# Sort locations: US first, then states alphabetically by name
us_loc <- "US"
state_locs <- setdiff(forecast_locations, us_loc)

# Get names for sorting
state_names <- location_names |>
  filter(location %in% state_locs) |>
  arrange(location_name)

sorted_locations <- c(us_loc[us_loc %in% forecast_locations], state_names$location)

# Create plots in batches of 6 per page
plots_per_page <- 6
n_pages <- ceiling(length(sorted_locations) / plots_per_page)

pdf(pdf_file, width = 14, height = 10)

for (page in 1:n_pages) {
  start_idx <- (page - 1) * plots_per_page + 1
  end_idx <- min(page * plots_per_page, length(sorted_locations))
  page_locs <- sorted_locations[start_idx:end_idx]

  cat(sprintf("  Page %d: %s\n", page, paste(page_locs, collapse = ", ")))

  plots <- lapply(page_locs, function(loc) {
    plot_location(loc, forecasts_wide, truth_data, historical_data, location_names, ref_date)
  })

  # Remove NULL plots
  plots <- plots[!sapply(plots, is.null)]

  if (length(plots) > 0) {
    # Arrange plots in grid
    grid_plot <- cowplot::plot_grid(plotlist = plots, ncol = 3, nrow = 2)

    # Add overall title
    title <- cowplot::ggdraw() +
      cowplot::draw_label(
        sprintf("UMass-%s Forecasts | Reference Date: %s", model_name, ref_date),
        fontface = 'bold', size = 16
      )

    final_plot <- cowplot::plot_grid(title, grid_plot, ncol = 1, rel_heights = c(0.05, 0.95))
    print(final_plot)
  }
}

dev.off()
cat(sprintf("\nSaved: %s\n", pdf_file))

# Also create individual PNGs for key locations
cat("\nGenerating individual PNGs for key locations...\n")
key_locations <- c("US", "06", "48", "36", "12", "42")  # US, CA, TX, NY, FL, PA
key_locations <- key_locations[key_locations %in% forecast_locations]

for (loc in key_locations) {
  p <- plot_location(loc, forecasts_wide, truth_data, historical_data, location_names, ref_date)
  if (!is.null(p)) {
    loc_name <- location_names |>
      filter(location == loc) |>
      pull(location_name) |>
      first()
    if (is.na(loc_name)) loc_name <- loc

    png_file <- file.path(output_dir, sprintf("%s-%s-%s.png", ref_date, model_name, loc))
    ggsave(png_file, p, width = 8, height = 6, dpi = 150)
    cat(sprintf("  Saved: %s\n", png_file))
  }
}

cat("\nDone!\n")
