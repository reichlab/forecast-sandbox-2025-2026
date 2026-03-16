# MEM Utilities for Flusion Spatial MEM Ensemble
#
# This file contains helper functions for:
# - Fetching historical NHSN hospitalization data
# - Formatting data for the MEM package
# - Calculating MEM thresholds
# - Determining epidemic phase and intensity
# - Phase-aware weight selection

library(dplyr)
library(tidyr)
library(lubridate)
library(MMWRweek)
library(mem)

# Load fiphde if available, otherwise use fallback
FIPHDE_AVAILABLE <- requireNamespace("fiphde", quietly = TRUE)

if (FIPHDE_AVAILABLE) {
  message("fiphde package loaded successfully")
} else {
  message("WARNING: fiphde package not available. Install with:")
  message("  remotes::install_github('signaturescience/fiphde')")
  message("Using fallback data methods.")
}


#' Get historical NHSN hospitalization data for MEM calculation
#'
#' Fetches historical weekly flu hospitalization data from NHSN for a given location.
#' Uses the fiphde package which provides access to extended NHSN time series data.
#'
#' @param location FIPS code for location (e.g., "US", "25" for Massachusetts)
#' @param ref_date Reference date for forecast (Date object)
#' @param num_seasons Number of historical seasons to use (default 5)
#' @param min_weeks_per_season Minimum weeks required per season (default 30)
#'
#' @return Data frame with columns: location, date, year_week, hosp, season
#'   Returns NULL if insufficient data available
#'
#' @examples
#' hist_data <- get_historical_hosp_data("US", as.Date("2024-12-07"), num_seasons = 5)
get_historical_hosp_data <- function(location, ref_date, num_seasons = 5, min_weeks_per_season = 30) {

  # Determine current season to exclude it (MEM needs complete past seasons only)
  current_month <- lubridate::month(ref_date)
  current_year <- lubridate::year(ref_date)
  current_season <- ifelse(current_month >= 10, current_year, current_year - 1)

  # Calculate end date: end of the most recent COMPLETE season
  # Complete season ends in September (week 39)
  # So we want data up to Sep 30 of the year the current season started
  end_date <- as.Date(paste0(current_season, "-09-30"))

  # Calculate start date: go back num_seasons + 1 years to ensure we get enough complete seasons
  # Add buffer to capture full seasons
  start_date <- end_date - lubridate::years(num_seasons + 1)

  message(sprintf("Fetching historical data for location %s from %s to %s (excluding current incomplete season %d)",
                  location, start_date, end_date, current_season))

  if (!FIPHDE_AVAILABLE) {
    message("fiphde not available - cannot fetch historical data")
    return(NULL)
  }

  tryCatch({
    # Fetch NHSN weekly data using fiphde
    # The fiphde package provides get_nhsn_weekly() which retrieves all available data
    # We then filter to the date range we need
    nhsn_data <- fiphde::get_nhsn_weekly()

    # Convert week_end to date if needed and filter to date range
    if ("week_end" %in% names(nhsn_data)) {
      nhsn_data <- nhsn_data %>%
        mutate(date = as.Date(week_end)) %>%
        filter(date >= start_date & date <= end_date)
    }

    # Filter to the specified location
    # fiphde uses 'abbreviation' for location codes
    if (location == "US") {
      loc_data <- nhsn_data %>%
        filter(abbreviation == "USA")
    } else {
      # Convert FIPS to state abbreviation if needed
      # For now, try both the location code and as-is
      loc_data <- nhsn_data %>%
        filter(abbreviation == !!location)
    }

    # Rename columns to match expected format
    if ("flu.admits" %in% names(loc_data)) {
      loc_data <- loc_data %>%
        rename(hosp = flu.admits,
               location = abbreviation)
    }

    if (nrow(loc_data) == 0) {
      message(sprintf("No data available for location %s", location))
      return(NULL)
    }

    # Add MMWR week and season information
    loc_data <- loc_data %>%
      mutate(
        mmwr_week = MMWRweek::MMWRweek(date)$MMWRweek,
        mmwr_year = MMWRweek::MMWRweek(date)$MMWRyear,
        # Flu season: Oct-May, labeled by year when it starts
        # Season 2023 = Oct 2023 - May 2024
        season = ifelse(
          lubridate::month(date) >= 10,
          lubridate::year(date),
          lubridate::year(date) - 1
        ),
        year_week = sprintf("%04d-W%02d", mmwr_year, mmwr_week)
      ) %>%
      arrange(date)

    # Check we have enough seasons
    seasons_available <- loc_data %>%
      group_by(season) %>%
      summarize(n_weeks = n(), .groups = "drop") %>%
      filter(n_weeks >= min_weeks_per_season)

    if (nrow(seasons_available) < num_seasons) {
      message(sprintf(
        "Insufficient seasons available for location %s: found %d, need %d",
        location, nrow(seasons_available), num_seasons
      ))
      return(NULL)
    }

    # Keep only complete seasons
    complete_seasons <- seasons_available$season
    loc_data <- loc_data %>%
      filter(season %in% complete_seasons)

    message(sprintf("Successfully fetched %d weeks across %d seasons",
                    nrow(loc_data), length(complete_seasons)))

    return(loc_data)

  }, error = function(e) {
    message(sprintf("Error fetching data for location %s: %s", location, e$message))
    return(NULL)
  })
}


#' Prepare historical data for MEM in required data frame format
#'
#' MEM expects data as a data frame where:
#' - Each row is a week number (1-52)
#' - Each column is a season
#' - Values are hospitalization rates or counts
#'
#' @param hosp_data Historical hospitalization data from get_historical_hosp_data()
#' @param value_col Column name containing hospitalization values (default "hosp")
#' @param max_weeks Maximum number of weeks per season (default 52)
#'
#' @return Data frame in MEM format (weeks x seasons), or NULL if data insufficient
#'
#' @examples
#' hist_data <- get_historical_hosp_data("US", as.Date("2024-12-07"))
#' mem_df <- prepare_mem_data(hist_data)
prepare_mem_data <- function(hosp_data, value_col = "hosp", max_weeks = 52) {

  if (is.null(hosp_data) || nrow(hosp_data) == 0) {
    return(NULL)
  }

  # Assign week-within-season for each observation
  # Week 1 starts in October (first week of flu season)
  hosp_data <- hosp_data %>%
    group_by(season) %>%
    arrange(date) %>%
    mutate(
      season_week = row_number()
    ) %>%
    ungroup()

  # Check that we don't exceed max_weeks
  if (max(hosp_data$season_week) > max_weeks) {
    message(sprintf("Warning: some seasons have >%d weeks, truncating", max_weeks))
    hosp_data <- hosp_data %>%
      filter(season_week <= max_weeks)
  }

  # Reshape to wide format: weeks as rows, seasons as columns
  # This is the format required by mem::memmodel()
  mem_df <- hosp_data %>%
    select(season, season_week, value = !!sym(value_col)) %>%
    pivot_wider(
      names_from = season,
      values_from = value,
      names_prefix = "season_"
    ) %>%
    arrange(season_week) %>%
    select(-season_week)  # Remove week column, keep only season columns

  # Replace NA with 0 (if any weeks are missing)
  mem_df[is.na(mem_df)] <- 0

  message(sprintf("Prepared MEM data frame: %d weeks x %d seasons",
                  nrow(mem_df), ncol(mem_df)))

  return(mem_df)
}


#' Calculate MEM thresholds for a location
#'
#' Calculates epidemic and intensity thresholds using the Moving Epidemic Method.
#' Returns threshold values and the full MEM model object.
#'
#' @param location FIPS code
#' @param ref_date Reference date for forecast
#' @param num_seasons Number of historical seasons to use (default 5)
#' @param i.type Intensity threshold type (default "geometric")
#' @param i.level Intensity levels confidence (default c(0.40, 0.90, 0.975))
#'
#' @return List with components:
#'   - use_adaptive: Boolean, TRUE if MEM was successfully calculated
#'   - epidemic_threshold: Numeric, epidemic start threshold
#'   - intensity_thresholds: Named vector of intensity level thresholds
#'   - mem_model: Full mem model object (if successful)
#'   - phase: Current phase (if determinable)
#'
#' @examples
#' thresholds <- calculate_mem_thresholds("US", as.Date("2024-12-07"))
calculate_mem_thresholds <- function(location, ref_date, num_seasons = 5,
                                     i.type.intensity = 6,
                                     i.level.intensity = c(0.40, 0.90, 0.975)) {

  message(sprintf("\nCalculating MEM thresholds for location %s", location))

  # Get historical data - try with requested num_seasons first
  hist_data <- get_historical_hosp_data(location, ref_date, num_seasons)

  # If we don't have enough seasons, try with fewer (minimum 4)
  if (is.null(hist_data) && num_seasons > 4) {
    message(sprintf("Trying with %d seasons instead...", num_seasons - 1))
    hist_data <- get_historical_hosp_data(location, ref_date, num_seasons - 1)
  }

  if (is.null(hist_data)) {
    message("No historical data available - using equal weights (no adaptation)")
    return(list(
      use_adaptive = FALSE,
      phase = "unknown",
      reason = "no_historical_data"
    ))
  }

  # Log if using fewer seasons than requested
  actual_seasons <- length(unique(hist_data$season))
  if (actual_seasons < num_seasons) {
    message(sprintf("Note: Using %d seasons (requested %d, limited by data availability)",
                    actual_seasons, num_seasons))
  }

  # Prepare data for MEM
  mem_data <- prepare_mem_data(hist_data)

  if (is.null(mem_data)) {
    message("Could not prepare data for MEM - using equal weights")
    return(list(
      use_adaptive = FALSE,
      phase = "unknown",
      reason = "data_preparation_failed"
    ))
  }

  # Calculate MEM model
  tryCatch({
    message("Calculating MEM model...")

    # Run memmodel with specified parameters
    # Use defaults for most parameters for simplicity
    mem_model <- memmodel(
      i.data = mem_data,
      i.level.intensity = i.level.intensity
    )

    # Extract thresholds
    epidemic_threshold <- mem_model$epidemic.thresholds[1]
    intensity_thresholds <- mem_model$intensity.thresholds

    message(sprintf("MEM thresholds calculated successfully:"))
    message(sprintf("  Epidemic threshold: %.2f", epidemic_threshold))
    message(sprintf("  Low intensity: %.2f", intensity_thresholds[1]))
    message(sprintf("  Medium intensity: %.2f", intensity_thresholds[2]))
    message(sprintf("  High intensity: %.2f", intensity_thresholds[3]))

    return(list(
      use_adaptive = TRUE,
      epidemic_threshold = epidemic_threshold,
      intensity_thresholds = intensity_thresholds,
      mem_model = mem_model,
      historical_data = hist_data
    ))

  }, error = function(e) {
    message(sprintf("Error calculating MEM: %s", e$message))
    message("Falling back to equal weights")
    return(list(
      use_adaptive = FALSE,
      phase = "unknown",
      reason = paste0("mem_calculation_failed: ", e$message)
    ))
  })
}


#' Determine current epidemic phase based on recent data and MEM thresholds
#'
#' Compares recent hospitalization activity against MEM thresholds to classify
#' the current epidemic phase and intensity level.
#'
#' @param location FIPS code
#' @param ref_date Reference date for forecast
#' @param mem_thresholds MEM threshold object from calculate_mem_thresholds()
#' @param lookback_weeks Number of recent weeks to examine (default 4)
#'
#' @return Character string: "baseline", "low", "medium", "high", "very_high", or "unknown"
#'
#' @examples
#' thresholds <- calculate_mem_thresholds("US", as.Date("2024-12-07"))
#' phase <- determine_epidemic_phase("US", as.Date("2024-12-07"), thresholds)
determine_epidemic_phase <- function(location, ref_date, mem_thresholds, lookback_weeks = 4) {

  if (!mem_thresholds$use_adaptive) {
    return("unknown")
  }

  message(sprintf("\nDetermining epidemic phase for location %s", location))

  # Get recent hospitalization data
  recent_start <- ref_date - lubridate::weeks(lookback_weeks)

  tryCatch({
    # Fetch recent data
    if (!FIPHDE_AVAILABLE) {
      message("fiphde not available - cannot determine phase")
      return("unknown")
    }

    # Get all NHSN data and filter to recent period
    recent_data <- fiphde::get_nhsn_weekly()

    if ("week_end" %in% names(recent_data)) {
      recent_data <- recent_data %>%
        mutate(date = as.Date(week_end)) %>%
        filter(date >= recent_start & date <= ref_date)
    }

    # Filter to location (use "USA" for US, state abbreviation otherwise)
    if (location == "US") {
      recent_data <- recent_data %>%
        filter(abbreviation == "USA")
    } else {
      recent_data <- recent_data %>%
        filter(abbreviation == !!location)
    }

    # Rename flu.admits to hosp
    if ("flu.admits" %in% names(recent_data)) {
      recent_data <- recent_data %>%
        rename(hosp = flu.admits) %>%
        arrange(date)
    }

    if (nrow(recent_data) == 0) {
      message("No recent data available")
      return("unknown")
    }

    # Use the most recent week's hospitalization count/rate
    current_value <- recent_data %>%
      slice_tail(n = 1) %>%
      pull(hosp)

    # Alternatively, use average over last 2 weeks for stability
    current_value_avg <- recent_data %>%
      slice_tail(n = 2) %>%
      summarize(mean_hosp = mean(hosp, na.rm = TRUE)) %>%
      pull(mean_hosp)

    # Use the average for more stability
    current_value <- current_value_avg

    message(sprintf("Current hospitalization value: %.2f", current_value))

    # Compare against thresholds
    epidemic_thresh <- mem_thresholds$epidemic_threshold
    intensity_thresh <- mem_thresholds$intensity_thresholds

    if (current_value < epidemic_thresh) {
      phase <- "baseline"
    } else if (current_value < intensity_thresh[1]) {
      phase <- "low"
    } else if (current_value < intensity_thresh[2]) {
      phase <- "medium"
    } else if (current_value < intensity_thresh[3]) {
      phase <- "high"
    } else {
      phase <- "very_high"
    }

    message(sprintf("Determined phase: %s", phase))
    message(sprintf("  (value: %.2f, epidemic threshold: %.2f)",
                    current_value, epidemic_thresh))

    return(phase)

  }, error = function(e) {
    message(sprintf("Error determining phase: %s", e$message))
    return("unknown")
  })
}


#' Get ensemble weights based on epidemic phase
#'
#' Returns model-specific weights for ensemble combination based on the
#' current epidemic phase. Weights are based on empirical performance
#' patterns observed in different epidemic phases.
#'
#' @param phase Epidemic phase: "baseline", "low", "medium", "high", "very_high", or "unknown"
#' @param location Location code (for potential location-specific weights in future)
#'
#' @return Named vector of weights: c(ar6 = weight, gbqr = weight)
#'   Weights sum to 1.0
#'
#' @examples
#' weights <- get_phase_weights("medium", "US")
#' # Returns: c(ar6 = 0.4, gbqr = 0.6)
get_phase_weights <- function(phase, location = NULL) {

  # Define phase-specific weights
  # These are initial estimates based on expected model performance characteristics:
  # - AR models perform better during stable/baseline periods
  # - GBQR models capture rapid changes and spatial patterns better during epidemics

  weights <- list(
    baseline   = c(ar6 = 0.50, gbqr = 0.50),  # Equal for stable baseline
    low        = c(ar6 = 0.45, gbqr = 0.55),  # Slight GBQR preference as epidemic starts
    medium     = c(ar6 = 0.40, gbqr = 0.60),  # More GBQR during growth
    high       = c(ar6 = 0.35, gbqr = 0.65),  # Strong GBQR preference at high intensity
    very_high  = c(ar6 = 0.30, gbqr = 0.70),  # Maximum GBQR weight at peak
    unknown    = c(ar6 = 0.50, gbqr = 0.50)   # Equal weights as safe fallback
  )

  # Return weights for the given phase
  if (phase %in% names(weights)) {
    selected_weights <- weights[[phase]]
  } else {
    message(sprintf("Unknown phase '%s', using equal weights", phase))
    selected_weights <- weights$unknown
  }

  message(sprintf("Phase '%s' weights: AR6=%.2f, GBQR=%.2f",
                  phase, selected_weights["ar6"], selected_weights["gbqr"]))

  return(selected_weights)
}


#' Create weighted linear pool ensemble
#'
#' Combines component model forecasts using a weighted linear pool (weighted average).
#' Applies weights at the quantile level.
#'
#' @param data Component forecast data with columns: model_id, value, and task columns
#' @param weights Named vector of weights (e.g., c(ar6 = 0.4, gbqr = 0.6))
#' @param ar6_model_id Model ID for AR6 component (default "UMass-AR6_pooled")
#' @param output_model_id Output model ID for ensemble (default "UMass-flusion_spatial_mem")
#'
#' @return Ensemble forecast data frame
#'
#' @examples
#' weights <- c(ar6 = 0.4, gbqr = 0.6)
#' ensemble <- weighted_linear_pool(forecasts, weights)
weighted_linear_pool <- function(data, weights,
                                 ar6_model_id = "UMass-AR6_pooled",
                                 output_model_id = "UMass-flusion_spatial_mem") {

  message(sprintf("Creating weighted ensemble with AR6 weight=%.2f, GBQR weight=%.2f",
                  weights["ar6"], weights["gbqr"]))

  # Group by forecast task and calculate weighted quantiles
  ensemble <- data %>%
    group_by(reference_date, target, horizon, location, target_end_date,
             output_type, output_type_id) %>%
    summarize(
      value = {
        # Get values for each model type
        ar6_val <- value[model_id == ar6_model_id]
        gbqr_val <- value[model_id != ar6_model_id]

        # Ensure we have exactly one value per model
        if (length(ar6_val) != 1 || length(gbqr_val) != 1) {
          # Fall back to median if we have unexpected model counts
          median(value, na.rm = TRUE)
        } else {
          # Weighted combination
          weights["ar6"] * ar6_val + weights["gbqr"] * gbqr_val
        }
      },
      .groups = "drop"
    ) %>%
    mutate(model_id = output_model_id)

  message(sprintf("Ensemble created: %d forecast rows", nrow(ensemble)))

  return(ensemble)
}
