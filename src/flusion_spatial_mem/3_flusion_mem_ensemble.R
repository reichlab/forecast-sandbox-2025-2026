library(dplyr)
library(hubData)
library(hubEnsembles)

# Source MEM utility functions
source("mem_utils.R")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
ref_date <- as.Date(args[1])

# ============================================================================
# Main Ensemble Creation with MEM-based Adaptive Weighting
# ============================================================================

message(paste0(rep("=", 70), collapse = ""))
message(sprintf("Starting MEM-based adaptive ensemble for %s", ref_date))
message(paste0(rep("=", 70), collapse = ""))

# Connect to intermediate model outputs
hub_con <- hubData::connect_model_output(
  "intermediate-output/model-output",
  schema = arrow::schema(location = arrow::utf8())
)

# Define component models (same as flusion_spatial2_prod)
state_models_to_blend <- c("UMass-gbqr_3src_spatial2", "UMass-AR6_pooled")
us_models_to_blend <- c("UMass-gbqr_3src", "UMass-AR6_pooled")

# Load state-level forecasts
state_dat <- hub_con |>
  filter(
    reference_date == ref_date,
    model_id %in% state_models_to_blend,
    location != "US",
    horizon >= 0
  ) |>
  collect_hub()

# Load US-level forecasts
us_dat <- hub_con |>
  filter(
    reference_date == ref_date,
    model_id %in% us_models_to_blend,
    location == "US",
    horizon >= 0
  ) |>
  collect_hub()

message("\nComponent forecasts loaded:")
message(sprintf("  State locations: %d", length(unique(state_dat$location))))
message(sprintf("  US location: %s", ifelse(nrow(us_dat) > 0, "YES", "NO")))
message(sprintf("  Total forecast rows: %d", nrow(state_dat) + nrow(us_dat)))

# ============================================================================
# MEM-Based Adaptive Weighting
# ============================================================================

message(paste0("\n", paste0(rep("=", 70), collapse = "")))
message("Calculating MEM thresholds and determining epidemic phase")
message(paste0(rep("=", 70), collapse = ""))

# Calculate MEM thresholds for US (representative for all locations)
# Future enhancement: could calculate location-specific thresholds
mem_thresholds_us <- calculate_mem_thresholds("US", ref_date, num_seasons = 5)

# Determine current epidemic phase
if (mem_thresholds_us$use_adaptive) {
  current_phase <- determine_epidemic_phase("US", ref_date, mem_thresholds_us)
} else {
  current_phase <- "unknown"
  message("\nMEM thresholds not available, using equal weights")
}

# Get ensemble weights for this phase
ensemble_weights <- get_phase_weights(current_phase, "US")

# ============================================================================
# Create Ensemble Forecast
# ============================================================================

message(paste0("\n", paste0(rep("=", 70), collapse = "")))
message("Creating weighted ensemble forecast")
message(paste0(rep("=", 70), collapse = ""))

# Apply weighted linear pool ensemble to state and US forecasts
ens_state <- weighted_linear_pool(
  state_dat,
  ensemble_weights,
  ar6_model_id = "UMass-AR6_pooled",
  output_model_id = "UMass-flusion_spatial_mem"
)

ens_us <- weighted_linear_pool(
  us_dat,
  ensemble_weights,
  ar6_model_id = "UMass-AR6_pooled",
  output_model_id = "UMass-flusion_spatial_mem"
)

# Combine state and US forecasts
ens_model <- bind_rows(ens_state, ens_us)

message(sprintf("\nEnsemble forecast created:"))
message(sprintf("  Total rows: %d", nrow(ens_model)))
message(sprintf("  Locations: %d", length(unique(ens_model$location))))
message(sprintf("  Horizons: %s", paste(sort(unique(ens_model$horizon)), collapse = ", ")))

# ============================================================================
# Save Outputs
# ============================================================================

# Save main ensemble forecast
output_dir <- "../../model-output/UMass-flusion_spatial_mem"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message(sprintf("Created output directory: %s", output_dir))
}

output_file <- file.path(
  output_dir,
  paste0(ref_date, "-UMass-flusion_spatial_mem.csv")
)

utils::write.csv(
  ens_model |> dplyr::select(-model_id),
  file = output_file,
  row.names = FALSE
)

message(sprintf("\nMain forecast saved to: %s", output_file))

# Save diagnostic information
diagnostic_file <- file.path(
  output_dir,
  paste0(ref_date, "-diagnostic-info.csv")
)

diagnostic_info <- data.frame(
  reference_date = ref_date,
  epidemic_phase = current_phase,
  ar6_weight = ensemble_weights["ar6"],
  gbqr_weight = ensemble_weights["gbqr"],
  mem_available = mem_thresholds_us$use_adaptive,
  epidemic_threshold = ifelse(
    mem_thresholds_us$use_adaptive,
    mem_thresholds_us$epidemic_threshold,
    NA_real_
  ),
  intensity_low = ifelse(
    mem_thresholds_us$use_adaptive && length(mem_thresholds_us$intensity_thresholds) >= 1,
    mem_thresholds_us$intensity_thresholds[1],
    NA_real_
  ),
  intensity_medium = ifelse(
    mem_thresholds_us$use_adaptive && length(mem_thresholds_us$intensity_thresholds) >= 2,
    mem_thresholds_us$intensity_thresholds[2],
    NA_real_
  ),
  intensity_high = ifelse(
    mem_thresholds_us$use_adaptive && length(mem_thresholds_us$intensity_thresholds) >= 3,
    mem_thresholds_us$intensity_thresholds[3],
    NA_real_
  ),
  timestamp = Sys.time()
)

utils::write.csv(
  diagnostic_info,
  file = diagnostic_file,
  row.names = FALSE
)

message(sprintf("Diagnostic info saved to: %s", diagnostic_file))

# ============================================================================
# Summary
# ============================================================================

message(paste0("\n", paste0(rep("=", 70), collapse = "")))
message("ENSEMBLE CREATION COMPLETE")
message(paste0(rep("=", 70), collapse = ""))
message(sprintf("Reference date: %s", ref_date))
message(sprintf("Epidemic phase: %s", current_phase))
message(sprintf("Ensemble weights: AR6=%.2f, GBQR=%.2f",
                ensemble_weights["ar6"], ensemble_weights["gbqr"]))
message(sprintf("Forecast file: %s", output_file))
message(sprintf("Diagnostic file: %s", diagnostic_file))
message(paste0(rep("=", 70), collapse = ""))
