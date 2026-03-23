#!/usr/bin/env Rscript
# Test script to verify MEM setup and data availability
#
# This script tests:
# 1. R package installation (mem, fiphde, etc.)
# 2. Data fetching from fiphde
# 3. MEM threshold calculation
# 4. Phase detection
#
# Usage: Rscript test_mem_setup.R

message(paste0(rep("=", 70), collapse = ""))
message("MEM Setup Verification Script")
message(paste0(rep("=", 70), collapse = ""))

# Test 1: Check required packages
message("\n[Test 1] Checking required R packages...")

required_packages <- c("dplyr", "tidyr", "lubridate", "MMWRweek", "mem", "fiphde",
                       "hubData", "hubEnsembles")

missing_packages <- c()
for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("  ✓ %s", pkg))
  } else {
    message(sprintf("  ✗ %s (MISSING)", pkg))
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) > 0) {
  message("\n⚠ Missing packages detected:")
  message("Run the following to install missing packages:")
  message("  renv::restore()  # For CRAN packages")
  if ("fiphde" %in% missing_packages) {
    message("  remotes::install_github('signaturescience/fiphde')  # For fiphde")
  }
  stop("Please install missing packages before proceeding")
}

message("\n✓ All required packages are installed")

# Test 2: Source MEM utilities
message("\n[Test 2] Loading MEM utility functions...")

tryCatch({
  source("mem_utils.R")
  message("✓ mem_utils.R loaded successfully")
}, error = function(e) {
  stop(sprintf("✗ Failed to load mem_utils.R: %s", e$message))
})

# Test 3: Test data fetching
message("\n[Test 3] Testing NHSN data fetching with fiphde...")

test_date <- as.Date("2024-12-07")
test_location <- "US"

tryCatch({
  hist_data <- get_historical_hosp_data(
    location = test_location,
    ref_date = test_date,
    num_seasons = 3  # Use fewer seasons for quick test
  )

  if (is.null(hist_data)) {
    message("⚠ No data returned - this may be expected if fiphde is not configured")
    message("   or if historical data is not available for the test date")
  } else {
    message(sprintf("✓ Successfully fetched %d rows of historical data", nrow(hist_data)))
    message(sprintf("  Date range: %s to %s", min(hist_data$date), max(hist_data$date)))
    message(sprintf("  Seasons: %s", paste(unique(hist_data$season), collapse = ", ")))
  }
}, error = function(e) {
  message(sprintf("⚠ Data fetch failed: %s", e$message))
  message("   This may be expected if fiphde is not yet configured")
  hist_data <- NULL
})

# Test 4: Test MEM data preparation
if (!is.null(hist_data)) {
  message("\n[Test 4] Testing MEM data preparation...")

  tryCatch({
    mem_matrix <- prepare_mem_data(hist_data)

    if (is.null(mem_matrix)) {
      message("⚠ MEM data preparation returned NULL")
    } else {
      message(sprintf("✓ MEM matrix created: %d seasons x %d weeks",
                      nrow(mem_matrix), ncol(mem_matrix)))
      message("  First few values:")
      print(head(mem_matrix[, 1:min(5, ncol(mem_matrix))]))
    }
  }, error = function(e) {
    message(sprintf("⚠ MEM data preparation failed: %s", e$message))
    mem_matrix <- NULL
  })
} else {
  message("\n[Test 4] Skipping MEM data preparation (no historical data)")
  mem_matrix <- NULL
}

# Test 5: Test MEM threshold calculation
if (!is.null(mem_matrix)) {
  message("\n[Test 5] Testing MEM threshold calculation...")

  tryCatch({
    mem_thresholds <- calculate_mem_thresholds(
      location = test_location,
      ref_date = test_date,
      num_seasons = 3
    )

    if (mem_thresholds$use_adaptive) {
      message("✓ MEM thresholds calculated successfully:")
      message(sprintf("  Epidemic threshold: %.2f", mem_thresholds$epidemic_threshold))
      message(sprintf("  Intensity thresholds:"))
      message(sprintf("    Low: %.2f", mem_thresholds$intensity_thresholds[1]))
      message(sprintf("    Medium: %.2f", mem_thresholds$intensity_thresholds[2]))
      message(sprintf("    High: %.2f", mem_thresholds$intensity_thresholds[3]))
    } else {
      message("⚠ MEM thresholds calculation failed - adaptive weighting not available")
      message(sprintf("  Reason: %s", mem_thresholds$reason))
    }
  }, error = function(e) {
    message(sprintf("⚠ MEM threshold calculation failed: %s", e$message))
    mem_thresholds <- list(use_adaptive = FALSE)
  })
} else {
  message("\n[Test 5] Skipping MEM threshold calculation (no MEM matrix)")
  mem_thresholds <- list(use_adaptive = FALSE)
}

# Test 6: Test phase determination
if (mem_thresholds$use_adaptive) {
  message("\n[Test 6] Testing epidemic phase determination...")

  tryCatch({
    current_phase <- determine_epidemic_phase(
      location = test_location,
      ref_date = test_date,
      mem_thresholds = mem_thresholds
    )

    message(sprintf("✓ Epidemic phase determined: %s", current_phase))

    # Get weights for this phase
    weights <- get_phase_weights(current_phase, test_location)
    message(sprintf("  Ensemble weights: AR6=%.2f, GBQR=%.2f",
                    weights["ar6"], weights["gbqr"]))

  }, error = function(e) {
    message(sprintf("⚠ Phase determination failed: %s", e$message))
  })
} else {
  message("\n[Test 6] Skipping phase determination (MEM thresholds not available)")
}

# Summary
message(paste0("\n", paste0(rep("=", 70), collapse = "")))
message("SETUP VERIFICATION SUMMARY")
message(paste0(rep("=", 70), collapse = ""))

if (length(missing_packages) == 0) {
  message("✓ All packages installed")
} else {
  message("✗ Some packages missing")
}

if (!is.null(hist_data)) {
  message("✓ Data fetching works")
} else {
  message("⚠ Data fetching needs configuration")
}

if (!is.null(mem_matrix)) {
  message("✓ MEM data preparation works")
} else {
  message("⚠ MEM data preparation not tested")
}

if (exists("mem_thresholds") && mem_thresholds$use_adaptive) {
  message("✓ MEM threshold calculation works")
  message("✓ Ready for adaptive ensemble weighting!")
} else {
  message("⚠ MEM threshold calculation needs data/configuration")
  message("  Model will fall back to equal weights until configured")
}

message("\nNext steps:")
message("1. If fiphde data fetch failed, install: remotes::install_github('signaturescience/fiphde')")
message("2. Run a test forecast: python main.py --today_date=2024-12-07 --short_run")
message("3. Check diagnostic output in: ../../model-output/UMass-flusion_spatial_mem/")
message(paste0(rep("=", 70), collapse = ""))
