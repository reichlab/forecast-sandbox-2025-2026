#!/usr/bin/env Rscript
# Install required R packages for flusion_spatial ensemble
# Run this once on Unity before submitting jobs

cat("Installing required R packages...\n")
cat("This may take 10-20 minutes on first run.\n\n")

# Set up user library if needed
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") {
  user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
}
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE)
  cat(sprintf("Created user library: %s\n", user_lib))
}
.libPaths(c(user_lib, .libPaths()))
cat(sprintf("Installing to: %s\n\n", user_lib))

# Function to install a package if not already installed
install_if_missing <- function(pkg, source = "CRAN", repo = NULL) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing %s from %s...\n", pkg, source))
    if (source == "CRAN") {
      install.packages(pkg, lib = user_lib, repos = "https://cloud.r-project.org")
    } else if (source == "GitHub") {
      if (!require("remotes", quietly = TRUE)) {
        install.packages("remotes", lib = user_lib, repos = "https://cloud.r-project.org")
      }
      remotes::install_github(repo, lib = user_lib, upgrade = "never")
    }
  } else {
    cat(sprintf("✓ %s already installed\n", pkg))
  }
}

# Install CRAN packages
cat("\n=== Installing CRAN packages ===\n")
cran_packages <- c("dplyr", "readr", "remotes")
for (pkg in cran_packages) {
  install_if_missing(pkg, source = "CRAN")
}

# Install Hubverse packages from GitHub
cat("\n=== Installing Hubverse packages ===\n")
github_packages <- list(
  list(pkg = "hubData", repo = "hubverse-org/hubData"),
  list(pkg = "hubEnsembles", repo = "hubverse-org/hubEnsembles")
)

for (item in github_packages) {
  install_if_missing(item$pkg, source = "GitHub", repo = item$repo)
}

# Install Reich Lab packages from GitHub
cat("\n=== Installing Reich Lab packages ===\n")
reichlab_packages <- list(
  list(pkg = "idforecastutils", repo = "reichlab/idforecastutils")
)

for (item in reichlab_packages) {
  install_if_missing(item$pkg, source = "GitHub", repo = item$repo)
}

# Verify all packages are installed
cat("\n=== Verifying installation ===\n")
required_packages <- c("dplyr", "hubData", "hubEnsembles", "readr", "idforecastutils")
all_installed <- TRUE

for (pkg in required_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("✓ %s\n", pkg))
  } else {
    cat(sprintf("✗ %s - FAILED TO INSTALL\n", pkg))
    all_installed <- FALSE
  }
}

if (all_installed) {
  cat("\n========================================\n")
  cat("SUCCESS: All R packages installed!\n")
  cat("You can now run the model.\n")
  cat("========================================\n")
} else {
  cat("\n========================================\n")
  cat("ERROR: Some packages failed to install.\n")
  cat("Please check the error messages above.\n")
  cat("========================================\n")
  quit(status = 1)
}
