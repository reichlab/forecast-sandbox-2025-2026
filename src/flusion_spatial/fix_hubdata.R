#!/usr/bin/env Rscript
# Fix hubData installation by installing dependencies first

cat("Fixing hubData installation...\n\n")

# Set up user library
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") {
  user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
}
.libPaths(c(user_lib, .libPaths()))

# Install remotes if needed
if (!require("remotes", quietly = TRUE)) {
  install.packages("remotes", lib = user_lib, repos = "https://cloud.r-project.org")
}

cat("\n=== Installing hubData dependencies ===\n")

# Common dependencies that hubData needs
dependencies <- c(
  "dplyr", "tidyr", "purrr", "stringr", "tibble",
  "arrow", "jsonlite", "yaml", "cli", "rlang",
  "fs", "checkmate", "lifecycle"
)

for (pkg in dependencies) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, lib = user_lib, repos = "https://cloud.r-project.org")
  } else {
    cat(sprintf("✓ %s\n", pkg))
  }
}

cat("\n=== Installing hubUtils (required by hubData) ===\n")
if (!require("hubUtils", quietly = TRUE)) {
  remotes::install_github("hubverse-org/hubUtils", lib = user_lib, upgrade = "never")
} else {
  cat("✓ hubUtils already installed\n")
}

cat("\n=== Installing hubData ===\n")
remotes::install_github("hubverse-org/hubData", lib = user_lib, upgrade = "never", force = TRUE)

# Verify installation
cat("\n=== Verifying hubData installation ===\n")
if (require("hubData", quietly = TRUE)) {
  cat("✓ SUCCESS: hubData installed!\n")
  cat(sprintf("  Version: %s\n", packageVersion("hubData")))
} else {
  cat("✗ FAILED: hubData still not installed\n")
  cat("\nTry manually:\n")
  cat("  R\n")
  cat("  > install.packages('arrow')  # This is often the problematic dependency\n")
  cat("  > remotes::install_github('hubverse-org/hubData')\n")
  quit(status = 1)
}
