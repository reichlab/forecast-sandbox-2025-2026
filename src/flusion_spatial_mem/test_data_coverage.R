#!/usr/bin/env Rscript
# Test script to examine NHSN data coverage for MEM

library(dplyr)
library(lubridate)
library(fiphde)

cat("=== Testing NHSN Data Coverage ===\n\n")

# Fetch all available NHSN data
cat("Fetching NHSN data...\n")
nhsn_data <- get_nhsn_weekly()

cat(sprintf("Total rows fetched: %d\n", nrow(nhsn_data)))
cat(sprintf("Date range: %s to %s\n",
            min(as.Date(nhsn_data$week_end)),
            max(as.Date(nhsn_data$week_end))))
cat("\n")

# Filter to US
us_data <- nhsn_data %>%
  filter(abbreviation == "USA") %>%
  mutate(
    date = as.Date(week_end),
    season = ifelse(
      month(date) >= 10,
      year(date),
      year(date) - 1
    )
  )

cat(sprintf("US data rows: %d\n", nrow(us_data)))
cat(sprintf("US date range: %s to %s\n", min(us_data$date), max(us_data$date)))
cat("\n")

# Count weeks per season
season_summary <- us_data %>%
  group_by(season) %>%
  summarize(
    n_weeks = n(),
    start_date = min(date),
    end_date = max(date),
    .groups = "drop"
  ) %>%
  arrange(season)

cat("Weeks per season:\n")
print(season_summary, n = Inf)
cat("\n")

# Check which seasons have >= 30 weeks
complete_seasons <- season_summary %>%
  filter(n_weeks >= 30)

cat(sprintf("Seasons with >= 30 weeks: %d\n", nrow(complete_seasons)))
cat("Complete seasons:\n")
print(complete_seasons, n = Inf)
cat("\n")

# Check for gaps in reporting
gaps <- us_data %>%
  arrange(date) %>%
  mutate(
    days_since_prev = as.numeric(date - lag(date))
  ) %>%
  filter(!is.na(days_since_prev) & days_since_prev > 14) %>%
  select(date, days_since_prev)

if (nrow(gaps) > 0) {
  cat("Reporting gaps (>14 days between reports):\n")
  print(gaps, n = Inf)
} else {
  cat("No major reporting gaps detected\n")
}
cat("\n")

# Check the specific period we're interested in: 2019-2024
ref_date <- as.Date("2024-12-07")
lookback_start <- ref_date - years(5) - days(60)

period_data <- us_data %>%
  filter(date >= lookback_start & date <= ref_date)

cat(sprintf("Data for period %s to %s:\n", lookback_start, ref_date))
cat(sprintf("  Rows: %d\n", nrow(period_data)))

period_seasons <- period_data %>%
  group_by(season) %>%
  summarize(n_weeks = n(), .groups = "drop") %>%
  arrange(season)

cat("  Seasons in this period:\n")
print(period_seasons, n = Inf)
