# Shared helpers for the peak-model EDA documents (peak-timing-size-eda.qmd, relative-size-eda.qmd).
# Descriptive analysis and plotting only; model code lives in Python (idmodels, relative_size_*.py).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(arrow)
  library(scales)
  library(knitr)
})

# ---------------------------------------------------------------------------------------------------------------
# constants

WINDOW <- c(10, 43)        # season weeks over which the peak is defined (week 1 = MMWR week 31)
REPLAY_START <- 5
LOG_EPS <- c(nhsn = 0.01, flusurvnet = 0.01, ilinet = 0.001)
DROP_ILINET_LOCATIONS <- c("72", "78")  # Puerto Rico, US Virgin Islands: all-zero ILINet x % positive series

SOURCE_GROUPS <- c("ILINet states", "ILINet national + regions", "FluSurv-NET")
SOURCE_COLORS <- c("ILINet states" = "#2a78d6", "ILINet national + regions" = "#eb6834", "FluSurv-NET" = "#1baf7a")
ACCENT <- "#eb6834"
INK <- "#0b0b0b"
INK2 <- "#52514e"
GRID <- "#e4e3df"
BLUE <- "#2a78d6"
Z_LIM <- c(-0.3, 9.5)

FIPS <- c(
  "01" = "AL", "02" = "AK", "04" = "AZ", "05" = "AR", "06" = "CA", "08" = "CO", "09" = "CT", "10" = "DE",
  "11" = "DC", "12" = "FL", "13" = "GA", "15" = "HI", "16" = "ID", "17" = "IL", "18" = "IN", "19" = "IA",
  "20" = "KS", "21" = "KY", "22" = "LA", "23" = "ME", "24" = "MD", "25" = "MA", "26" = "MI", "27" = "MN",
  "28" = "MS", "29" = "MO", "30" = "MT", "31" = "NE", "32" = "NV", "33" = "NH", "34" = "NJ", "35" = "NM",
  "36" = "NY", "37" = "NC", "38" = "ND", "39" = "OH", "40" = "OK", "41" = "OR", "42" = "PA", "44" = "RI",
  "45" = "SC", "46" = "SD", "47" = "TN", "48" = "TX", "49" = "UT", "50" = "VT", "51" = "VA", "53" = "WA",
  "54" = "WV", "55" = "WI", "56" = "WY", "72" = "PR", "78" = "VI", "US" = "US"
)
loc_label <- function(location) {
  out <- unname(FIPS[location])
  ifelse(is.na(out), sub("Region ", "R", location), out)
}

source_group <- function(source, agg_level) {
  case_when(
    source == "ilinet" & agg_level == "state" ~ "ILINet states",
    source == "ilinet" ~ "ILINet national + regions",
    source == "flusurvnet" ~ "FluSurv-NET",
    TRUE ~ "other"
  )
}

# ---------------------------------------------------------------------------------------------------------------
# season weeks and dates

# Saturday ending season week w of a season, e.g. season_week_date("2017/18", 10); week 1 is MMWR week 31
season_week_date <- function(season, w) {
  y0 <- as.integer(substr(season, 1, 4))
  # MMWR week 31 ends on the 31st Saturday counted from the Saturday ending MMWR week 1 (first Saturday >= Jan 4)
  j4 <- as.Date(paste0(y0, "-01-04"))
  w1 <- j4 + ((6 - as.integer(format(j4, "%u"))) %% 7)
  w1 + 7 * (30 + w - 1)
}

# approximate month labels for season-week axes (dates shift by a few days between seasons)
week_breaks <- c(5, 10, 15, 20, 25, 30, 35, 40)
week_labels <- function(breaks = week_breaks) {
  paste0(breaks, "\n", format(season_week_date("2017/18", breaks), "%b"))
}
scale_x_season_week <- function(breaks = week_breaks, limits = c(REPLAY_START - 0.5, WINDOW[2] + 0.5), ...) {
  scale_x_continuous("Season week (approx. month)", breaks = breaks, labels = week_labels(breaks),
                     limits = limits, expand = expansion(0), ...)
}
pre_window <- function(log_y = FALSE) {
  # shade the weeks before the window opens (use log_y = TRUE on log-scale y axes, where -Inf is not allowed)
  annotate("rect", xmin = -Inf, xmax = WINDOW[1] - 0.5, ymin = if (log_y) 0 else -Inf, ymax = Inf, fill = "#f1f0ec")
}
plain_log_labels <- function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)

# ---------------------------------------------------------------------------------------------------------------
# theme

theme_eda <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = GRID, linewidth = 0.4),
      axis.text = element_text(color = INK2),
      axis.title = element_text(color = INK2),
      plot.title = element_text(face = "bold", size = rel(1.05), color = INK),
      plot.title.position = "plot",
      strip.text = element_text(face = "bold", hjust = 0, color = INK),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}
theme_set(theme_eda())

pct <- function(x, acc = 1) scales::percent(x, accuracy = acc)
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)

# ---------------------------------------------------------------------------------------------------------------
# data

# Long-format ILINet (x % positive) and FluSurv-NET data, final values, as cached by validate_ilinet.py.
load_surveillance <- function(path, drop_ilinet = DROP_ILINET_LOCATIONS) {
  read_parquet(path) |>
    mutate(wk_end_date = as.Date(wk_end_date, tz = "UTC")) |>
    filter(!(source == "ilinet" & location %in% drop_ilinet)) |>
    mutate(group = factor(source_group(source, agg_level), levels = SOURCE_GROUPS))
}
