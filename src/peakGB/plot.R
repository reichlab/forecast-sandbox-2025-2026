# PDF of one UMass-peakGB forecast: one panel per location (6 per page), drawn with plot_peak_forecast() from
# analysis/peak-models/R/slides_helpers.R in the forecast sandbox (copied here so this folder is self-contained).
#
# Usage (from this directory, after main.py): Rscript plot.R <reference_date>
# Reads output/model-output/UMass-peakGB/<reference_date>-UMass-peakGB.csv and the NHSN release of the Wednesday
# before the reference date; writes output/plots/<reference_date>-UMass-peakGB.pdf.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

MODEL <- "UMass-peakGB"
MODEL_COLORS <- c("UMass-peakGB" = "#2a78d6")
INK <- "#0b0b0b"
INK2 <- "#52514e"
PER_PAGE <- 6  # 3 columns x 2 rows

theme_set(
  theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "#e4e2dc", linewidth = 0.3),
          plot.title = element_text(face = "bold", size = 10), plot.title.position = "plot",
          plot.background = element_rect(fill = "white", color = NA))
)

season_of <- function(date) {
  d <- as.Date(date)
  y <- as.integer(format(d, "%Y")) - (as.integer(format(d, "%m")) < 7)
  paste0(y, "/", substr(y + 1, 3, 4))
}

# earlier seasons shifted onto the calendar of `season_` (whole weeks), for context lines
past_seasons <- function(weekly, season_) {
  y0 <- as.integer(substr(season_, 1, 4))
  weekly |> mutate(season = season_of(date)) |>
    filter(season < season_) |>
    mutate(date = date + round((y0 - as.integer(substr(season, 1, 4))) * 52.1775) * 7) |>
    select(season, date, value)
}

# One forecast: data as published on the forecast date (solid), final data (dashed, optional), past seasons (thin
# grey, on this season's calendar), the peak-size forecast as horizontal bands (95% light, 50% darker, median line)
# at each model's most likely peak week, and below it one peak-week pmf panel per model on the same date axis.
#   prelim: date, value (as published on ref_date) · final: date, value · past: season, date, value
#   pmf: model, week, value · size_q: model, tau, value (at least 0.025, 0.25, 0.5, 0.75, 0.975)
# Copied from analysis/peak-models/R/slides_helpers.R (forecast sandbox), with a `title` argument added.
plot_peak_forecast <- function(prelim, pmf, ref_date, final = NULL, past = NULL, size_q = NULL,
                               models = unique(pmf$model), xlim = NULL, truth = NULL,
                               ylab = "Weekly admissions\n(log scale)", heights = NULL, band_width = 7, title = NULL) {
  models <- intersect(models, unique(as.character(pmf$model)))
  if (is.null(xlim)) xlim <- range(pmf$week) + c(-4, 4)
  inx <- function(d) if (is.null(d)) NULL else filter(d, date >= xlim[1], date <= xlim[2], value > 0)
  prelim <- inx(prelim); final <- inx(final); past <- inx(past)
  # gridlines on the 1st of each month, except within two weeks of the right edge (no room for that month's label)
  month_breaks <- function(lim) {
    b <- seq(as.Date(format(lim[1], "%Y-%m-01")), lim[2], by = "1 month")
    b[b >= lim[1] & b <= lim[2] - 14]
  }
  xs <- scale_x_date(NULL, limits = xlim, breaks = month_breaks, date_labels = "%b", expand = expansion(0))
  top <- ggplot() + geom_vline(xintercept = ref_date, color = INK2, linewidth = 0.6)
  if (!is.null(past)) top <- top + geom_line(data = past, aes(date, value, group = season), color = "#c9ced6",
                                             linewidth = 0.5)
  if (!is.null(size_q)) {
    sq <- size_q |> filter(model %in% models) |>
      mutate(tau = round(as.numeric(tau), 3)) |> filter(tau %in% c(0.025, 0.25, 0.5, 0.75, 0.975)) |>
      select(model, tau, value) |> pivot_wider(names_from = tau, values_from = value, names_prefix = "q") |>
      left_join(pmf |> group_by(model) |> slice_max(value, n = 1, with_ties = FALSE) |>
                  transmute(model, mode = as.Date(week)), by = "model") |>
      mutate(model_ord = match(model, models)) |> arrange(model_ord) |>
      group_by(mode) |> mutate(lane = row_number(), n = n()) |> ungroup()
    # each band sits at its model's predicted peak week (pmf mode); models sharing a mode sit side by side around it
    sq <- sq |> mutate(x0 = mode + (lane - 1 - n / 2) * band_width, x1 = x0 + band_width - 1)
    top <- top +
      geom_rect(data = sq, aes(xmin = x0, xmax = x1, ymin = q0.025, ymax = q0.975, fill = model), alpha = 0.18) +
      geom_rect(data = sq, aes(xmin = x0, xmax = x1, ymin = q0.25, ymax = q0.75, fill = model), alpha = 0.4) +
      geom_segment(data = sq, aes(x = x0, xend = x1, y = q0.5, yend = q0.5, color = model), linewidth = 1)
  }
  if (!is.null(final)) top <- top + geom_line(data = final, aes(date, value), color = INK, linewidth = 0.8,
                                              linetype = "22")
  top <- top + geom_line(data = prelim, aes(date, value), color = INK, linewidth = 1.2) +
    geom_point(data = slice_max(prelim, date, n = 1), aes(date, value), color = INK, size = 2.5)
  if (!is.null(truth)) top <- top + geom_vline(xintercept = truth, color = "#d9531e", linewidth = 0.8, linetype = "22")
  top <- top + xs + scale_fill_manual(values = MODEL_COLORS, guide = "none") +
    scale_color_manual(values = MODEL_COLORS, guide = "none") +
    scale_y_log10(ylab, labels = label_comma()) + theme(axis.text.x = element_blank())
  if (!is.null(title)) top <- top + labs(title = title)
  # weeks where the cumulative peak-week probability first reaches 25%, 50% and 75% (a dot above that week's bar;
  # dots for quartiles falling in the same week are stacked)
  qd <- pmf |> filter(model %in% models) |> mutate(model = factor(model, levels = models)) |>
    group_by(model) |> arrange(week, .by_group = TRUE) |> mutate(cum = cumsum(value), top = max(value)) |>
    reframe(p = c(0.25, 0.5, 0.75), week = week[vapply(p, function(x) which(cum >= x - 1e-9)[1], 1L)],
            value = value[vapply(p, function(x) which(cum >= x - 1e-9)[1], 1L)], top = top[1]) |>
    group_by(model, week) |> mutate(y = value + top * (0.1 + 0.12 * (row_number() - 1))) |> ungroup()
  bottom <- pmf |> filter(model %in% models) |> mutate(model = factor(model, levels = models)) |>
    ggplot(aes(week, value, fill = model)) + geom_col(width = 5.5) +
    geom_point(data = qd, aes(week, y), inherit.aes = FALSE, color = INK, size = 1.4) +
    geom_vline(xintercept = ref_date, color = INK2, linewidth = 0.6) +
    facet_wrap(~model, ncol = 1, strip.position = "right") + xs +
    scale_fill_manual(values = MODEL_COLORS, guide = "none") +
    scale_y_continuous("P(peak in week)", labels = label_percent(), n.breaks = 3) +
    theme(strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
          # gridlines mark the 1st of each month; start each month's name at its gridline
          axis.text.x = element_text(hjust = 0, margin = margin(t = 2)))
  if (!is.null(truth)) bottom <- bottom + geom_vline(xintercept = truth, color = "#d9531e", linewidth = 0.8,
                                                     linetype = "22")
  if (is.null(heights)) heights <- c(1.6, 0.55 * length(models))
  top / bottom + plot_layout(heights = heights)
}

args <- commandArgs(trailingOnly = TRUE)
ref_date <- as.Date(args[1])
data_date <- ref_date - 3  # the NHSN release (Wednesday) used by the forecast

fc <- read_csv(paste0("output/model-output/", MODEL, "/", ref_date, "-", MODEL, ".csv"),
               col_types = cols(.default = col_character()))
locations <- read_csv(
  "https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/refs/heads/main/auxiliary-data/locations.csv",
  col_types = cols(.default = col_character())
)
nhsn <- read_csv(paste0("https://infectious-disease-data.s3.amazonaws.com/data-raw/influenza-nhsn/nhsn-", data_date,
                        ".csv"), show_col_types = FALSE) |>
  transmute(date = as.Date(`Week Ending Date`), abbreviation = `Geographic aggregation`,
            value = as.numeric(`Total Influenza Admissions`)) |>
  mutate(abbreviation = ifelse(abbreviation == "USA", "US", abbreviation)) |>
  inner_join(locations |> select(abbreviation, location), by = "abbreviation") |>
  filter(!is.na(value))
loc_name <- setNames(locations$location_name, locations$location)
loc_name["US"] <- "United States"
out_pdf <- paste0("output/plots/", ref_date, "-", MODEL, ".pdf")

season <- season_of(ref_date)
fc <- fc |> mutate(model = MODEL, value = as.numeric(value))
pmf_all <- fc |> filter(target == "peak week inc flu hosp") |> transmute(location, model, week = as.Date(output_type_id), value)
size_all <- fc |> filter(target == "peak inc flu hosp") |> transmute(location, model, tau = output_type_id, value)

# US first, then alphabetical by name
locs <- unique(fc$location)
locs <- locs[order(locs != "US", loc_name[locs])]

# x axis: from 9 weeks before the first peak-week Saturday (shows the season so far) to the end of the window
xlim <- c(min(pmf_all$week) - 9 * 7, max(pmf_all$week) + 4)

panels <- lapply(locs, function(loc) {
  weekly <- nhsn |> filter(location == loc) |> select(date, value)
  plot_peak_forecast(
    prelim = weekly |> filter(season_of(date) == season),
    pmf = pmf_all |> filter(location == loc) |> select(-location),
    ref_date = ref_date,
    past = past_seasons(weekly, season) |> filter(season >= "2022/23"),
    size_q = size_all |> filter(location == loc) |> select(-location),
    models = MODEL, xlim = xlim, ylab = "Weekly admissions\n(log scale)",
    title = unname(ifelse(is.na(loc_name[loc]), loc, loc_name[loc]))
  ) & theme(strip.text.y = element_blank())
})

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
pages <- split(panels, ceiling(seq_along(panels) / PER_PAGE))
pdf(out_pdf, width = 13, height = 8.5)
for (p in seq_along(pages)) {
  print(wrap_plots(pages[[p]], ncol = 3, nrow = 2) +
          plot_annotation(
            title = paste0(MODEL, ", forecast date ", format(ref_date, "%b %d %Y"), " (page ", p, " of ",
                           length(pages), ")"),
            subtitle = paste("Black: NHSN data released", format(data_date, "%b %d %Y"), "(grey vertical line: forecast date);",
                             "grey: earlier seasons (since 2022/23) on this season's calendar.",
                             "\nBand: peak-size median, 50% and 95% intervals, at the most likely peak week.",
                             "Bars: probability that each week is the peak week;",
                             "dots: weeks where the cumulative probability reaches 25%, 50% and 75%."),
            theme = theme(plot.title = element_text(face = "bold", size = 13))))
}
invisible(dev.off())
cat("wrote", out_pdf, "\n")
