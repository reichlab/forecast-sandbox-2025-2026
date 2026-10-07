# Data loaders and plotting helpers for peak-models-slides.qmd (figures only; no model fitting).
# Paths are relative to analysis/peak-models/.

source("R/eda_helpers.R")
suppressPackageStartupMessages({
  library(readr)
  library(jsonlite)
  library(purrr)
  library(stringr)
})

ROOT <- "../.."
MODEL_COLORS <- c("GBQR SB / core" = "#2a78d6", "GBQR core / core+hol" = "#eb6834", "Baseline" = "#6b7280",
                  "FluSight-ensemble" = "#1baf7a", "other FluSight models" = "#b7bdc6")
MODEL_IDS <- c("UMass-peak_gbqr_sb" = "GBQR SB / core", "UMass-peak_gbqr_core_hol" = "GBQR core / core+hol",
               "UMass-peak_baseline" = "Baseline")
# the two GBQR models differ in both parts, so name each by the feature set of the target being shown
TIMING_LABELS <- c("GBQR SB / core" = "GBQR timing: core", "GBQR core / core+hol" = "GBQR timing: core + holiday",
                   "Baseline" = "Baseline")
SIZE_LABELS <- c("GBQR SB / core" = "GBQR size: SB", "GBQR core / core+hol" = "GBQR size: core", "Baseline" = "Baseline")
Q_LEVELS <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
              0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

theme_slides <- function(base_size = 15) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = GRID, linewidth = 0.4),
          plot.title = element_text(face = "bold", size = base_size + 1), plot.title.position = "plot",
          plot.subtitle = element_text(color = INK2), plot.caption = element_text(color = INK2, size = base_size - 4),
          strip.text = element_text(face = "bold", hjust = 0), legend.position = "bottom",
          plot.background = element_rect(fill = "white", color = NA))
}
theme_set(theme_slides())

season_of <- function(date) {
  d <- as.Date(date)
  y <- as.integer(format(d, "%Y")) - (as.integer(format(d, "%m")) < 7)
  paste0(y, "/", substr(y + 1, 3, 4))
}

# ---------------------------------------------------------------------------------------------------------------
# data

peak_windows <- function() {
  t <- fromJSON(file.path(ROOT, "hub-config/tasks.json"), simplifyVector = FALSE)
  out <- list()
  for (mt in t$rounds[[1]]$model_tasks) {
    if (identical(unlist(mt$task_ids$target$required), "peak week inc flu hosp")) {
      dates <- as.Date(unlist(mt$output_type$pmf$output_type_id$required))
      out[[season_of(dates[1])]] <- dates
    }
  }
  out
}

load_weekly_nhsn <- function() {
  read_csv(file.path(ROOT, "target-data/oracle-output.csv"), show_col_types = FALSE,
           col_types = cols(location = col_character(), .default = col_guess())) |>
    filter(target == "wk inc flu hosp") |>
    distinct(target_end_date, location, .keep_all = TRUE) |>
    transmute(date = as.Date(target_end_date), location, value = oracle_value)
}

load_peak_oracle <- function() {
  read_csv("peak-oracle.csv", show_col_types = FALSE, col_types = cols(location = col_character())) |>
    mutate(peak_week = as.Date(peak_week))
}

# forecasts of the three final models (model-output), long format
load_model_output <- function(models = names(MODEL_IDS)) {
  map_dfr(models, function(m) {
    files <- list.files(file.path(ROOT, "model-output", m), pattern = "\\.csv$", full.names = TRUE)
    map_dfr(files, ~ read_csv(.x, show_col_types = FALSE,
                              col_types = cols(location = col_character(), output_type_id = col_character(),
                                               .default = col_guess()))) |>
      mutate(model = MODEL_IDS[[m]])
  }) |>
    mutate(reference_date = as.Date(reference_date), season = season_of(reference_date))
}

# NHSN as published on each forecast date and final, from the forecast viewer's data file
load_viewer_data <- function() {
  f <- list.files("viz", pattern = "^data-[0-9a-f]+\\.json$", full.names = TRUE)[1]
  fromJSON(f, simplifyVector = TRUE)
}
asof_series <- function(V, season, loc) {
  S <- V$seasons[[season]]
  dates <- as.Date(S$dates)
  map_dfr(names(S$asof), function(r) {
    v <- S$asof[[r]][[loc]]
    if (is.null(v)) return(NULL)
    tibble(reference_date = as.Date(r), date = dates, value = as.numeric(v))
  })
}

load_scores <- function(path = "peak-scores.csv") {
  read_csv(path, show_col_types = FALSE, col_types = cols(location = col_character())) |>
    mutate(label = ifelse(model %in% names(MODEL_IDS), MODEL_IDS[model], model))
}

# ---------------------------------------------------------------------------------------------------------------
# plots

model_scale_color <- function(...) scale_color_manual(values = MODEL_COLORS, ...)
model_scale_fill <- function(...) scale_fill_manual(values = MODEL_COLORS, ...)

# peak-week pmf heat strip: rows = reference dates, columns = window weeks, one panel per model
pmf_strip <- function(fc, loc, season_, truth = NULL, models = unname(MODEL_IDS), refs = NULL, xlim = NULL) {
  d <- fc |> filter(location == loc, season == season_, target == "peak week inc flu hosp", model %in% models) |>
    mutate(week = as.Date(output_type_id), model = factor(model, levels = models))
  if (!is.null(refs)) d <- d |> filter(reference_date %in% refs)
  p <- ggplot(d, aes(week, reference_date, fill = value)) +
    geom_tile(height = 6.2, width = 6.2) +
    scale_fill_gradientn("P(peak in week)", colours = c("#f4f6f9", "#cde2fb", "#6da7ec", "#256abf", "#0d366b"),
                         values = scales::rescale(c(0, 0.05, 0.15, 0.35, 1)), limits = c(0, 1)) +
    geom_point(aes(x = reference_date, y = reference_date), shape = 124, size = 3, color = INK, inherit.aes = FALSE,
               data = distinct(d, model, reference_date)) +
    facet_wrap(~model, nrow = 1) +
    scale_x_date(NULL, date_breaks = "1 month", date_labels = "%b", limits = xlim) +
    scale_y_date("Forecast date", date_labels = "%b %d") +
    theme(legend.position = "right", legend.key.height = unit(1.4, "cm"))
  if (!is.null(truth)) p <- p + geom_vline(xintercept = truth, color = "#d9531e", linewidth = 0.8, linetype = "22")
  p
}

# peak-size intervals over reference dates, one panel per model, with the observed peak (horizontal) and the
# observed peak week (vertical)
size_ribbons <- function(fc, loc, season_, peak = NULL, models = unname(MODEL_IDS),
                         ylab = "Peak weekly admissions (log scale)", peak_week = NULL) {
  d <- fc |> filter(location == loc, season == season_, target == "peak inc flu hosp", model %in% models) |>
    mutate(tau = as.numeric(output_type_id)) |>
    filter(tau %in% c(0.025, 0.25, 0.5, 0.75, 0.975)) |>
    select(model, reference_date, tau, value) |>
    pivot_wider(names_from = tau, values_from = value, names_prefix = "q") |>
    mutate(model = factor(model, levels = models))
  p <- ggplot(d, aes(reference_date)) +
    geom_ribbon(aes(ymin = q0.025, ymax = q0.975, fill = model), alpha = 0.2) +
    geom_ribbon(aes(ymin = q0.25, ymax = q0.75, fill = model), alpha = 0.45) +
    geom_line(aes(y = q0.5, color = model), linewidth = 1) +
    facet_wrap(~model, nrow = 1) + model_scale_color(guide = "none") + model_scale_fill(guide = "none") +
    scale_y_log10(ylab, labels = label_comma()) +
    scale_x_date("Forecast date", date_breaks = "1 month", date_labels = "%b")
  if (!is.null(peak)) p <- p + geom_hline(yintercept = peak, color = "#d9531e", linewidth = 0.8)
  if (!is.null(peak_week)) p <- p + geom_vline(xintercept = peak_week, color = "#d9531e", linewidth = 0.8,
                                               linetype = "22")
  p
}

# small final-data panel to go under a pmf_strip, one copy per model column so the x axes line up
season_data_strip <- function(final, models, xlim, peak_week = NULL) {
  d <- tidyr::crossing(final, model = factor(models, levels = models))
  p <- ggplot(d, aes(date, value)) + geom_line(color = INK, linewidth = 0.8) + facet_wrap(~model, nrow = 1) +
    scale_x_date(NULL, date_breaks = "1 month", date_labels = "%b", limits = xlim) +
    scale_y_continuous("weekly\nadmissions", labels = label_comma(), n.breaks = 3) +
    theme(strip.text = element_blank(), axis.title.y = element_text(size = 12))
  if (!is.null(peak_week)) p <- p + geom_vline(xintercept = peak_week, color = "#d9531e", linewidth = 0.8,
                                               linetype = "22")
  p
}

# pmf_strip with the season's final data below each model's panel
pmf_strip_data <- function(fc, loc, season_, final, truth = NULL, models = unname(MODEL_IDS), heights = c(4, 1),
                           mark_refs = NULL) {
  w <- windows[[season_]]
  xlim <- c(min(w) - 4, max(w) + 4)
  top <- pmf_strip(fc, loc, season_, truth = truth, models = models, xlim = xlim) +
    theme(axis.text.x = element_blank())
  if (!is.null(mark_refs)) top <- top + geom_hline(yintercept = mark_refs, color = INK2, linetype = "dotted")
  bottom <- season_data_strip(final |> filter(date >= xlim[1], date <= xlim[2]), models, xlim, truth)
  top / bottom + plot_layout(heights = heights, guides = "collect") & theme(legend.position = "right")
}

# One forecast: data as published on the forecast date (solid), final data (dashed), past seasons (thin grey, on this
# season's calendar), the peak-size forecast as horizontal bands after the forecast date (95% light, 50% darker,
# median line), each at its model's predicted peak week (pmf mode), and below it one peak-week pmf panel per model on the same date axis.
#   prelim: date, value (as published on ref_date) · final: date, value · past: season, date, value
#   pmf: model, week, value · size_q: model, tau, value (at least 0.025, 0.25, 0.5, 0.75, 0.975)
plot_peak_forecast <- function(prelim, pmf, ref_date, final = NULL, past = NULL, size_q = NULL,
                               models = unique(pmf$model), xlim = NULL, truth = NULL,
                               ylab = "Weekly admissions\n(log scale)", heights = NULL, band_width = 7) {
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
  top <- top + xs + model_scale_fill(guide = "none") + model_scale_color(guide = "none") +
    scale_y_log10(ylab, labels = label_comma()) + theme(axis.text.x = element_blank())
  # weeks where the cumulative peak-week probability first reaches 25%, 50% and 75% (a dot above that week's bar;
  # dots for quartiles falling in the same week are stacked)
  qd <- pmf |> filter(model %in% models) |> mutate(model = factor(model, levels = models)) |>
    group_by(model) |> arrange(week, .by_group = TRUE) |> mutate(cum = cumsum(value), top = max(value)) |>
    reframe(p = c(0.25, 0.5, 0.75), week = week[vapply(p, function(x) which(cum >= x - 1e-9)[1], 1L)],
            value = value[vapply(p, function(x) which(cum >= x - 1e-9)[1], 1L)], top = top[1]) |>
    group_by(model, week) |> mutate(y = value + top * (0.1 + 0.12 * (row_number() - 1))) |> ungroup()
  bottom <- pmf |> filter(model %in% models) |> mutate(model = factor(model, levels = models)) |>
    ggplot(aes(week, value, fill = model)) + geom_col(width = 5.5) +
    geom_point(data = qd, aes(week, y), inherit.aes = FALSE, color = INK, size = 2) +
    geom_vline(xintercept = ref_date, color = INK2, linewidth = 0.6) +
    facet_wrap(~model, ncol = 1, strip.position = "right") + xs + model_scale_fill(guide = "none") +
    scale_y_continuous("P(peak in week)", labels = label_percent(), n.breaks = 3) +
    theme(strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
          # gridlines mark the 1st of each month; start each month's name at its gridline
          axis.text.x = element_text(hjust = 0, margin = margin(t = 2)))
  if (!is.null(truth)) bottom <- bottom + geom_vline(xintercept = truth, color = "#d9531e", linewidth = 0.8,
                                                     linetype = "22")
  if (is.null(heights)) heights <- c(1.6, 0.55 * length(models))
  top / bottom + plot_layout(heights = heights)
}

# series as published on a date (from the viewer data) and past NHSN seasons shifted onto a season's calendar
prelim_on <- function(V, season, loc, ref_date) {
  asof_series(V, season, loc) |> filter(reference_date == ref_date, !is.na(value)) |> select(date, value)
}
past_seasons <- function(weekly, loc, season_) {
  y0 <- as.integer(substr(season_, 1, 4))
  weekly |> filter(location == loc) |> mutate(season = season_of(date)) |>
    filter(season < season_) |>
    mutate(date = date + round((y0 - as.integer(substr(season, 1, 4))) * 52.1775) * 7) |>
    select(season, date, value)
}
