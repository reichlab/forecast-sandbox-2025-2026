# Helpers for relative-size-eda.qmd: summaries and plots of the replay rows and of the cross-validation scores
# written by relative_size_eda.py / relative_size_predictors.py. Descriptive only.

source("R/eda_helpers.R")

BINS <- list(c(5, 11), c(12, 16), c(17, 21), c(22, 26), c(27, 31), c(32, 43))
CORR_BINS <- list(c(12, 16), c(17, 21), c(22, 26), c(27, 31))
bin_label <- function(b) paste0(b[1], "–", b[2])
week_bin <- function(w, bins = BINS) {
  out <- rep(NA_character_, length(w))
  for (b in bins) out[w >= b[1] & w <= b[2]] <- bin_label(b)
  factor(out, levels = sapply(bins, bin_label))
}

read_rows <- function(path) {
  read_parquet(path) |> mutate(group = factor(group, levels = SOURCE_GROUPS))
}

# quantiles of z by grouping variables
ztab <- function(df, ..., col = z) {
  df |> group_by(...) |>
    summarise(q10 = quantile({{ col }}, 0.1), q25 = quantile({{ col }}, 0.25), q50 = median({{ col }}),
              q75 = quantile({{ col }}, 0.75), q90 = quantile({{ col }}, 0.9), p0 = mean(at_zero), n = n(),
              .groups = "drop")
}

# ---------------------------------------------------------------------------------------------------------------
# example series: top = series, running max and peak (log scale); bottom = eventual z

example_panel <- function(arr, rows, src, loc, ssn, title, ylab, t_arrow = 18, show_legend = FALSE) {
  eps <- LOG_EPS[[src]]
  a <- arr |> filter(source == src, location == loc, season == ssn, season_week >= REPLAY_START, season_week <= WINDOW[2])
  r <- rows |> filter(source == src, location == loc, season == ssn)
  pk_w <- r$peak_week[1]
  pk <- a$y[a$season_week == pk_w]
  ma <- a$M[a$season_week == t_arrow]
  za <- log(pk + eps) - log(ma + eps)
  top <- ggplot(a, aes(season_week)) +
    pre_window(log_y = TRUE) +
    geom_hline(yintercept = pk + eps, linetype = "dashed", color = BLUE, linewidth = 0.4) +
    geom_line(aes(y = y + eps, color = "weekly value"), linewidth = 0.5) +
    geom_point(aes(y = y + eps, color = "weekly value"), size = 0.7) +
    geom_step(aes(y = M + eps, color = "running max M[t]"), linewidth = 0.9) +
    annotate("point", x = pk_w, y = pk + eps, shape = 8, color = BLUE, size = 2.5) +
    annotate("text", x = pk_w + 0.7, y = (pk + eps) * 1.6, label = paste("peak, week", pk_w), hjust = 0, size = 3) +
    annotate("segment", x = t_arrow, xend = t_arrow, y = ma + eps, yend = pk + eps, color = BLUE,
             arrow = arrow(ends = "both", length = unit(0.12, "cm"))) +
    annotate("text", x = t_arrow - 0.6, y = sqrt((ma + eps) * (pk + eps)), label = paste0("z[", t_arrow, "] == ", f2(za, 1)),
             parse = TRUE, hjust = 1, color = BLUE, size = 3) +
    scale_y_log10(paste(ylab, "+ eps (log)"), limits = c(eps * 0.6, (pk + eps) * 3), labels = plain_log_labels) +
    scale_color_manual(values = c("weekly value" = INK2, "running max M[t]" = ACCENT),
                       labels = c("weekly value" = "weekly value", "running max M[t]" = expression(running ~ max ~ M[t])),
                       breaks = c("weekly value", "running max M[t]")) +
    scale_x_season_week(breaks = c(5, 10, 20, 30, 40)) +
    labs(title = title, x = NULL) +
    theme(legend.position = if (show_legend) "top" else "none", axis.text.x = element_blank(),
          axis.title.x = element_blank())
  bottom <- ggplot(r, aes(season_week, z)) +
    pre_window() +
    geom_vline(xintercept = pk_w, linetype = "dashed", color = BLUE, linewidth = 0.4) +
    geom_line(color = BLUE, linewidth = 0.8) + geom_point(color = BLUE, size = 0.9) +
    scale_x_season_week(breaks = c(5, 10, 20, 30, 40)) +
    scale_y_continuous("eventual z", limits = Z_LIM)
  top / bottom + plot_layout(heights = c(1.25, 1))
}

example_row <- function(arr, rows, specs) {
  panels <- lapply(seq_along(specs), function(i) {
    s <- specs[[i]]
    example_panel(arr, rows, s$src, s$loc, s$season, s$title, s$ylab, show_legend = i == 1)
  })
  wrap_plots(panels, nrow = 1)
}

# ---------------------------------------------------------------------------------------------------------------
# quantile bars (median, 25-75%, 10-90%) with share at z = 0 printed above

quantile_bars <- function(df, xvar, xlab, color = BLUE) {
  s <- df |> filter(!is.na({{ xvar }})) |> group_by(week_lab, x = {{ xvar }}) |>
    summarise(q10 = quantile(z, 0.1), q25 = quantile(z, 0.25), q50 = median(z), q75 = quantile(z, 0.75),
              q90 = quantile(z, 0.9), p0 = mean(at_zero), n = n(), .groups = "drop") |>
    filter(n >= 10)
  ggplot(s, aes(x)) +
    geom_linerange(aes(ymin = q10, ymax = q90), color = color, alpha = 0.5, linewidth = 0.6) +
    geom_linerange(aes(ymin = q25, ymax = q75), color = color, alpha = 0.85, linewidth = 3.5) +
    geom_errorbar(aes(ymin = q50, ymax = q50), width = 0.45, color = INK, linewidth = 0.8) +
    geom_text(aes(y = Z_LIM[2] * 0.97, label = paste0(percent(p0, 1), "\nn=", n)), size = 2.4, color = INK2,
              vjust = 1, lineheight = 0.9) +
    facet_wrap(~week_lab, nrow = 1) +
    scale_y_continuous("eventual z", limits = Z_LIM) +
    labs(x = xlab) +
    theme(panel.grid.major.x = element_blank())
}

# ---------------------------------------------------------------------------------------------------------------
# within-week-bin Spearman correlations

spearman_bins <- function(df, feats, target, bins = CORR_BINS) {
  df$bin <- week_bin(df$season_week, bins)
  expand_grid(feature = feats, bin = levels(df$bin)) |>
    rowwise() |>
    mutate(rho = {
      x <- df[df$bin %in% bin, c(feature, target)]
      suppressWarnings(cor(x[[1]], as.numeric(x[[2]]), method = "spearman", use = "complete.obs"))
    }) |>
    ungroup()
}

corr_heatmap <- function(df, groups, first, first_label = "reference") {
  feats <- c(first, unlist(groups, use.names = FALSE))
  grp <- c(rep(first_label, length(first)), rep(names(groups), lengths(groups)))
  d <- bind_rows(spearman_bins(df, feats, "z") |> mutate(target = "Spearman with z"),
                 spearman_bins(df, feats, "pos") |> mutate(target = "Spearman with 1{z > 0}")) |>
    mutate(feature = factor(feature, levels = rev(feats)), group = factor(grp[match(feature, feats)], levels = unique(grp)),
           target = factor(target, levels = c("Spearman with z", "Spearman with 1{z > 0}")),
           bin = paste("wk", bin))
  p <- ggplot(d, aes(bin, feature, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = f2(rho), color = abs(rho) > 0.6), size = 2.4) +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
    scale_fill_distiller("Spearman correlation within the week bin (same scale in both panels)", palette = "RdBu",
                         limits = c(-1, 1)) +
    facet_grid(group ~ target, scales = "free_y", space = "free_y") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme(panel.grid = element_blank(), legend.position = "bottom", strip.text.y = element_text(angle = 0, hjust = 0),
          strip.text.x = element_text(hjust = 0)) +
    guides(fill = guide_colorbar(barwidth = 18, barheight = 0.6, title.position = "top"))
  list(plot = p, data = d)
}

# ---------------------------------------------------------------------------------------------------------------
# binned relationship (deciles of a feature within week bins)

binned_plot <- function(df, feats, bins = list(c(12, 16), c(17, 21), c(22, 26)), current = character()) {
  shades <- c("#6baed6", "#2171b5", "#08306b")
  d <- bind_rows(lapply(feats, function(f) {
    bind_rows(lapply(seq_along(bins), function(i) {
      b <- bins[[i]]
      x <- df |> filter(season_week >= b[1], season_week <= b[2], !is.na(.data[[f]]))
      if (n_distinct(x[[f]]) < 5) return(NULL)
      x |> mutate(dec = ntile(.data[[f]], 10)) |> group_by(dec) |>
        summarise(mid = median(.data[[f]]), q25 = quantile(z, 0.25), q50 = median(z), q75 = quantile(z, 0.75),
                  .groups = "drop") |>
        mutate(feature = f, weeks = paste("weeks", bin_label(b)))
    }))
  })) |> mutate(feature = factor(ifelse(feature %in% current, paste(feature, "(current feature)"), feature),
                                 levels = ifelse(feats %in% current, paste(feats, "(current feature)"), feats)))
  ggplot(d, aes(mid, q50, color = weeks, fill = weeks)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) + geom_point(size = 1) +
    facet_wrap(~feature, scales = "free_x", ncol = 3) +
    scale_color_manual(values = shades) + scale_fill_manual(values = shades) +
    scale_y_continuous("eventual z (median, 25–75%)", limits = Z_LIM) +
    labs(x = NULL)
}

# ---------------------------------------------------------------------------------------------------------------
# cross-validation summaries

cv_rel <- function(cv, ref, bins = BINS) {
  cv <- cv |> mutate(bin = week_bin(season_week, bins))
  by_bin <- cv |> group_by(set, bin) |> summarise(pinball = mean(pinball), logloss = mean(logloss), .groups = "drop")
  all <- cv |> group_by(set) |> summarise(pinball = mean(pinball), logloss = mean(logloss), .groups = "drop") |>
    mutate(bin = "all")
  x <- bind_rows(by_bin |> mutate(bin = as.character(bin)), all)
  r <- x |> filter(set == ref) |> select(bin, p_ref = pinball, l_ref = logloss)
  x |> left_join(r, by = "bin") |>
    mutate(rel_pinball = pinball / p_ref, d_logloss = logloss - l_ref,
           bin = factor(bin, levels = c(sapply(bins, bin_label), "all")))
}

cv_table <- function(rel, value, order, digits = 3) {
  rel |> select(set, bin, v = {{ value }}) |>
    pivot_wider(names_from = bin, values_from = v) |>
    mutate(set = factor(set, levels = order)) |> arrange(set) |>
    mutate(across(-set, ~ formatC(.x, format = "f", digits = digits))) |>
    rename(`feature set` = set)
}

# display names for the feature sets (internal names in the saved CV outputs are "current" = core, "SB" = SB);
# used only for tables and figures, never for filtering
vlabel <- function(x) {
  x <- as.character(x)
  x <- sub("^current features$", "core", x)
  x <- sub("^current(?=$| )", "core", x, perl = TRUE)
  x <- gsub("(?<![A-Za-z-])SB(?![A-Za-z])", "SB", x, perl = TRUE)
  sub("^core \\+ synchrony \\+ burden$", "SB (core + synchrony + burden)", x)
}
relabel_df <- function(df) {
  df <- as.data.frame(df, check.names = FALSE)
  for (j in seq_along(df)) if (is.character(df[[j]]) || is.factor(df[[j]])) df[[j]] <- vlabel(df[[j]])
  names(df) <- vlabel(names(df))
  df
}

cv_dotplot <- function(rel, order, ref_label) {
  d <- rel |> filter(set %in% order) |> mutate(set = factor(vlabel(set), levels = vlabel(rev(order))))
  shades <- setNames(c(colorRampPalette(c("#9ecae1", "#08306b"))(nlevels(d$bin) - 1), INK), levels(d$bin))
  labs_b <- setNames(ifelse(levels(d$bin) == "all", "all weeks", paste("weeks", levels(d$bin))), levels(d$bin))
  pd <- position_dodge(width = 0.7)
  p1 <- ggplot(d, aes(rel_pinball, set, color = bin, shape = bin == "all")) +
    geom_vline(xintercept = 1, color = INK2) +
    geom_point(position = pd, size = 2) +
    labs(x = paste("pinball loss relative to", vlabel(ref_label)), y = NULL, title = "Quantile forecasts of z (lower is better)")
  p2 <- ggplot(d, aes(d_logloss, set, color = bin, shape = bin == "all")) +
    geom_vline(xintercept = 0, color = INK2) +
    geom_point(position = pd, size = 2) +
    labs(x = paste("change in log loss for z > 0 vs", vlabel(ref_label)), y = NULL,
         title = "Probability that the peak is still ahead (lower is better)") +
    theme(axis.text.y = element_blank())
  (p1 | p2) + plot_layout(guides = "collect") &
    scale_color_manual(values = shades, labels = labs_b) &
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") &
    theme(legend.position = "bottom", panel.grid.major.y = element_blank())
}

boot_get <- function(boot, fam, s, ref, metric, what = "boot_win") {
  boot |> filter(family == fam, set == s, reference == ref, metric == !!metric) |> pull(all_of(what))
}

# season-cluster bootstrap mean and 95% interval
cluster_ci <- function(v, cl, n_boot = 2000, seed = 0) {
  g <- tibble(v = v, c = cl) |> group_by(c) |> summarise(s = sum(v), n = n())
  set.seed(seed)
  bs <- replicate(n_boot, {
    i <- sample.int(nrow(g), replace = TRUE)
    sum(g$s[i]) / sum(g$n[i])
  })
  c(mean = mean(v), lo = unname(quantile(bs, 0.025)), hi = unname(quantile(bs, 0.975)), n = length(v))
}
