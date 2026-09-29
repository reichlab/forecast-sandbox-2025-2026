# Influenza virologic surveillance (WHO/NREVSS) by type and subtype, from CDC FluView Interactive.
# Descriptive use only. Files are cached under eda/strain/ (gitignored).
#
# Source: the FluView Interactive portal (https://gis.cdc.gov/grasp/fluview/fluportaldashboard.html) posts to
# https://gis.cdc.gov/flu2/PostPhase02DataDownload and returns a zip with three CSVs per geography:
#   ICL_NREVSS_Combined_prior_to_2015_16.csv  weekly, clinical + public health labs combined, 1997/98-2014/15
#                                             (states from 2010/11): A (2009 H1N1), A (H1), A (H3), A (subtyping not
#                                             performed), A (unable to subtype), B, H3N2v, A (H5)
#   ICL_NREVSS_Clinical_Labs.csv              weekly, 2015/16 on: total specimens, total A, total B
#   ICL_NREVSS_Public_Health_Labs.csv         2015/16 on: A (2009 H1N1), A (H3), A (subtyping not performed), B, BVic,
#                                             BYam, H3N2v, A (H5); weekly for national and HHS regions, but only
#                                             season totals for states

STRAIN_DIR <- "eda/strain"

download_strain <- function(dir = STRAIN_DIR, season_ids = 37:65) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  geos <- list(national = list(rt = 3, subs = list(list(ID = 0, Name = ""))),
               hhs = list(rt = 1, subs = lapply(1:10, function(i) list(ID = i, Name = as.character(i)))),
               state = list(rt = 5, subs = lapply(1:59, function(i) list(ID = i, Name = as.character(i)))))
  for (g in names(geos)) {
    body <- list(AppVersion = "Public", DatasourceDT = list(list(ID = 1, Name = "WHO_NREVSS")),
                 RegionTypeId = geos[[g]]$rt, SubRegionsDT = geos[[g]]$subs,
                 SeasonsDT = lapply(season_ids, function(i) list(ID = i, Name = as.character(i))))
    resp <- httr2::request("https://gis.cdc.gov/flu2/PostPhase02DataDownload") |>
      httr2::req_body_json(body, auto_unbox = TRUE) |> httr2::req_timeout(180) |> httr2::req_perform()
    zf <- tempfile(fileext = ".zip")
    writeBin(httr2::resp_body_raw(resp), zf)
    unzip(zf, exdir = file.path(dir, paste0("cdc_", g)))
  }
  invisible(dir)
}

read_cdc <- function(path) {
  x <- readr::read_csv(path, skip = 1, na = c("X", "", "NA"), show_col_types = FALSE, guess_max = 1e5)
  names(x) <- gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(names(x))))
  x
}

# Saturday ending MMWR week `week` of `year` (week 1 ends on the first Saturday on or after Jan 4)
mmwr_end <- function(year, week) {
  j4 <- as.Date(paste0(year, "-01-04"))
  j4 + ((6 - as.integer(format(j4, "%u"))) %% 7) + 7 * (week - 1)
}

add_season <- function(x) {
  x |> mutate(date = mmwr_end(year, week),
              y0 = ifelse(week >= 31, year, year - 1),
              season = paste0(y0, "/", substr(y0 + 1, 3, 4)),
              season_week = as.integer(round(as.numeric(date - season_week_date(season, 1)) / 7)) + 1L) |>
    select(-y0)
}

location_code <- function(region_type, region) {
  abb <- setNames(c(state.abb, "DC"), c(state.name, "District of Columbia"))
  case_when(region_type == "National" ~ "US",
            region_type == "HHS Regions" ~ region,
            TRUE ~ unname(abb[region]))
}

# Weekly A and B positives (and H1 / H3 / B lineage where reported) by location
load_strain_weekly <- function(dir = STRAIN_DIR) {
  out <- list()
  for (g in c("national", "hhs", "state")) {
    pre <- read_cdc(file.path(dir, paste0("cdc_", g), "ICL_NREVSS_Combined_prior_to_2015_16.csv")) |>
      transmute(region_type, region, year, week, specimens = total_specimens,
                A = rowSums(across(c(a_2009_h1n1, a_h1, a_h3, a_subtyping_not_performed, a_unable_to_subtype, h3n2v)), na.rm = TRUE),
                B = b, H1 = a_2009_h1n1 + a_h1, H3 = a_h3 + coalesce(h3n2v, 0), BVic = NA_real_, BYam = NA_real_)
    cl <- read_cdc(file.path(dir, paste0("cdc_", g), "ICL_NREVSS_Clinical_Labs.csv")) |>
      transmute(region_type, region, year, week, specimens = total_specimens, A = total_a, B = total_b)
    post <- cl
    if (g != "state") {
      ph <- read_cdc(file.path(dir, paste0("cdc_", g), "ICL_NREVSS_Public_Health_Labs.csv")) |>
        transmute(region_type, region, year, week, H1 = a_2009_h1n1, H3 = a_h3 + coalesce(h3n2v, 0), BVic = bvic, BYam = byam)
      post <- cl |> left_join(ph, by = c("region_type", "region", "year", "week"))
    } else {
      post <- post |> mutate(H1 = NA_real_, H3 = NA_real_, BVic = NA_real_, BYam = NA_real_)
    }
    out[[g]] <- bind_rows(pre, post) |> mutate(level = g)
  }
  bind_rows(out) |> mutate(location = location_code(region_type, region)) |>
    filter(!is.na(location)) |> add_season()
}

# Season totals of subtypes from public health labs, states, 2015/16 on
load_strain_state_phl <- function(dir = STRAIN_DIR) {
  read_cdc(file.path(dir, "cdc_state", "ICL_NREVSS_Public_Health_Labs.csv")) |>
    transmute(location = location_code(region_type, region),
              season = sub("Season (\\d{4})-(\\d{2}).*", "\\1/\\2", season_description),
              H1 = a_2009_h1n1, H3 = a_h3 + coalesce(h3n2v, 0), B = b + coalesce(bvic, 0) + coalesce(byam, 0),
              BVic = bvic, BYam = byam) |>
    filter(!is.na(location))
}
