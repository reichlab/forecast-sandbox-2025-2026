"""
Search for better real-time predictors of the relative size target z (companion to relative_size_eda.py).

Model-side computations only; all figures and descriptive tables are made in R by relative-size-eda.qmd.
Builds the season-replay rows exactly as the models do (via relative_size_eda.load_rows, which calls
idmodels.peak.series.build_replay_rows; NHSN excluded, Puerto Rico / Virgin Islands ILINet dropped), adds candidate
features computed only from data available at season week t, runs leakage checks, and runs leave-one-season-out
LightGBM quantile regression of z and classification of z > 0 for many feature sets. Writes to
analysis/peak-models/eda/:
  predictor_rows.parquet     replay rows with all candidate features
  cv_results.parquet         per-row CV scores (pinball, log loss) for the section 9 feature sets; cv_sets.csv
  cv_holiday.parquet, cv_holiday_adjust.parquet, holiday_excess.parquet, holiday_calendar.csv   (section 10)
  cv_lit.parquet, reflection_check.csv                                                          (section 11)
  cv_sb_groups.parquet       section 9 groups and holiday added to SB (so every group has an SB reference)
  cv_bootstrap.csv           season-bootstrap win rates and per-season comparisons for every feature set
  importance.csv             LightGBM gain importances
  predictors_numbers.json    leakage / consistency checks, holiday adjustment factors, runtimes

Usage (from the repository root; about 10 minutes from scratch, 2-3 minutes with all --reuse_* flags):
    OMP_NUM_THREADS=1 DYLD_FALLBACK_LIBRARY_PATH=<venv>/lib/python3.12/site-packages/sklearn/.dylibs \
        python analysis/peak-models/relative_size_predictors.py [--reuse_cv --reuse_holiday_cv --reuse_lit_cv]
"""
import argparse
import json
import sys
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from idmodels.peak.base import season_week_to_date  # noqa: E402
from idmodels.peak.series import build_replay_rows, build_season_arrays, running_max, season_peaks  # noqa: E402
from relative_size_eda import MIN_OBS, OUT, REPLAY_START, W0, W1, load_rows  # noqa: E402

KEYS = ["source", "agg_level", "location", "season", "season_week"]
CURRENT = ["season_week", "rel_max", "wks_since_max", "g1", "g2", "g3", "rm3", "cum_rel", "hist_rel", "nat_rel_max",
           "nat_wks_since_max", "nat_g3", "src_code"]
GROUPS_NEW = {
    "trend": ["tay1_w4_lvl", "tay1_w4_slope", "tay1_w6_lvl", "tay1_w6_slope", "tay2_w6_lvl", "tay2_w6_slope",
              "tay2_w6_curv", "tay2_w8_lvl", "tay2_w8_slope", "tay2_w8_curv", "rm2", "rm4", "lag1_rel", "lag2_rel",
              "lag3_rel", "lag4_rel"],
    "timing": ["t_minus_hist_pk", "t_minus_hist_pk_med", "t_minus_clim_pk", "hist_pk_sd", "hist_z_same_week"],
    "synchrony": ["sync_frac_past2", "sync_frac_half", "sync_med_rel_max", "sync_med_g3"],
    "onset": ["wks_since_onset2", "wks_since_onset4", "max_over_base", "cur_over_base"],
    "burden": ["cum_vs_hist_total", "cum_vs_hist_same_week"],
    "level_hist": ["cur_vs_hist_peak", "lvl_vs_hist_same_week", "max_vs_hist_same_week"],
}
ALL_NEW = [f for g in GROUPS_NEW.values() for f in g]
# holiday features are evaluated separately (section 10), so ALL_NEW above is unchanged
HOLIDAY = ["hol_now", "max_in_holiday", "holiday_excess_max", "weeks_since_holiday_end", "rel_max_nohol",
           "wks_since_max_nohol"]
HOL_FLAGS = ["hol_now", "max_in_holiday", "weeks_since_holiday_end"]
HOL_ADJ = ["holiday_excess_max", "rel_max_nohol", "wks_since_max_nohol"]
# section 11: candidates suggested by analogous problems in other fields
GROUPS_LIT = {
    "reflection": ["sigma_hat", "snr", "p_exc0", "p_exc", "ez_rw", "p_exc0_alt", "ez_rw_alt"],
    "severity": ["sev_peaked", "sev_gap", "n_peaked"],
    "chain_ladder": ["cl_z", "p0_hist", "z_bf"],
    "deceleration": ["r_max", "gfrac", "wks_since_gmax", "g_decel", "z_par"],
    "recession": ["rec_rate", "rec_consist", "rec_anom"],
    "records": ["n_rec4", "rec_frac"],
    "bass": ["bass_rel"],
    "cross_source": ["xs_rel_max", "xs_wsm", "xs_g3"],
}
LIT_LABELS = {"reflection": "reflection principle", "severity": "severity in peaked locations",
              "chain_ladder": "chain ladder", "deceleration": "deceleration", "recession": "recession",
              "records": "record counts", "bass": "Bass implied peak", "cross_source": "cross-source"}
ALL_LIT = [f for g in GROUPS_LIT.values() for f in g]
# section 12: regional / neighbour synchrony (state series only; FIPS -> centroid lat, lon and HHS region)
STATE_CENTROIDS = {"01": (32.59, -86.75), "02": (64.0, -150.0), "04": (34.22, -111.62), "05": (34.73, -92.30), "06": (36.53, -119.77), "08": (38.68, -105.51), "09": (41.59, -72.36), "10": (38.68, -74.98), "11": (38.90, -77.03), "12": (27.87, -81.69), "13": (32.33, -83.37), "15": (20.8, -156.3), "16": (43.56, -113.93), "17": (40.05, -89.38), "18": (40.05, -86.08), "19": (41.94, -93.37), "20": (38.42, -98.12), "21": (37.39, -84.77), "22": (30.62, -92.27), "23": (45.62, -68.98), "24": (39.28, -76.65), "25": (42.36, -71.58), "26": (43.14, -84.69), "27": (46.39, -94.60), "28": (32.68, -89.81), "29": (38.33, -92.51), "30": (46.82, -109.32), "31": (41.34, -99.59), "32": (39.11, -116.85), "33": (43.39, -71.39), "34": (39.96, -74.23), "35": (34.48, -105.94), "36": (43.14, -75.14), "37": (35.42, -78.47), "38": (47.25, -100.10), "39": (40.22, -82.60), "40": (35.51, -97.12), "41": (43.91, -120.07), "42": (40.91, -77.45), "44": (41.59, -71.12), "45": (33.62, -80.51), "46": (44.34, -99.72), "47": (35.68, -86.46), "48": (31.39, -98.79), "49": (39.11, -111.33), "50": (44.25, -72.55), "51": (37.56, -78.20), "53": (47.42, -119.75), "54": (38.42, -80.67), "55": (44.59, -89.99), "56": (43.05, -107.26)}  # noqa: E501
HHS_REGION = {f: r for r, fs in {
    1: ["09", "23", "25", "33", "44", "50"], 2: ["34", "36", "72", "78"], 3: ["10", "11", "24", "42", "51", "54"],
    4: ["01", "12", "13", "21", "28", "37", "45", "47"], 5: ["17", "18", "26", "27", "39", "55"],
    6: ["05", "22", "35", "40", "48"], 7: ["19", "20", "29", "31"], 8: ["08", "30", "38", "46", "49", "56"],
    9: ["04", "06", "15", "32"], 10: ["02", "16", "41", "53"]}.items() for f in fs}
GROUPS_REG = {
    "nbr5": ["nbr5_frac_past2", "nbr5_med_rel_max", "nbr5_med_wsm", "nbr5_med_g3"],
    "hhs": ["hhs_frac_past2", "hhs_med_rel_max"],
    "wave": ["wave_lead"],
    "latlon": ["lat", "lon"],
}
REG_LABELS = {"nbr5": "5 nearest states", "hhs": "HHS region", "wave": "wave position", "latlon": "latitude/longitude"}
ALL_REG = [f for g in GROUPS_REG.values() for f in g]
# section 13: influenza type / subtype (WHO/NREVSS, final values from the iddata S3 file cached in eda/strain/)
STRAIN_PATH = OUT / "strain" / "who-nrevss.csv"
STRAIN_URL = "https://infectious-disease-data.s3.amazonaws.com/data-raw/influenza-who-nrevss/who-nrevss.csv"
STATE_FIPS = {"Alabama": "01", "Alaska": "02", "Arizona": "04", "Arkansas": "05", "California": "06", "Colorado": "08",
              "Connecticut": "09", "Delaware": "10", "District of Columbia": "11", "Florida": "12", "Georgia": "13",
              "Hawaii": "15", "Idaho": "16", "Illinois": "17", "Indiana": "18", "Iowa": "19", "Kansas": "20",
              "Kentucky": "21", "Louisiana": "22", "Maine": "23", "Maryland": "24", "Massachusetts": "25",
              "Michigan": "26", "Minnesota": "27", "Mississippi": "28", "Missouri": "29", "Montana": "30",
              "Nebraska": "31", "Nevada": "32", "New Hampshire": "33", "New Jersey": "34", "New Mexico": "35",
              "New York": "36", "North Carolina": "37", "North Dakota": "38", "Ohio": "39", "Oklahoma": "40",
              "Oregon": "41", "Pennsylvania": "42", "Rhode Island": "44", "South Carolina": "45", "South Dakota": "46",
              "Tennessee": "47", "Texas": "48", "Utah": "49", "Vermont": "50", "Virginia": "51", "Washington": "53",
              "West Virginia": "54", "Wisconsin": "55", "Wyoming": "56"}
MIN_POS = 20  # minimum positives (A + B, or subtyped A) in a window for a share to be computed
GROUPS_TYPE = {
    "B": ["b_share_3wk", "b_share_cum", "b_share_trend", "b_rising", "a_rising", "b_minus_a_growth", "b_frac_of_peak"],
    "H3": ["h3_share_cum", "h3_share_3wk", "h3_share_nat"],
}
ALL_TYPE = [f for g in GROUPS_TYPE.values() for f in g]
N_WORKERS = 3  # overridden by --workers
SB = CURRENT + ["sync_frac_past2", "sync_frac_half", "sync_med_rel_max", "sync_med_g3", "cum_vs_hist_total",
                "cum_vs_hist_same_week"]
GROUP_LABELS = {"trend": "trend (Taylor, rolling means, lags)", "timing": "timing vs history",
                "synchrony": "synchrony", "onset": "onset", "burden": "burden to date",
                "level_hist": "level vs history"}
QS = [0.1, 0.25, 0.5, 0.75, 0.9]
BINS = [(5, 11), (12, 16), (17, 21), (22, 26), (27, 31), (32, 43)]
CORR_BINS = [(12, 16), (17, 21), (22, 26), (27, 31)]
LGB_PARAMS = dict(n_estimators=200, learning_rate=0.05, min_child_samples=50, num_leaves=15, verbose=-1, n_jobs=3)


def bin_label(lo, hi):
    return f"{lo}–{hi}"


def week_bin(w: pd.Series, bins=BINS) -> pd.Series:
    out = pd.Series(pd.NA, index=w.index, dtype="object")
    for lo, hi in bins:
        out[(w >= lo) & (w <= hi)] = bin_label(lo, hi)
    return out


# ---------------------------------------------------------------------------------------------------------------
# candidate features

def _ffill(y):
    idx = np.where(np.isnan(y), 0, np.arange(y.shape[1]))
    np.maximum.accumulate(idx, axis=1, out=idx)
    return y[np.arange(y.shape[0])[:, None], idx]


def _taylor_pinv(w, deg):
    """Least-squares map from a trailing window of w values (lags -w+1..0) to Taylor coefficients (level, slope,
    curvature) at the last week, as in timeseriesutils.featurize.windowed_taylor_coefs (trailing)."""
    lags = np.arange(-w + 1, 1, dtype=float)
    X = np.column_stack([np.ones(w)] + [lags ** d / np.prod(np.arange(1, d + 1)) for d in range(1, deg + 1)])
    return np.linalg.pinv(X)


def history_index(keys: pd.DataFrame) -> list[np.ndarray]:
    """For each series, the row indices of strictly earlier seasons of the same source and location."""
    out = []
    grouped = {k: g for k, g in keys.reset_index().groupby(["source", "location"])}
    for src, loc, ssn in zip(keys["source"], keys["location"], keys["season"]):
        g = grouped[(src, loc)]
        out.append(g.loc[g["season"] < ssn, "index"].to_numpy())
    return out


def _hist_median(vals: np.ndarray, hist: list[np.ndarray]) -> np.ndarray:
    out = np.full((len(hist),) + vals.shape[1:], np.nan)
    for i, h in enumerate(hist):
        if len(h):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                out[i] = np.nanmedian(vals[h], axis=0)
    return out


def robust_sd_diffs(Lraw: np.ndarray, t: int, width=8, min_obs=4) -> np.ndarray:
    """1.4826 * MAD of the weekly log changes over the last `width` weeks ending at t (raw, unfilled values)."""
    lo = max(t - width, 1)
    d = Lraw[:, lo:t] - Lraw[:, lo - 1:t - 1]
    ok = np.sum(~np.isnan(d), axis=1) >= min_obs
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        med = np.nanmedian(d, axis=1)
        mad = np.nanmedian(np.abs(d - med[:, None]), axis=1)
    return np.where(ok, 1.4826 * mad, np.nan)


def bass_implied(yt: np.ndarray, t: int, eps: np.ndarray, m: np.ndarray) -> np.ndarray:
    """
    Quadratic "Bass" fit x_s = a + b C_{s-1} + c C_{s-1}^2 (C = cumulative sum from week 1) over observed weeks from
    onset (first week above 2 x (early-season minimum + eps)) through t, at least 5 points. If c < 0 the implied peak
    is a - b^2 / (4c); returns log((max(x*, M_t) + eps) / (M_t + eps)) clipped to [0, 5], else NaN.
    """
    n = yt.shape[0]
    out = np.full(n, np.nan)
    C = np.concatenate([np.zeros((n, 1)), np.nancumsum(yt[:, :t], axis=1)], axis=1)  # C[:, s] = sum of weeks < s+1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        base = np.nanmin(yt[:, :min(t, 8)], axis=1)
    for i in range(n):
        x = yt[i, :t]
        if np.isnan(base[i]):
            continue
        above = np.flatnonzero(x > 2 * (base[i] + eps[i]))
        if len(above) == 0:
            continue
        s = np.arange(above[0], t)
        s = s[~np.isnan(x[s])]
        if len(s) < 5:
            continue
        c = C[i, s]
        scale = max(c.max(), 1e-12)
        cs = c / scale
        X = np.column_stack([np.ones(len(s)), cs, cs ** 2])
        a, b, cc = np.linalg.lstsq(X, x[s], rcond=None)[0]
        if cc < 0:
            xs = a - b ** 2 / (4 * cc)
            out[i] = np.clip(np.log((max(xs, m[i]) + eps[i]) / (m[i] + eps[i])), 0, 5)
    return out


def _hist_mean(vals: np.ndarray, hist: list[np.ndarray]) -> np.ndarray:
    """Mean of vals (n_series,) or (n_series, T) over each series' history rows; NaN if no history."""
    out = np.full((len(hist),) + vals.shape[1:], np.nan)
    for i, h in enumerate(hist):
        if len(h):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                out[i] = np.nanmean(vals[h], axis=0)
    return out


# ---------------------------------------------------------------------------------------------------------------
# reflection-principle quantities for a Gaussian random walk with drift mu and step SD sigma, started a distance
# D >= 0 below its running maximum, over tau further steps

BG = 0.5826  # -zeta(1/2) / sqrt(2 pi)


def refl_p_exceed(D, mu, sigma, tau):
    """P(max of Brownian motion with drift mu, variance sigma^2 per week, over tau weeks >= D)."""
    from scipy.stats import norm

    s = sigma * np.sqrt(tau)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        t1 = norm.sf((D - mu * tau) / s)
        lt2 = 2 * mu * D / sigma ** 2 + norm.logsf((D + mu * tau) / s)
        return np.clip(t1 + np.exp(np.minimum(lt2, 50)), 0, 1)


def refl_ez0(D, sigma, tau):
    """E[(max_{u <= tau} W_u - D)^+] for driftless Brownian motion (max ~ |N(0, sigma^2 tau)|)."""
    from scipy.stats import norm

    s = sigma * np.sqrt(tau)
    a = D / s
    return 2 * (s * norm.pdf(a) - D * norm.sf(a))


def verify_reflection(n_sim=40000, seed=1) -> pd.DataFrame:
    """Compare the closed forms with simulated Gaussian random walks: weekly steps (what the data are) and fine
    steps (100 per week, approximating Brownian motion)."""
    rng = np.random.default_rng(seed)
    rows = []
    for D, mu, sigma, tau in [(0.5, 0.0, 0.3, 4), (1.0, 0.0, 0.3, 10), (0.3, 0.1, 0.25, 6), (1.0, -0.1, 0.3, 15),
                              (0.2, 0.2, 0.2, 3), (2.0, 0.05, 0.4, 20)]:
        out = {"D": D, "mu": mu, "sigma": sigma, "tau": tau,
               "P formula": float(refl_p_exceed(D, mu, sigma, tau)),
               # Broadie-Glasserman continuity correction for weekly (discrete) monitoring: shift D by 0.5826 sigma
               "P formula BG": float(refl_p_exceed(D + BG * sigma, mu, sigma, tau))}
        for label, k in [("weekly", 1), ("fine", 100)]:
            steps = rng.normal(mu / k, sigma / np.sqrt(k), size=(n_sim, tau * k)).astype(np.float32)
            mx = np.maximum(np.cumsum(steps, axis=1).max(axis=1), 0)
            out[f"P sim {label}"] = float((mx >= D).mean())
            if mu == 0:
                out[f"E sim {label}"] = float(np.maximum(mx - D, 0).mean())
        if mu == 0:
            out["E formula"] = float(refl_ez0(D, sigma, tau))
            out["E formula BG"] = float(refl_ez0(D + BG * sigma, sigma, tau))
        rows.append(out)
    return pd.DataFrame(rows)


def holiday_weeks(season: str) -> list[int]:
    """Season weeks whose Saturday week-ending date falls in Dec 22 - Jan 7."""
    y0 = int(season[:4])
    lo, hi = pd.Timestamp(y0, 12, 22), pd.Timestamp(y0 + 1, 1, 7)
    return [w for w in range(15, 30) if lo <= pd.Timestamp(season_week_to_date(season, w)) <= hi]


def mmwr_week(d) -> int:
    """MMWR week of a Saturday week-ending date (week 1 ends on the first Saturday on or after Jan 4)."""
    d = pd.Timestamp(d)

    def first_end(year):
        j4 = pd.Timestamp(year, 1, 4)
        return j4 + pd.Timedelta(days=(5 - j4.weekday()) % 7)

    fe = first_end(d.year)
    if d < fe:
        fe = first_end(d.year - 1)
    return int((d - fe).days // 7 + 1)


def holiday_matrix(keys: pd.DataFrame) -> np.ndarray:
    """H[i, j] is True when season week j + 1 of series i is a holiday week."""
    cal = {s: holiday_weeks(s) for s in keys["season"].unique()}
    H = np.zeros((len(keys), 53), dtype=bool)
    for i, s in enumerate(keys["season"]):
        H[i, np.array(cal[s]) - 1] = True
    return H


def candidate_features(arrays) -> pd.DataFrame:
    """
    One row per (series, season week t in replay_start..W1) with candidate features. Current-season quantities are
    computed from a copy of the series with every week after t set to missing, so no future data can enter.
    Historical quantities use only strictly earlier seasons of the same source and location (the climatological
    peak week uses all series of strictly earlier seasons).
    """
    n_obs = np.sum(~np.isnan(arrays.y[:, W0 - 1:W1]), axis=1)
    arrays = arrays.subset(n_obs >= MIN_OBS)  # the same series as build_replay_rows
    keys, y = arrays.keys, arrays.y
    eps = arrays.eps[:, None]
    n = len(keys)

    # completed-season quantities, used only through the history of later seasons
    peak, peak_week = season_peaks(y, W0, W1)
    yf_full = _ffill(y)
    L_full = np.log(yf_full + eps)
    M_full = np.column_stack([running_max(y, t, W0)[0] for t in range(1, W1 + 1)])
    lM_full = np.log(M_full + eps)
    z_full = np.log(peak + eps[:, 0])[:, None] - lM_full
    C_full = np.column_stack([np.nansum(y[:, 4:t], axis=1) if t > 4 else np.nansum(y[:, :t], axis=1)
                              for t in range(1, W1 + 1)])
    total = np.nansum(y[:, 4:W1], axis=1)

    hist = history_index(keys)
    hist_pk_mean = _hist_mean(peak_week, hist)
    hist_pk_med = np.array([np.median(peak_week[h]) if len(h) else np.nan for h in hist])
    hist_pk_sd = np.array([np.std(peak_week[h]) if len(h) >= 2 else np.nan for h in hist])
    hist_lpk = _hist_mean(np.log(peak + eps[:, 0]), hist)
    hist_total = _hist_mean(total, hist)
    hist_C = _hist_mean(C_full, hist)
    hist_L = _hist_mean(L_full[:, :W1], hist)
    hist_lM = _hist_mean(lM_full, hist)
    hist_z = _hist_mean(z_full, hist)
    seasons = keys["season"].to_numpy()
    clim_pk = np.array([np.nanmean(peak_week[seasons < s]) if (seasons < s).any() else np.nan for s in seasons])

    tay = {(1, 4): _taylor_pinv(4, 1), (1, 6): _taylor_pinv(6, 1), (2, 6): _taylor_pinv(6, 2),
           (2, 8): _taylor_pinv(8, 2)}
    sync_group = (keys["source"] + "|" + keys["season"]).to_numpy()
    not_nat = (keys["agg_level"] != "national").to_numpy()
    H = holiday_matrix(keys)
    last_hol = H.shape[1] - np.argmax(H[:, ::-1], axis=1)  # last holiday season week of each series
    first_hol = np.argmax(H, axis=1) + 1

    # section 11 history: weekly-change SD, chain-ladder z and share peaked, post-peak decline rates
    Lraw_full = np.log(y + eps)
    sig_full = np.column_stack([robust_sd_diffs(Lraw_full, t) for t in range(1, W1 + 1)])
    hist_sig = _hist_mean(sig_full, hist)
    hist_pk_q90 = np.array([np.quantile(peak_week[h], 0.9) if len(h) else np.nan for h in hist])
    peaked_by = (peak_week[:, None] <= np.arange(1, W1 + 1)[None, :]).astype(float)
    hist_clz = _hist_median(z_full, hist)
    hist_p0 = _hist_mean(peaked_by, hist)
    n_hist = np.array([len(h) for h in hist])
    src_arr = keys["source"].to_numpy()
    pooled_clz = np.full_like(hist_clz, np.nan)
    pooled_p0 = np.full_like(hist_p0, np.nan)
    for (src, ssn), idx in keys.groupby(["source", "season"]).groups.items():
        prev = np.flatnonzero((src_arr == src) & (seasons < ssn))
        if len(prev):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                pooled_clz[idx] = np.nanmedian(z_full[prev], axis=0)
                pooled_p0[idx] = np.nanmean(peaked_by[prev], axis=0)
    use_loc = (n_hist >= 3)[:, None]
    clz = np.where(use_loc, hist_clz, pooled_clz)
    p0h = np.where(use_loc, hist_p0, pooled_p0)
    pk_idx = np.nan_to_num(peak_week, nan=W0).astype(int) - 1
    pk_l = L_full[np.arange(n), pk_idx]
    decl_full = np.full((n, 21), np.nan)
    for k in range(1, 21):
        j = np.minimum(pk_idx + k, L_full.shape[1] - 1)
        decl_full[:, k] = np.where(pk_idx + k < L_full.shape[1], (pk_l - L_full[np.arange(n), j]) / k, np.nan)
    hist_decl = _hist_median(decl_full, hist)
    # cross-source partner: the same location and season in the other of ILINet / FluSurv-NET
    kk = keys.reset_index()
    part = kk.merge(kk.assign(source=kk["source"].map({"ilinet": "flusurvnet", "flusurvnet": "ilinet"})),
                    on=["source", "location", "season"], how="left", suffixes=("", "_p"))
    partner = part["index_p"].fillna(-1).astype(int).to_numpy()
    sev_groups = [np.flatnonzero((sync_group == g) & not_nat) for g in np.unique(sync_group)]
    st_lat, st_lon, is_state, nbr_idx, hhs_idx, inv_w = regional_structure(keys)

    frames = []
    for t in range(REPLAY_START, W1 + 1):
        yt = y.copy()
        yt[:, t:] = np.nan  # nothing after week t
        yf = _ffill(yt[:, :t])
        L = np.log(yf + eps)
        lx = L[:, t - 1]
        m, m_week = running_max(yt, t, W0)
        lm = np.log(m + eps[:, 0])
        f = {"season_week": np.full(n, float(t)), "rel_max_chk": lx - lm, "wsm_chk": t - m_week}
        # a/b. local Taylor fits on the log scale, rolling means and lags, relative to the running max
        for (deg, w), P in tay.items():
            if t >= w:
                B = L[:, t - w:t] @ P.T  # (n, deg + 1)
            else:
                B = np.full((n, deg + 1), np.nan)
            f[f"tay{deg}_w{w}_lvl"] = B[:, 0] - lm
            f[f"tay{deg}_w{w}_slope"] = B[:, 1]
            if deg == 2:
                f[f"tay{deg}_w{w}_curv"] = B[:, 2]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            for k in (2, 4):
                f[f"rm{k}"] = np.log(np.nanmean(yf[:, max(t - k, 0):t], axis=1) + eps[:, 0]) - lm
        for lag in range(1, 5):
            f[f"lag{lag}_rel"] = (L[:, t - 1 - lag] - lm) if t - 1 - lag >= 0 else np.full(n, np.nan)
        # c. timing relative to history
        f["t_minus_hist_pk"] = t - hist_pk_mean
        f["t_minus_hist_pk_med"] = t - hist_pk_med
        f["t_minus_clim_pk"] = t - clim_pk
        f["hist_pk_sd"] = hist_pk_sd
        f["hist_z_same_week"] = hist_z[:, t - 1]
        # e. onset: baseline = mean of weeks 1..min(t, 8); onset = first week >= 9 at 2x / 4x the baseline
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            lb = np.log(np.nanmean(yt[:, :min(t, 8)], axis=1) + eps[:, 0])
        rel_base = np.log(yt[:, :t] + eps) - lb[:, None]
        for mult in (2, 4):
            above = np.zeros_like(rel_base, dtype=bool)
            above[:, 8:] = rel_base[:, 8:] >= np.log(mult)  # NaN compares False
            first = np.where(above.any(axis=1), above.argmax(axis=1) + 1, 0)
            f[f"wks_since_onset{mult}"] = np.where(first > 0, t - first, -1.0)
        f["max_over_base"] = lm - lb
        f["cur_over_base"] = lx - lb
        # f. burden to date vs history
        cum = np.nansum(yt[:, 4:t], axis=1) if t > 4 else np.nansum(yt[:, :t], axis=1)
        f["cum_vs_hist_total"] = np.log(cum + eps[:, 0]) - np.log(hist_total + eps[:, 0])
        f["cum_vs_hist_same_week"] = np.log(cum + eps[:, 0]) - np.log(hist_C[:, t - 1] + eps[:, 0])
        # g. level vs history
        f["cur_vs_hist_peak"] = lx - hist_lpk
        f["lvl_vs_hist_same_week"] = lx - hist_L[:, t - 1]
        f["max_vs_hist_same_week"] = lm - hist_lM[:, t - 1]
        # d. synchrony across the other non-national series of the same source and season, as observed at t
        obs = ~np.isnan(yt[:, t - 1])
        rel = lx - lm
        wsm = t - m_week
        df_s = pd.DataFrame({"g": sync_group, "use": obs & not_nat, "past2": (wsm >= 2) & (t >= W0),
                             "half": rel >= np.log(0.5), "rel": rel,
                             "g3": lx - L[:, t - 4] if t >= 4 else np.nan})
        agg = df_s[df_s["use"]].groupby("g").agg(sync_frac_past2=("past2", "mean"), sync_frac_half=("half", "mean"),
                                                 sync_med_rel_max=("rel", "median"), sync_med_g3=("g3", "median"))
        for c in agg.columns:
            f[c] = pd.Series(sync_group).map(agg[c]).to_numpy()
        # h. holiday weeks (calendar is known in advance; the non-holiday max uses in-window weeks <= t only)
        f["hol_now"] = H[:, t - 1].astype(float)
        mw = m_week.astype(int)
        f["max_in_holiday"] = (H[np.arange(n), mw - 1] & (t >= W0)).astype(float)
        f["weeks_since_holiday_end"] = np.where(t < first_hol, -1.0,
                                                np.where(t <= last_hol, 0.0, np.minimum(t - last_hol, 10)))
        if t >= W0:
            w = yt[:, W0 - 1:t].copy()
            w[H[:, W0 - 1:t]] = np.nan
            has = ~np.all(np.isnan(w), axis=1)
            m_nh = np.full(n, np.nan)
            wk_nh = np.full(n, np.nan)
            m_nh[has] = np.nanmax(w[has], axis=1)
            wk_nh[has] = W0 + np.nanargmax(w[has], axis=1)
        else:
            m_nh, wk_nh = m, m_week
        lm_nh = np.log(m_nh + eps[:, 0])
        f["holiday_excess_max"] = lm - lm_nh
        f["rel_max_nohol"] = lx - lm_nh
        f["wks_since_max_nohol"] = t - wk_nh
        # section 11 -------------------------------------------------------------------------------------------
        Lraw = np.log(yt[:, :t] + eps)
        D = lm - lx
        r = (lx - L[:, t - 4]) / 3 if t >= 4 else np.full(n, np.nan)
        s_own = robust_sd_diffs(Lraw, t)
        s_hist = hist_sig[:, t - 1]
        sig = np.where(np.isnan(s_hist), s_own, np.where(np.isnan(s_own), s_hist, 0.5 * (s_own + s_hist)))
        sig = np.maximum(sig, 0.05)
        tau = max(W1 - t, 1)
        tau_alt = np.maximum(hist_pk_q90 - t, 1)
        Dc = D + BG * sig  # continuity correction for weekly data
        f["sigma_hat"] = sig
        f["snr"] = r / sig
        f["p_exc0"] = refl_p_exceed(Dc, 0.0, sig, tau)
        f["p_exc"] = refl_p_exceed(Dc, np.nan_to_num(r), sig, tau)
        f["ez_rw"] = refl_ez0(Dc, sig, tau)
        f["p_exc0_alt"] = refl_p_exceed(Dc, 0.0, sig, tau_alt)
        f["ez_rw_alt"] = refl_ez0(Dc, sig, tau_alt)
        # severity realized so far in other locations of the same source and season that have visibly peaked
        hr = lm - hist_lpk
        elig = obs & (wsm >= 2) & (t >= W0) & ~np.isnan(hr)
        sev = np.full(n, np.nan)
        npk = np.full(n, np.nan)
        for g in sev_groups:
            e = elig[g]
            for ii, i in enumerate(g):
                others = e.copy()
                others[ii] = False
                npk[i] = others.sum() / max(len(g) - 1, 1)
                if others.any():
                    sev[i] = np.median(hr[g][others])
        f["sev_peaked"] = sev
        f["sev_gap"] = sev - hr
        f["n_peaked"] = npk
        # chain ladder
        f["cl_z"] = clz[:, t - 1]
        f["p0_hist"] = p0h[:, t - 1]
        f["z_bf"] = np.log1p((1 - np.exp(-clz[:, t - 1])) * np.exp(-hr))
        # deceleration of the 3-week growth rate r_s = (l_s - l_{s-3}) / 3, s = 5..t
        if t >= REPLAY_START:
            R = (L[:, 4:t] - L[:, 1:t - 3]) / 3
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                rmax = np.nanmax(R, axis=1)
            allnan = np.all(np.isnan(R), axis=1)
            tg = np.where(allnan, np.nan, REPLAY_START + np.nanargmax(np.where(np.isnan(R), -np.inf, R), axis=1))
        else:
            rmax, tg = np.full(n, np.nan), np.full(n, np.nan)
        dec = (rmax - r) / np.maximum(t - tg, 1)
        f["r_max"] = rmax
        with np.errstate(divide="ignore", invalid="ignore"):
            f["gfrac"] = np.where(rmax > 0, r / rmax, np.nan)
        f["wks_since_gmax"] = t - tg
        f["g_decel"] = dec
        with np.errstate(divide="ignore", invalid="ignore"):
            f["z_par"] = np.where((r > 0) & (dec > 0), np.minimum(r ** 2 / (2 * dec), 5), 0.0)
        # recession since the running max
        f["rec_rate"] = D / np.maximum(wsm, 1)
        neg = np.zeros((n, t + 1))
        neg[:, 2:t + 1] = np.cumsum((L[:, 1:t] - L[:, :t - 1]) < 0, axis=1) if t >= 2 else 0
        mwi = m_week.astype(int)
        with np.errstate(divide="ignore", invalid="ignore"):
            f["rec_consist"] = np.where(wsm > 0, (neg[:, t] - neg[np.arange(n), mwi]) / wsm, np.nan)
        kd = np.clip(wsm, 1, 20).astype(int)
        f["rec_anom"] = np.where(wsm > 0, f["rec_rate"] - hist_decl[np.arange(n), kd], np.nan)
        # records: weeks (from week 10) that set a new in-window maximum
        if t >= W0:
            w = yt[:, W0 - 1:t]
            prev = np.fmax.accumulate(np.concatenate([np.full((n, 1), np.nan), w[:, :-1]], axis=1), axis=1)
            recs = ~np.isnan(w) & (np.isnan(prev) | (w > prev))
            f["n_rec4"] = recs[:, -4:].sum(axis=1).astype(float)
            f["rec_frac"] = recs.sum(axis=1) / (t - W0 + 1)
        else:
            f["n_rec4"] = np.full(n, np.nan)
            f["rec_frac"] = np.full(n, np.nan)
        f["bass_rel"] = bass_implied(yt, t, eps[:, 0], m)
        # cross-source discrepancy (partner observed at t)
        hp = (partner >= 0) & obs[np.maximum(partner, 0)]
        pi = np.maximum(partner, 0)
        g3 = lx - L[:, t - 4] if t >= 4 else np.full(n, np.nan)
        f["xs_rel_max"] = np.where(hp, rel[pi] - rel, np.nan)
        f["xs_wsm"] = np.where(hp, wsm[pi] - wsm, np.nan)
        f["xs_g3"] = np.where(hp, g3[pi] - g3, np.nan)
        # section 12: neighbour / regional synchrony among state series of the same source and season observed at t
        past2 = ((wsm >= 2) & (t >= W0)).astype(float)
        g3v = lx - L[:, t - 4] if t >= 4 else np.full(n, np.nan)
        cols = {c: np.full(n, np.nan) for c in ALL_REG}
        state_share = {}
        for i in np.flatnonzero(is_state):
            nb = nbr_idx[i][obs[nbr_idx[i]]]
            if len(nb):
                cols["nbr5_frac_past2"][i] = past2[nb].mean()
                cols["nbr5_med_rel_max"][i] = np.median(rel[nb])
                cols["nbr5_med_wsm"][i] = np.median(wsm[nb])
                cols["nbr5_med_g3"][i] = np.nanmedian(g3v[nb]) if np.any(~np.isnan(g3v[nb])) else np.nan
            rg = hhs_idx[i][obs[hhs_idx[i]]]
            if len(rg):
                cols["hhs_frac_past2"][i] = past2[rg].mean()
                cols["hhs_med_rel_max"][i] = np.median(rel[rg])
            oi, w = inv_w[i]
            ok = obs[oi]
            if ok.any():
                # distance-weighted share of other states past their max, minus the unweighted share of all of them
                cols["wave_lead"][i] = np.sum(w[ok] * past2[oi[ok]]) / np.sum(w[ok]) - past2[oi[ok]].mean()
            cols["lat"][i] = st_lat[i]
            cols["lon"][i] = st_lon[i]
        f.update(cols)
        fr = pd.concat([keys, pd.DataFrame(f)], axis=1)
        frames.append(fr.loc[obs])
    return pd.concat(frames, ignore_index=True)


def load_rows_adjusted(hol_adjust: dict | str | None = None):
    """
    As relative_size_eda.load_rows, with the ILINet values in holiday weeks adjusted before anything is computed:
    hol_adjust = {k: v} divides the k-th holiday week (k = 0, 1, 2) of each season by exp(v); hol_adjust = "interp"
    replaces the holiday weeks by log-linear interpolation between the weeks just before and just after the block.
    """
    from relative_size_eda import DROP_ILINET_LOCATIONS, ILI_PATH, group_of
    from idmodels.peak.series import LOG_EPS

    data = pd.read_parquet(ILI_PATH)
    data = data.loc[~((data["source"] == "ilinet") & data["location"].isin(DROP_ILINET_LOCATIONS))]
    arrays = build_season_arrays(data)
    if hol_adjust == "interp":
        H = holiday_matrix(arrays.keys)
        eps = arrays.eps
        for i in np.flatnonzero((arrays.keys["source"] == "ilinet").to_numpy()):
            hw = np.flatnonzero(H[i])
            a, b = hw[0] - 1, hw[-1] + 1
            ya, yb = arrays.y[i, a], arrays.y[i, b]
            if np.isnan(ya) or np.isnan(yb):
                continue
            la, lb = np.log(ya + eps[i]), np.log(yb + eps[i])
            frac = (hw - a) / (b - a)
            arrays.y[i, hw] = np.exp(la + frac * (lb - la)) - eps[i]
    elif hol_adjust:
        H = holiday_matrix(arrays.keys)
        ili = (arrays.keys["source"] == "ilinet").to_numpy()
        pos = np.cumsum(H, axis=1) - 1  # position within the holiday block
        for k, v in hol_adjust.items():
            sel = H & (pos == k) & ili[:, None]
            arrays.y[sel] = arrays.y[sel] * np.exp(-v)
    rows = build_replay_rows(arrays, W0, W1, REPLAY_START, MIN_OBS)
    rows["group"] = group_of(rows["source"], rows["agg_level"])
    rows["eps"] = rows["source"].map(LOG_EPS)
    rows["M"] = np.exp(rows["lm"]) - rows["eps"]
    rows["at_zero"] = rows["z"].abs() < 1e-9
    rows["log_peak"] = rows["z"] + rows["lm"]
    return arrays, rows


def load_strain():
    """Weekly A, B, H1, H3 positives by geography ('US', 'Region k', state FIPS) and season, season weeks 1..53."""
    if not STRAIN_PATH.exists():
        STRAIN_PATH.parent.mkdir(parents=True, exist_ok=True)
        import urllib.request
        urllib.request.urlretrieve(STRAIN_URL, STRAIN_PATH)
    x = pd.read_csv(STRAIN_PATH)
    x["geo"] = np.where(x["region_type"] == "National", "US",
                        np.where(x["region_type"] == "HHS Regions", x["region"], x["region"].map(STATE_FIPS)))
    x = x.dropna(subset=["geo"])
    out = {}
    for (geo, season), g in x.groupby(["geo", "season"]):
        arr = np.full((4, N_SEASON_WEEKS_T), np.nan)
        wk = g["season_week"].to_numpy().astype(int) - 1
        for r, c in enumerate(["a", "b", "a_h1", "a_h3"]):
            arr[r, wk] = g[c].to_numpy(dtype=float)
        out[(geo, season)] = arr
    return out


N_SEASON_WEEKS_T = 53


def type_arrays(keys: pd.DataFrame, strain: dict) -> dict:
    """For each series: A and B arrays of its own geography (state / HHS region / nation; FluSurv-NET sites use their
    state), and H1, H3 arrays of its HHS region (the nation for national series), plus the national H1, H3."""
    n = len(keys)
    out = {k: np.full((n, N_SEASON_WEEKS_T), np.nan) for k in ["A", "B", "H1", "H3", "H1n", "H3n"]}
    for i, (loc, agg, ssn) in enumerate(zip(keys["location"], keys["agg_level"], keys["season"])):
        own = "US" if loc == "US" else (loc if agg == "hhs region" else loc)
        reg = "US" if loc == "US" else (loc if agg == "hhs region" else
                                        (f"Region {HHS_REGION[loc]}" if loc in HHS_REGION else None))
        if (own, ssn) in strain:
            out["A"][i], out["B"][i] = strain[(own, ssn)][0], strain[(own, ssn)][1]
        if reg is not None and (reg, ssn) in strain:
            out["H1"][i], out["H3"][i] = strain[(reg, ssn)][2], strain[(reg, ssn)][3]
        if ("US", ssn) in strain:
            out["H1n"][i], out["H3n"][i] = strain[("US", ssn)][2], strain[("US", ssn)][3]
    return out


def type_features(keys: pd.DataFrame, T: dict) -> pd.DataFrame:
    """Type / subtype features for every series and season week t = REPLAY_START..W1, from weeks <= t only."""
    n = len(keys)
    frames = []

    def wsum(x, lo, hi):  # sum over season weeks lo..hi (1-based, inclusive), NaN if all missing
        lo = max(lo, 1)
        if hi < lo:
            return np.full(n, np.nan)
        v = x[:, lo - 1:hi]
        return np.where(np.all(np.isnan(v), axis=1), np.nan, np.nansum(v, axis=1))

    def share(num, den):
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(den >= MIN_POS, num / den, np.nan)

    for t in range(REPLAY_START, W1 + 1):
        Tt = {k: v.copy() for k, v in T.items()}
        for v in Tt.values():
            v[:, t:] = np.nan  # nothing after week t
        A, B = Tt["A"], Tt["B"]
        A3, B3 = wsum(A, t - 2, t), wsum(B, t - 2, t)
        A3p, B3p = wsum(A, t - 5, t - 3), wsum(B, t - 5, t - 3)
        f = {"season_week": np.full(n, float(t))}
        f["b_share_3wk"] = share(B3, A3 + B3)
        f["b_share_cum"] = share(wsum(B, 5, t), wsum(A, 5, t) + wsum(B, 5, t))
        f["b_share_trend"] = f["b_share_3wk"] - share(B3p, A3p + B3p)
        f["b_rising"] = np.log((B3 + 1) / (B3p + 1))
        f["a_rising"] = np.log((A3 + 1) / (A3p + 1))
        f["b_minus_a_growth"] = f["b_rising"] - f["a_rising"]
        roll = np.column_stack([wsum(B, s - 2, s) for s in range(3, t + 1)]) if t >= 3 else np.full((n, 1), np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            bmax = np.nanmax(roll, axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            f["b_frac_of_peak"] = np.where(bmax > 0, B3 / bmax, np.nan)
        H1c, H3c = wsum(Tt["H1"], 5, t), wsum(Tt["H3"], 5, t)
        f["h3_share_cum"] = share(H3c, H1c + H3c)
        f["h3_share_3wk"] = share(wsum(Tt["H3"], t - 2, t), wsum(Tt["H1"], t - 2, t) + wsum(Tt["H3"], t - 2, t))
        H1n, H3n = wsum(Tt["H1n"], 5, t), wsum(Tt["H3n"], 5, t)
        f["h3_share_nat"] = share(H3n, H1n + H3n)
        frames.append(pd.concat([keys, pd.DataFrame(f)], axis=1))
    return pd.concat(frames, ignore_index=True)


def type_leakage_check(arrays, weeks=(14, 22, 30)) -> dict:
    """Corrupt the type / subtype series after week t (x 100 + 5): type features at t must not change."""
    n_obs = np.sum(~np.isnan(arrays.y[:, W0 - 1:W1]), axis=1)
    keys = arrays.subset(n_obs >= MIN_OBS).keys
    T = type_arrays(keys, load_strain())
    base = type_features(keys, T).set_index(KEYS)
    out = {}
    for t in weeks:
        T2 = {k: v.copy() for k, v in T.items()}
        for v in T2.values():
            v[:, t:] = v[:, t:] * 100 + 5
        f2 = type_features(keys, T2)
        f2 = f2[f2["season_week"] == t].set_index(KEYS)
        b = base.loc[f2.index, ALL_TYPE]
        out[f"type_future_weeks_t{t}"] = float(np.nanmax(np.abs((b - f2[ALL_TYPE]).to_numpy())))
    return out


def build_dataset(hol_adjust=None):
    arrays, rows = load_rows_adjusted(hol_adjust) if hol_adjust else load_rows(keep_extra=True)
    feats = candidate_features(arrays)
    # newer idmodels versions compute some of these features themselves: keep ours, and record the agreement
    overlap = [c for c in feats.columns if c in rows.columns and c not in KEYS]
    rows = rows.rename(columns={c: f"{c}_idm" for c in overlap})
    d = rows.merge(feats, on=KEYS, how="left", validate="one_to_one")
    idm_diff = {c: float(np.nanmax(np.abs(d[c] - d[f"{c}_idm"]))) for c in overlap}
    idm_nan_mismatch = {c: int((d[c].isna() != d[f"{c}_idm"].isna()).sum()) for c in overlap}
    d = d.drop(columns=[f"{c}_idm" for c in overlap])
    # consistency with the model code: our recomputed rel_max and wks_since_max must match exactly
    chk_rel = np.nanmax(np.abs(d["rel_max"] - d["rel_max_chk"]))
    chk_wsm = np.nanmax(np.abs(d["wks_since_max"] - d["wsm_chk"]))
    assert chk_rel < 1e-9 and chk_wsm < 1e-9, (chk_rel, chk_wsm)
    assert d[ALL_NEW[0]].notna().any() and len(d) == len(rows)
    n_obs = np.sum(~np.isnan(arrays.y[:, W0 - 1:W1]), axis=1)
    keys_c = arrays.subset(n_obs >= MIN_OBS).keys
    d = d.merge(type_features(keys_c, type_arrays(keys_c, load_strain())), on=KEYS, how="left", validate="one_to_one")
    d["pos"] = (~d["at_zero"]).astype(int)
    return d.drop(columns=["rel_max_chk", "wsm_chk"]), {"max_abs_diff_rel_max": float(chk_rel),
                                                        "idmodels_feature_max_abs_diff": idm_diff,
                                                        "idmodels_feature_nan_mismatch": idm_nan_mismatch}


def leakage_check(d: pd.DataFrame, arrays, weeks=(14, 22, 30), season="2015/16") -> dict:
    """
    Two perturbation checks, returning the largest absolute change in the features that must not move:
      1. corrupt every week after t in every series (x100 + 5): current-season features at t must not change;
      2. corrupt every series of `season` and later seasons entirely: historical features of `season` rows must not
         change (history uses strictly earlier seasons only).
    """
    out = {}
    cur_only = (GROUPS_NEW["trend"] + GROUPS_NEW["onset"] + GROUPS_NEW["synchrony"] + HOLIDAY + ALL_REG
                + ["r_max", "gfrac", "wks_since_gmax", "g_decel", "z_par", "rec_rate", "rec_consist", "n_rec4",
                   "rec_frac", "bass_rel", "xs_rel_max", "xs_wsm", "xs_g3"])
    base = d.set_index(KEYS)
    for t in weeks:
        a2 = type(arrays)(keys=arrays.keys, y=arrays.y.copy())
        a2.y[:, t:] = a2.y[:, t:] * 100 + 5
        f2 = candidate_features(a2)
        f2 = f2[f2["season_week"] == t].set_index(KEYS)
        idx = f2.index.intersection(base.index)
        out[f"future_weeks_t{t}"] = float(np.nanmax((base.loc[idx, cur_only] - f2.loc[idx, cur_only]).abs().to_numpy()))
    a2 = type(arrays)(keys=arrays.keys, y=arrays.y.copy())
    later = (arrays.keys["season"] >= season).to_numpy()
    a2.y[later] = a2.y[later] * 100 + 5
    f2 = candidate_features(a2)
    f2 = f2[f2["season"] == season].set_index(KEYS)
    b = base[base.index.get_level_values("season") == season]
    idx = f2.index.intersection(b.index)
    pure_hist = GROUPS_NEW["timing"] + ["cl_z", "p0_hist"]  # purely historical features
    out["later_seasons_hist_features"] = float(np.nanmax((b.loc[idx, pure_hist] - f2.loc[idx, pure_hist]).abs().to_numpy()))
    # 3. corrupt only the weeks after t of one season: every feature of that season at t must be unchanged
    allf = ALL_NEW + HOLIDAY + ALL_LIT + ALL_REG
    for t in weeks:
        a2 = type(arrays)(keys=arrays.keys, y=arrays.y.copy())
        rows_s = (arrays.keys["season"] == season).to_numpy()
        a2.y[np.ix_(rows_s, np.arange(t, a2.y.shape[1]))] = a2.y[np.ix_(rows_s, np.arange(t, a2.y.shape[1]))] * 100 + 5
        f2 = candidate_features(a2)
        f2 = f2[(f2["season_week"] == t) & (f2["season"] == season)].set_index(KEYS)
        idx = f2.index.intersection(base.index)
        diff = (base.loc[idx, allf] - f2.loc[idx, allf]).abs().to_numpy()
        out[f"season_{season}_future_weeks_t{t}_all_features"] = float(np.nanmax(diff))
    return out


# ---------------------------------------------------------------------------------------------------------------
# evaluation

def pinball(y, q, tau):
    d = y - q
    return np.maximum(tau * d, (tau - 1) * d)


def regional_structure(keys: pd.DataFrame):
    """For each state series: indices of its 5 nearest other state series (same source and season), of the other
    state series in its HHS region, and inverse-squared-distance weights over all other state series of the source-season."""
    n = len(keys)
    lat = np.array([STATE_CENTROIDS.get(l, (np.nan, np.nan))[0] for l in keys["location"]])
    lon = np.array([STATE_CENTROIDS.get(l, (np.nan, np.nan))[1] for l in keys["location"]])
    is_state = (keys["agg_level"] == "state").to_numpy() & ~np.isnan(lat)
    hhs = np.array([HHS_REGION.get(l, -1) for l in keys["location"]])
    grp = (keys["source"] + "|" + keys["season"]).to_numpy()
    nbr, reg, inv = [np.array([], int)] * n, [np.array([], int)] * n, [(np.array([], int), np.array([]))] * n
    r = np.pi / 180
    for g in np.unique(grp[is_state]):
        idx = np.flatnonzero((grp == g) & is_state)
        la, lo = lat[idx] * r, lon[idx] * r
        D = 6371 * np.arccos(np.clip(np.sin(la)[:, None] * np.sin(la)[None, :] +
                                     np.cos(la)[:, None] * np.cos(la)[None, :] * np.cos(lo[:, None] - lo[None, :]), -1, 1))
        for a, i in enumerate(idx):
            order = [b for b in np.argsort(D[a]) if b != a]
            nbr[i] = idx[order[:5]]
            reg[i] = idx[[b for b in range(len(idx)) if b != a and hhs[idx[b]] == hhs[i]]]
            others = np.array([b for b in range(len(idx)) if b != a], dtype=int)
            inv[i] = (idx[others], 1.0 / np.maximum(D[a, others], 50.0) ** 2)  # inverse squared distance (km)
    return lat, lon, is_state, nbr, reg, inv


def _fit_fold(tr_X, tr_z, tr_pos, te_X):
    import lightgbm as lgb

    params = dict(LGB_PARAMS, n_jobs=1)
    out = {}
    for tau in QS:
        m = lgb.LGBMRegressor(objective="quantile", alpha=tau, **params)
        m.fit(tr_X, tr_z)
        out[f"q{tau}"] = m.predict(te_X)
    c = lgb.LGBMClassifier(objective="binary", **params)
    c.fit(tr_X, tr_pos)
    out["p_pos"] = np.clip(c.predict_proba(te_X)[:, 1], 1e-6, 1 - 1e-6)
    return out


def cv_feature_sets(d: pd.DataFrame, sets: dict, quick=False, n_jobs=None) -> pd.DataFrame:
    """Leave-one-season-out: each season (all sources and locations) is held out in turn."""
    from joblib import Parallel, delayed

    seasons = sorted(d["season"].unique())
    if quick:
        seasons = seasons[::4]
    jobs = [(name, s) for name in sets for s in seasons]
    t0 = time.time()
    res = Parallel(n_jobs=n_jobs or N_WORKERS)(
        delayed(_fit_fold)(d.loc[d["season"] != s, sets[name]], d.loc[d["season"] != s, "z"],
                           d.loc[d["season"] != s, "pos"], d.loc[d["season"] == s, sets[name]])
        for name, s in jobs)
    print(f"cv fits: {time.time() - t0:.0f}s", flush=True)
    out = []
    for (name, s), r in zip(jobs, res):
        p = d.loc[d["season"] == s, KEYS + ["z", "pos"]].copy()
        for k, v in r.items():
            p[k] = v
        p["pinball"] = np.mean([pinball(p["z"].to_numpy(), p[f"q{tau}"].to_numpy(), tau) for tau in QS], axis=0)
        p["logloss"] = -(p["pos"] * np.log(p["p_pos"]) + (1 - p["pos"]) * np.log(1 - p["p_pos"]))
        p["set"] = name
        out.append(p[KEYS + ["z", "pos", "pinball", "logloss", "set"]])
    return pd.concat(out, ignore_index=True)


def season_bootstrap(cv, a, b, metric="pinball", n_boot=2000, seed=0):
    """Share of season-level bootstrap resamples in which set a beats set b (lower mean loss), rows weeks 12-31."""
    x = cv[(cv["season_week"] >= 12) & (cv["season_week"] <= 31)]
    s = x.pivot_table(index=[*KEYS], columns="set", values=metric)[[a, b]].dropna()
    per = s.groupby(level="season").agg(["sum", "count"])
    seasons = per.index.to_numpy()
    rng = np.random.default_rng(seed)
    wins = 0
    for _ in range(n_boot):
        pick = rng.choice(seasons, size=len(seasons), replace=True)
        pp = per.loc[pick]
        wins += (pp[(a, "sum")].sum() / pp[(a, "count")].sum()) < (pp[(b, "sum")].sum() / pp[(b, "count")].sum())
    return wins / n_boot


def per_season_rel(cv, a, b, metric="pinball"):
    """Per held-out season (weeks 12-31): mean loss of set a divided by that of set b."""
    x = cv[(cv["season_week"] >= 12) & (cv["season_week"] <= 31)]
    g = x.groupby(["set", "season"])[metric].mean().unstack("set")
    return (g[a] / g[b])



# ---------------------------------------------------------------------------------------------------------------
# holiday weeks (section 10)

def holiday_excess(arrays) -> pd.DataFrame:
    """
    Per series: log value in each holiday week minus the mean log value of the 2 weeks before and the 2 weeks after
    the season's holiday block (observed values only), plus the same quantity for "placebo" blocks shifted by
    -6 and +6 weeks.
    """
    from relative_size_eda import group_of

    n_obs = np.sum(~np.isnan(arrays.y[:, W0 - 1:W1]), axis=1)
    a = arrays.subset(n_obs >= MIN_OBS)
    H = holiday_matrix(a.keys)
    L = np.log(a.y + a.eps[:, None])
    out = []
    for i in range(len(a.keys)):
        hw = np.flatnonzero(H[i])
        for shift, lab in [(0, "holiday"), (-6, "placebo −6 wk"), (6, "placebo +6 wk")]:
            b = hw + shift
            ref_idx = np.r_[b[0] - 2, b[0] - 1, b[-1] + 1, b[-1] + 2]
            ref = L[i, ref_idx]
            if np.isnan(ref).any() or np.isnan(L[i, b]).any():
                continue
            ex = L[i, b] - ref.mean()
            out.append({**a.keys.iloc[i].to_dict(), "block": lab, "excess": ex.mean(), "n_weeks": len(b),
                        **{f"excess_k{k}": ex[k] if k < len(ex) else np.nan for k in range(3)}})
    df = pd.DataFrame(out)
    df["group"] = group_of(df["source"], df["agg_level"])
    return df


def cv_fsn_target(dtrain: pd.DataFrame, sets: dict, n_jobs=None) -> pd.DataFrame:
    """Leave-one-season-out, but score only the FluSurv-NET rows of the held-out season (targets unaffected by any
    ILINet adjustment)."""
    from joblib import Parallel, delayed

    seasons = sorted(dtrain.loc[dtrain["source"] == "flusurvnet", "season"].unique())
    jobs = [(name, s) for name in sets for s in seasons]
    te_mask = {s: (dtrain["season"] == s) & (dtrain["source"] == "flusurvnet") for s in seasons}
    res = Parallel(n_jobs=n_jobs or N_WORKERS)(
        delayed(_fit_fold)(dtrain.loc[dtrain["season"] != s, sets[name]], dtrain.loc[dtrain["season"] != s, "z"],
                           dtrain.loc[dtrain["season"] != s, "pos"], dtrain.loc[te_mask[s], sets[name]])
        for name, s in jobs)
    out = []
    for (name, s), r in zip(jobs, res):
        p = dtrain.loc[te_mask[s], KEYS + ["z", "pos"]].copy()
        for k, v in r.items():
            p[k] = v
        p["pinball"] = np.mean([pinball(p["z"].to_numpy(), p[f"q{tau}"].to_numpy(), tau) for tau in QS], axis=0)
        p["logloss"] = -(p["pos"] * np.log(p["p_pos"]) + (1 - p["pos"]) * np.log(1 - p["p_pos"]))
        p["set"] = name
        out.append(p[KEYS + ["z", "pos", "pinball", "logloss", "set"]])
    return pd.concat(out, ignore_index=True)


def holiday_outputs(d, dcv, cv_main, nums, reuse=False):
    """Calendar, per-series holiday excess, holiday-feature CV and the ILINet-adjustment CV."""
    hn = {}
    cal = []
    for s in sorted(d["season"].unique()):
        for k, w in enumerate(holiday_weeks(s)):
            dt = season_week_to_date(s, w)
            cal.append({"season": s, "position": k, "season_week": w, "date": str(dt), "mmwr_week": mmwr_week(dt)})
    pd.DataFrame(cal).to_csv(OUT / "holiday_calendar.csv", index=False)
    arrays, _ = load_rows()
    ex = holiday_excess(arrays)
    ex.to_parquet(OUT / "holiday_excess.parquet")

    sets = {"current": CURRENT, "+ holiday": CURRENT + HOLIDAY, "+ holiday flags": CURRENT + HOL_FLAGS,
            "+ holiday-adjusted max": CURRENT + HOL_ADJ,
            "+ synchrony + burden": CURRENT + GROUPS_NEW["synchrony"] + GROUPS_NEW["burden"],
            "+ synchrony + burden + holiday": CURRENT + GROUPS_NEW["synchrony"] + GROUPS_NEW["burden"] + HOLIDAY}
    path = OUT / "cv_holiday.parquet"
    if reuse and path.exists():
        cvh = pd.read_parquet(path)
    else:
        new = {k: v for k, v in sets.items() if k not in cv_main["set"].unique()}
        cvh = pd.concat([cv_main[cv_main["set"].isin(list(sets))], cv_feature_sets(dcv, new)], ignore_index=True)
        cvh.to_parquet(path)
    boot = []
    for a_, b_ in [("+ holiday", "current"), ("+ holiday flags", "current"), ("+ holiday-adjusted max", "current"),
                   ("+ synchrony + burden + holiday", "+ synchrony + burden")]:
        for metric in ["pinball", "logloss"]:
            boot.append({"family": "holiday", "set": a_, "reference": b_, "metric": metric,
                         "boot_win": season_bootstrap(cvh, a_, b_, metric),
                         "seasons_better": int((per_season_rel(cvh, a_, b_, metric) < 1).sum())})

    # adjusting the ILINet training series: score FluSurv-NET rows only. Scale = mean paired ILINet - FluSurv-NET
    # holiday excess by position in the holiday block (0 if negative)
    h = ex[ex["block"] == "holiday"]
    pair = h[h["source"] == "ilinet"].merge(h[h["source"] == "flusurvnet"], on=["location", "season"],
                                           suffixes=("_ili", "_fsn"))
    shift = {}
    for k in range(3):
        dk = (pair[f"excess_k{k}_ili"] - pair[f"excess_k{k}_fsn"]).dropna()
        if len(dk) >= 10:
            shift[k] = max(float(dk.mean()), 0.0)
    hn["adjust_shift"] = {str(k): v for k, v in shift.items()}
    path = OUT / "cv_holiday_adjust.parquet"
    if not (reuse and path.exists()):
        base_sets = {"current": CURRENT, "+ synchrony + burden": CURRENT + GROUPS_NEW["synchrony"] + GROUPS_NEW["burden"]}
        parts = []
        r = cv_fsn_target(dcv, dict(base_sets, **{"+ holiday": CURRENT + HOLIDAY}))
        r["data"] = "ILINet as reported"
        parts.append(r)
        for lab, adj in [("ILINet holiday weeks scaled down", shift), ("ILINet holiday weeks interpolated", "interp")]:
            da, _ = build_dataset(adj)
            dca = da[da["season_week"] % 2 == 0].reset_index(drop=True)
            r = cv_fsn_target(dca, base_sets)
            r["data"] = lab
            parts.append(r)
        pd.concat(parts, ignore_index=True).to_parquet(path)
    nums["holiday"] = hn
    return boot


# ---------------------------------------------------------------------------------------------------------------
# ideas from other fields (section 11)

def literature_outputs(d, dcv, cv_main, nums, reuse=False):
    verify_reflection().to_csv(OUT / "reflection_check.csv", index=False)
    sets = {"SB": SB}
    for g, fs in GROUPS_LIT.items():
        sets[f"SB + {LIT_LABELS[g]}"] = SB + fs
    sets["SB + all"] = SB + ALL_LIT
    sets["SB + recession + records"] = SB + GROUPS_LIT["recession"] + GROUPS_LIT["records"]
    path = OUT / "cv_lit.parquet"
    if reuse and path.exists():
        cvl = pd.read_parquet(path)
        missing = {k: v for k, v in sets.items() if k not in cvl["set"].unique()}
        if missing:
            cvl = pd.concat([cvl, cv_feature_sets(dcv, missing)], ignore_index=True)
            cvl.to_parquet(path)
    else:
        base = cv_main[cv_main["set"] == "+ synchrony + burden"].assign(set="SB")
        cvl = pd.concat([base, cv_feature_sets(dcv, {k: v for k, v in sets.items() if k != "SB"})],
                        ignore_index=True)
        cvl.to_parquet(path)
    boot = []
    for st in list(sets)[1:]:
        for metric in ["pinball", "logloss"]:
            boot.append({"family": "literature", "set": st, "reference": "SB", "metric": metric,
                         "boot_win": season_bootstrap(cvl, st, "SB", metric),
                         "seasons_better": int((per_season_rel(cvl, st, "SB", metric) < 1).sum())})
    return boot


def regional_outputs(d, dcv, cv_main, reuse=False):
    """Section 12: neighbour / HHS-region synchrony and lat/lon, added to SB."""
    sets = {"SB": SB}
    for g, fs in GROUPS_REG.items():
        sets[f"SB + {REG_LABELS[g]}"] = SB + fs
    sets["SB + all regional"] = SB + ALL_REG
    path = OUT / "cv_regional.parquet"
    if reuse and path.exists():
        cvr = pd.read_parquet(path)
        missing = {k: v for k, v in sets.items() if k not in cvr["set"].unique()}
        if missing:
            cvr = pd.concat([cvr, cv_feature_sets(dcv, missing)], ignore_index=True)
            cvr.to_parquet(path)
    else:
        base = cv_main[cv_main["set"] == "+ synchrony + burden"].assign(set="SB")
        cvr = pd.concat([base, cv_feature_sets(dcv, {k: v for k, v in sets.items() if k != "SB"})], ignore_index=True)
        cvr.to_parquet(path)
    ili_state = cvr[(cvr["source"] == "ilinet") & (cvr["agg_level"] == "state")]
    boot = []
    for st in list(sets)[1:]:
        for fam, x in [("regional", cvr), ("regional, ILINet states", ili_state)]:
            for metric in ["pinball", "logloss"]:
                boot.append({"family": fam, "set": st, "reference": "SB", "metric": metric,
                             "boot_win": season_bootstrap(x, st, "SB", metric),
                             "seasons_better": int((per_season_rel(x, st, "SB", metric) < 1).sum())})
    return boot


def type_outputs(d, dcv, cv_main, reuse=False):
    """Section 13: B share / B rising and H3 share features added to SB."""
    sets = {"SB": SB, "SB + B share / B rising": SB + GROUPS_TYPE["B"], "SB + H3 share": SB + GROUPS_TYPE["H3"],
            "SB + both": SB + ALL_TYPE}
    path = OUT / "cv_types.parquet"
    if reuse and path.exists():
        cvt = pd.read_parquet(path)
    else:
        base = cv_main[cv_main["set"] == "+ synchrony + burden"].assign(set="SB")
        cvt = pd.concat([base, cv_feature_sets(dcv, {k: v for k, v in sets.items() if k != "SB"})], ignore_index=True)
        cvt.to_parquet(path)
    subsets = [("types", cvt), ("types, ILINet states", cvt[(cvt["source"] == "ilinet") & (cvt["agg_level"] == "state")]),
               ("types, FluSurv-NET", cvt[cvt["source"] == "flusurvnet"])]
    boot = []
    for st in list(sets)[1:]:
        for fam, x in subsets:
            for metric in ["pinball", "logloss"]:
                boot.append({"family": fam, "set": st, "reference": "SB", "metric": metric,
                             "boot_win": season_bootstrap(x, st, "SB", metric),
                             "seasons_better": int((per_season_rel(x, st, "SB", metric) < 1).sum())})
    return boot


def sb_group_outputs(d, dcv, cv_main, reuse=False):
    """Section 9 groups (and the holiday features) added to SB, so that every group is compared with SB."""
    sets = {"SB": SB}
    for g in ["trend", "timing", "onset", "level_hist"]:
        sets[f"SB + {GROUP_LABELS[g]}"] = SB + GROUPS_NEW[g]
    sets["SB + holiday"] = SB + HOLIDAY
    path = OUT / "cv_sb_groups.parquet"
    # SB and SB + level vs history were already fit in section 9 (same features, same folds)
    done = {"SB": "+ synchrony + burden", "SB + level vs history": "+ synchrony + burden + level"}
    if reuse and path.exists():
        cvs = pd.read_parquet(path)
    else:
        cvs = pd.concat([cv_main[cv_main["set"] == v].assign(set=k) for k, v in done.items()], ignore_index=True)
    missing = {k: v for k, v in sets.items() if k not in cvs["set"].unique()}
    if missing:
        cvs = pd.concat([cvs, cv_feature_sets(dcv, missing)], ignore_index=True)
        cvs.to_parquet(path)
    boot = []
    for st in list(sets)[1:]:
        for metric in ["pinball", "logloss"]:
            boot.append({"family": "SB groups", "set": st, "reference": "SB", "metric": metric,
                         "boot_win": season_bootstrap(cvs, st, "SB", metric),
                         "seasons_better": int((per_season_rel(cvs, st, "SB", metric) < 1).sum())})
    return boot


def gain_importance(dcv, feats, label):
    import lightgbm as lgb

    out = []
    for kind, model in [("median quantile model", lgb.LGBMRegressor(objective="quantile", alpha=0.5,
                                                                       importance_type="gain", **dict(LGB_PARAMS, n_jobs=N_WORKERS))),
                        ("z > 0 classifier", lgb.LGBMClassifier(objective="binary", importance_type="gain",
                                                                 **dict(LGB_PARAMS, n_jobs=N_WORKERS)))]:
        model.fit(dcv[feats], dcv["z"] if kind.startswith("median") else dcv["pos"])
        imp = pd.Series(model.feature_importances_, index=feats)
        out.append(pd.DataFrame({"feature set": label, "model": kind, "feature": feats,
                                 "gain_share": (imp / imp.sum()).to_numpy()}))
    return pd.concat(out, ignore_index=True)


# ---------------------------------------------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true", help="every 4th held-out season only, for testing")
    parser.add_argument("--reuse_cv", action="store_true", help="reuse eda/cv_results.parquet")
    parser.add_argument("--reuse_holiday_cv", action="store_true", help="reuse eda/cv_holiday*.parquet")
    parser.add_argument("--reuse_lit_cv", action="store_true", help="reuse eda/cv_lit.parquet")
    parser.add_argument("--reuse_regional_cv", action="store_true", help="reuse eda/cv_regional.parquet")
    parser.add_argument("--reuse_types_cv", action="store_true", help="reuse eda/cv_types.parquet")
    parser.add_argument("--reuse_sb_groups_cv", action="store_true", help="reuse eda/cv_sb_groups.parquet")
    parser.add_argument("--sb_groups_only", action="store_true",
                        help="with --reuse_cv: only add the SB-referenced group sets (cv_sb_groups.parquet), then stop")
    parser.add_argument("--workers", type=int, default=3, help="parallel single-threaded LightGBM fits")
    args = parser.parse_args()
    global N_WORKERS
    N_WORKERS = args.workers
    t_start = time.time()
    OUT.mkdir(exist_ok=True)

    d, chk = build_dataset()
    nums = {"consistency": chk, "n_rows": len(d)}
    arrays, _ = load_rows()
    nums["leakage_check_max_abs_change"] = leakage_check(d, arrays)
    nums["leakage_check_max_abs_change"].update(type_leakage_check(arrays))
    d.to_parquet(OUT / "predictor_rows.parquet")
    t_feat = time.time() - t_start
    groups = ([("current", "current", f) for f in CURRENT]
              + [("section 9", GROUP_LABELS[g], f) for g, fs in GROUPS_NEW.items() for f in fs]
              + [("section 10", "holiday", f) for f in HOLIDAY]
              + [("section 11", LIT_LABELS[g], f) for g, fs in GROUPS_LIT.items() for f in fs]
              + [("section 12", REG_LABELS[g], f) for g, fs in GROUPS_REG.items() for f in fs]
              + [("section 13", g, f) for g, fs in GROUPS_TYPE.items() for f in fs])
    pd.DataFrame(groups, columns=["section", "group", "feature"]).to_csv(OUT / "feature_groups.csv", index=False)

    # CV on every other origin week (even weeks) to keep the runtime down
    dcv = d[d["season_week"] % 2 == 0].reset_index(drop=True)
    sets = {"current": CURRENT}
    for g, fs in GROUPS_NEW.items():
        sets[f"+ {GROUP_LABELS[g]}"] = CURRENT + fs
    sets["+ all new"] = CURRENT + ALL_NEW
    for g, fs in GROUPS_NEW.items():
        sets[f"all − {GROUP_LABELS[g]}"] = CURRENT + [f for f in ALL_NEW if f not in fs]
    sets["+ synchrony + burden"] = CURRENT + GROUPS_NEW["synchrony"] + GROUPS_NEW["burden"]
    sets["+ synchrony + burden + level"] = (CURRENT + GROUPS_NEW["synchrony"] + GROUPS_NEW["burden"]
                                            + GROUPS_NEW["level_hist"])
    sets["current − hist_rel"] = [f for f in CURRENT if f != "hist_rel"]
    pd.DataFrame([{"set": k, "order": i, "n_features": len(v), "features": " ".join(v)}
                  for i, (k, v) in enumerate(sets.items())]).to_csv(OUT / "cv_sets.csv", index=False)
    cv_path = OUT / "cv_results.parquet"
    t0 = time.time()
    if args.reuse_cv and cv_path.exists():
        cv = pd.read_parquet(cv_path)
    else:
        cv = cv_feature_sets(dcv, sets, quick=args.quick)
        cv.to_parquet(cv_path)
    boot = []
    for a in [k for k in sets if k != "current"]:
        for metric in ["pinball", "logloss"]:
            boot.append({"family": "main", "set": a, "reference": "current", "metric": metric,
                         "boot_win": season_bootstrap(cv, a, "current", metric),
                         "seasons_better": int((per_season_rel(cv, a, "current", metric) < 1).sum())})
    t_cv = time.time() - t0
    if args.sb_groups_only:
        sb_group_outputs(d, dcv, cv, reuse=args.reuse_sb_groups_cv)
        return
    boot += holiday_outputs(d, dcv, cv, nums, reuse=args.reuse_holiday_cv)
    boot += literature_outputs(d, dcv, cv, nums, reuse=args.reuse_lit_cv)
    t_reg = time.time()
    boot += regional_outputs(d, dcv, cv, reuse=args.reuse_regional_cv)
    nums["runtime_regional_cv_s"] = round(time.time() - t_reg)
    t_typ = time.time()
    boot += type_outputs(d, dcv, cv, reuse=args.reuse_types_cv)
    nums["runtime_types_cv_s"] = round(time.time() - t_typ)
    boot += sb_group_outputs(d, dcv, cv, reuse=args.reuse_sb_groups_cv)
    pd.DataFrame(boot).to_csv(OUT / "cv_bootstrap.csv", index=False)
    pd.concat([gain_importance(dcv, CURRENT + ALL_NEW, "current + all new"),
               gain_importance(dcv, SB + ALL_LIT, "SB + all literature"),
               gain_importance(dcv, SB + ALL_TYPE, "SB + types")], ignore_index=True
              ).to_csv(OUT / "importance.csv", index=False)
    nums["runtime_s"] = {"features": round(t_feat), "cv": round(t_cv), "total": round(time.time() - t_start)}
    (OUT / "predictors_numbers.json").write_text(json.dumps(nums, indent=1, default=str))
    print(json.dumps({k: nums[k] for k in ["consistency", "leakage_check_max_abs_change", "runtime_s"]},
                     indent=1, default=str))


if __name__ == "__main__":
    main()
