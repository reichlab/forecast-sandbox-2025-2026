"""
Search for better predictors of the relative size target z (companion to relative_size_eda.py).

Builds the season-replay rows exactly as the models do (via relative_size_eda.load_rows, which calls
idmodels.peak.series.build_replay_rows; NHSN excluded, Puerto Rico / Virgin Islands ILINet dropped), adds candidate
features computed only from data available at season week t, and evaluates them by
  - within-week-bin Spearman correlations with z and with 1{z > 0}, and
  - leave-one-season-out LightGBM quantile regression of z and classification of z > 0, comparing feature sets.
Outputs (figures, markdown tables, predictors_numbers.json, cv_results.parquet) go to analysis/peak-models/eda/.

Usage (from the repository root):
    OMP_NUM_THREADS=4 DYLD_FALLBACK_LIBRARY_PATH=<venv>/lib/python3.12/site-packages/sklearn/.dylibs \
        python analysis/peak-models/relative_size_predictors.py [--quick]
"""
import argparse
import json
import sys
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from idmodels.peak.base import season_week_to_date  # noqa: E402
from idmodels.peak.series import build_replay_rows, build_season_arrays, running_max, season_peaks  # noqa: E402
from relative_size_eda import (INK, INK2, MIN_OBS, OUT, REPLAY_START, W0, W1, Z_YLIM, load_rows,  # noqa: E402
                               md_table, save)

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


def build_dataset(hol_adjust=None):
    arrays, rows = load_rows_adjusted(hol_adjust) if hol_adjust else load_rows()
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
    cur_only = (GROUPS_NEW["trend"] + GROUPS_NEW["onset"] + GROUPS_NEW["synchrony"] + HOLIDAY
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
    allf = ALL_NEW + HOLIDAY + ALL_LIT
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


def cv_feature_sets(d: pd.DataFrame, sets: dict, quick=False, n_jobs=3) -> pd.DataFrame:
    """Leave-one-season-out: each season (all sources and locations) is held out in turn."""
    from joblib import Parallel, delayed

    seasons = sorted(d["season"].unique())
    if quick:
        seasons = seasons[::4]
    jobs = [(name, s) for name in sets for s in seasons]
    t0 = time.time()
    res = Parallel(n_jobs=n_jobs)(
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


def summarize_cv(cv: pd.DataFrame, ref="current") -> tuple[pd.DataFrame, pd.DataFrame]:
    cv = cv.assign(bin=week_bin(cv["season_week"]))
    by = cv.groupby(["set", "bin"])[["pinball", "logloss"]].mean().unstack("bin")
    tot = cv.groupby("set")[["pinball", "logloss"]].mean()
    for m in ["pinball", "logloss"]:
        by[(m, "all")] = tot[m]
    rel_pin = by["pinball"].div(by["pinball"].loc[ref], axis=1)
    d_ll = by["logloss"].sub(by["logloss"].loc[ref], axis=1)
    order = [bin_label(*b) for b in BINS] + ["all"]
    return rel_pin[order], d_ll[order], by


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


def per_season_rel(cv, a, b):
    x = cv[(cv["season_week"] >= 12) & (cv["season_week"] <= 31)]
    g = x.groupby(["set", "season"])["pinball"].mean().unstack("set")
    return (g[a] / g[b])


# ---------------------------------------------------------------------------------------------------------------
# figures

def fig_corr_heatmap(d, groups=None, name="fig11_predictor_correlations", first=None):
    groups = GROUPS_NEW if groups is None else groups
    first = CURRENT[1:-1] if first is None else first  # drop season_week and src_code
    feats = first + [f for g in groups.values() for f in g]
    d = d.assign(bin=week_bin(d["season_week"], CORR_BINS))
    labels = [bin_label(*b) for b in CORR_BINS]
    res = {}
    for target in ["z", "pos"]:
        mat = pd.DataFrame(index=feats, columns=labels, dtype=float)
        for lab in labels:
            x = d[d["bin"] == lab]
            for f in feats:
                mat.loc[f, lab] = x[[f, target]].corr(method="spearman").iloc[0, 1]
        res[target] = mat
    fig, axes = plt.subplots(1, 2, figsize=(10, 0.28 * len(feats) + 1.5), sharey=True)
    cmap = matplotlib.colormaps["RdBu_r"]
    for ax, (target, title) in zip(axes, [("z", "Spearman with z"), ("pos", "Spearman with 1{z > 0}")]):
        mat = res[target]
        im = ax.imshow(mat.to_numpy(dtype=float), cmap=cmap, vmin=-1, vmax=1, aspect="auto")
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat.iat[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                            color="white" if abs(v) > 0.6 else INK)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels([f"wk {lab}" for lab in labels], fontsize=8.5)
        ax.xaxis.tick_top()
        ax.set_title(title, loc="left", pad=22)
        ax.grid(False)
        # separators between feature groups
        edges = np.cumsum([len(first)] + [len(v) for v in groups.values()])[:-1]
        for e in edges:
            ax.axhline(e - 0.5, color="white", lw=3)
    axes[0].set_yticks(range(len(feats)))
    axes[0].set_yticklabels(feats, fontsize=8)
    # group labels on the right
    starts = np.concatenate([[0], np.cumsum([len(first)] + [len(v) for v in groups.values()])])
    labels_g = {**GROUP_LABELS, **LIT_LABELS}
    names = ["current" if first == CURRENT[1:-1] else "reference"] + [labels_g[g] for g in groups]
    for s0, s1, nm in zip(starts[:-1], starts[1:], names):
        axes[1].text(len(labels) - 0.35, (s0 + s1 - 1) / 2, nm, ha="left", va="center", fontsize=8.5, color=INK2)
    cb = fig.colorbar(im, ax=axes, orientation="horizontal", fraction=0.02, pad=0.02, aspect=50)
    cb.set_label("Spearman correlation within the week bin (same scale in both panels)")
    save(fig, name)
    return res


def fig_cv(rel_pin, d_ll, sets_order, name="fig12_cv_feature_sets", ref_label="current features"):
    labels = [c for c in rel_pin.columns]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True)
    y = np.arange(len(sets_order))
    blues = matplotlib.colormaps["Blues"]
    cols = {lab: blues(0.35 + 0.6 * i / (len(labels) - 2)) for i, lab in enumerate(labels[:-1])}
    cols["all"] = INK
    for ax, (tab, xlab, ref) in zip(axes, [(rel_pin, f"pinball loss relative to {ref_label}", 1.0),
                                          (d_ll, f"change in log loss for z > 0 vs {ref_label}", 0.0)]):
        ax.axvline(ref, color=INK2, lw=1)
        for k, lab in enumerate(labels):
            off = (k - (len(labels) - 1) / 2) * 0.1
            ax.scatter(tab.loc[sets_order, lab], y + off, s=22 if lab != "all" else 46, color=cols[lab],
                       marker="o" if lab != "all" else "D", label=f"weeks {lab}" if lab != "all" else "all weeks",
                       zorder=3, edgecolor="white", linewidth=0.6)
        ax.set_xlabel(xlab)
        ax.grid(axis="y", visible=False)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(sets_order)
    axes[0].invert_yaxis()
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=len(lab), fontsize=8.5, bbox_to_anchor=(0.55, -0.04))
    axes[0].set_title("Quantile forecasts of z (lower is better)", loc="left")
    axes[1].set_title("Probability that the peak is still ahead (lower is better)", loc="left")
    fig.tight_layout()
    save(fig, name)


def fig_importance(imp):
    top = imp.head(25)[::-1]
    fig, ax = plt.subplots(figsize=(7, 6.5))
    cols = ["#b7b6b0" if f in CURRENT else "#2a78d6" for f in top.index]
    ax.barh(range(len(top)), top.values, color=cols, height=0.7)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top.index, fontsize=8.5)
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    ax.set_xlabel("share of total gain (median quantile model, all features)")
    ax.grid(axis="y", visible=False)
    ax.text(0.98, 0.04, "grey: current features\nblue: new candidates", transform=ax.transAxes, ha="right",
            fontsize=8.5, color=INK2)
    fig.tight_layout()
    save(fig, "fig13_importance")


def fig_binned(d, feats, name, title_extra=None):
    bins = [(12, 16), (17, 21), (22, 26)]
    shades = [matplotlib.colormaps["Blues"](v) for v in (0.45, 0.7, 0.95)]
    ncol = 3
    nrow = int(np.ceil(len(feats) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(11.5, 3.3 * nrow), sharey=True)
    axes = np.atleast_2d(axes)
    for ax, f in zip(axes.flat, feats):
        for (lo, hi), c in zip(bins, shades):
            x = d[(d["season_week"] >= lo) & (d["season_week"] <= hi)].dropna(subset=[f])
            if x[f].nunique() < 5:
                continue
            q = pd.qcut(x[f].rank(method="first"), 10, labels=False)
            g = x.groupby(q)
            mid, med = g[f].median(), g["z"].median()
            ax.fill_between(mid, g["z"].quantile(0.25), g["z"].quantile(0.75), color=c, alpha=0.18, lw=0)
            ax.plot(mid, med, color=c, lw=2, marker="o", ms=3.5, label=f"weeks {lo}–{hi}")
        ax.set_title(f + (" (current feature)" if f in CURRENT else ""), loc="left")
        ax.set_ylim(*Z_YLIM)
    for ax in axes.flat[len(feats):]:
        ax.set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel("eventual z (median, 25–75%)")
    axes.flat[0].legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    save(fig, name)


# ---------------------------------------------------------------------------------------------------------------

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


def cluster_ci(v: pd.Series, clusters: pd.Series, n_boot=2000, seed=0):
    """Mean and 95% bootstrap interval, resampling seasons."""
    g = pd.DataFrame({"v": v.to_numpy(), "c": clusters.to_numpy()}).groupby("c")["v"].agg(["sum", "count"])
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(g), size=(n_boot, len(g)))
    bs = g["sum"].to_numpy()[idx].sum(axis=1) / g["count"].to_numpy()[idx].sum(axis=1)
    return float(v.mean()), float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))


def fig_holiday_excess(ex):
    from relative_size_eda import COLORS, GROUPS

    blocks = ["placebo −6 wk", "holiday", "placebo +6 wk"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4), gridspec_kw={"width_ratios": [1.3, 1]})
    ax = axes[0]
    res = {}
    for gi, g in enumerate(GROUPS):
        for bi, b in enumerate(blocks):
            v = ex[(ex["group"] == g) & (ex["block"] == b)]
            mu, lo, hi = cluster_ci(v["excess"], v["season"])
            res[(g, b)] = (mu, lo, hi, len(v))
            x = bi + (gi - 1) * 0.22
            ax.plot([x, x], [lo, hi], color=COLORS[g], lw=2)
            ax.plot([x], [mu], marker="o", ms=7, color=COLORS[g], mec="white",
                    label=g if bi == 0 else None)
    ax.axhline(0, color=INK2, lw=1)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["6 weeks earlier\n(placebo)", "holiday weeks\n(Dec 22 – Jan 7)", "6 weeks later\n(placebo)"])
    ax.set_ylabel("mean log excess vs 2 weeks before + 2 after")
    ax.set_title("Holiday excess by source (mean, 95% season-bootstrap CI)", loc="left")
    ax.legend(fontsize=8.5, loc="upper left")
    ax.grid(axis="x", visible=False)
    # paired: ILINet minus FluSurv-NET, same location and season
    ax = axes[1]
    h = ex[ex["block"] == "holiday"]
    pair = h[h["source"] == "ilinet"].merge(h[h["source"] == "flusurvnet"], on=["location", "season"],
                                           suffixes=("_ili", "_fsn"))
    diffs = {}
    for k in range(3):
        dk = (pair[f"excess_k{k}_ili"] - pair[f"excess_k{k}_fsn"]).dropna()
        if len(dk) < 10:
            continue
        mu, lo, hi = cluster_ci(dk, pair.loc[dk.index, "season"])
        diffs[k] = (mu, lo, hi, len(dk))
        ax.plot([k, k], [lo, hi], color=INK, lw=2)
        ax.plot([k], [mu], marker="D", ms=7, color=INK, mec="white")
    dall = pair["excess_ili"] - pair["excess_fsn"]
    mu, lo, hi = cluster_ci(dall, pair["season"])
    diffs["block"] = (mu, lo, hi, len(dall))
    ax.plot([3, 3], [lo, hi], color=COLORS["ILINet states"], lw=2.5)
    ax.plot([3], [mu], marker="D", ms=8, color=COLORS["ILINet states"], mec="white")
    ax.axhline(0, color=INK2, lw=1)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["1st holiday\nweek", "2nd", "3rd", "block\nmean"])
    ax.set_ylabel("ILINet excess − FluSurv-NET excess")
    ax.set_title(f"Same location & season ({len(pair)} pairs)", loc="left")
    ax.grid(axis="x", visible=False)
    fig.tight_layout()
    save(fig, "fig15_holiday_excess")
    return res, diffs


def fig_holiday_max(d):
    """P(eventual z = 0) given the lag since the running max was set, split by whether that max week was a holiday
    week; rows at season weeks 20-30."""
    from relative_size_eda import COLORS, GROUPS

    x = d[(d["season_week"] >= 20) & (d["season_week"] <= 30) & (d["wks_since_max"] >= 1) &
          (d["wks_since_max"] <= 5)]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), sharey=True)
    res = {}
    for ax, g in zip(axes, GROUPS):
        v = x[x["group"] == g]
        for flag, ls, lab in [(1.0, "-", "max set in a holiday week"), (0.0, "--", "max set in another week")]:
            u = v[v["max_in_holiday"] == flag]
            p = u.groupby("wks_since_max")["at_zero"].agg(["mean", "size"])
            ax.plot(p.index, p["mean"], color=COLORS[g], lw=2, ls=ls, marker="o", ms=4, label=lab)
            for w, r in p.iterrows():
                res[(g, int(flag), int(w))] = (round(float(r["mean"]), 3), int(r["size"]))
        ax.set_title(g, loc="left")
        ax.set_xticks(range(1, 6))
        ax.set_xlabel("weeks since the running max was set")
        ax.set_ylim(0, 1)
        ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    axes[0].set_ylabel("share whose max was the final peak\n(eventual z = 0)")
    axes[0].legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    save(fig, "fig16_holiday_max_final")
    return res


def post_holiday_table(d):
    """Rows 1-4 weeks after the last holiday week: how often the running max is a holiday week and what happens."""
    from relative_size_eda import GROUPS

    x = d[(d["weeks_since_holiday_end"] >= 1) & (d["weeks_since_holiday_end"] <= 4)]
    rows = []
    for g in GROUPS:
        for j in range(1, 5):
            v = x[(x["group"] == g) & (x["weeks_since_holiday_end"] == j)]
            hmax = v[v["max_in_holiday"] == 1]
            other = v[v["max_in_holiday"] == 0]
            rows.append({"source": g, "weeks after holidays": j, "rows": len(v),
                         "max in holiday week": f"{len(hmax) / max(len(v), 1):.0%}",
                         "P(z=0) holiday max": f"{hmax['at_zero'].mean():.0%}" if len(hmax) else "–",
                         "P(z=0) other max": f"{other['at_zero'].mean():.0%}" if len(other) else "–",
                         "median z holiday max": round(float(hmax["z"].median()), 2) if len(hmax) else np.nan,
                         "median z other max": round(float(other["z"].median()), 2) if len(other) else np.nan,
                         "q90 z holiday max": round(float(hmax["z"].quantile(0.9)), 2) if len(hmax) else np.nan,
                         "q90 z other max": round(float(other["z"].quantile(0.9)), 2) if len(other) else np.nan})
    tab = pd.DataFrame(rows)
    md_table(tab, OUT / "table_holiday_post.md")
    return tab


def hol_bin(w_after: pd.Series, hol_now: pd.Series) -> pd.Series:
    out = pd.Series("other weeks", index=w_after.index, dtype="object")
    out[hol_now == 1] = "holiday weeks"
    out[(w_after >= 1) & (w_after <= 3)] = "1–3 wk after"
    out[(w_after >= 4) & (w_after <= 6)] = "4–6 wk after"
    return out


def summarize_holiday_cv(cv, dmeta, sets, ref="current"):
    x = cv.merge(dmeta, on=KEYS, how="left")
    x["hbin"] = hol_bin(x["weeks_since_holiday_end"], x["hol_now"])
    x["src"] = np.where(x["source"] == "ilinet", "ILINet", "FluSurv-NET")
    out = []
    for src in ["all", "ILINet", "FluSurv-NET"]:
        v = x if src == "all" else x[x["src"] == src]
        for bin_name, sel in [("all weeks", slice(None)), ("weeks 17–31", (v["season_week"] >= 17) & (v["season_week"] <= 31))] + \
                [(b, v["hbin"] == b) for b in ["holiday weeks", "1–3 wk after", "4–6 wk after"]]:
            u = v.loc[sel] if not isinstance(sel, slice) else v
            g = u.groupby("set")[["pinball", "logloss"]].mean()
            for s in sets:
                if s == ref or s not in g.index:
                    continue
                out.append({"rows": src, "weeks": bin_name, "feature set": s,
                            "pinball ratio": g.loc[s, "pinball"] / g.loc[ref, "pinball"],
                            "Δ log loss": g.loc[s, "logloss"] - g.loc[ref, "logloss"],
                            "n": int((u["set"] == ref).sum())})
    return pd.DataFrame(out)


def cv_fsn_target(dtrain: pd.DataFrame, sets: dict, n_jobs=3) -> pd.DataFrame:
    """Leave-one-season-out, but score only the FluSurv-NET rows of the held-out season (targets unaffected by any
    ILINet adjustment)."""
    from joblib import Parallel, delayed

    seasons = sorted(dtrain.loc[dtrain["source"] == "flusurvnet", "season"].unique())
    jobs = [(name, s) for name in sets for s in seasons]
    te_mask = {s: (dtrain["season"] == s) & (dtrain["source"] == "flusurvnet") for s in seasons}
    res = Parallel(n_jobs=n_jobs)(
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


def holiday_section(d, dcv, cv_main, nums, reuse=False):
    """Section 10: holiday weeks."""
    from relative_size_eda import GROUPS

    hn = {}
    # calendar
    seasons = sorted(d["season"].unique())
    cal = []
    for s in seasons:
        hw = holiday_weeks(s)
        cal.append({"season": s, "season weeks": ", ".join(str(w) for w in hw),
                    "week-ending dates": ", ".join(f"{season_week_to_date(s, w):%b %d}" for w in hw),
                    "MMWR weeks": ", ".join(str(mmwr_week(season_week_to_date(s, w))) for w in hw)})
    cal = pd.DataFrame(cal)
    md_table(cal, OUT / "table_holiday_calendar.md")
    hn["calendar"] = cal.set_index("season").to_dict(orient="index")

    arrays, _ = load_rows()
    ex = holiday_excess(arrays)
    res, diffs = fig_holiday_excess(ex)
    hn["excess"] = {f"{g} | {b}": [round(v, 3) for v in r[:3]] + [r[3]] for (g, b), r in res.items()}
    hn["excess_paired_ili_minus_fsn"] = {str(k): [round(v, 3) for v in r[:3]] + [r[3]] for k, r in diffs.items()}
    h = ex[ex["block"] == "holiday"]
    hn["excess_by_position"] = {g: {k: round(float(h.loc[h["group"] == g, f"excess_k{k}"].mean()), 3)
                                    for k in range(3)} for g in GROUPS}
    # how often the peak itself falls in a holiday week, per week, vs the 3 weeks either side
    ser = d.drop_duplicates(["source", "agg_level", "location", "season"])
    pk_share = {}
    for g in GROUPS:
        v = ser[ser["group"] == g]
        hol_rate, near_rate = [], []
        for _, r in v.iterrows():
            hw = holiday_weeks(r["season"])
            hol_rate.append(r["peak_week"] in hw)
            near = list(range(hw[0] - 3, hw[0])) + list(range(hw[-1] + 1, hw[-1] + 4))
            near_rate.append(r["peak_week"] in near)
        nh = np.mean([len(holiday_weeks(s)) for s in v["season"]])
        pk_share[g] = {"share_peaks_in_holiday": round(float(np.mean(hol_rate)), 3),
                       "per_week_holiday": round(float(np.mean(hol_rate) / nh), 3),
                       "per_week_adjacent6": round(float(np.mean(near_rate) / 6), 3), "n_series": len(v)}
    hn["peak_in_holiday"] = pk_share
    hn["holiday_max_final"] = {f"{g} | {'holiday' if f else 'other'} | wsm {w}": v
                               for (g, f, w), v in fig_holiday_max(d).items()}
    post = post_holiday_table(d)
    hn["post_holiday_table"] = post.to_dict(orient="records")

    # correlations of the holiday features with z in weeks 20-30, by source group
    x = d[(d["season_week"] >= 20) & (d["season_week"] <= 30)]
    hn["spearman_20_30"] = {g: {f: round(float(x.loc[x["group"] == g, [f, "z"]].corr(method="spearman").iloc[0, 1]), 2)
                                for f in HOLIDAY} for g in GROUPS}

    # CV: holiday features
    sets = {"current": CURRENT, "+ holiday": CURRENT + HOLIDAY, "+ holiday flags": CURRENT + HOL_FLAGS,
            "+ holiday-adjusted max": CURRENT + HOL_ADJ,
            "+ synchrony + burden": CURRENT + GROUPS_NEW["synchrony"] + GROUPS_NEW["burden"],
            "+ synchrony + burden + holiday": CURRENT + GROUPS_NEW["synchrony"] + GROUPS_NEW["burden"] + HOLIDAY}
    path = OUT / "cv_holiday.parquet"
    t0 = time.time()
    if reuse and path.exists():
        cvh = pd.read_parquet(path)
    else:
        new = {k: v for k, v in sets.items() if k not in cv_main["set"].unique()}
        cvh = pd.concat([cv_main[cv_main["set"].isin(list(sets))], cv_feature_sets(dcv, new)], ignore_index=True)
        cvh.to_parquet(path)
    hn["runtime_cv_s"] = round(time.time() - t0)
    meta = dcv[KEYS + ["weeks_since_holiday_end", "hol_now"]]
    tab = summarize_holiday_cv(cvh, meta, list(sets))
    hn["cv"] = tab.round(4).to_dict(orient="records")
    wide = tab.pivot_table(index=["rows", "feature set"], columns="weeks", values="pinball ratio")
    wide = wide[["all weeks", "weeks 17–31", "holiday weeks", "1–3 wk after", "4–6 wk after"]].reset_index()
    md_table(wide, OUT / "table_holiday_cv_pinball.md", floatfmt=3)
    wide = tab.pivot_table(index=["rows", "feature set"], columns="weeks", values="Δ log loss")
    wide = wide[["all weeks", "weeks 17–31", "holiday weeks", "1–3 wk after", "4–6 wk after"]].reset_index()
    md_table(wide, OUT / "table_holiday_cv_logloss.md", floatfmt=4)
    fig_holiday_cv(tab)
    comp = {}
    for a_, b_ in [("+ holiday", "current"), ("+ holiday flags", "current"), ("+ holiday-adjusted max", "current"),
                   ("+ synchrony + burden + holiday", "+ synchrony + burden")]:
        for metric in ["pinball", "logloss"]:
            comp[f"{a_} vs {b_} | {metric}"] = {"boot_win": round(season_bootstrap(cvh, a_, b_, metric), 3),
                                                "seasons_better": int((cvh[cvh["season_week"].between(12, 31)]
                                                                       .groupby(["set", "season"])[metric].mean()
                                                                       .unstack("set").eval(f"`{a_}` < `{b_}`")).sum())}
    hn["cv_season_comparisons"] = comp

    # adjusting the ILINet training series: score FluSurv-NET rows only
    shift = {k: max(diffs[k][0], 0.0) for k in range(3) if k in diffs}
    hn["adjust_shift"] = {str(k): round(v, 3) for k, v in shift.items()}
    path = OUT / "cv_holiday_adjust.parquet"
    t0 = time.time()
    if reuse and path.exists():
        cva = pd.read_parquet(path)
    else:
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
        cva = pd.concat(parts, ignore_index=True)
        cva.to_parquet(path)
    hn["runtime_adjust_s"] = round(time.time() - t0)
    cva = cva.merge(meta, on=KEYS, how="left")
    cva["hbin"] = hol_bin(cva["weeks_since_holiday_end"], cva["hol_now"])
    ref = cva[(cva["data"] == "ILINet as reported") & (cva["set"] == "current")]
    rows = []
    for (dat, st), v in cva.groupby(["data", "set"]):
        r = {"training data": dat, "feature set": st}
        for b, sel_v, sel_r in [("all weeks", slice(None), slice(None))] + \
                [(hb, v["hbin"] == hb, ref["hbin"] == hb) for hb in ["holiday weeks", "1–3 wk after", "4–6 wk after"]]:
            vv = v if isinstance(sel_v, slice) else v[sel_v]
            rr = ref if isinstance(sel_r, slice) else ref[sel_r]
            r[b] = vv["pinball"].mean() / rr["pinball"].mean()
        r["Δ log loss (all)"] = v["logloss"].mean() - ref["logloss"].mean()
        rows.append(r)
    adj_tab = pd.DataFrame(rows)
    md_table(adj_tab, OUT / "table_holiday_adjust.md", floatfmt=3)
    hn["adjust"] = adj_tab.round(4).to_dict(orient="records")
    nums["holiday"] = hn


def fig_holiday_cv(tab):
    bins = ["all weeks", "holiday weeks", "1–3 wk after", "4–6 wk after"]
    sets = ["+ holiday", "+ holiday flags", "+ holiday-adjusted max", "+ synchrony + burden",
            "+ synchrony + burden + holiday"]
    shades = [INK] + [matplotlib.colormaps["Blues"](v) for v in (0.45, 0.7, 0.95)]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True, sharex=True)
    for ax, src in zip(axes, ["ILINet", "FluSurv-NET"]):
        v = tab[tab["rows"] == src]
        for k, (b, c) in enumerate(zip(bins, shades)):
            u = v[v["weeks"] == b].set_index("feature set").reindex(sets)
            ax.scatter(u["pinball ratio"], np.arange(len(sets)) + (k - 1.5) * 0.14, color=c, s=30 if k else 44,
                       marker="D" if k == 0 else "o", edgecolor="white", linewidth=0.6, label=b, zorder=3)
        ax.axvline(1, color=INK2, lw=1)
        ax.set_title(f"{src} rows", loc="left")
        ax.set_xlabel("pinball loss relative to current features")
        ax.grid(axis="y", visible=False)
    axes[0].set_yticks(range(len(sets)))
    axes[0].set_yticklabels(sets)
    axes[0].invert_yaxis()
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=4, fontsize=8.5, bbox_to_anchor=(0.55, -0.05))
    fig.tight_layout()
    save(fig, "fig17_holiday_cv")


def lit_section(d, dcv, cv_main, nums, reuse=False):
    """Section 11: features suggested by analogous problems (reflection principle, chain ladder, ...)."""
    from relative_size_eda import GROUPS

    ln = {}
    ver = verify_reflection()
    md_table(ver.round(3), OUT / "table_reflection_check.md", floatfmt=3)
    ln["reflection_check"] = ver.round(4).to_dict(orient="records")
    corr = fig_corr_heatmap(d, GROUPS_LIT, "fig18_lit_correlations", first=["hist_rel", "cum_vs_hist_total", "rel_max"])
    ln["spearman_z"] = corr["z"].round(2).to_dict()
    ln["spearman_pos"] = corr["pos"].round(2).to_dict()
    ln["coverage"] = {f: round(float(d[f].notna().mean()), 3) for f in ALL_LIT}
    ln["bf_check"] = {"min": float(np.nanmin(d["z_bf"])), "max": float(np.nanmax(d["z_bf"])),
                      "spearman_with_z": round(float(d[["z_bf", "z"]].corr(method="spearman").iloc[0, 1]), 3),
                      "spearman_cl_with_z": round(float(d[["cl_z", "z"]].corr(method="spearman").iloc[0, 1]), 3)}

    sets = {"SB": SB}
    for g, fs in GROUPS_LIT.items():
        sets[f"SB + {LIT_LABELS[g]}"] = SB + fs
    sets["SB + all"] = SB + ALL_LIT
    sets["SB + recession + records"] = SB + GROUPS_LIT["recession"] + GROUPS_LIT["records"]
    path = OUT / "cv_lit.parquet"
    t0 = time.time()
    if reuse and path.exists():
        cvl = pd.read_parquet(path)
        missing = {k: v for k, v in sets.items() if k not in cvl["set"].unique()}
        if missing:
            cvl = pd.concat([cvl, cv_feature_sets(dcv, missing, n_jobs=3)], ignore_index=True)
            cvl.to_parquet(path)
    else:
        base = cv_main[cv_main["set"] == "+ synchrony + burden"].assign(set="SB")
        cvl = pd.concat([base, cv_feature_sets(dcv, {k: v for k, v in sets.items() if k != "SB"}, n_jobs=3)],
                        ignore_index=True)
        cvl.to_parquet(path)
    ln["runtime_cv_s"] = round(time.time() - t0)
    rel_pin, d_ll, raw = summarize_cv(cvl, ref="SB")
    order = list(sets)
    fig_cv(rel_pin, d_ll, order, name="fig19_lit_cv", ref_label="SB")
    for tabname, tab, fmt in [("table_lit_cv_pinball.md", rel_pin, 3), ("table_lit_cv_logloss.md", d_ll, 4)]:
        t = tab.loc[order].copy()
        t.insert(0, "feature set", t.index)
        md_table(t.reset_index(drop=True), OUT / tabname, floatfmt=fmt)
    ln["rel_pinball"] = rel_pin.round(4).to_dict(orient="index")
    ln["delta_logloss"] = d_ll.round(4).to_dict(orient="index")
    ln["abs_SB"] = {f"{k[0]} {k[1]}": v for k, v in raw.loc["SB"].round(4).to_dict().items()}
    # by source of the held-out rows
    x = cvl.merge(d[KEYS + ["group"]], on=KEYS, how="left")
    src_rows = []
    for g in GROUPS:
        rp, dl, _ = summarize_cv(x[x["group"] == g].drop(columns="group"), ref="SB")
        for st in order[1:]:
            src_rows.append({"rows": g, "feature set": st, "pinball all": rp.loc[st, "all"],
                             "pinball 17–21": rp.loc[st, "17–21"], "pinball 22–26": rp.loc[st, "22–26"],
                             "Δ log loss all": dl.loc[st, "all"]})
    # rows that have a cross-source partner
    haspart = d.loc[d["xs_rel_max"].notna(), KEYS]
    rp, dl, _ = summarize_cv(x.merge(haspart, on=KEYS).drop(columns="group"), ref="SB")
    for st in order[1:]:
        src_rows.append({"rows": "rows with a partner source", "feature set": st, "pinball all": rp.loc[st, "all"],
                         "pinball 17–21": rp.loc[st, "17–21"], "pinball 22–26": rp.loc[st, "22–26"],
                         "Δ log loss all": dl.loc[st, "all"]})
    src_tab = pd.DataFrame(src_rows)
    md_table(src_tab, OUT / "table_lit_cv_by_source.md", floatfmt=3)
    ln["by_source"] = src_tab.round(4).to_dict(orient="records")
    ln["n_partner_rows_cv"] = int(len(x.merge(haspart, on=KEYS)) // len(order))
    boot = {}
    for st in order[1:]:
        boot[st] = {"pinball_win": round(season_bootstrap(cvl, st, "SB"), 3),
                    "logloss_win": round(season_bootstrap(cvl, st, "SB", metric="logloss"), 3),
                    "seasons_better_pinball": int((per_season_rel(cvl, st, "SB") < 1).sum())}
    ln["bootstrap_vs_SB"] = boot
    bt = pd.DataFrame(boot).T.reset_index().rename(columns={"index": "feature set"})
    md_table(bt, OUT / "table_lit_bootstrap.md", floatfmt=3)

    import lightgbm as lgb

    m = lgb.LGBMRegressor(objective="quantile", alpha=0.5, importance_type="gain", **dict(LGB_PARAMS, n_jobs=3))
    m.fit(dcv[SB + ALL_LIT], dcv["z"])
    imp = pd.Series(m.feature_importances_, index=SB + ALL_LIT).sort_values(ascending=False)
    imp = imp / imp.sum()
    ln["importance_SB_all"] = imp.head(20).round(3).to_dict()
    top = [f for f in imp.index if f in ALL_LIT][:6]
    ln["top_lit"] = top
    fig_binned(d, top, "fig20_lit_top_features")
    nums["literature"] = ln


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true", help="every 4th held-out season only, for testing")
    parser.add_argument("--reuse_cv", action="store_true", help="reuse eda/cv_results.parquet")
    parser.add_argument("--reuse_holiday_cv", action="store_true", help="reuse eda/cv_holiday*.parquet")
    parser.add_argument("--reuse_lit_cv", action="store_true", help="reuse eda/cv_lit.parquet")
    args = parser.parse_args()
    t_start = time.time()
    OUT.mkdir(exist_ok=True)

    d, chk = build_dataset()
    nums = {"consistency": chk, "n_rows": len(d)}
    t_feat = time.time() - t_start
    arrays, _ = load_rows()
    nums["leakage_check_max_abs_change"] = leakage_check(d, arrays)
    d.to_parquet(OUT / "predictor_rows.parquet")

    # decomposition z = (log peak - hist mean log peak) - hist_rel: the first term does not depend on t
    ser = d.dropna(subset=["hist_rel"]).drop_duplicates(["source", "agg_level", "location", "season"])
    anom = ser["log_peak"] - (ser["lm"] - ser["hist_rel"])
    nums["size_anomaly_sd"] = {g: round(float(v), 2) for g, v in anom.groupby(ser["group"]).std().items()}
    nums["size_anomaly_sd_all"] = round(float(anom.std()), 2)
    nums["size_anomaly_q10_q90"] = [round(float(anom.quantile(0.1)), 2), round(float(anom.quantile(0.9)), 2)]
    nums["share_rows_hist_rel_missing"] = round(float(d["hist_rel"].isna().mean()), 3)
    zb = d[(d["season_week"] >= 17) & (d["season_week"] <= 26)]
    nums["sd_z_weeks17_26"] = round(float(zb["z"].std()), 2)
    corr = fig_corr_heatmap(d)
    nums["spearman_z"] = corr["z"].round(2).to_dict()
    nums["spearman_pos"] = corr["pos"].round(2).to_dict()

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
    cv_path = OUT / "cv_results.parquet"
    t0 = time.time()
    if args.reuse_cv and cv_path.exists():
        cv = pd.read_parquet(cv_path)
    else:
        cv = cv_feature_sets(dcv, sets, quick=args.quick)
        cv.to_parquet(cv_path)
    t_cv = time.time() - t0
    rel_pin, d_ll, raw = summarize_cv(cv)
    order = [k for k in sets if k in rel_pin.index] if args.reuse_cv else list(sets)
    fig_cv(rel_pin, d_ll, [k for k in order if k != "current − hist_rel"])  # hist_rel ablation is off scale
    tab = rel_pin.loc[order].copy()
    tab.insert(0, "feature set", tab.index)
    md_table(tab.reset_index(drop=True), OUT / "table_cv_pinball.md", floatfmt=3)
    tab = d_ll.loc[order].copy()
    tab.insert(0, "feature set", tab.index)
    md_table(tab.reset_index(drop=True), OUT / "table_cv_logloss.md", floatfmt=3)
    nums["rel_pinball"] = rel_pin.round(3).to_dict(orient="index")
    nums["delta_logloss"] = d_ll.round(4).to_dict(orient="index")
    nums["abs_pinball_current"] = raw["pinball"].loc["current"].round(4).to_dict()
    nums["abs_logloss_current"] = raw["logloss"].loc["current"].round(4).to_dict()
    for a in ["+ all new", "+ synchrony + burden", "+ synchrony + burden + level", "current − hist_rel",
              "all − timing vs history"] + [f"+ {GROUP_LABELS[g]}" for g in GROUPS_NEW]:
        nums.setdefault("boot_win_vs_current", {})[a] = round(season_bootstrap(cv, a, "current"), 3)
        nums.setdefault("seasons_better_than_current", {})[a] = int((per_season_rel(cv, a, "current") < 1).sum())
    nums["n_seasons"] = int(cv["season"].nunique())

    # gain importance for the all-features median model fitted to all rows
    import lightgbm as lgb

    m = lgb.LGBMRegressor(objective="quantile", alpha=0.5, importance_type="gain", **LGB_PARAMS)
    m.fit(dcv[CURRENT + ALL_NEW], dcv["z"])
    imp = pd.Series(m.feature_importances_, index=CURRENT + ALL_NEW).sort_values(ascending=False)
    imp = imp / imp.sum()
    c = lgb.LGBMClassifier(objective="binary", importance_type="gain", **LGB_PARAMS)
    c.fit(dcv[CURRENT + ALL_NEW], dcv["pos"])
    imp_c = pd.Series(c.feature_importances_, index=CURRENT + ALL_NEW).sort_values(ascending=False)
    imp_c = imp_c / imp_c.sum()
    fig_importance(imp)
    nums["importance_median_model"] = imp.head(20).round(3).to_dict()
    nums["importance_classifier"] = imp_c.head(20).round(3).to_dict()
    top_new = [f for f in imp.index if f in ALL_NEW][:5]
    nums["top_new_features"] = top_new
    fig_binned(d, ["hist_rel"] + top_new, "fig14_top_new_features")
    holiday_section(d, dcv, cv, nums, reuse=args.reuse_holiday_cv)
    lit_section(d, dcv, cv, nums, reuse=args.reuse_lit_cv)
    nums["runtime_s"] = {"features": round(t_feat), "cv": round(t_cv), "total": round(time.time() - t_start)}
    (OUT / "predictors_numbers.json").write_text(json.dumps(nums, indent=1, default=str))
    print(json.dumps({k: nums[k] for k in ["rel_pinball", "delta_logloss", "boot_win_vs_current",
                                           "seasons_better_than_current", "top_new_features", "runtime_s",
                                           "leakage_check_max_abs_change"]}, indent=1, default=str))
    print(json.dumps({k: v for k, v in nums["literature"].items() if k not in ("spearman_z", "spearman_pos")},
                     indent=1, default=str))


if __name__ == "__main__":
    main()
