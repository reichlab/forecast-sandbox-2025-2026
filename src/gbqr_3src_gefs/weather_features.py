"""
Build weekly weather anomaly features from the GEFS analysis / forecast tables
written by gefs_data.py, and merge them onto surveillance rows.

For a surveillance row with week-ending date W (a Saturday), the forecast
"available" for that row is the GEFS 00Z run on the following Wednesday
(W + 4 days), which is when a FluSight submission using week W as the last
observed week would be made.

Features (anomalies relative to a 2000-2019 smoothed daily climatology):
- observed:  wx_{temp,pw}_obs0 (week W), wx_{temp,pw}_obs1 (week W-1)
- forecast:  wx_{temp,pw}_fc{1,2,3} for weeks ending W+7, W+14, W+21
    * fc1 combines Sun-Tue analysis with Wed-Sat forecast
    * forecasts are bias-corrected by (location, week-ahead, variable) using
      inits before BIAS_CUTOFF
    * rows without an archived forecast (before Oct 2020, or summer) get a
      synthetic forecast: observed anomaly + the full error vector (all
      locations, weeks ahead and variables) from one randomly chosen archived
      init, so that synthetic errors keep the real lead-dependence and spatial
      correlation of GEFS errors
- oracle:    wx_{temp,pw}_true{1,2,3}, the observed anomalies for those weeks
             (an upper bound on what perfect weather foresight could add)
"""

import zlib
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from gefs_data import add_aggregate_locations, load_cell_weights, week_end

DATA_DIR = Path(__file__).parent / "data"
CLIM_YEARS = (2000, 2019)
BIAS_CUTOFF = pd.Timestamp("2023-07-01")
VARS = ["temp", "pw"]  # pw = precipitable water, the humidity variable (2m RH is missing in GEFS analysis pre-2020)
N_HARMONICS = 3


@lru_cache(maxsize=1)
def load_daily_obs() -> pd.DataFrame:
    d = pd.read_csv(DATA_DIR / "gefs_analysis_daily.csv", dtype={"location": str}, parse_dates=["day"])
    # days straddling a pull block boundary appear twice; average the partial-day means
    d = d.groupby(["location", "day"], as_index=False)[VARS].mean()
    d = add_aggregate_locations(d, VARS, ["day"], load_cell_weights())
    return d.sort_values(["location", "day"]).reset_index(drop=True)


def _harmonic_design(doy: np.ndarray) -> np.ndarray:
    w = 2 * np.pi * doy / 365.25
    cols = [np.ones_like(w)]
    for h in range(1, N_HARMONICS + 1):
        cols += [np.sin(h * w), np.cos(h * w)]
    return np.column_stack(cols)


@lru_cache(maxsize=1)
def daily_anomalies() -> pd.DataFrame:
    """Daily obs with smoothed-climatology columns and anomalies."""
    d = load_daily_obs().copy()
    doy = d["day"].dt.dayofyear.to_numpy().astype(float)
    X = _harmonic_design(doy)
    in_clim = d["day"].dt.year.between(*CLIM_YEARS).to_numpy()
    for v in VARS:
        d[f"{v}_clim"] = np.nan
        y = d[v].to_numpy()
        clim = np.full(len(d), np.nan)
        for loc, idx in d.groupby("location").indices.items():
            ok = idx[in_clim[idx] & ~np.isnan(y[idx])]
            beta, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None)
            clim[idx] = X[idx] @ beta
        d[f"{v}_clim"] = clim
        d[f"{v}_anom"] = d[v] - d[f"{v}_clim"]
    return d


@lru_cache(maxsize=1)
def weekly_obs() -> pd.DataFrame:
    """Weekly (Sun-Sat) observed values, climatology and anomalies."""
    d = daily_anomalies()
    d = d.assign(wk_end_date=week_end(d["day"]))
    cols = [c for c in d.columns if c.startswith(tuple(VARS))]
    w = d.groupby(["location", "wk_end_date"])[cols].mean()
    w["n_days"] = d.groupby(["location", "wk_end_date"]).size()
    return w[w["n_days"] == 7].drop(columns="n_days").reset_index()


@lru_cache(maxsize=1)
def forecast_anomalies() -> pd.DataFrame:
    """
    One row per (init_date, location, k) with bias-corrected forecast anomalies
    {v}_fc and matching observed anomalies {v}_obs for the target week.
    """
    f = pd.read_csv(DATA_DIR / "gefs_forecast_weekly.csv", dtype={"location": str},
                    parse_dates=["init_date", "wk_end_date"])
    f = f.pivot_table(index=["init_date", "location", "wk_end_date", "k"], columns="var", values="mean").reset_index()
    f = add_aggregate_locations(f, VARS, ["init_date", "wk_end_date", "k"], load_cell_weights())
    daily = daily_anomalies()
    clim_cols = [f"{v}_clim" for v in VARS]
    # weekly climatology for each target week
    dd = daily.assign(wk_end_date=week_end(daily["day"]))
    wclim = dd.groupby(["location", "wk_end_date"])[clim_cols].mean().reset_index()
    f = f.merge(wclim, on=["location", "wk_end_date"], how="left")
    # week 1: forecast covers Wed-Sat only; blend in analysis for Sun-Tue
    sun_tue = dd[dd["day"].dt.weekday.isin([6, 0, 1])].groupby(["location", "wk_end_date"])[VARS].mean()
    sun_tue = sun_tue.rename(columns=lambda c: f"{c}_suntue").reset_index()
    f = f.merge(sun_tue, on=["location", "wk_end_date"], how="left")
    k1 = f["k"] == 1
    for v in VARS:
        f.loc[k1, v] = (3 * f.loc[k1, f"{v}_suntue"] + 4 * f.loc[k1, v]) / 7
        f[f"{v}_fc"] = f[v] - f[f"{v}_clim"]
    obs = weekly_obs()[["location", "wk_end_date"] + [f"{v}_anom" for v in VARS]]
    f = f.merge(obs.rename(columns={f"{v}_anom": f"{v}_obs" for v in VARS}), on=["location", "wk_end_date"], how="left")
    # bias correction by (location, k, var), estimated on early inits only
    early = f[f["init_date"] < BIAS_CUTOFF]
    for v in VARS:
        bias = (early[f"{v}_fc"] - early[f"{v}_obs"]).groupby([early["location"], early["k"]]).mean().rename(f"{v}_bias")
        f = f.merge(bias.reset_index(), on=["location", "k"], how="left")
        f[f"{v}_fc"] = f[f"{v}_fc"] - f[f"{v}_bias"].fillna(0)
    return f[["init_date", "location", "wk_end_date", "k"] + [f"{v}_fc" for v in VARS] + [f"{v}_obs" for v in VARS]]


@lru_cache(maxsize=1)
def error_pool() -> pd.DataFrame:
    """Forecast-minus-observed anomaly errors, indexed by init, from inits before BIAS_CUTOFF."""
    f = forecast_anomalies()
    f = f[f["init_date"] < BIAS_CUTOFF].dropna()
    e = f[["init_date", "location", "k"]].copy()
    for v in VARS:
        e[f"{v}_err"] = f[f"{v}_fc"] - f[f"{v}_obs"]
    return e


def weather_feature_table(locations, wk_end_dates) -> pd.DataFrame:
    """All weather features for each (location, W) pair."""
    obs = weekly_obs().set_index(["location", "wk_end_date"])
    keys = pd.MultiIndex.from_product([sorted(set(locations)), sorted(set(pd.to_datetime(wk_end_dates)))],
                                      names=["location", "wk_end_date"])
    out = pd.DataFrame(index=keys).reset_index()
    W = out["wk_end_date"]

    def obs_anom(v, shift_weeks):
        idx = pd.MultiIndex.from_arrays([out["location"], W + pd.Timedelta(weeks=shift_weeks)])
        return obs[f"{v}_anom"].reindex(idx).to_numpy()

    for v in VARS:
        out[f"wx_{v}_obs0"] = obs_anom(v, 0)
        out[f"wx_{v}_obs1"] = obs_anom(v, -1)
        for k in (1, 2, 3):
            out[f"wx_{v}_true{k}"] = obs_anom(v, k)

    # real archived forecasts where available
    fc = forecast_anomalies()
    fc = fc.assign(W=fc["init_date"] - pd.Timedelta(days=4))
    fcw = fc.pivot_table(index=["location", "W"], columns="k", values=[f"{v}_fc" for v in VARS])
    idx = pd.MultiIndex.from_arrays([out["location"], W])
    has_fc = np.zeros(len(out), dtype=bool)
    for v in VARS:
        for k in (1, 2, 3):
            vals = fcw[(f"{v}_fc", k)].reindex(idx).to_numpy() if (f"{v}_fc", k) in fcw.columns else np.full(len(out), np.nan)
            out[f"wx_{v}_fc{k}"] = vals
            has_fc |= ~np.isnan(vals)

    # synthetic forecasts elsewhere: true anomaly + one archived init's error vector per W
    pool = error_pool()
    pool_inits = np.array(sorted(pool["init_date"].unique()))
    perr = pool.set_index(["init_date", "location", "k"])
    for w_date, rows in out[~has_fc].groupby("wk_end_date").groups.items():
        rng = np.random.default_rng(zlib.crc32(str(w_date.date()).encode()))
        init = pool_inits[rng.integers(len(pool_inits))]
        for v in VARS:
            for k in (1, 2, 3):
                ei = pd.MultiIndex.from_arrays([np.repeat(init, len(rows)), out.loc[rows, "location"], np.repeat(k, len(rows))])
                err = perr[f"{v}_err"].reindex(ei).to_numpy()
                out.loc[rows, f"wx_{v}_fc{k}"] = out.loc[rows, f"wx_{v}_true{k}"].to_numpy() + err
    out["wx_fc_is_real"] = has_fc
    return out


FEATURE_SETS = {
    "none": [],
    "obs": [f"wx_{v}_obs{i}" for v in VARS for i in (0, 1)],
    "obs_fc": [f"wx_{v}_obs{i}" for v in VARS for i in (0, 1)] + [f"wx_{v}_fc{k}" for v in VARS for k in (1, 2, 3)],
    "obs_oracle": [f"wx_{v}_obs{i}" for v in VARS for i in (0, 1)] + [f"wx_{v}_true{k}" for v in VARS for k in (1, 2, 3)],
}


def add_weather_features(df: pd.DataFrame, variant: str, ref_date=None) -> tuple:
    """Merge the variant's weather features onto surveillance rows."""
    feats = FEATURE_SETS[variant]
    if not feats:
        return df, []
    tab = weather_feature_table(df["location"].unique(), df["wk_end_date"].unique())
    if ref_date is not None and variant == "obs_fc":
        # live-use guard: the test rows must have a real forecast, never a synthetic one
        last = tab[tab["wk_end_date"] == pd.Timestamp(df["wk_end_date"].max())]
        n_syn = (~last["wx_fc_is_real"]).sum()
        if n_syn:
            print(f"Warning: {n_syn} test-week rows have no archived GEFS forecast (synthetic used)")
    df = df.merge(tab[["location", "wk_end_date"] + feats], on=["location", "wk_end_date"], how="left")
    return df, feats
