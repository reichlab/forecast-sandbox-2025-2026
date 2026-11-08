"""WISAR6_fourthroot_adaptive_t CLI entrypoint: generate a FluSight-format
quantile forecast for one reference date, for all FluSight locations.

Usage (matches the hub's other models' CLI convention):

    python main.py --today_date=2024-01-06

See model.py for the model itself (pooled AR(6) + Fourier seasonality,
fit via either a standard Gaussian likelihood or a WIS-based Gibbs
posterior loss -- see model.py's module docstring), and README.md for
environment setup.
"""

import datetime
import os
from pathlib import Path

import click
import numpy as np
import pandas as pd
from dateutil import relativedelta
from iddata.enums import Disease
from iddata.sources.nhsn import NHSNDataSource

import model as wisar6

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

# Data source. NHSN admissions come from reichlab's `iddata` package, the
# same source and as_of convention (as_of = reference date) the other UMass
# components use through idmodels. iddata picks the newest NHSN release
# dated on or before as_of; releases are named by their (weekly, Wednesday)
# release date, so this is point-in-time correct for historical runs and
# has no skipped weeks. (The hub's target-data-archive, used previously,
# skips some weeks, which left the model a week short of data on 9 of the
# 80 backtest dates.) Before 2024-11-15 iddata serves the HHS-era series.
# Only locations.csv is still read from the hub clone.
HUB_ROOT = Path(__file__).resolve().parents[4] / "FluSight-forecast-hub"
LOCATIONS_PATH = HUB_ROOT / "auxiliary-data" / "locations.csv"

MODEL_ABBR = "UMass-WISAR6_fourthroot_adaptive_t"
OUTPUT_ROOT = Path(__file__).resolve().parents[1] / "intermediate-output" / "model-output" / MODEL_ABBR

# Model configuration. Pooling settings are left as an explicit, visible
# choice here (not buried in model.py) since which setting to use
# operationally is still being decided based on
# reading-literature/pymc-tutorial/AR6_FOURIER_POOLED_RESULTS.md and its
# all-53-location follow-up -- update these three lines once that's settled.
OBJECTIVE = "wis"          # "wis" or "likelihood"
THETA_POOLING = "shared"   # "shared" or "none"
FOURIER_POOLING = "shared"  # "shared" or "none"

MAX_HORIZON = 4  # horizons 0, 1, 2, 3
N_PARAM_DRAWS = 200
N_INNOV_PER_PARAM = 5
MCMC_DRAWS = 2000
MCMC_TUNE = 2000
# Override with the WISAR6_MCMC_CHAINS env var (e.g. 4 on a 4-core SLURM
# allocation, one chain per core). Changing the chain count changes the
# posterior sample, so outputs won't be bit-identical across chain counts.
MCMC_CHAINS = int(os.environ.get("WISAR6_MCMC_CHAINS", "2"))
# Defaults to the safe, validated cores=1 (see model.fit_model's docstring
# for the macOS EOFError this avoids). Override with the WISAR6_MCMC_CORES
# env var on Linux/Unity, where that crash has not been observed -- e.g.
# `WISAR6_MCMC_CORES=2 uv run python main.py --today_date=...` runs the 2
# chains in parallel instead of sequentially. submit-unity-parallel.sh sets
# this to 2.
MCMC_CORES = int(os.environ.get("WISAR6_MCMC_CORES", "1"))
TARGET_ACCEPT = 0.99

# Learned-prior-scale hyperparameters (see model.fit_model's docstring for
# the full rationale): each of sigma, phi, and the Fourier coefficients gets
# its own HalfCauchy(*_PRIOR_SCALE) hyperprior on its SD, mirroring
# production SARIX's own defaults exactly (all three default to 1.0 there
# too). If divergences reappear for the US series under this looser
# specification, tighten THETA_SD_PRIOR_SCALE first (see fit_model's
# docstring) rather than reintroducing a fixed phi scale.
SIGMA_PRIOR_SCALE = 1.0
THETA_SD_PRIOR_SCALE = 1.0
FOURIER_BETA_SD_PRIOR_SCALE = 1.0

# Innovation distribution under test in this sibling (see model.py's module
# docstring and model-output/UMass-WISAR6_fourthroot/PERFORMANCE_NOTES.md's
# "Attempt 3" section for the full motivation): a heavier-tailed Student-t
# in place of Gaussian, tested in isolation (no differencing) after the
# combined differencing+Student-t test (Attempt 2) could not separate the
# two changes' individual effects. T_DOF is NOT sampled under
# objective="wis" -- see model.fit_model's docstring for why.
INNOVATION_DIST = "studentt"  # "studentt" or "gaussian"
T_DOF = 5.0

# How much training history to use. iddata's NHSN series is complete across
# every location from 2022-08-06 (the HHS-era series used for as_of before
# 2024-11-15 starts 2022-09-03); this matches idmodels' WISARModel.
DEFAULT_TRAIN_START = "2022-08-06"

# Fixed seeds. All randomness in this model (MCMC sampling and Monte Carlo
# forecast simulation) is seeded from these two values plus a per-location
# offset, so re-running this script for the same --today_date reproduces
# bit-identical output. See reading-literature/pymc-tutorial/
# wis_vs_likelihood_ar.py's Section 4.1 for why Python's built-in hash() is
# never used for this (it is randomized per-process and silently breaks
# reproducibility).
SEED_FIT = 9001
SEED_FORECAST = 9002

# Optional progress log: a lightweight sanity check that a long-running batch
# job is still making progress. Set WISAR6_PROGRESS_LOG to a file path to
# enable it; WISAR6_PROGRESS_EVERY sets how many MCMC iterations (per chain)
# pass between lines. Logging never affects results (no randomness involved).
PROGRESS_LOG = os.environ.get("WISAR6_PROGRESS_LOG")
PROGRESS_EVERY = int(os.environ.get("WISAR6_PROGRESS_EVERY", "250"))


def log_progress(msg: str):
    if PROGRESS_LOG is None:
        return
    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(PROGRESS_LOG, "a") as f:
        f.write(f"[{stamp}] {msg}\n")


def make_mcmc_progress_callback():
    """PyMC sample() callback logging every PROGRESS_EVERY iterations per
    chain, with the running divergence count. Returns None (no callback) when
    progress logging is disabled."""
    if PROGRESS_LOG is None:
        return None
    divergences = {}

    def callback(trace, draw):
        stats = draw.stats[0] if isinstance(draw.stats, list) else draw.stats
        if not draw.tuning and stats.get("diverging", False):
            divergences[draw.chain] = divergences.get(draw.chain, 0) + 1
        i = draw.draw_idx + 1
        if i % PROGRESS_EVERY == 0 or draw.is_last:
            phase = "tune" if draw.tuning else "draw"
            done = i if draw.tuning else i - MCMC_TUNE
            total = MCMC_TUNE if draw.tuning else MCMC_DRAWS
            log_progress(f"chain {draw.chain}: {phase} {done}/{total} "
                         f"(divergences so far: {divergences.get(draw.chain, 0)})")

    return callback


# --------------------------------------------------------------------------
# Data loading
# --------------------------------------------------------------------------

def load_target_data(as_of: datetime.date) -> pd.DataFrame:
    """NHSN weekly flu admissions (counts) as available on as_of, via iddata.

    Returns hub-style columns (location, date, value). Values can be NaN for
    weeks a location didn't report (MA, MN, and WV in the 2024 off-season);
    generate_forecast interpolates those interior gaps. Raises if the newest
    week isn't the one ending as_of - 7 days, e.g. a live run before this
    week's NHSN release."""
    raw = NHSNDataSource(disease=Disease.FLU, rates=False).load(as_of=as_of)
    raw = raw.dropna(subset=["location"]).rename(columns={"wk_end_date": "date", "inc": "value"})
    raw = raw[["location", "date", "value"]].astype({"location": str})
    raw["date"] = pd.to_datetime(raw["date"])
    latest = raw["date"].max().date()
    expected = as_of - datetime.timedelta(days=7)
    click.echo(f"Target data: iddata NHSN release as of {as_of} (latest week {latest})")
    if latest < expected:
        raise RuntimeError(f"iddata NHSN data as of {as_of} ends {latest}, expected {expected}; "
                           "has this week's NHSN release been posted?")
    return raw


def load_locations() -> tuple[dict, dict]:
    loc_info = pd.read_csv(LOCATIONS_PATH)
    locations = dict(zip(loc_info["location"].astype(str), loc_info["location_name"]))
    population = dict(zip(loc_info["location"].astype(str), loc_info["population"]))
    return locations, population


def flu_season(d: pd.Timestamp) -> int:
    """FluSight seasons run Oct 1 (year Y) - Sep 30 (year Y+1), labeled by Y."""
    return d.year if d.month >= 10 else d.year - 1


# --------------------------------------------------------------------------
# Forecast generation for one reference date
# --------------------------------------------------------------------------

def generate_forecast(reference_date: datetime.date, train_start: str = DEFAULT_TRAIN_START):
    locations, population = load_locations()
    locs_ordered = list(locations.keys())
    n_loc = len(locs_ordered)

    raw = load_target_data(reference_date)
    data = raw[raw["location"].isin(locations)].sort_values(["location", "date"]).reset_index(drop=True)

    # Training data: everything strictly before the reference date, from
    # train_start onward. (The reference date's own week is the horizon-(-1)
    # target in FluSight's usual convention; this model does not forecast
    # horizon -1, matching AR6_fourierP_thetaP's own horizon range of 0..3.)
    train_start_ts = pd.Timestamp(train_start)
    ref_ts = pd.Timestamp(reference_date)
    train_df = data[(data["date"] >= train_start_ts) & (data["date"] < ref_ts)].reset_index(drop=True)

    if train_df.empty:
        raise ValueError(
            f"No training data found before reference_date={reference_date} "
            f"with train_start={train_start}. Check the iddata NHSN release used."
        )

    # Horizons are labeled relative to the reference date, but simulation
    # steps forward from the last training week, normally reference_date - 7
    # days. load_target_data enforces that, so n_gap_steps is 0; this stays
    # as a safeguard if the newest week is ever missing.
    last_train_date = train_df["date"].max()
    n_gap_steps = (ref_ts - last_train_date).days // 7 - 1
    if n_gap_steps > 0:
        click.echo(f"Note: latest data available before {reference_date} ends "
                   f"{last_train_date.date()}; simulating {n_gap_steps} extra "
                   "leading step(s) to reach horizon 0.")

    # Per-location interior-gap interpolation, matching production's
    # idmodels._interpolate_by_location (see model.py's module docstring
    # and reading-literature/pymc-tutorial/ar6_fourier_all_locations.py for
    # why: some locations have genuine NHSN reporting gaps, all in the
    # off-season trough, and this matches how production handles them).
    train_df = train_df.copy()
    train_df["value"] = train_df.groupby("location")["value"].transform(
        lambda s: s.interpolate(method="linear", limit_area="inside")
    )

    if train_df["value"].isna().any():
        missing = train_df.loc[train_df["value"].isna(), ["location", "date"]]
        raise ValueError("Training data has missing values that interpolation can't fill "
                         f"(leading or trailing gaps):\n{missing.to_string(index=False)}")

    # Transform each location's series.
    transform_params = {}
    z_train_by_loc, dates_by_loc = {}, {}
    for loc in locs_ordered:
        sub = train_df[train_df["location"] == loc]
        pop = population[loc]
        r_train = wisar6.fourth_root_transform(sub["value"].values, pop)
        q95, mean_scaled = wisar6.fit_scale_params(r_train)
        transform_params[loc] = (q95, mean_scaled)
        z_train_by_loc[loc] = wisar6.apply_center_scale(r_train, q95, mean_scaled)
        dates_by_loc[loc] = sub["date"].values

    week_counts = {loc: len(z) for loc, z in z_train_by_loc.items()}
    if len(set(week_counts.values())) != 1:
        raise ValueError(
            f"Locations have inconsistent training week counts: {week_counts}. "
            "Every location must share the same date range for the pooled design "
            "matrix construction below."
        )

    y_train, loc_idx_train, F_target_train, F_lags_train, z_lags_train = wisar6.build_full_design(
        z_train_by_loc, dates_by_loc, locs_ordered
    )

    log_progress(f"data prepared: {n_loc} locations x {len(y_train) // n_loc} training weeks")

    loss_scale = 1.0
    if OBJECTIVE == "wis":
        log_progress("calibrating WIS loss scale")
        loss_scale = wisar6.calibrate_loss_scale_pooled(
            y_train, loc_idx_train, F_target_train, F_lags_train, z_lags_train
        )

    log_progress(f"starting MCMC: {MCMC_CHAINS} chains x ({MCMC_TUNE} tune + {MCMC_DRAWS} draws), "
                 f"cores={MCMC_CORES} (first lines may lag while PyTensor compiles)")
    _, idata = wisar6.fit_model(
        y_train, loc_idx_train, F_target_train, F_lags_train, z_lags_train, n_loc,
        objective=OBJECTIVE, theta_pooling=THETA_POOLING, fourier_pooling=FOURIER_POOLING,
        loss_scale=loss_scale, draws=MCMC_DRAWS, tune=MCMC_TUNE, chains=MCMC_CHAINS,
        cores=MCMC_CORES, seed=SEED_FIT, target_accept=TARGET_ACCEPT,
        sigma_prior_scale=SIGMA_PRIOR_SCALE, theta_sd_prior_scale=THETA_SD_PRIOR_SCALE,
        fourier_beta_sd_prior_scale=FOURIER_BETA_SD_PRIOR_SCALE,
        innovation_dist=INNOVATION_DIST, t_dof=T_DOF,
        callback=make_mcmc_progress_callback(),
    )
    log_progress("MCMC finished; simulating forecasts")

    # Forecast each location at horizons 0..MAX_HORIZON-1.
    rows = []
    for li, loc in enumerate(locs_ordered):
        if li % 10 == 0:
            log_progress(f"forecast simulation: location {li + 1}/{n_loc}")
        z_hist = z_train_by_loc[loc]
        dates_hist = dates_by_loc[loc]
        pop = population[loc]
        q95, mean_scaled = transform_params[loc]

        sims = wisar6.simulate_forecast_draws(
            idata, THETA_POOLING, li, z_hist, dates_hist, h_max=n_gap_steps + MAX_HORIZON,
            n_param_draws=N_PARAM_DRAWS, n_innov_per_param=N_INNOV_PER_PARAM,
            seed=SEED_FORECAST + li, innovation_dist=INNOVATION_DIST, t_dof=T_DOF,
        )  # shape (n_sims, n_gap_steps + MAX_HORIZON)

        for h in range(MAX_HORIZON):
            draws_z = sims[:, n_gap_steps + h]
            draws_r = wisar6.invert_center_scale(draws_z, q95, mean_scaled)
            draws_natural = wisar6.invert_fourth_root_transform(draws_r, pop)
            target_end_date = reference_date + relativedelta.relativedelta(weeks=h)
            for lv in wisar6.FLUSIGHT_QUANTILE_LEVELS:
                value = float(np.quantile(draws_natural, lv))
                rows.append({
                    "location": loc,
                    "horizon": h,
                    "output_type_id": lv,
                    "value": value,
                    "target_end_date": target_end_date.isoformat(),
                    "reference_date": reference_date.isoformat(),
                    "output_type": "quantile",
                    "target": "wk inc flu hosp",
                })

    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

@click.command()
@click.option(
    "--today_date",
    type=str,
    required=False,
    help="Date to use as effective model run date (YYYY-MM-DD).",
)
@click.option(
    "--train-start",
    type=str,
    default=DEFAULT_TRAIN_START,
    show_default=True,
    help="Earliest date of training data to use.",
)
def main(today_date: str | None = None, train_start: str = DEFAULT_TRAIN_START):
    """Generate a WISAR6_fourthroot_adaptive_t flu hospitalization forecast for one reference date."""
    try:
        today = datetime.date.fromisoformat(today_date)
    except (TypeError, ValueError):
        today = datetime.date.today()
    reference_date = today + relativedelta.relativedelta(weekday=5)  # nearest Saturday

    click.echo(f"Generating WISAR6_fourthroot_adaptive_t forecast for reference_date={reference_date} "
               f"(objective={OBJECTIVE}, theta_pooling={THETA_POOLING}, "
               f"fourier_pooling={FOURIER_POOLING})")

    forecast_df = generate_forecast(reference_date, train_start=train_start)

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    out_path = OUTPUT_ROOT / f"{reference_date.isoformat()}-{MODEL_ABBR}.csv"
    forecast_df.to_csv(out_path, index=False)
    click.echo(f"Wrote {len(forecast_df)} rows to {out_path}")
    log_progress(f"done: wrote {len(forecast_df)} rows to {out_path}")


if __name__ == "__main__":
    main()
