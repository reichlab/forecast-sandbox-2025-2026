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

import model as wisar6

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

# Data sources. TARGET_DATA_PATH is the hub's own curated target-data file;
# for historical (--today_date in the past) runs we instead pull the
# point-in-time-correct snapshot from the hub's own dated archive, since the
# live file reflects revisions made after that date (see README.md's
# "Historical vs. live data" section for why this matters and why we do NOT
# use reichlab's `iddata` package for this -- the hub's own archive already
# gives point-in-time correctness with a plain CSV read, so no extra
# dependency is needed).
HUB_ROOT = Path(__file__).resolve().parents[3] / "FluSight-forecast-hub"
TARGET_DATA_PATH = HUB_ROOT / "target-data" / "target-hospital-admissions.csv"
TARGET_DATA_ARCHIVE_DIR = HUB_ROOT / "auxiliary-data" / "target-data-archive"
LOCATIONS_PATH = HUB_ROOT / "auxiliary-data" / "locations.csv"

MODEL_ABBR = "UMass-WISAR6_fourthroot_adaptive_t"
OUTPUT_ROOT = Path(__file__).resolve().parents[2] / "model-output" / MODEL_ABBR

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

# How much training history to use. The hub's target-data has a hard floor
# of 2022-02-05 across every location (verified directly); production's own
# idmodels/sarix pipeline hardcodes 2022-10-01 as its start date. We use the
# fuller available history instead (see
# reading-literature/pymc-tutorial/ar6_fourier_all_locations.py for the
# investigation behind this choice) -- override via --train-start if needed.
DEFAULT_TRAIN_START = "2022-02-05"

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

def _find_archive_snapshot(as_of: datetime.date) -> Path:
    """Return the most recent target-data-archive snapshot with a date <=
    as_of, for point-in-time-correct historical backtesting. Falls back to
    the live target-data file if as_of is on/after today (no archive needed)
    or if no earlier archive snapshot exists (e.g. the very first forecast
    dates in the hub's history)."""
    if not TARGET_DATA_ARCHIVE_DIR.is_dir():
        return TARGET_DATA_PATH
    candidates = sorted(TARGET_DATA_ARCHIVE_DIR.glob("target-hospital-admissions_*.csv"))
    usable = [
        p for p in candidates
        if p.stem.split("_")[-1] <= as_of.isoformat()
    ]
    if not usable:
        return TARGET_DATA_PATH
    return usable[-1]


def load_target_data(as_of: datetime.date) -> pd.DataFrame:
    path = _find_archive_snapshot(as_of)
    raw = pd.read_csv(path)
    raw["date"] = pd.to_datetime(raw["date"])
    raw["location"] = raw["location"].astype(str)
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
            f"with train_start={train_start}. Check the archive snapshot used."
        )

    # Per-location interior-gap interpolation, matching production's
    # idmodels._interpolate_by_location (see model.py's module docstring
    # and reading-literature/pymc-tutorial/ar6_fourier_all_locations.py for
    # why: some locations have genuine NHSN reporting gaps, all in the
    # off-season trough, and this matches how production handles them).
    train_df = train_df.copy()
    train_df["value"] = train_df.groupby("location")["value"].transform(
        lambda s: s.interpolate(method="linear", limit_area="inside")
    )

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
            idata, THETA_POOLING, li, z_hist, dates_hist, h_max=MAX_HORIZON,
            n_param_draws=N_PARAM_DRAWS, n_innov_per_param=N_INNOV_PER_PARAM,
            seed=SEED_FORECAST + li, innovation_dist=INNOVATION_DIST, t_dof=T_DOF,
        )  # shape (n_sims, MAX_HORIZON)

        for h in range(MAX_HORIZON):
            draws_z = sims[:, h]
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
