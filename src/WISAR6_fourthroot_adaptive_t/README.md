# WISAR6_fourthroot_adaptive_t: WIS-loss AR(6) with Fourier seasonality (fourth-root transform, adaptive priors, Student-t innovations)

An AR(6) model with Fourier seasonality, fit jointly across all FluSight
locations, using a **WIS-based Gibbs posterior loss** (Bissiri, Holmes, &
Walker 2016) as the estimation objective instead of a likelihood. See
`model-metadata/UMass-WISAR6_fourthroot_adaptive_t.yml` for the full
methods description.

Model components:

- **Data transform:** rate per 100,000, then `(rate + 0.01)^0.25`, then
  center/scale by each location's in-season 95th percentile. Preliminary
  held-out comparisons suggested that fitting on the log scale
  (`log(rate + 1)`, aligned with log-scale WIS) did not perform as well,
  so this model uses the fourth-root transform.
- **Seasonality:** K=2 Fourier harmonics with a 365.25-day period. AR and
  Fourier coefficients are shared across locations; the innovation scale
  `sigma` is location-specific.
- **Innovations:** Student-t with 5 degrees of freedom, scaled by `sigma`.
  The degrees of freedom are fixed rather than sampled, because the
  WIS-loss objective needs gradients through the Student-t quantile
  function, and those are not available with respect to the degrees of
  freedom.
- **Priors:** prior scales are learned from the data:
  `sigma ~ HalfCauchy(1)`, `phi ~ Normal(0, theta_sd)` with
  `theta_sd ~ HalfCauchy(1)`, and Fourier coefficients
  `~ Normal(0, fourier_beta_sd)` with `fourier_beta_sd ~ HalfCauchy(1)`.
- **Forecasts:** Monte Carlo simulation of the AR(6) recursion forward
  from posterior parameter draws, summarized to the FluSight 23-quantile
  grid.

This is a from-scratch PyMC/PyTensor implementation. It does **not** use
the `idmodels`/`sarix` packages that other AR(6) models on this hub use.

## Files

- `model.py`: the model itself: data transform, WIS loss, pooled AR(6) +
  Fourier design matrix, MCMC fitting, and Monte Carlo forecast
  simulation. Read the module docstring first.
- `main.py`: CLI entrypoint (`--today_date`). It connects `model.py` to
  the hub's data sources and writes hubverse-format output. Model settings
  (objective, pooling, priors, innovation distribution, MCMC and forecast
  sample sizes) are the constants near the top of this file.
- `pyproject.toml` / `uv.lock`: pinned dependencies, managed with
  [`uv`](https://docs.astral.sh/uv/).
- `submit-unity-season-{2023-24,2024-25,2025-26}.sh`: SLURM array jobs,
  one task per reference date in that season (29, 23, and 28 dates: the
  valid reference dates in `hub-config/tasks.json`). These produced the
  forecasts in `model-output/UMass-WISAR6_fourthroot_adaptive_t/`.
- `submit-unity-all-seasons.sh`: submits the three season arrays in
  sequence, each starting after the previous one finishes.
- `wisar6-task.sh`: the shared per-task body the season scripts source.
- `run-all-forecasts-2023-2025.sh`: runs all 80 dates sequentially, for
  local use.
- `submit-unity-parallel.sh`: an older single 80-task array job (2 chains,
  no progress log). The season scripts replace it.

## To run locally

```bash
cd src/WISAR6_fourthroot_adaptive_t
uv sync
uv run python main.py --today_date=2024-01-06
```

This writes a forecast CSV to
`../../model-output/UMass-WISAR6_fourthroot_adaptive_t/`. Locally the
model runs 2 MCMC chains sequentially by default (see the note on
`cores=1` below).

## Running on Unity (UMass's HPC cluster)

From a Unity login node:

```bash
cd src/WISAR6_fourthroot_adaptive_t
uv sync                          # build .venv once (must be built on Linux, not copied from a laptop)
mkdir -p logs
./submit-unity-all-seasons.sh    # or: sbatch submit-unity-season-2024-25.sh
```

Each array task runs one reference date with 4 MCMC chains on 4 CPUs
and 4G memory, with at most 10 tasks running at once per season. Each task:

- skips a date whose output CSV already exists, so resubmitting a season
  reruns only missing dates (set `WISAR6_OVERWRITE=1` to refit);
- writes a progress log (`logs/<job>-<id>_<task>.progress`) with
  MCMC draw counts, divergences, and a CPU/memory heartbeat;
- uses its own PyTensor compile cache and user cache directory on
  node-local disk, so simultaneous tasks don't collide.

Two environment variables control the chains: `WISAR6_MCMC_CHAINS`
(default 2) and `WISAR6_MCMC_CORES` (default 1). `wisar6-task.sh` sets
both to the task's CPU count. Changing the chain count changes the
posterior sample, so outputs are only reproducible for the same chain
count.

Time limits are 2 h (2023-24), 3 h (2024-25), and 4 h (2025-26). Run
time grows with the length of the training data. In the Unity run that
produced the committed forecasts, the slowest date per season took 42 min,
56 min, and 2 h 13 min, and no MCMC divergences occurred.

## Historical vs. live data

`main.py` reads target data from the FluSight Forecast Hub repo, assumed
to be checked out as a sibling directory (`../../../FluSight-forecast-hub`
relative to this file; adjust `HUB_ROOT` in `main.py` if your checkout is
laid out differently).

For a **historical** `--today_date`, `main.py` selects the most recent
snapshot in the hub's `auxiliary-data/target-data-archive/` that
predates the requested date, rather than the live
`target-data/target-hospital-admissions.csv`. The live file includes
revisions made after that date, which would leak future information into
a backtest.

## A note on `cores=1`

`model.fit_model()` defaults to `cores=1` (chains sampled sequentially in
one process). At this model's scale (53 locations fit jointly), PyMC's
multiprocessing crashed with `EOFError` on the macOS development machine.
That failure is specific to macOS's process-start behavior and has not
been seen on Linux, so the Unity scripts run chains in parallel via
`WISAR6_MCMC_CORES`.
