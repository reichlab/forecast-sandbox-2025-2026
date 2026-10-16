# GBQR 3-source spatial x WISAR6 ensemble, + gbqr_4src_nssp (production)

Production pipeline for two published models:

- **`UMass-gbqr_3src_spatial_wisar`** - an ensemble (target: `wk inc flu hosp`) of:
  - **GBQR 3-source spatial** (`gbqr/`) - state-level forecasts, with
    directional wave features. No US-level forecast (the wave features are
    state-level only).
  - **GBQR 3-source** (`gbqr_us/`) - the same GBQR model without spatial
    features, used only for its US-level forecast.
  - **WISAR6_fourthroot_adaptive_t** (`wisar/`) - a from-scratch PyMC AR(6) +
    Fourier model with a WIS-loss Gibbs posterior objective, covering both
    state and US levels.
- **`UMass-gbqr_3src_spatial_wisar_nssp`** - the above ensemble row-bound with:
  - **gbqr_4src_nssp** (`nssp/`) - an independent GBQR model (NSSP as main
    source) predicting a different target, `wk inc flu prop ed visits`. Not
    blended with anything; its forecast is simply added as extra rows. Only
    valid from the 2025/26 season onward.

`3_ensemble.R` builds both: first the `wk inc flu hosp` ensemble (equally
weighted quantile average of gbqr + wisar for states, gbqr_us + wisar for the
US), published on its own to `../../model-output/UMass-gbqr_3src_spatial_wisar/`;
then that same ensemble row-bound with `nssp/`'s output (when available for
the reference date) and published to
`../../model-output/UMass-gbqr_3src_spatial_wisar_nssp/`, which is the file
actually submitted to FluSight.

This directory is self-contained: each component's code is copied in rather
than shelling out to its sibling `src/` directory, so this pipeline keeps
working even if `gbqr_3src_spatial`, `gbqr_3src`, `WISAR6_fourthroot_adaptive_t`,
or `gbqr_4src_nssp` change independently elsewhere in the repo.

## Why three Python environments

- `gbqr/` and `gbqr_us/` depend on `idmodels` pinned to commit `bffe814c`,
  which pins `pandas~=2.0` (pandas 3.0 is a stated incompatibility in
  `idmodels`' own `pyproject.toml`). They share one venv, `.venv-gbqr/`.
- `wisar/` is a from-scratch PyMC/PyTensor implementation with no `idmodels`
  dependency, and its lockfile resolves `pandas>=3.0.6` alongside
  `pymc==5.28.5` / `pytensor==2.38.3` - incompatible with `idmodels`' pin, so
  it needs its own venv, `wisar/.venv/` (managed with `uv`).
- `nssp/` also depends on `idmodels`, but pinned to a later commit
  (`a8ad6840`, with GBQR/feature changes not in `bffe814c`). Reusing
  `.venv-gbqr` would silently run it against the older, untested pin, so it
  gets its own venv too, `.venv-nssp/`.

## To run locally

From this directory (`src/gbqr_3src_spatial_wisar`):

### Python setup (gbqr + gbqr_us, shared venv)

```bash
python -m venv .venv-gbqr
source .venv-gbqr/bin/activate
python -m pip install -r requirements-gbqr.txt
deactivate
```

On Mac, you may also need `brew install libomp` for LightGBM (used by GBQR).

### Python setup (wisar, separate venv via uv)

```bash
cd wisar
uv sync   # must be built on Linux for Unity; do not copy a macOS .venv there
cd ..
```

### Python setup (nssp, separate venv)

```bash
python -m venv .venv-nssp
source .venv-nssp/bin/activate
python -m pip install -r requirements-nssp.txt
deactivate
```

### R setup (ensemble step)

```bash
Rscript -e "renv::restore()"
```

### Running the full pipeline

`main.py` itself only needs `click` and `python-dateutil` (both already in
`.venv-gbqr`), so run it with that venv's Python; it shells out to each
component's own venv for that component's step:

```bash
.venv-gbqr/bin/python main.py --today_date=2024-12-07
```

For a shorter test run (gbqr/gbqr_us/nssp only; wisar/ has no short-run mode):

```bash
.venv-gbqr/bin/python main.py --today_date=2024-12-07 --short_run
```

This runs `gbqr/main.py`, `gbqr_us/main.py`, `wisar/main.py`, and
`nssp/main.py` **concurrently** (none of them depend on each other's output,
only the final ensemble step does), then once all have finished:

- `3_ensemble.R` builds `UMass-gbqr_3src_spatial_wisar` (gbqr + gbqr_us +
  wisar) and `UMass-gbqr_3src_spatial_wisar_nssp` (the same, combined with
  nssp)

`nssp` is skipped automatically for the three known-invalid NSSP reference
dates (see below), and a non-zero exit from it does not stop the pipeline -
only gbqr/gbqr_us/wisar failing does, since the combined output still makes
sense without nssp but not without any of the other three.

Since `gbqr`/`gbqr_us`/`nssp` all fit LightGBM models (which default to using
every available core), `main.py` splits the machine's CPU count three ways
and caps each of their `OMP_NUM_THREADS`/`MKL_NUM_THREADS`/
`OPENBLAS_NUM_THREADS` to that share, so they don't oversubscribe the machine
running at the same time; `wisar`'s own MCMC chain count is controlled
separately via `WISAR6_MCMC_CHAINS`/`WISAR6_MCMC_CORES`.

Each component writes to its own subdirectory under
`intermediate-output/model-output/` (gitignored scratch space); only the two
final outputs above are written to the top-level `../../model-output/`. Each
component's live output is also tailed to the console with a `[label]`
prefix and written to its own `logs/<component>-<reference_date>.log`, since
four processes' output interleaved on one stream would otherwise be
unreadable.

## Running on Unity

```bash
cd src/gbqr_3src_spatial_wisar
module load r/4.4.1-4oqw6bx
python -m venv .venv-gbqr && source .venv-gbqr/bin/activate && \
  python -m pip install -r requirements-gbqr.txt && deactivate
python -m venv .venv-nssp && source .venv-nssp/bin/activate && \
  python -m pip install -r requirements-nssp.txt && deactivate
cd wisar && uv sync && cd ..
Rscript -e "renv::restore()"
mkdir -p logs
sbatch submit-unity-season-2025-26.sh
```

`submit-unity-season-2025-26.sh` runs one reference date per SLURM array
task (this season's dates from `hub-config/tasks.json`), skipping dates whose
combined output (`UMass-gbqr_3src_spatial_wisar_nssp`) already exists (set
`ENSEMBLE_OVERWRITE=1` to force a rerun). Each task calls
`python main.py --today_date=<date>`, which runs all four components
concurrently (12 CPUs requested per task: 3-way split across
gbqr/gbqr_us/nssp, the rest for wisar's 2 MCMC chains) and then the R
ensemble step.

### NSSP's known-invalid dates

`gbqr_4src_nssp` cannot produce valid output for reference dates 2025-11-29,
2025-12-27, and 2026-01-03 (NSSP's holiday-week reporting lag). `main.py`
skips the nssp step automatically for these dates (see
`NSSP_KNOWN_INVALID_DATES`), and `3_ensemble.R` publishes the hosp-only
ensemble under the combined model name for them. If the nssp step fails for
any *other* date, the pipeline logs a warning and continues the same way
rather than failing the whole run - check `logs/` if a combined output file
is missing the ED-visits target for a date not in that known-invalid list.

## Files

- `main.py` - orchestrator; runs the four components and the R ensemble/rbind
  step for one reference date.
- `gbqr/main.py` - GBQR 3-source spatial model (copied from
  `src/gbqr_3src_spatial/main.py`).
- `gbqr_us/main.py` - GBQR 3-source model (copied from `src/gbqr_3src/main.py`).
- `wisar/main.py`, `wisar/model.py` - WISAR6_fourthroot_adaptive_t (copied
  from `src/WISAR6_fourthroot_adaptive_t/`).
- `nssp/main.py` - gbqr_4src_nssp (copied from `src/gbqr_4src_nssp/main.py`).
- `3_ensemble.R` - builds the hosp-only ensemble, then the combined
  (hosp + ED-visits) submission file.
- `requirements-gbqr.txt` - shared pip requirements for `gbqr/` + `gbqr_us/`.
- `requirements-nssp.txt` - pip requirements for `nssp/` (its own idmodels pin).
- `wisar/pyproject.toml`, `wisar/uv.lock` - `wisar/`'s own `uv`-managed deps.
- `renv.lock`, `.Rprofile`, `renv/` - R environment for the ensemble step
  (own copy of `dplyr`/`hubData`/`hubEnsembles`, not shared with other
  models' renv libraries).
- `submit-unity-season-2025-26.sh` - SLURM array job, one task per reference
  date this season.

## Provenance

`UMass-gbqr_3src_spatial_wisar` reproduces the ensemble that previously
existed only as a backfilled pairwise blend, built from the three
components' already-existing `model-output/` files by
`src/gbqr_ar_ensembles/ensemble_pair.R`. That backfill's values were verified
to match this pipeline's ensembling logic exactly (state = mean of gbqr +
wisar, US = mean of gbqr_us + wisar) before this directory was built, so
historical dates in `model-output/UMass-gbqr_3src_spatial_wisar/` need not be
regenerated.

`UMass-gbqr_3src_spatial_wisar_nssp` is new: the combined hosp + ED-visits
submission file this pipeline produces going forward. Note this is a
different model from the pre-existing `UMass-flusion` (an unrelated,
longer-running designated model - a plain GBQR + plain AR6 ensemble, hosp
target only); the two names are easy to confuse but are unrelated models.
