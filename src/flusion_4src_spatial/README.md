# Flusion 4-Source Spatial Model

This is a combine-only ensemble model: it does not fit anything itself. It takes a
per-quantile median (via `hubEnsembles::simple_ensemble`) of forecasts already
published in `../../model-output/` by other models in this repo:

- **State-level forecasts**: `UMass-gbqr_4src_spatial` + `UMass-AR6_pooled`
- **US-level forecasts**: `UMass-gbqr_4src` + `UMass-AR6_pooled`

`gbqr_4src_spatial` has no US-level output because its directional wave spatial
features only support a single aggregation level, so `gbqr_4src` (the non-spatial
4-source GBQR variant with the same source list) stands in for US-level forecasts.

This is the sibling of `flusion_4src_spatial_fourierP`, which uses
`AR6_fourierP_thetaP` instead of `AR6_pooled` as its AR(6) component.

Both component models must already have a published forecast for the reference
date before running this model -- `main.py` checks for this and exits with an
error naming any missing file rather than silently producing a partial ensemble.

## To run locally

Run the following with `flusion_4src_spatial` as your working directory.

### Python setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Alternative using `uv`:

```bash
uv venv
uv pip install -r requirements.txt
```

### R setup

Restore R packages from the lockfile:

```bash
Rscript -e "renv::restore()"
```

### Running the model

```bash
python main.py --today_date=2026-01-07
```

Output is saved to `model-output/UMass-flusion_4src_spatial/`.

## requirements.txt and renv.lock details

Python dependencies (`requirements.txt`):
- `click` - command line interface
- `python-dateutil` - reference date calculation

R dependencies (via `renv.lock`, copied from `flusion_spatial2_prod` since the
package set is a superset of what's needed):
- `dplyr` - data manipulation
- `readr` - reading component model CSVs directly
- `hubEnsembles` - creating ensemble forecasts
