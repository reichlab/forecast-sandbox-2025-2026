# Flusion Spatial2 Production Model

This is the production version of the flusion_spatial2 ensemble model for influenza forecasting. It combines:
- AR(6) pooled model
- GBQR 3-source spatial2 model (for state-level forecasts)
- GBQR 3-source model (for US-level forecasts)

## To run locally

Run the following with `flusion_spatial2_prod` as your working directory.

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

On Unity, you can log in to an RStudio Server, navigate to the `src/flusion_spatial2_prod` folder, set it as your working directory and type `renv::restore()` to activate the project and create a project library. 
It may require restarting the R session and re-running `renv::restore()`.

### Running the model

```bash
python main.py --today_date=2024-12-07
```

For a shorter test run:

```bash
python main.py --today_date=2024-12-07 --short_run
```

This runs all four steps:
1. AR(6) pooled model
2. GBQR 3-source spatial2 model
3. GBQR 3-source model
4. R ensemble combining the outputs

Output is saved to `model-output/UMass-flusion_spatial/`.

## requirements.txt and renv.lock details

Python dependencies (`requirements.txt`):
- `click` - command line interface
- `idmodels` - Reich Lab forecasting models (pinned to specific commit)

R dependencies (via `renv.lock`):
- `dplyr` - data manipulation
- `hubData` - accessing hub model outputs
- `hubEnsembles` - creating ensemble forecasts

Note: On Mac, you may need to install `libomp` for the GBQR models:

```bash
brew install libomp
```
