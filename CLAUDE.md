# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is the Reich Lab's forecasting testbed for experimental respiratory virus forecasts for the 2025/2026 season. It follows the Hubverse format (modeled after the FluSight Forecast Hub) and contains implementations of forecasting models from the [idmodels](https://github.com/reichlab/idmodels) library, with data pipelines from [iddata](https://github.com/reichlab/iddata/).

The repository structure follows a hub-based architecture with:
- **Model implementations** in `src/` (each model has its own directory)
- **Model outputs** stored in `model-output/` (organized by model name)
- **Model metadata** in `model-metadata/` (YAML files describing each model)
- **Hub configuration** in `hub-config/` (defines tasks, schema, and admin settings)

## Repository Structure

### Hub Configuration (`hub-config/`)
- `tasks.json`: Defines forecast targets (weekly incident flu hospitalizations), horizons (-1 to 3 weeks), locations (US + states), required quantile levels (23 quantiles from 0.01 to 0.99), and optional sample outputs
- `admin.json`: Hub metadata including maintainer info and repository details
- `model-metadata-schema.json`: Schema for model metadata files

### Model Implementations (`src/`)
Each model directory contains:
- `main.py`: Entry point for running the model
- `requirements.txt`: Python dependencies (usually includes `click` and `idmodels` from GitHub)
- `README.md`: Model-specific documentation

Models are named descriptively (e.g., `AR6_pooled`, `AR6_fourier_pooled_both`) and use variants of:
- **SARIX models**: Seasonal AutoRegressive Integrated with eXogenous variables
- **GBQR models**: Gradient Boosted Quantile Regression
- **Ensemble models**: Combinations of multiple base models

### Flusion Ensemble Model
The `flu_flusion` model is a multi-stage pipeline:
1. `0_ar6_pooled.py`: Runs AR(6) pooled model
2. `1_gbqr.py`: Runs GBQR model
3. `2_flusion_ensemble.R`: Combines outputs using `hubEnsembles::simple_ensemble()`
4. All three are orchestrated by `main.py`

Intermediate outputs are stored in `src/flu_flusion/intermediate-output/model-output/` before final ensemble creation.

## Common Commands

### Running a Model Locally

Models use Python virtual environments. From a model directory (e.g., `src/AR6_pooled/`):

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python main.py --today_date=2024-01-06
```

Most models accept:
- `--today_date` (optional): Effective model run date in YYYY-MM-DD format; defaults to today
- Models automatically calculate reference_date as the following Saturday using `relativedelta(weekday=5)`

The flusion model also accepts:
- `--short_run`: Flag for abbreviated runs during development

### Running Tests and Validation

Pre-commit hooks are configured for code quality:

```bash
pip install pre-commit
pre-commit install
pre-commit run --all-files
```

Hooks include:
- Trailing whitespace removal
- YAML validation (allows multiple documents)
- AWS credentials detection
- Private key detection
- Markdown linting (markdownlint with specific rules disabled: MD007, MD013, MD024, MD033, MD041)
- Spelling checks (codespell)

### GitHub Workflows

Three automated workflows validate hub submissions:
1. **validate-submission.yaml**: Validates model outputs using R's `hubValidations` package on PRs
2. **validate-config.yaml**: Validates hub configuration files
3. **cache-hubval-deps.yaml**: Caches R dependencies for faster validation

## Model Configuration Patterns

All models follow a consistent configuration pattern using `SimpleNamespace`:

### Model Config
- `model_class`: Type of model ("sarix", "gbqr")
- `model_name`: Unique identifier
- `sources`: Data sources (e.g., ["nhsn"], ["flusurvnet", "nhsn", "ilinet"])
- `fit_locations_separately`: Boolean for pooled vs. unpooled fitting
- `power_transform`: Transform applied to data (typically "4rt" for fourth root)
- Model-specific parameters (e.g., `p`, `P`, `d`, `D` for SARIX; `num_bags` for GBQR)
- Pooling parameters: `theta_pooling`, `sigma_pooling` (values: "shared", "none")
- `x`: List of covariates (can be empty)

### Run Config
- `disease`: Always "flu" in this hub
- `ref_date`: Saturday reference date for the forecast
- `output_root`: Path to output directory (usually `Path("../../model-output/")`)
- `locations`: List of all forecast locations (US + 52 states/territories, excluding AK, CT, HI)
- `q_levels` and `q_labels`: 23 quantile levels matching hub requirements
- `max_horizon`: Maximum forecast horizon (typically 4 for horizons -1 to 3)
- MCMC parameters for Bayesian models: `num_warmup`, `num_samples`, `num_chains`

## Data Transformations

Models typically apply these transformations in sequence:
1. Convert hospitalization counts to rate per 100,000 population
2. Apply power transform (fourth root)
3. Scale and center (divide by 95th percentile, subtract mean)

Predictions are inverted back to original scale before output.

## Output Format

Model outputs are CSV files following the Hubverse format:
- Filename: `{reference_date}-{team_abbr}-{model_abbr}.csv`
- Required columns: `reference_date`, `target`, `horizon`, `location`, `target_end_date`, `output_type`, `output_type_id`, `value`
- Output types: `quantile` (required, 23 levels) and optionally `sample` (100 samples)
- All quantile values must be non-negative

## Model Metadata

Each model requires a YAML metadata file in `model-metadata/` with:
- Team information (`team_name`, `team_abbr`)
- Model details (`model_name`, `model_abbr`, `model_version`)
- Contributors with affiliations and emails
- `license`: Typically "CC-BY-4.0"
- `designated_model`: Boolean indicating if model is designated for official submissions
- `methods`: Short description
- `data_inputs`: Data sources used
- `methods_long`: Detailed methodology including transformations
- `ensemble_of_models` and `ensemble_of_hub_models`: Booleans

## Date Handling

All dates use ISO format (YYYY-MM-DD). The reference_date is always a Saturday, calculated from today_date using:
```python
reference_date = today_date + relativedelta.relativedelta(weekday=5)
```

Forecast horizons:
- Horizon -1: Previous week (nowcast of past data)
- Horizon 0: Current week ending on reference_date
- Horizon 1-3: Future weeks

Target end dates are Saturdays corresponding to the forecast horizon.

## Dependencies

Python models depend on:
- `idmodels`: Reich Lab's forecasting models library (installed from GitHub)
- `click`: Command-line interface creation
- `dateutil`: Date manipulation

R scripts (ensemble models) use:
- `hubData`: Loading model outputs from Hubverse hubs
- `hubEnsembles`: Creating ensemble forecasts
- `idforecastutils`: Reich Lab utilities for forecast transformations

Both Python and R may use additional packages for data manipulation (`pandas`, `numpy`, `dplyr`, `readr`).
