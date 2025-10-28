# Flusion Fourier Ensemble Model for Influenza

This is an ensemble model that combines:
1. AR(6) pooled model with Fourier seasonality (using NHSN data)
2. Gradient Boosted Quantile Regression (using NHSN, FluSurv-NET, and ILINet data)

The two component models are combined using equal-weighted quantile averaging via `hubEnsembles::simple_ensemble()`.

## Key Difference from UMass-Flusion

The main difference from the standard `flusion` model is that the AR(6) component includes **Fourier seasonality terms** (K=2 harmonics) to explicitly model annual seasonal patterns in flu hospitalizations.

## Prerequisites

This model requires local development versions of `sarix` and `idmodels` packages to be available at:
- `../../../sarix` (relative to this directory)
- `../../../idmodels` (relative to this directory)

The Fourier seasonality functionality may not yet be available in released versions.

## Running Locally

To test this out locally, run the following with this directory (`flusion_fourier`) as your working directory.

### Full Run

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --today_date=2024-01-06
```

### Short Run (for testing)

For faster testing with reduced MCMC samples and fewer quantiles:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --today_date=2024-01-06 --short_run
```

The `--short_run` flag will:
- Reduce MCMC warmup and samples from 2000 to 100 for AR6
- Reduce quantile levels from 23 to 7
- Reduce GBQR bagging from 100 to 10 bags

## Output

The model will create outputs in:
- `intermediate-output/model-output/` - Individual component model outputs
- `../../model-output/UMass-flusion_fourier/` - Final ensemble output

## Model Components

### 0_ar6_pooled.py
Runs an AR(6) model with:
- Pooled AR parameters (shared across locations)
- Location-specific variance
- Fourier seasonality with K=2 harmonics
- NHSN data only
- Fourth-root power transformation

### 1_gbqr.py
Runs gradient boosted quantile regression with:
- Multiple data sources (NHSN, FluSurv-NET, ILINet)
- 100 bootstrap bags (10 in short run)
- Level features included
- Fourth-root power transformation

### 2_flusion_ensemble.R
Creates an equal-weighted quantile average of the two component models using the hubEnsembles package.
