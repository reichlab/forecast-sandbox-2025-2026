# Flusion3 - Three-Component Ensemble Model for Influenza

This is an ensemble model that combines three components:
1. AR(6) pooled model (using NHSN data)
2. AR(6) pooled model with Fourier seasonality (K=2 harmonics, using NHSN data)
3. Gradient Boosted Quantile Regression (using NHSN, FluSurv-NET, and ILINet data)

The three component models are combined using equal-weighted quantile averaging via `hubEnsembles::simple_ensemble()`.

## Rationale

This model explores whether combining both AR6 variants (with and without Fourier) alongside GBQR can improve forecast performance by:
- Capturing both short-term AR dynamics and longer-term seasonal patterns
- Diversifying the ensemble with complementary modeling approaches
- Reducing variance through increased ensemble size

## Prerequisites

This model requires:
- Python 3.11+ (for reichlab/sarix with Fourier support)
- Local development versions of `sarix` and `idmodels` packages at:
  - `../../../sarix` (relative to this directory)
  - `../../../idmodels` (relative to this directory)

## Running Locally

To test this out locally, run the following with this directory (`flusion3`) as your working directory.

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
- Reduce MCMC warmup and samples from 2000 to 100 for both AR6 models
- Reduce quantile levels from 23 to 7
- Reduce GBQR bagging from 100 to 10 bags

## Output

The model will create outputs in:
- `intermediate-output/model-output/` - Individual component model outputs
- `../../model-output/UMass-flusion3/` - Final ensemble output

## Model Components

### 0_ar6_pooled.py
Runs an AR(6) model with:
- Pooled AR parameters (shared across locations)
- Location-specific variance
- **NO Fourier seasonality**
- NHSN data only
- Fourth-root power transformation

### 1_ar6_pooled_fourier.py
Runs an AR(6) model with:
- Pooled AR parameters (shared across locations)
- Location-specific variance
- **Fourier seasonality with K=2 harmonics**
- NHSN data only
- Fourth-root power transformation

### 2_gbqr.py
Runs gradient boosted quantile regression with:
- Multiple data sources (NHSN, FluSurv-NET, ILINet)
- 100 bootstrap bags (10 in short run)
- Level features included
- Fourth-root power transformation

### 3_flusion3_ensemble.R
Creates an equal-weighted quantile average of all three component models using the hubEnsembles package.
