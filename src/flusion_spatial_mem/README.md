# Flusion Spatial MEM Model

This is an experimental version of the flusion_spatial2 ensemble model that incorporates **Moving Epidemic Method (MEM)** based adaptive ensemble weighting. It combines the same component models as flusion_spatial2_prod but uses epidemic phase detection to adaptively weight the ensemble components.

**Status:** v0.2 - MEM-based adaptive weighting is **fully implemented and active**!

## Model Components

The ensemble combines:
- **AR(6) pooled model** - Autoregressive baseline model (SARIX)
- **GBQR 3-source spatial2 model** - Gradient Boosted Quantile Regression with spatial wave features (for state-level forecasts)
- **GBQR 3-source model** - GBQR without spatial features (for US-level forecasts)

## Key Innovation: MEM-Based Adaptive Weighting

Instead of using simple median ensemble (as in flusion_spatial2_prod), this model:

1. **Calculates MEM thresholds** from historical NHSN hospitalization data
2. **Detects current epidemic phase** and intensity level at forecast time
3. **Applies phase-specific ensemble weights** based on which models perform better in different epidemic phases

### Epidemic Phases

The MEM approach classifies current epidemic status into:
- **Baseline**: Activity below epidemic threshold
- **Low intensity**: Above epidemic threshold but below medium intensity
- **Medium intensity**: Moderate epidemic activity
- **High intensity**: High epidemic activity
- **Very high intensity**: Extreme epidemic activity

### Adaptive Weighting Strategy

The ensemble weights are adjusted based on epidemic phase:
- **Baseline/Low**: Higher weight to AR(6) model for stability
- **Medium/High**: Balanced weights across models
- **Very High/Peak**: Higher weight to GBQR spatial models for rapid changes

## To Run Locally

Run the following with `flusion_spatial_mem` as your working directory.

### Python Setup

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

### R Setup

Restore R packages from the lockfile:

```bash
Rscript -e "renv::restore()"
```

This will install the required R packages including:
- `dplyr` - data manipulation
- `hubData` - accessing hub model outputs
- `hubEnsembles` - ensemble utilities
- `mem` - Moving Epidemic Method implementation

### Verifying Setup

Before running the model, verify that all packages are installed correctly:

```bash
Rscript test_mem_setup.R
```

This will test:
- R package installation (mem, fiphde, etc.)
- Data fetching from NHSN via fiphde
- MEM threshold calculation
- Phase detection

If fiphde is not installed, install it with:

```r
remotes::install_github('signaturescience/fiphde')
```

### Running the Model

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
4. **MEM-based adaptive weighted ensemble** (phase detection + adaptive weights)

Output is saved to `model-output/UMass-flusion_spatial_mem/`:
- `YYYY-MM-DD-UMass-flusion_spatial_mem.csv` - Main forecast
- `YYYY-MM-DD-diagnostic-info.csv` - Phase, weights, and threshold info

## Retrospective Evaluation

To generate forecasts for the past two flu seasons (2023-2024 and 2024-2025):

```bash
# Easy: Use the shell script (85 forecasts)
./run_retrospective_evaluation.sh

# Or for testing (faster)
./run_retrospective_evaluation.sh --short-run

# Or use Python directly
python evaluate_mem_ensemble.py --start_date=2023-10-07 --end_date=2025-05-17

# Dry run to see dates without running
python evaluate_mem_ensemble.py --start_date=2023-10-07 --end_date=2025-05-17 --dry_run
```

**Important:** MEM thresholds are calculated adaptively for each forecast date using only data available *at that time*. This ensures truly retrospective evaluation with no look-ahead bias. Early forecasts may have fewer historical seasons available, and the system adapts automatically.

## Dependencies

### Python (`requirements.txt`)
- `click` - command line interface
- `idmodels` - Reich Lab forecasting models (pinned to specific commit)

### R (via `renv.lock`)
- `dplyr` - data manipulation
- `hubData` - accessing hub model outputs
- `hubEnsembles` - creating ensemble forecasts
- `mem` - Moving Epidemic Method for threshold calculation and phase detection

Note: On Mac, you may need to install `libomp` for the GBQR models:

```bash
brew install libomp
```

## Implementation Details

### Files
- `main.py` - Orchestrates the four-step forecasting pipeline
- `0_ar6_pooled.py`, `1_gbqr_3src_spatial2.py`, `2_gbqr_3src.py` - Component models
- `3_flusion_mem_ensemble.R` - MEM-based adaptive ensemble script
- `mem_utils.R` - Helper functions for data fetching, MEM calculation, and phase detection
- `test_mem_setup.R` - Setup verification script
- `evaluate_mem_ensemble.py` - Retrospective evaluation script
- `IMPLEMENTATION_NOTES.md` - Detailed technical documentation

### How It Works

1. **Fetch Historical Data**: Uses fiphde package to get 5 years of NHSN hospitalization data
2. **Calculate MEM Thresholds**: Formats data and applies MEM to establish epidemic/intensity thresholds
3. **Detect Current Phase**: Compares recent hospitalization levels to thresholds
4. **Apply Adaptive Weights**: Uses phase-specific weights (baseline: 50/50, high: 35/65 AR6/GBQR)
5. **Create Ensemble**: Weighted linear pool of component forecasts
6. **Save Diagnostics**: Tracks phase, weights, and thresholds for analysis

### Fallback Behavior

If MEM thresholds cannot be calculated (e.g., insufficient historical data), the model automatically falls back to equal weights (50/50), equivalent to median ensemble.

## Experimental Status

This model is experimental and designed to test whether MEM-based adaptive weighting improves forecast accuracy compared to simple median ensemble (flusion_spatial2_prod). Results will be compared to assess:
- Overall forecast accuracy across all horizons
- Performance during different epidemic phases
- Stability and robustness of adaptive weighting
- Calibration and sharpness by phase

## References

- Vega T, Lozano JE, Meerhoff T, et al. Influenza surveillance in Europe: establishing epidemic thresholds by the moving epidemic method. Influenza Other Respir Viruses. 2013;7(4):546-558.
- Moving Epidemic Method R Package: https://cran.r-project.org/package=mem
