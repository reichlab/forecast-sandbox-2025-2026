# Flusion Spatial Ensemble Model

This is an ensemble model inspired by the original flusion model, but incorporating spatial features through the gbqr_3src_spatial component. The ensemble combines forecasts from two component models:

1. **AR(6) Pooled Model**: A seasonal autoregressive model with 6 lags using NHSN hospitalization data
2. **GBQR 3-Source Spatial Model**: A gradient boosted quantile regression model using FluSurvNet, NHSN, and ILINet data sources with directional wave features for spatial modeling

The ensemble is created using the `simple_ensemble()` function from the `hubEnsembles` package, which computes the median of the quantile forecasts from the two component models.

## Model Components

### 0_ar6_pooled.py
- Model: AR(6) with pooled estimation across locations
- Data source: NHSN hospitalization data
- Transformation: Fourth root power transform
- Parameter pooling: Shared theta, separate sigma by location
- MCMC sampling: 2000 warmup + 2000 samples per chain

### 1_gbqr_3src_spatial.py
- Model: Gradient Boosted Quantile Regression with spatial features
- Data sources: FluSurvNet, NHSN, ILINet
- Transformation: Fourth root power transform
- Spatial features: Directional wave features (8 directions, 2 temporal lags, 1500km max distance)
- Bagging: 100 bags with 70% sample fraction

### 2_flusion_spatial_ensemble.R
- Combines component model outputs using median ensemble
- Outputs quantile forecasts for horizons 0-3 weeks ahead
- Saves results to `model-output/UMass-flusion_spatial/`

## Running Locally

To test this model locally, navigate to this directory and run:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

# Run with default (today's) date
python main.py

# Run with specific date
python main.py --today_date=2024-01-06

# Run with short_run flag for faster testing
python main.py --today_date=2024-01-06 --short_run
```

The `--short_run` flag reduces computation time by:
- Using fewer quantile levels (7 instead of 23)
- Reducing MCMC iterations (100 warmup + 100 samples)
- Using only 10 bags for GBQR (instead of 100)

## Running on Unity (SLURM)

This model includes a SLURM submission script for parallel execution on the Unity cluster.

### Setup

1. Create and activate a virtual environment:
```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

2. Ensure the `logs/` directory exists:
```bash
mkdir -p logs
```

### Submit Jobs

The `submit-unity-parallel.sh` script runs the model for multiple forecast dates in parallel:

```bash
sbatch submit-unity-parallel.sh
```

This will submit 55 array jobs (one for each forecast date) to the Unity cluster. Each job:
- Uses 8 CPU cores
- Requests 32GB memory
- Has a 2-hour time limit
- Runs the full ensemble pipeline (AR6 + GBQR + Ensemble)

### Monitor Jobs

Check job status:
```bash
squeue -u $USER
```

View output logs:
```bash
ls -lh logs/
tail logs/slurm-<JOB_ID>_<ARRAY_ID>.out
```

### Job Configuration

The dates array in `submit-unity-parallel.sh` covers:
- 2023-24 season: October 2023 through May 2024
- 2024-25 season: November 2024 through May 2025

To modify the dates or add new forecast dates, edit the `dates` array in the script.

## Model Metadata

Model metadata is stored in `model-metadata/UMass-flusion_spatial.yaml` following the Hubverse schema.

## Output Format

Model outputs follow the Hubverse format:
- Filename: `{reference_date}-UMass-flusion_spatial.csv`
- Columns: `reference_date`, `target`, `horizon`, `location`, `target_end_date`, `output_type`, `output_type_id`, `value`
- Output type: `quantile` (23 quantile levels from 0.01 to 0.99)
- Horizons: 0, 1, 2, 3 weeks ahead
- Locations: US and all states/territories (excluding AK, CT, HI)
