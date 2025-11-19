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

### One-Time Setup

**Step 1: Clone and navigate to model**
```bash
cd ~/forecast-sandbox-2025-2026
git checkout flusion-spatial-ensemble
git pull origin flusion-spatial-ensemble
cd src/flusion_spatial
```

**Step 2: Create logs directory**
```bash
mkdir -p logs
```

**Step 3: Set up Python virtual environment**
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

**Step 4: Install R packages** (takes 10-20 minutes first time)
```bash
module load r/4.4.0
Rscript install_r_packages.R
```

This installs required R packages to your user library (`~/R/library`):
- CRAN packages: dplyr, readr, remotes
- Hubverse packages: hubData, hubEnsembles (from GitHub)
- Reich Lab packages: idforecastutils (from GitHub)

**Note**: If `hubData` fails to install, you may need to install it manually:
```bash
R
```
Then in R console:
```r
install.packages("remotes")
remotes::install_github("hubverse-org/hubData")
quit()
```

**Step 5: Verify setup**
```bash
./test_unity_env.sh
```

All 5 checks should pass with ✓:
1. ✓ Virtual environment found and activated
2. ✓ Python packages installed
3. ✓ R module loaded (auto-detects r/4.4.0 or R/4.3.2-gfbf-2023a)
4. ✓ R packages installed
5. ✓ Python can find Rscript

If any checks fail, see `SETUP.md` for troubleshooting.

**Step 6: Test a single forecast** (optional but recommended)
```bash
source .venv/bin/activate  # IMPORTANT: Must activate venv first!
module load r/4.4.0
python main.py --today_date=2024-12-28 --short_run
```

**Note**: Always activate the virtual environment with `source .venv/bin/activate` before running Python scripts. You should see `(.venv)` at the start of your prompt when activated.

Should complete in ~3-5 minutes with the `--short_run` flag.

### Submit Jobs

Once setup is complete:

```bash
sbatch submit-unity-parallel.sh
```

This will submit 55 array jobs (one for each forecast date) to the Unity cluster. Each job:
- Uses 8 CPU cores
- Requests 32GB memory
- Has a 2-hour time limit
- Runs the full ensemble pipeline (AR6 + GBQR + Ensemble)
- Expected runtime: ~15-20 minutes per forecast date

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

Check for errors:
```bash
grep -r "Error\|error" logs/ | head -20
```

### Job Configuration

The dates array in `submit-unity-parallel.sh` covers:
- 2023-24 season: October 2023 through May 2024
- 2024-25 season: November 2024 through May 2025

To modify the dates or add new forecast dates, edit the `dates` array in the script.

### Troubleshooting

**Common Issues:**

1. **ModuleNotFoundError: No module named 'dateutil'**
   - Solution: Activate the virtual environment first: `source .venv/bin/activate`
   - You should see `(.venv)` at the start of your prompt

2. **hubData installation fails**
   - Solution: Install manually in R:
     ```r
     install.packages("remotes")
     remotes::install_github("hubverse-org/hubData")
     ```

See `SETUP.md` for more detailed troubleshooting steps, including:
- R package installation issues
- Module loading problems
- Memory/timeout issues
- Performance tuning

## Model Metadata

Model metadata is stored in `model-metadata/UMass-flusion_spatial.yaml` following the Hubverse schema.

## Output Format

Model outputs follow the Hubverse format:
- Filename: `{reference_date}-UMass-flusion_spatial.csv`
- Columns: `reference_date`, `target`, `horizon`, `location`, `target_end_date`, `output_type`, `output_type_id`, `value`
- Output type: `quantile` (23 quantile levels from 0.01 to 0.99)
- Horizons: 0, 1, 2, 3 weeks ahead
- Locations: US and all states/territories (excluding AK, CT, HI)
