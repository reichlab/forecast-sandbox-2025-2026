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

### R setup on Unity

The renv lockfile requires **R 4.5.1**. On Unity, we use the Rocker ML Verse container via Apptainer to get this version.

#### About the R module

The module `r-rocker-ml-verse/4.5.1_cuda12.8.1+apptainer` is:
- **r-rocker-ml-verse**: A container image from the [Rocker Project](https://rocker-project.org/) that includes R, tidyverse, and machine learning libraries
- **4.5.1**: R version 4.5.1
- **cuda12.8.1**: CUDA toolkit for GPU computing (not used by this model, but included in the container)
- **apptainer**: The container runtime (formerly Singularity) that runs the R environment

#### Quick installation check

Before setting up renv, verify the module is available and working:

```bash
# Check if the module exists
module spider r-rocker-ml-verse/4.5.1

# Load the module
module load r-rocker-ml-verse/4.5.1_cuda12.8.1+apptainer

# Verify R version and that it can see your project files
Rscript -e "print(R.version.string); print(.libPaths())"
```

If successful, you should see `R version 4.5.1` in the output.

#### Option 1: Interactive setup via command line

```bash
# Load the required R module (containerized R 4.5.1)
module load r-rocker-ml-verse/4.5.1_cuda12.8.1+apptainer

# Navigate to the project directory
cd /path/to/forecast-sandbox-2025-2026/src/flusion_spatial2_prod

# Restore R packages
Rscript -e "renv::restore()"
```

#### Option 2: Interactive setup via RStudio Server

1. Log in to Unity's RStudio Server (https://unity.rc.umass.edu)
2. In the terminal pane, load the R module:
   ```bash
   module load r-rocker-ml-verse/4.5.1_cuda12.8.1+apptainer
   ```
3. Restart R (Session > Restart R) to pick up the new R version
4. Navigate to `src/flusion_spatial2_prod` and set it as your working directory
5. Run `renv::restore()` to install packages into the project library

**Important**: If renv was previously initialized with a different R version, you may need to:
1. Delete the `renv/library/` directory
2. Run `renv::restore()` again with R 4.5.1 loaded

#### Verifying your setup

Before submitting jobs, verify R is correctly configured:

```bash
module load r-rocker-ml-verse/4.5.1_cuda12.8.1+apptainer
Rscript --version  # Should show R version 4.5.1
cd src/flusion_spatial2_prod
Rscript -e "renv::status()"  # Should show all packages are installed
```

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
