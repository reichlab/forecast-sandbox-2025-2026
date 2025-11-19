# Unity Setup Instructions

Follow these steps to set up the flusion_spatial model on Unity.

## One-Time Setup

### 1. Clone and Navigate to Model

```bash
cd ~/forecast-sandbox-2025-2026
git checkout flusion-spatial-ensemble
git pull origin flusion-spatial-ensemble
cd src/flusion_spatial
```

### 2. Create Logs Directory

```bash
mkdir -p logs
```

### 3. Set Up Python Virtual Environment

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### 4. Install R Packages

This will take 10-20 minutes the first time:

```bash
module load r/4.4.0  # or whatever module test script found
Rscript install_r_packages.R
```

The script will:
- Create `~/R/library` for your user-specific packages
- Install CRAN packages: dplyr, readr, remotes
- Install Hubverse packages: hubData, hubEnsembles
- Install Reich Lab packages: idforecastutils
- Verify all packages installed correctly

### 5. Verify Setup

```bash
./test_unity_env.sh
```

All checks should pass with ✓.

## Running Forecasts

### Test Single Forecast

Before submitting all jobs, test with one date:

```bash
source .venv/bin/activate
module load r/4.4.0
python main.py --today_date=2024-12-28 --short_run
```

This runs with:
- 7 quantile levels (instead of 23)
- 100 MCMC iterations (instead of 2000+2000)
- 10 GBQR bags (instead of 100)

Should complete in ~3-5 minutes.

### Submit All Forecasts

```bash
sbatch submit-unity-parallel.sh
```

This submits 55 array jobs (one per forecast date).

## Monitoring Jobs

### Check Job Status

```bash
squeue -u $USER
```

### View Job Output

```bash
# List log files
ls -lh logs/

# View specific job output
tail -f logs/slurm-JOBID_ARRAYID.out

# Check for errors
grep -r "Error\|error" logs/ | head -20
```

### Cancel Jobs

```bash
# Cancel all your jobs
scancel -u $USER

# Cancel specific job
scancel JOBID

# Cancel specific array task
scancel JOBID_ARRAYID
```

## Output Files

Model outputs are saved to:
```
../../model-output/UMass-flusion_spatial/
```

Each forecast creates a file:
```
YYYY-MM-DD-UMass-flusion_spatial.csv
```

## Troubleshooting

### R Packages Not Found

If jobs fail with R package errors:

```bash
module load r/4.4.0
Rscript -e "library(hubData); library(hubEnsembles); library(idforecastutils)"
```

If any fail, re-run:
```bash
Rscript install_r_packages.R
```

### Python Packages Not Found

```bash
source .venv/bin/activate
pip install -r requirements.txt
```

### Jobs Timeout

Default time limit is 2 hours. If jobs timeout:

1. Check logs to see which component is slow
2. Consider increasing time limit in `submit-unity-parallel.sh`:
   ```bash
   #SBATCH -t 03:00:00  # 3 hours instead of 2
   ```

### Memory Issues

Default memory is 32GB. If jobs run out of memory, increase in `submit-unity-parallel.sh`:
```bash
#SBATCH --mem=64G  # 64GB instead of 32GB
```

## Performance

Expected runtimes per forecast (full run, 8 cores):
- AR(6) model: ~10-15 seconds
- GBQR spatial model: ~10-15 minutes
- R ensemble: ~5-10 seconds
- **Total**: ~15-20 minutes per forecast date

With 55 dates running in parallel, all forecasts should complete in ~15-20 minutes.
