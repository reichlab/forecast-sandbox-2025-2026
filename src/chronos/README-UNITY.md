# Running Chronos on Unity Cluster

This document describes how to run the Chronos time series foundation model on the Unity HPC cluster for flu forecasting.

## Prerequisites

1. **Unity Account**: Ensure you have an active Unity cluster account
2. **Python Environment**: Set up a virtual environment with all dependencies
3. **File Access**: Ensure the repository is accessible from Unity (clone or mount)

## Setup

### 1. Clone repository and navigate to chronos directory

```bash
cd ~
git clone https://github.com/reichlab/forecast-sandbox-2025-2026.git
cd forecast-sandbox-2025-2026/src/chronos
```

### 2. Create logs directory

```bash
mkdir -p logs
```

### 3. Set up Python virtual environment

```bash
# Load Python module (check available versions with: module avail python)
module load python/3.11  # or appropriate version

# Create virtual environment
python -m venv .venv

# Activate it
source .venv/bin/activate

# Upgrade pip
python -m pip install --upgrade pip

# Install dependencies
python -m pip install -r requirements.txt

# Test installation
python -c "from chronos import ChronosPipeline; print('Chronos installed successfully!')"
```

**Note**: The first forecast run will download pre-trained model weights (~200MB for small, ~800MB for base). These are cached for subsequent runs.

### 4. Configure the submission scripts

Edit `run_chronos.slurm` or `run_chronos_array.slurm` and update:

- `--mail-user=trobacker@umass.edu` - Your email address for notifications
- Adjust resource requests based on model size:
  - `--cpus-per-task=8` - Number of CPU cores
  - `--mem=16G` - Memory per job (12G for small, 16G for base, 32G for large)
  - `--time=00:30:00` - Time limit (15 min for small, 20 min for base, 45 min for large)
- Verify `--partition=cpu` matches your allocation

## Submission Scripts

### run_chronos.slurm (Single Forecast)

Runs a single forecast for a specific date.

**Features**:
- Configurable via environment variables
- Supports all 5 model sizes (tiny/mini/small/base/large)
- Adjustable context length and sampling parameters

**Submit with:**
```bash
# Default settings (today's date, small model)
sbatch run_chronos.slurm

# Specific date and model size
sbatch --export=TODAY_DATE=2025-01-04,MODEL_SIZE=base run_chronos.slurm

# Full configuration
sbatch --export=TODAY_DATE=2025-01-11,MODEL_SIZE=small,CONTEXT_LEN=104,NUM_SAMPLES=200 run_chronos.slurm
```

### run_chronos_array.slurm (Multiple Dates in Parallel)

Runs multiple forecast dates as parallel array jobs.

**Features**:
- Predefined list of forecast dates (edit DATES array in script)
- Each date runs as independent array task
- Automatic array index to date mapping

**Configure and submit:**
```bash
# 1. Edit run_chronos_array.slurm and update DATES array
# 2. Update --array parameter to match number of dates (e.g., --array=0-9 for 10 dates)
# 3. Submit
sbatch run_chronos_array.slurm

# Or with throttling (max 10 concurrent jobs)
# Edit script to change: --array=0-9%10
sbatch run_chronos_array.slurm
```

**Example dates configuration** (weekly forecasts):
```bash
DATES=(
    "2024-12-07"
    "2024-12-14"
    "2024-12-21"
    "2024-12-28"
    "2025-01-04"
    "2025-01-11"
    "2025-01-18"
    "2025-01-25"
    "2025-02-01"
    "2025-02-08"
)
```

## Monitoring Jobs

### Check job status
```bash
squeue -u $USER
```

### Check specific job details
```bash
scontrol show job JOBID
```

### View output logs
```bash
# List all logs
ls logs/

# View specific log (single job)
cat logs/chronos_JOBID.out

# View array job log
cat logs/chronos_array_ARRAYJOBID_TASKID.out

# View errors
cat logs/chronos_JOBID.err

# Tail running job
tail -f logs/chronos_JOBID.out
```

### Check array job summary
```bash
sacct -j JOBID --format=JobID,JobName,Partition,State,ExitCode,Elapsed
```

### Check forecast outputs
```bash
# List generated forecasts
ls ../../model-output/UMass-chronos_small/

# Count successful forecasts
ls ../../model-output/UMass-chronos_small/ | wc -l
```

## Managing Jobs

### Cancel job
```bash
# Cancel single job
scancel JOBID

# Cancel entire array job
scancel ARRAYJOBID

# Cancel specific array task
scancel ARRAYJOBID_TASKID
# Example: scancel 123456_5

# Cancel all your jobs
scancel -u $USER
```

### Adjust concurrent job limit (after submission)
```bash
scontrol update jobid=JOBID arraytaskthrottle=5
```

## Model Size Selection

Choose based on your accuracy vs. speed requirements:

| Model Size | Parameters | CPUs | Memory | Time Limit | Typical Runtime |
|-----------|-----------|------|--------|------------|-----------------|
| tiny      | 8M        | 4    | 8GB    | 10 min     | ~2 min          |
| mini      | 20M       | 8    | 10GB   | 10 min     | ~3 min          |
| **small** | **46M**   | **8**| **12GB**| **15 min** | **~5 min**      |
| base      | 200M      | 8    | 16GB   | 20 min     | ~8 min          |
| large     | 710M      | 16   | 32GB   | 45 min     | ~25 min         |

**Recommended**: Start with **small** model - best balance of accuracy and speed.

## Measured Performance

Based on testing with **Chronos-T5-Base** on Unity CPU nodes:
- **Average runtime**: 8 min 10 sec (std: 5 sec)
- **Configuration**: 8 CPUs, 16GB RAM
- **Test dates**: 2024-12-28, 2025-01-04, 2025-01-11
- **Output**: 6,095 forecast rows (23 quantiles × 5 horizons × 53 locations)

## GPU Acceleration (Optional)

For faster inference, use GPU partitions:

### 1. Modify SLURM script header:
```bash
#SBATCH --partition=gpu          # or gpu-long, gpu-preempt
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
```

### 2. Update main.py (line 89):
```python
device_map="cuda",  # Change from "cpu"
```

**Expected speedup**: 3-5x faster than CPU

## Troubleshooting

### Jobs fail with memory errors
- Increase `--mem` in SBATCH header
- Use smaller model size (small instead of base)
- Reduce `--num_samples` parameter
- Check logs for "Out of Memory" or "OOM" messages

### Jobs timeout
- Increase `--time` limit in SBATCH header
- Use smaller model size
- Reduce context length: `--export=CONTEXT_LEN=52`

### Module not found errors
```bash
# Check available Python modules
module avail python

# Load appropriate version
module load python/3.11

# Verify virtual environment
source .venv/bin/activate
python -c "import chronos"
```

### First run takes extra long
- Normal! First run downloads model weights (~200-800MB depending on size)
- Subsequent runs use cached weights and are much faster
- Downloaded to `~/.cache/huggingface/`

### Virtual environment not found
```bash
cd ~/forecast-sandbox-2025-2026/src/chronos
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Some array jobs succeed, others fail
```bash
# Check which jobs failed
sacct -j ARRAYJOBID --format=JobID,State,ExitCode | grep FAILED

# Check individual error logs
cat logs/chronos_array_*_TASKID.err

# Resubmit specific failed indices
# Edit run_chronos_array.slurm to change --array=5,12,33
sbatch run_chronos_array.slurm
```

## Example Workflows

### Workflow 1: Single forecast (testing)
```bash
# 1. Set up environment (one-time)
cd ~/forecast-sandbox-2025-2026/src/chronos
mkdir -p logs
module load python/3.11
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# 2. Run test forecast
sbatch --export=TODAY_DATE=2025-01-04,MODEL_SIZE=small run_chronos.slurm

# 3. Monitor
squeue -u $USER
tail -f logs/chronos_*.out

# 4. Check output
ls ../../model-output/UMass-chronos_small/
```

### Workflow 2: Multiple weekly forecasts (production)
```bash
# 1. Edit run_chronos_array.slurm
nano run_chronos_array.slurm
# - Update DATES array with your forecast dates
# - Update --array parameter (e.g., --array=0-9 for 10 dates)
# - Add throttling if needed (e.g., --array=0-9%10 for max 10 concurrent)
# - Update --mail-user to your email

# 2. Submit array job
sbatch run_chronos_array.slurm

# 3. Monitor progress
watch -n 30 'squeue -u $USER'

# 4. Check results
sacct -j JOBID --format=JobID,State,ExitCode,Elapsed
grep "completed successfully" logs/chronos_array_*.out | wc -l

# 5. Check for failures
grep -l "Error\|Traceback" logs/chronos_array_*.err
```

### Workflow 3: Batch processing multiple dates sequentially
```bash
# Submit multiple single jobs with slight delays
for date in 2024-12-28 2025-01-04 2025-01-11 2025-01-18; do
    sbatch --export=TODAY_DATE=$date,MODEL_SIZE=small run_chronos.slurm
    sleep 2  # Avoid overloading scheduler
done
```

## Resource Optimization Tips

1. **Model size**: Start with `small` (46M) - nearly as accurate as `base` (200M) but 40% faster
2. **Context length**: Default 104 weeks (2 years) is usually sufficient; more context ≠ better forecasts
3. **Sampling**: 200 samples gives smooth quantiles; increasing to 500+ has diminishing returns
4. **Parallelization**: Use array jobs for multiple dates to maximize cluster utilization
5. **Throttling**: Add `%10` to array directive (e.g., `--array=0-9%10`) to limit concurrent jobs

## Output Location

Forecasts are saved to:
```
../../model-output/UMass-chronos_<model_size>/YYYY-MM-DD-UMass-chronos_<model_size>.csv
```

Example:
```
model-output/UMass-chronos_small/2025-01-04-UMass-chronos_small.csv
```

Each file contains 6,095 rows (23 quantiles × 5 horizons × 53 locations).

## Additional Resources

- [Unity Documentation](https://docs.unity.rc.umass.edu/)
- [Slurm Array Jobs](https://docs.unity.rc.umass.edu/documentation/jobs/sbatch/arrays/)
- [Unity Support](mailto:hpc@umass.edu)
- [Chronos Documentation](https://github.com/amazon-science/chronos-forecasting)
- [Model Details](../../docs/chronos-foundation-model.md)
- [SLURM_README.md](./SLURM_README.md) - Detailed SLURM usage guide
