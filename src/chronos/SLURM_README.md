# Running Chronos on Unity HPC Cluster

This directory contains SLURM scripts for running Chronos forecasts on the UMass Unity cluster.

## Initial Setup on Unity

### 1. Clone Repository and Navigate to Chronos Directory

```bash
# SSH to Unity
ssh your_username@unity.rc.umass.edu

# Clone repo (if not already done)
cd ~
git clone https://github.com/reichlab/forecast-sandbox-2025-2026.git
cd forecast-sandbox-2025-2026/src/chronos
```

### 2. Set Up Python Environment

```bash
# Load Python module (check available versions with: module avail python)
module load python/3.11  # or whatever version is available

# Create virtual environment
python -m venv .venv

# Activate environment
source .venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install dependencies
pip install -r requirements.txt

# Test installation
python -c "from chronos import ChronosPipeline; print('Success!')"
```

**Note**: The first time you run a forecast, Chronos will download the pre-trained model weights (~200MB for small, ~800MB for base). This is cached for future runs.

### 3. Create Logs Directory

```bash
mkdir -p logs
```

## Running Forecasts

### Single Forecast Job

Run a single forecast for a specific date:

```bash
# Default (today, small model)
sbatch run_chronos.slurm

# Specific date and model size
sbatch --export=TODAY_DATE=2025-01-04,MODEL_SIZE=base run_chronos.slurm

# With all options
sbatch --export=TODAY_DATE=2025-01-11,MODEL_SIZE=small,CONTEXT_LEN=104,NUM_SAMPLES=200 run_chronos.slurm
```

### Array Job (Multiple Dates in Parallel)

To run multiple forecast dates simultaneously:

1. **Edit** `run_chronos_array.slurm` and update the `DATES` array with your desired dates
2. **Update** the `--array` parameter to match the number of dates (e.g., `--array=0-9` for 10 dates)
3. **Submit**:
   ```bash
   sbatch run_chronos_array.slurm
   ```

Example for 10 weekly forecasts:
```bash
# Edit DATES array in run_chronos_array.slurm:
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

# Then submit:
sbatch run_chronos_array.slurm
```

## Model Size Selection

Choose based on your accuracy vs. speed requirements:

| Model Size | Parameters | Typical Runtime | Memory | Recommended Use |
|------------|-----------|----------------|--------|-----------------|
| tiny       | 8M        | ~2 min         | 8GB    | Quick tests, prototyping |
| mini       | 20M       | ~3 min         | 10GB   | Fast production runs |
| small      | 46M       | ~5 min         | 12GB   | **Default - best balance** |
| base       | 200M      | ~10 min        | 16GB   | Higher accuracy needed |
| large      | 710M      | ~25 min        | 32GB   | Maximum accuracy, research |

## GPU Acceleration (Optional)

To use GPUs for faster inference, modify the SLURM script:

```bash
#SBATCH --partition=gpu          # or gpu-long, gpu-preempt
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
```

Then edit `main.py` line 89 to change:
```python
device_map="cpu",  # Change to "cuda"
```

GPU runtime is typically 3-5x faster than CPU.

## Monitoring Jobs

```bash
# Check job status
squeue -u $USER

# Check specific job
squeue -j <job_id>

# View output (while running or after completion)
tail -f logs/chronos_<job_id>.out

# View all array job outputs
tail -f logs/chronos_array_<array_job_id>_*.out

# Cancel a job
scancel <job_id>

# Cancel all your jobs
scancel -u $USER

# Cancel array job
scancel <array_job_id>
```

## Output Location

Forecasts are saved to:
```
../../model-output/UMass-chronos_<model_size>/YYYY-MM-DD-UMass-chronos_<model_size>.csv
```

Example:
```
model-output/UMass-chronos_small/2025-01-04-UMass-chronos_small.csv
```

## Troubleshooting

### Virtual Environment Not Found
```bash
cd ~/forecast-sandbox-2025-2026/src/chronos
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Out of Memory
- Reduce `--mem` if job fails to start
- Use smaller model size (tiny/mini instead of base/large)
- Reduce `--num_samples` parameter

### Job Time Limit Exceeded
- Increase `--time` in SLURM script
- Use smaller model size
- Check if you have access to longer time limit partitions

### Module Not Found Errors
Check available Python modules:
```bash
module avail python
module load python/3.11  # or appropriate version
```

### First Run Takes Longer
The first run downloads model weights. Subsequent runs use cached weights and are much faster.

## Resource Requirements Summary

### Recommended Settings by Model Size

**Small Model (Default)**:
```bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=00:15:00
```

**Base Model**:
```bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:20:00
```

**Large Model**:
```bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:45:00
```

## Advanced Usage

### Running with Different Context Lengths

```bash
# Use 3 years of history instead of 2
sbatch --export=TODAY_DATE=2025-01-04,CONTEXT_LEN=156 run_chronos.slurm
```

### Running with More Samples (Smoother Quantiles)

```bash
# Increase from 200 to 500 samples
sbatch --export=TODAY_DATE=2025-01-04,NUM_SAMPLES=500 run_chronos.slurm
```

### Batch Processing Multiple Dates Sequentially

```bash
for date in 2024-12-28 2025-01-04 2025-01-11 2025-01-18; do
    sbatch --export=TODAY_DATE=$date run_chronos.slurm
    sleep 1  # Avoid overloading scheduler
done
```

## Getting Help

- **Unity Documentation**: https://docs.unity.rc.umass.edu/
- **Unity Support**: unity-support@umass.edu
- **Chronos Issues**: https://github.com/amazon-science/chronos-forecasting
- **Model Issues**: Thomas Robacker (trobacker@umass.edu)

## Partition Selection Guide

Unity likely has partitions like:
- `cpu`: Standard CPU jobs (default)
- `cpu-long`: Longer time limits
- `gpu`, `gpu-long`, `gpu-preempt`: GPU jobs

Check available partitions:
```bash
sinfo
scontrol show partition
```

Adjust the `--partition` directive in the SLURM scripts based on your needs and available allocations.
