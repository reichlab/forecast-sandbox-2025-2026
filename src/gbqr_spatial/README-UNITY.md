# Running gbqr_spatial on Unity Cluster

This document describes how to run the gbqr_spatial forecasting model on the Unity HPC cluster using parallel array jobs.

## Prerequisites

1. **Unity Account**: Ensure you have an active Unity cluster account
2. **Python Environment**: Set up a virtual environment with all dependencies
3. **File Access**: Ensure the repository is accessible from Unity (clone or mount)

## Setup

### 1. Create logs directory

```bash
cd src/gbqr_spatial
mkdir -p logs
```

### 2. Set up Python virtual environment

```bash
# Create virtual environment
python -m venv .venv

# Activate it
source .venv/bin/activate

# Install dependencies
python -m pip install -r requirements.txt
```

### 3. Configure the submission script

Edit either `submit-unity-parallel.sh` or `submit-unity-throttled.sh` and update:

- `--mail-user=YOUR_EMAIL@umass.edu` - Your email address for notifications
- Adjust resource requests if needed:
  - `-c 4` - Number of CPU cores (increase if model needs more)
  - `--mem=16G` - Memory per job (increase if jobs run out of memory)
  - `-t 04:00:00` - Time limit (increase if jobs timeout)
- Adjust Python environment activation section (lines 72-81)

## Submission Scripts

### submit-unity-parallel.sh

Runs all 56 forecast dates as parallel array jobs with **no concurrency limit**.

- All 56 jobs will start simultaneously (subject to cluster availability)
- Fastest option if cluster resources are available
- Use this for urgent runs or when the cluster is not busy

**Submit with:**
```bash
sbatch submit-unity-parallel.sh
```

### submit-unity-throttled.sh

Runs all 56 forecast dates with a **maximum of 10 concurrent jobs** (`--array=0-55%10`).

- Only 10 jobs run at once; new jobs start as others complete
- More considerate of shared cluster resources
- Recommended for routine runs

**Submit with:**
```bash
sbatch submit-unity-throttled.sh
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

# View specific log
cat logs/slurm-JOBID_ARRAYINDEX.out

# View errors
cat logs/slurm-JOBID_ARRAYINDEX.err

# Tail running job
tail -f logs/slurm-JOBID_ARRAYINDEX.out
```

### Check array job summary
```bash
sacct -j JOBID --format=JobID,JobName,Partition,State,ExitCode,Elapsed
```

## Managing Jobs

### Cancel entire array job
```bash
scancel JOBID
```

### Cancel specific array task
```bash
scancel JOBID_ARRAYINDEX
# Example: scancel 123456_5
```

### Adjust concurrent job limit (after submission)
```bash
scontrol update jobid=JOBID arraytaskthrottle=5
```

## Array Job Details

- **Total dates**: 56 (2023-10-18 through 2025-05-14)
- **Array indices**: 0-55
- **Date selection**: Uses `$SLURM_ARRAY_TASK_ID` to index into dates array
- **Output files**: `logs/slurm-JOBID_ARRAYINDEX.out` (one per date)
- **Error files**: `logs/slurm-JOBID_ARRAYINDEX.err` (one per date)

## Resource Recommendations

Based on typical GBQR model runs:

| Resource | Recommended | Notes |
|----------|-------------|-------|
| Cores | 4-8 | Depends on `num_bags` and parallelization |
| Memory | 16-32G | Increase if model uses large datasets |
| Time | 2-4 hours | Increase for first runs or large models |
| Concurrent jobs | 10-20 | Balance speed vs. cluster fairness |

## Troubleshooting

### Jobs fail with memory errors
- Increase `--mem` in SBATCH header
- Check logs for "Out of Memory" or "OOM" messages

### Jobs timeout
- Increase `-t` time limit
- Consider using `--mail-type=TIME_LIMIT_80` to get warning at 80% time

### Environment activation fails
- Verify `.venv` exists in `src/gbqr_spatial/`
- Check that all dependencies are installed
- Try absolute paths in activation command

### Some jobs succeed, others fail
- Check individual error logs: `logs/slurm-*_ARRAYINDEX.err`
- Identify which dates failed: `grep -l "Error" logs/*.err`
- Resubmit specific failed indices: `sbatch --array=5,12,33 submit-unity-parallel.sh`

## Example Workflow

```bash
# 1. Set up environment (one-time setup)
cd src/gbqr_spatial
mkdir -p logs
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# 2. Edit submission script
nano submit-unity-throttled.sh
# Update email address and resource requests

# 3. Submit jobs
sbatch submit-unity-throttled.sh

# 4. Monitor progress
watch -n 30 'squeue -u $USER'

# 5. Check results
ls logs/
grep "Successfully completed" logs/*.out | wc -l

# 6. Check for failures
grep -l "Error" logs/*.err
```

## Additional Resources

- [Unity Documentation](https://docs.unity.rc.umass.edu/)
- [Slurm Array Jobs](https://docs.unity.rc.umass.edu/documentation/jobs/sbatch/arrays/)
- [Unity Support](mailto:hpc@umass.edu)
