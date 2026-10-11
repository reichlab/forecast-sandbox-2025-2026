# gbqr_4src_spatial

GBQR flu model variant with NHSN as the main source, supplemented with NSSP, ILINet, and
FluSurv-NET as additional training data, with directional wave spatial features (all 8
directions) added.

- Main source: NHSN
- Supplementary sources: FluSurvNet, ILINet, NSSP

# To run locally without Docker

To test this out locally, run the following with this directory (`gbqr_4src_spatial`) as your
working directory.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --today_date=2025-11-19 --short_run
```

This should result in a model output file under `../../model-output/UMass-gbqr_4src_spatial/`.

# Generating forecasts for all reference dates

`run-all-forecasts.sh` runs this model locally, sequentially, once per `reference_date` in the
2025-26 season from `../../hub-config/tasks.json`. The dates in the script are run dates (the
Wednesday before each Saturday reference date, matching the production run schedule); `main.py`
rolls each one forward to the following Saturday internally. Dates start at 2025-11-19 because
NSSP, a supplementary source for this model, only supports as_of >= 2025-09-17, and that's also
when the 2025-26 season's reference dates begin. `submit-unity-parallel.sh` runs the same set of
dates in parallel as a Slurm job array on Unity, throttled to a handful of concurrent tasks since
the shared Unity "cpu" partition nodes see heavy multi-tenant contention that makes per-bag
training time highly variable.

# requirements.txt

`requirements.txt` was generated according to [README.md](../README.md).
