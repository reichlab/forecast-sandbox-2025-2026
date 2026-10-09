# gbqr_4src_spatial_nssp

GBQR flu model variant for predicting ED visit proportion. Main source is NSSP, supplemented with
all other standard surveillance sources (NHSN, ILINet, FluSurv-NET) as additional training data,
with directional wave spatial features (all 8 directions) added.

# To run locally without Docker

To test this out locally, run the following with this directory (`gbqr_4src_spatial_nssp`) as your
working directory.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --today_date=2025-11-19 --short_run
```

This should result in a model output file under `../../model-output/UMass-gbqr_4src_spatial_nssp/`.

# Generating forecasts for all reference dates

`run-all-forecasts.sh` runs this model locally, sequentially, once per `reference_date` in the
2025-26 season. The dates in the script are run dates (the Wednesday before each Saturday
reference date, matching the production run schedule); `main.py` rolls each one forward to the
following Saturday internally. `submit-unity-parallel.sh` runs the same set of dates in parallel
as a Slurm job array on Unity.

# requirements.txt

`requirements.txt` was generated according to [README.md](../README.md).
