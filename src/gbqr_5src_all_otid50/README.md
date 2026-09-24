# gbqr_5src_all_otid50

GBQR flu model variant exploring supplementary data sources and SMH scenario-hub trajectory
filtering.

- Main source: NHSN
- Supplementary sources: NSSP, ILINet, FluSurvNet, SMH
- SMH model filter: None (all available SMH models, sampled independently per model)
- SMH output_type_id sample size (`smh_num_otid`): 50 per SMH model
- SMH otid sampling seed (`smh_otid_seed`): 42

# To run locally without Docker

To test this out locally, run the following with this directory (`gbqr_5src_all_otid50`) as your
working directory.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --today_date=2024-01-06 --short_run
```

This should result in a model output file under `../../model-output/UMass-gbqr_5src_all_otid50/`.

# Generating forecasts for all reference dates

`run-all-forecasts.sh` runs this model locally, sequentially, once per `reference_date` in
`../../hub-config/tasks.json`. The dates in the script are run dates (the Wednesday before each
Saturday reference date, matching the production run schedule); `main.py` rolls each one forward
to the following Saturday internally. `submit-unity-parallel.sh` runs the same set of dates in
parallel as a Slurm job array on Unity.

# requirements.txt

`requirements.txt` was generated according to [README.md](../README.md).
