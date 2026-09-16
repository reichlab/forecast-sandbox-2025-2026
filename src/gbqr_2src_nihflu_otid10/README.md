# gbqr_2src_nihflu_otid10

GBQR flu model variant exploring supplementary data sources and SMH scenario-hub trajectory
filtering.

- Main source: NHSN
- Supplementary sources: SMH
- SMH model filter: NIH-FLU_TS
- SMH output_type_id sample size (`smh_num_otid`): 10
- SMH otid sampling seed (`smh_otid_seed`): 42

# To run locally without Docker

To test this out locally, run the following with this directory (`gbqr_2src_nihflu_otid10`) as your
working directory.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --today_date=2024-01-06 --short_run
```

This should result in a model output file under `../../model-output/UMass-gbqr_2src_nihflu_otid10/`.

# Generating forecasts for all reference dates

`run-all-forecasts.sh` runs this model locally, sequentially, once per `reference_date` in
`../../hub-config/tasks.json`. The dates in the script are run dates (the Wednesday before each
Saturday reference date, matching the production run schedule); `main.py` rolls each one forward
to the following Saturday internally. `submit-unity-parallel.sh` runs the same set of dates in
parallel as a Slurm job array on Unity.

# requirements.txt

`requirements.txt` was generated according to [README.md](../README.md).
