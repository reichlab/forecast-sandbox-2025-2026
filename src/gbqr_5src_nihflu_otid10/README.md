# gbqr_5src_nihflu_otid10

GBQR flu model variant exploring supplementary data sources and SMH scenario-hub trajectory
filtering.

- Main source: NHSN
- Supplementary sources: NSSP, ILINet, FluSurvNet, SMH
- SMH model filter: NIH-Flu_TS
- SMH output_type_id sample size (`smh_num_otid`): 10
- SMH otid sampling seed (`smh_otid_seed`): 42

# To run locally without Docker

To test this out locally, run the following with this directory (`gbqr_5src_nihflu_otid10`) as your
working directory.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --today_date=2024-01-06 --short_run
```

This should result in a model output file under `../../model-output/UMass-gbqr_5src_nihflu_otid10/`.

# Generating forecasts for all reference dates

`run-all-forecasts.sh` runs this model locally, sequentially, once per `reference_date` in
`../../hub-config/tasks.json`. The dates in the script are run dates (the Wednesday before each
Saturday reference date, matching the production run schedule); `main.py` rolls each one forward
to the following Saturday internally. `submit-unity-parallel.sh` runs the same set of dates in
parallel as a Slurm job array on Unity.

`run-all-forecasts-2023-2025.sh` / `submit-unity-backfill-2023-2025.sh` run the same model over
the 2023-24 and 2024-25 seasons, now that NSSP as a supplementary source supports `as_of` dates
before 2025-09-17. NSSP has no versioned data that far back, so for those dates it falls back to
static/current NSSP data filtered to the historical date rather than true real-time vintage. SMH
trajectory data availability also varies by `as_of`: round 4 only before 2024-08-11 (covers the
2023-24 dates), rounds 4+5 from 2024-08-11 through 2025-08-10 (covers the 2024-25 dates) — both
lighter than the current season's 3-round tier, so the existing 48G/3h budget should be a
conservative upper bound.

# requirements.txt

`requirements.txt` was generated according to [README.md](../README.md).
