# sync_flusight

Keeps the data this sandbox borrows from [cdcepi/FluSight-forecast-hub](https://github.com/cdcepi/FluSight-forecast-hub) current during an ongoing season:

- `target-data/oracle-output.csv`: copied as-is from FluSight. The dashboard's eval page scores against this file, so new weeks get scored only after a refresh.
- `model-output/FluSight-baseline/`: FluSight-baseline forecasts filtered to what this hub accepts (quantile output, horizons 0–3, `wk inc flu hosp` and `wk inc flu prop ed visits`, and the hub's locations). The eval page uses this as the relative-skill baseline, so every round you want relative WIS for needs a baseline file.

Run it weekly from the hub root after FluSight's Wednesday target-data update, then open a PR with the changes:

```bash
python3 src/sync_flusight/sync_flusight.py
```

The script uses only the standard library. By default it fetches baseline files for reference dates from `2026-10-10` up to today that are listed in `hub-config/tasks.json` and missing locally. Use `--since`, `--overwrite`, `--skip-oracle` or `--skip-baseline` to change that.

The fan chart's observed-data JSON files live on the dashboard repo's `ptc/data` branch and are synced separately; see that repo's `src/sync_ptc_targets.py`.
