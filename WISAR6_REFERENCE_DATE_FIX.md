# WISAR6: fix look-ahead data leak (reference-date snapshot selection)

Instructions for Claude on Unity. Apply the fix below to **every** WISAR6 variant's
`main.py`: `src/WISAR6`, `src/WISAR6_fourthroot`, `src/WISAR6_fourthroot_adaptive`,
`src/WISAR6_fourthroot_adaptive_diff`, `src/WISAR6_fourthroot_adaptive_t`. Ask before
regenerating any outputs.

## The bug

`_find_archive_snapshot(as_of)` picks the newest
`FluSight-forecast-hub/auxiliary-data/target-data-archive/target-hospital-admissions_<DATE>.csv`
with `<DATE> <= reference_date`. But the hub's `archive_target_data.R` names each file
after the **last week of data it contains** (`max(target_data$date)`), not after the date
it was published. The file named for week W is committed the Wednesday *after* W. So for
reference date R, the code loaded the file named R, which was published about 4 days
after R and contains revisions that weren't available at R.

Examples (value of the last training week):

| ref date | location | correct (file named R-7) | what was used (file named R) |
|---|---|---|---|
| 2025-04-05 | NH (33), week 3/29 | 124 | 34 |
| 2026-03-21 | ID (16), week 3/14 | 16 | 75 |

70 of the 80 `UMass-WISAR6_fourthroot_adaptive_t` reference dates were affected.

## The fix (3 parts, same in each `main.py`)

1. **`_find_archive_snapshot`**: use only snapshots named strictly before R
   (`p.stem.split("_")[-1] < as_of.isoformat()`). Use the live `TARGET_DATA_PATH` only
   when `as_of >= datetime.date.today()` (an operational run). For a historical run with
   no usable snapshot, `raise FileNotFoundError` instead of silently falling back to the
   live (revised) file. Update the docstring to explain the naming convention.

2. **Handle archive gaps.** The archive skips some weeks, so the newest usable file can
   end at R-14 or earlier. `simulate_forecast_draws` steps forward from the last training
   date, so without a correction, horizon 0 would be mislabeled. Right after the
   `if train_df.empty: raise ...` block, add:
   ```python
   last_train_date = train_df["date"].max()
   n_gap_steps = (ref_ts - last_train_date).days // 7 - 1
   if n_gap_steps > 0:
       click.echo(f"Note: latest data available before {reference_date} ends "
                  f"{last_train_date.date()}; simulating {n_gap_steps} extra "
                  "leading step(s) to reach horizon 0.")
   ```
3. In the forecast loop, change `h_max=MAX_HORIZON` to `h_max=n_gap_steps + MAX_HORIZON`
   in the `simulate_forecast_draws(...)` call, and change `draws_z = sims[:, h]` to
   `draws_z = sims[:, n_gap_steps + h]`. In `_adaptive_diff`, keep the
   `undifference_forecast` step on the full simulated path *before* indexing. It is a
   cumulative sum from `z_last`, so the leading gap steps have to be included.

Leave the training filter (`data["date"] < ref_ts`) and everything else unchanged.

## Verify (no MCMC needed)

First `git -C ../FluSight-forecast-hub pull` so the Unity hub clone's archive is current.
Then, from a variant's directory:
```bash
uv run python -c "
import datetime as dt, pandas as pd, main
for ref, loc in [('2025-04-05','33'),('2026-03-21','16'),('2026-03-14','16')]:
    r=dt.date.fromisoformat(ref); p=main._find_archive_snapshot(r); d=main.load_target_data(r)
    t=d[(d.location==loc)&(d.date<pd.Timestamp(r))].sort_values('date')
    print(ref, p.name, t.date.max().date(), t.value.iloc[-1], (pd.Timestamp(r)-t.date.max()).days//7-1)
"
```
Expected output:
- `2025-04-05 ..._2025-03-29.csv 2025-03-29 124.0 0`
- `2026-03-21 ..._2026-03-14.csv 2026-03-14 16.0 0`
- `2026-03-14 ..._2026-02-28.csv 2026-02-28 91.0 1` (archive gap → 1 extra step)

## After the fix

- **Regenerate all WISAR6 model-output files.** Every existing WISAR6 output was produced
  with the leak.
- **Rescore.** Earlier WISAR6 scores and writeups (for example `ENSEMBLE_EVALUATION.md`,
  `LOG_WIS_EVALUATION_FINDING.md`, and the `scores_*.csv` files) are optimistic.
- Reference dates where the archive gap applies (expect the "Note:" line):
  2024-11-30, 2024-12-14, 2025-01-11, 2025-11-29, 2025-12-27, 2026-02-07, 2026-03-14,
  2026-05-02 and 2026-05-30.
