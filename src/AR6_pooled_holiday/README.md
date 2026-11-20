# AR6_pooled_holiday

AR(6) SARIX models for weekly flu hospital admissions (NHSN), the same as `AR6_pooled` (fourth root, scaled and
centred by location, AR coefficients shared across locations), with multiplicative holiday effects from sarix's
calendar effects ([sarix#17](https://github.com/reichlab/sarix/pull/17), not yet merged). One `main.py` produces
five models, chosen with `--variant`:

| variant   | model                          | holiday weeks                          | pooling                                        |
|-----------|--------------------------------|----------------------------------------|------------------------------------------------|
| `none`    | `UMass-AR6_pooled_cal_none`    | none: baseline, same code versions as the holiday variants | |
| `shared`  | `UMass-AR6_pooled_cal_shared`  | Christmas, New Year                    | one effect per week, shared by all locations   |
| `partial` | `UMass-AR6_pooled_cal_partial` | Christmas, New Year                    | per location, shrunk to a common effect        |
| `xm1`     | `UMass-AR6_pooled_cal_xm1`     | week before Christmas, Christmas, New Year | shared                                     |
| `xm1_ny1` | `UMass-AR6_pooled_cal_xm1_ny1` | week before Christmas, Christmas, New Year, week after New Year | shared                |

Holiday weeks are the Sunday-Saturday weeks containing Christmas and New Year's Day, and the weeks before and after
them, each with its own effect (see `VARIANTS` in `main.py` and `idmodels.calendar`). The effects multiply the rate
per 100k plus the floors and offsets that idmodels adds before the fourth root.

`cal_none` is there because `AR6_pooled`'s existing outputs were made with older idmodels and sarix versions, so
comparing the holiday models to it would mix holiday effects with version changes.

**Held-out seasons.** For now only the 2023/24 season is run (`run-all-forecasts-2023-2024.sh`). Do not run or look at
forecasts or data for 2024/25 or later until we decide to unblind them; see the README of the `sarix-flu-dev` repo.

**2023/24 is not out-of-sample for the choice of holiday weeks.** The holiday weeks in the `xm1` and `xm1_ny1`
variants were chosen partly from NHSN data for 2023/24 (see the comment above `HOLIDAYS` in `main.py`). Treat their
2023/24 scores as development results, and compare variants on the held-out seasons.

# Running locally

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --variant=partial --today_date=2023-12-20 --short_run
```

This writes `../../model-output/UMass-AR6_pooled_cal_partial/2023-12-23-UMass-AR6_pooled_cal_partial.csv`. A short
run's file has fewer quantiles than the hub requires, so don't commit it.

`run-all-forecasts-2023-2024.sh [variant ...]` runs all 2023/24 reference dates sequentially;
`submit-unity-parallel.sh` runs them as a Slurm job array on Unity (5 variants x 29 dates).
