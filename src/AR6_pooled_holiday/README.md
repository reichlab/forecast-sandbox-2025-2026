# AR6_pooled_holiday

AR(6) SARIX models for weekly flu hospital admissions (NHSN), the same as `AR6_pooled` (fourth root, scaled and
centred by location, AR coefficients shared across locations), with multiplicative holiday effects from sarix's
calendar effects ([sarix#17](https://github.com/reichlab/sarix/pull/17), not yet merged). One `main.py` produces
seven models, chosen with `--variant`:

| variant   | model                          | holiday weeks                          | pooling                                        |
|-----------|--------------------------------|----------------------------------------|------------------------------------------------|
| `none`    | `UMass-AR6_pooled_cal_none`    | none: baseline, same code versions as the holiday variants | |
| `shared`  | `UMass-AR6_pooled_cal_shared`  | Christmas, New Year                    | one effect per week, shared by all locations   |
| `partial` | `UMass-AR6_pooled_cal_partial` | Christmas, New Year                    | per location, shrunk to a common effect        |
| `xm1`     | `UMass-AR6_pooled_cal_xm1`     | week before Christmas, Christmas, New Year | shared                                     |
| `xm1_ny1` | `UMass-AR6_pooled_cal_xm1_ny1` | week before Christmas, Christmas, New Year, week after New Year | shared                |
| `fp_none` | `UMass-AR6_fourierP_thetaP_cal_none` | none: baseline for `fp_xm1` | |
| `fp_xm1`  | `UMass-AR6_fourierP_thetaP_cal_xm1`  | as `xm1` | shared |

The `fp_` variants add `AR6_fourierP_thetaP`'s Fourier seasonality (K = 2 harmonics, shared across locations). They
were run after `xm1` was chosen, on all three seasons, so their scores are exploratory: no season is held out for them.

Holiday weeks are the Sunday-Saturday weeks containing Christmas and New Year's Day, and the weeks before and after
them, each with its own effect (see `VARIANTS` in `main.py` and `idmodels.calendar`). The effects multiply the rate
per 100k plus the floors and offsets that idmodels adds before the fourth root.

`cal_none` is there because `AR6_pooled`'s existing outputs were made with older idmodels and sarix versions, so
comparing the holiday models to it would mix holiday effects with version changes.

**Selected variant: `xm1`** (2026-11-21), after scoring 2023/24 and the held-out seasons 2024/25 and 2025/26. Pooled
over the three seasons, its relative WIS against `cal_none` is 0.849 (counts) and 0.932 (log scale), against 0.864
and 0.947 for `shared`. The scores and the reasons are in notebook `05-scores-by-season.ipynb` of the `sarix-flu-dev`
repo.

**Held-out seasons.** The 2023/24 season (`run-all-forecasts-2023-2024.sh`) and, since 2026-11-20, 2024/25
(`run-all-forecasts-2024-2025.sh`, variants `none`, `shared`, `xm1` and `xm1_ny1`) and 2025/26
(`run-all-forecasts-2025-2026.sh`, variants `none`, `shared` and `xm1`) are run. Do not run or look at forecasts or data
for 2026/27 until we decide to unblind it; see the README of the `sarix-flu-dev` repo.

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

`run-all-forecasts-2023-2024.sh [variant ...]` runs all 2023/24 reference dates sequentially, and
`run-all-forecasts-2024-2025.sh [variant ...]` the 2024/25 ones; `submit-unity-parallel.sh` runs the 2023/24 dates as
a Slurm job array on Unity (5 variants x 29 dates).
