# gbqr_3src_gefs: GBQR with GEFS weather features (prototype)

An alternative to `../gbqr_3src_weather`. It uses real archived NOAA GEFS ensemble
forecasts instead of "observed weather plus noise". It also uses population-weighted
state averages instead of single city centroids.

## Data (`gefs_data.py`)

| Source | Use |
|---|---|
| dynamical.org GEFS analysis zarr (`noaa/gefs/analysis`, 2000–present, 0.25°, 3-hourly) | observed weather and the 2000–2019 climatology |
| dynamical.org GEFS 35-day forecast zarr (`noaa/gefs/forecast-35-day`, 00Z inits from 2020-10-01, 31 members) | weekly forecasts from Wednesday inits |
| Census 2020 county population centroids (`data/CenPop2020_Mean_CO*.txt`) | weights for grid cells to states; regions and US are population-weighted averages of states |

```bash
python gefs_data.py analysis    # -> data/gefs_analysis_daily.csv  (~1–2 h)
python gefs_data.py forecasts   # -> data/gefs_forecast_weekly.csv (212 inits, ~1 h)
```

**Variables.** The variables are 2 m temperature and precipitable water (PW).
PW stands in for absolute humidity because GEFS analysis has no RH before 2020,
so absolute humidity cannot be computed earlier than that. PW's weekly winter
anomalies correlate 0.88 with absolute humidity anomalies (per-state median 0.89).

## Features (`weather_features.py`)

All features are weekly anomalies against a 3-harmonic daily climatology.
Forecasts are bias-corrected per (location, lead, var) using inits before 2023-07-01.

- `obs`: the last two observed weeks (`wx_{temp,pw}_obs0/1`). 4 features.
- `obs_fc`: `obs` plus GEFS ensemble-mean anomalies for target weeks 1–3, taken from
  the Wednesday init before the reference date. 10 features. Training seasons before
  Oct 2020 have no archived forecast. For those, the forecast is the true anomaly plus
  a resampled real GEFS error vector, which keeps the correlation across locations,
  leads and variables. The synthetic error SD matches the real one (temp k=1/2/3:
  0.35/1.33/1.93 synthetic vs 0.33/1.29/1.96 real). Test weeks always use real forecasts.
- `obs_oracle`: `obs` plus the *true* future anomalies. This is an upper bound on what a
  perfect weather forecast could add.

Forecast skill after the bias-correction cutoff (Nov–Mar, anomaly correlation, k=1/2/3):
temperature 0.99 / 0.92 / 0.63, PW 0.99 / 0.82 / 0.35.

## Model and evaluation

`main.py` runs the same model as `gbqr_3src` (3 sources, 4th-root transform, states only)
with the weather columns added to `curr_feat_names`. It requires the pre-refactor
idmodels/iddata pins in `requirements.txt`.

```bash
uv venv && uv pip install -r requirements.txt   # then reinstall iddata with --no-deps
python main.py --today_date 2024-01-10 --variant obs_fc --short_q
VARIANTS="none obs obs_fc obs_oracle" ./run_eval.sh OUT 20 3
python score_eval.py OUT
```

## Results

The ablation used 26 reference dates (every other date of 2023/24 and 2024/25), 20 bags
and 7 quantiles. It covers 52 states/territories over horizons 0–3, with 5,404 matched
forecasts per variant. Relative WIS is shown against `none`, with a bootstrap 95% CI
over reference dates.

| variant | rel WIS | 95% CI | 2023/24 | 2024/25 |
|---|---|---|---|---|
| obs | 0.994 | 0.980–1.005 | 0.981 | 0.999 |
| obs_fc | 0.999 | 0.977–1.015 | 0.978 | 1.006 |
| obs_oracle | 1.001 | 0.978–1.020 | 0.980 | 1.008 |

By horizon, the only CI that excludes 1 is `obs` at horizon 1 (0.975, 0.949–0.997).
Coverage is unchanged (about 41% at 50%, about 88% at 95%).

**Interpretation.** Weather adds nothing detectable to GBQR at the state level for
1–4 week flu hospitalization forecasts. Even *perfect* knowledge of future weather does
not help, so better weather forecasts cannot fix this. Weather does get used: in a test
run it took about 25% of split importance. That is consistent with it fitting noise. The
contemporaneous incidence trend already carries most of the climate signal, and weekly
state-level anomalies are a weak driver next to immunity, strain and holiday mixing.

Possible directions if this is pursued further:
- Test season-onset timing specifically, at low incidence early in the season, where
  humidity effects are best documented, rather than overall WIS.
- Use longer-window features (2–4 week cumulative cold or dry spells), or interact
  weather with `season_week`.
- Use mechanistic models (SIRS forced by absolute humidity), where the covariate enters
  through transmission rather than as a tree split.
- Run 100 bags and 23 quantiles to tighten CIs. This is unlikely to change the conclusion
  given the oracle result.
