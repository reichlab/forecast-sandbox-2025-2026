# AR6_fourierP_holiday_sl

`AR6_fourierP_thetaP_cal_xm1` (variant `fp_xm1` of `src/AR6_pooled_holiday`: AR(6), Fourier seasonality with K = 2
harmonics, and holiday effects for the week before Christmas week and the Christmas and New Year's weeks, all shared
across locations) with an innovation standard deviation that depends on the previous week's level. Each location's
innovation sd at week t is sigma * softplus(a + coef * x) / softplus(a), where x = y[t-1] minus the location's mean of
y, and y is the series as modeled (fourth root, centered and scaled; `inc_trans_cs`). coef has sarix's data-based
default prior and is shared by all locations; the offset a is learned and shared by all locations. The sd grows at
most linearly in x: at high levels it is about proportional to the distance of x from -a / coef. Two models, chosen
with `--variant`:

| variant           | model                                 | innovations                           |
|-------------------|---------------------------------------|---------------------------------------|
| `slsp_shared`     | `UMass-AR6_fP_tP_cal_xm1_slsp_shared` | normal, as in `fp_xm1`                |
| `slsp_tdf_shared` | `UMass-AR6_fP_tP_cal_xm1_slsp_tdf`    | Student's t, one df for all locations |

Compare them with `UMass-AR6_fourierP_thetaP_cal_xm1`, which differs only in the innovation sd (and, for
`slsp_tdf_shared`, its distribution), and with each other. The motivation is the residual check in `sarix-flu-dev`
(`sarix/holiday-effects/06-residuals.ipynb`): the residual spread of `fp_xm1` rises with the level on the modeled
scale. `slsp_tdf_shared` follows a second check (`07-residuals-slsp.ipynb`), where the residuals of `slsp_shared` fit
a t with about 7 df; in the fits, though, the df went to its lower bound of 2, and its scores are worse than
`slsp_shared`'s.

An earlier version used sd = sigma * exp(coef * x). The feedback of simulated values on the sd made some forecast
paths overflow near the 2024/25 peak (absurd upper quantiles in small locations, and no forecast at all for
2025-01-18), so it was dropped, and sarix no longer offers it.

**Exploratory.** These models were run on 2023/24, 2024/25 and 2025/26 after those seasons' data had been seen, so no
season is held out for them. Do not run or look at forecasts or data for 2026/27 until we decide to unblind it; see the
README of the `sarix-flu-dev` repo.

**Development versions.** The level-dependent sd is in a sarix branch (`elr/sigma-level`, reichlab/sarix PR 19, stacked on PRs 18
and 17) and an idmodels branch (`sigma-level`, reichlab/idmodels PR 49, stacked on PRs 48 and 47), neither merged
yet. `requirements.txt` pins the idmodels branch by commit, and that commit pins the sarix branch.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --variant=slsp_shared --today_date=2023-12-20 --short_run
```

`run-all-forecasts.sh [variant ...]` runs all 80 reference dates of the three seasons sequentially.
