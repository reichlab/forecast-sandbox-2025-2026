# AR6_fourierP_holiday_t

`AR6_fourierP_thetaP_cal_xm1` (variant `fp_xm1` of `src/AR6_pooled_holiday`: AR(6), Fourier seasonality with K = 2
harmonics, and holiday effects for the week before Christmas week and the Christmas and New Year's weeks, all shared
across locations) with Student's t innovations instead of normal ones. The t is standardized so that sigma is still
the innovation standard deviation, and its degrees of freedom are above 2, with sarix's default prior (Gamma(2, 0.1)
truncated to values above 2). Two models, chosen with `--variant`:

| variant       | model                                 | degrees of freedom                                       |
|---------------|---------------------------------------|----------------------------------------------------------|
| `tdf_shared`  | `UMass-AR6_fP_tP_cal_xm1_tdf_shared`  | one value shared by all locations                        |
| `tdf_partial` | `UMass-AR6_fP_tP_cal_xm1_tdf_partial` | one per location, partially pooled toward a common value |

Compare them with `UMass-AR6_fourierP_thetaP_cal_xm1`, which differs only in the innovation distribution.

**Exploratory.** These models were run on 2023/24, 2024/25 and 2025/26 after those seasons' data had been seen, so
no season is held out for them. Do not run or look at forecasts or data for 2026/27 until we decide to unblind it;
see the README of the `sarix-flu-dev` repo.

**Development versions.** t innovations are in a sarix branch (`elr/tdist-innovations`, reichlab/sarix PR 18, stacked on PR 17) and an
idmodels branch (`t-innovations`, reichlab/idmodels PR 48, stacked on PR 47), neither merged yet. `requirements.txt`
pins the idmodels branch by commit, and that commit pins the sarix branch.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --variant=tdf_shared --today_date=2023-12-20 --short_run
```

`run-all-forecasts.sh [variant ...]` runs all 80 reference dates of the three seasons sequentially.
