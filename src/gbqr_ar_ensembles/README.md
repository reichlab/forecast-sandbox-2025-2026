# GBQR x AR ensembles

Equally weighted quantile averages (`hubEnsembles::simple_ensemble`) of one GBQR model and one AR model, built from the
components' existing files in `model-output/`. No model fitting happens here; building all pairs takes about a minute.

| Ensemble | GBQR component | AR component | Seasons |
|---|---|---|---|
| UMass-gbqr_4src_ar6p | UMass-gbqr_4src | UMass-AR6_pooled | 2025/26 |
| UMass-gbqr_4src_ar6fp | UMass-gbqr_4src | UMass-AR6_fourierP_thetaP | 2025/26 |
| UMass-gbqr_4src_wisar | UMass-gbqr_4src | UMass-WISAR6_fourthroot_adaptive_t | 2025/26 |
| UMass-gbqr_nih10_ar6p | UMass-gbqr_5src_nihflu_otid10 | UMass-AR6_pooled | 2025/26 |
| UMass-gbqr_nih10_ar6fp | UMass-gbqr_5src_nihflu_otid10 | UMass-AR6_fourierP_thetaP | 2025/26 |
| UMass-gbqr_nih10_wisar | UMass-gbqr_5src_nihflu_otid10 | UMass-WISAR6_fourthroot_adaptive_t | 2025/26 |
| UMass-gbqr_3src_ar6p | UMass-gbqr_3src | UMass-AR6_pooled | 2023/24 – 2025/26 |
| UMass-gbqr_3src_ar6fp | UMass-gbqr_3src | UMass-AR6_fourierP_thetaP | 2023/24 – 2025/26 |
| UMass-gbqr_3src_wisar | UMass-gbqr_3src | UMass-WISAR6_fourthroot_adaptive_t | 2023/24 – 2025/26 |
| UMass-gbqr_3src_spatial_wisar | UMass-gbqr_3src_spatial (US: UMass-gbqr_3src) | UMass-WISAR6_fourthroot_adaptive_t | 2023/24 – 2025/26 |

The two other spatial GBQR x AR pairs already exist as UMass-flusion_3src_spatial (with UMass-AR6_pooled) and
UMass-flusion_3src_spatial_fourierP (with UMass-AR6_fourierP_thetaP); they were built elsewhere, but their files equal
`ensemble_pair.R <model> UMass-gbqr_3src_spatial <ar model> UMass-gbqr_3src` exactly (checked with validate.py).
`gbqr_3src_spatial` has no US forecasts, since its directional wave features are state-level only.

An ensemble covers every reference date its two components share: the `gbqr_4src` and `nih10` GBQR components use
NSSP or NIH flu scenario data available only for 2025/26, while `gbqr_3src` and the AR models cover all three seasons.

## Files

- `ensemble_pair.R <output_model_id> <gbqr_model_id> <ar_model_id> [us_gbqr_model_id]`: builds one ensemble for every
  shared reference date (horizons 0–3); the optional fourth argument uses a different GBQR model for the US forecast.
- `run_all.sh`: builds all the pairs above (written for the Unity cluster: `module load r/4.4.0` and the renv library of
  `src/flusion_spatial2_prod`; elsewhere, run the `Rscript ensemble_pair.R ...` lines with dplyr, hubUtils and
  hubEnsembles installed).
- `validate.py "<ensemble> <gbqr> <ar> [us_gbqr]" ...`: checks each ensemble file against its components (value = mean of
  the two), for missing or negative values, monotone quantiles, and reference dates / locations / quantile levels in
  `hub-config/tasks.json`.

Run all three from this directory.
