# AR_GBQR_ensemble_nssp

A single parameterized R script (`main.R`) that builds a simple (median) two-model ensemble from
one already-run AR/SARIX model and one already-run GBQR model, both NSSP-based. It does not fit
any new model: it reads the two component models' existing `model-output/UMass-<model>/` CSVs for
a given reference date and combines their quantiles with `hubEnsembles::simple_ensemble()`.

`run-all-forecasts.sh` runs this for all 15 pairs formed by crossing the top 3 AR models with the
top 5 GBQR models (by relative WIS against `FluSight-baseline`):

- AR: `AR4_fourierP_thetaP_nssp`, `AR2_fourierP_thetaP_nssp`, `AR6_fourierP_thetaP_nssp`
- GBQR: `gbqr_2src_flusurvnet_nssp`, `gbqr_4src_nssp`, `gbqr_4src_spatial_nssp`,
  `gbqr_2src_ilinet_nssp`, `gbqr_2src_nhsn_nssp`

Each pair gets its own ensemble model name and its own `model-output/UMass-<name>/` directory,
named by joining the two component names (minus their shared `_nssp` suffix) with `_ens_nssp`,
e.g. `AR4_fourierP_thetaP_gbqr_2src_flusurvnet_ens_nssp`.

# To run locally

This requires both component models' forecasts to already exist for the reference date in
question. With this directory (`AR_GBQR_ensemble_nssp`) as your working directory:

```bash
Rscript main.R --ar_model=AR2_nssp --gbqr_model=gbqr_1src_nssp --today_date=2025-11-19
```

This should result in a model output file under
`../../model-output/UMass-AR2_gbqr_1src_ens_nssp/`.

# Generating forecasts for all reference dates and pairs

`run-all-forecasts.sh` loops over all 15 (AR, GBQR) pairs and every `reference_date` in the
2025-26 season (375 runs total). Since this only combines already-computed CSVs (no model
fitting), it's fast and intended to run locally; there's no Unity submission script for it.

# R package dependencies

`dplyr`, `readr`, `lubridate`, `hubEnsembles` (all already installed globally on this machine; see
the root `CLAUDE.md`/environment notes).
