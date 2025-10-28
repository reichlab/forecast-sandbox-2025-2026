# Ensemble of baseline models with trends for influenza

The model outputs for the trends ensemble live in a directory called UMass-trends_ensemble. These contain quantile forecasts only.

# requirements.txt and renv.lock details

`requirements.txt` and `renv.lock` were generated according to [README.md](../README.md). For `renv.lock`, we installed these specific libraries:

```bash
Rscript -e "renv::install(c('readr', 'fs', 'lubridate'))"
Rscript -e "renv::install('arrow')"
Rscript -e "renv::install('hubverse-org/hubData@*release')"
Rscript -e "renv::install('hubverse-org/hubVis@*release')"
Rscript -e "renv::install('reichlab/trendsEnsemble')"
Rscript -e "renv::install('reichlab/idforecastutils')"  # NB: installs dev versions of above
```
