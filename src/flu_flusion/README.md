# Flusion model for influenza

# To run locally without Docker

To test this out locally, run the following with `flu_flusion` as your working directory.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --today_date=2024-01-06 --short_run
```

Alternative instructions using `uv` in your working directory:

```bash
uv venv
uv pip install -r requirements.txt
uv run main.py --today_date=2025-04-05 --short_run
```

This should result in a model output file `model-output/UMass-flusion/`.

# requirements.txt and renv.lock details

`requirements.txt` and `renv.lock` were generated according to [README.md](../README.md). For `renv.lock`, we installed these specific libraries:

```bash
Rscript -e "renv::install(c('dplyr'))"
Rscript -e "renv::install('arrow')"
Rscript -e "renv::install('reichlab/zoltr')"
Rscript -e "renv::install('hubverse-org/hubData@*release')"
Rscript -e "renv::install('hubverse-org/hubVis@*release')"
Rscript -e "renv::install('hubverse-org/hubEnsembles@*release')"
Rscript -e "renv::install('reichlab/covidData')"
Rscript -e "renv::install('reichlab/idforecastutils')"  # NB: installs dev versions of above
```

Note: On Mac, one issue we ran into was a OS library dependency in running `01_gbqr.py` in which `libomp` was not found. 
If you encouter this, you may need to install `libomp`, one such way is:

```bash
brew install libomp
brew info libomp
```