# peakGB

Operational candidate for the FluSight seasonal targets `peak week inc flu hosp` (pmf over the 34 window Saturdays)
and `peak inc flu hosp` (23 quantiles), submitted as **UMass-peakGB**.

GBQR direct peak model (`idmodels.peak.PeakGBQRModel`) with separate feature sets for its two parts:

- peak size: **core + SB** (synchrony + burden) features, as in `peak_gbqr_sb`
- peak timing: **core + holiday** features, as in `peak_gbqr_core_hol`

Everything else is the `PeakGBQRModelConfig` default: training on NHSN, ILINet and FluSurv-NET seasons before the
forecast season, 25 bags over seasons, NHSN revision simulation (200 draws), peak-week probabilities floored at 1e-4.
Size and timing are fitted independently, so the hindcast scores of this combination equal those of `peak_gbqr_sb`
for peak size and `peak_gbqr_core_hol` for peak timing (`model-output/`, 2023/24–2025/26). Methods:
`analysis/peak-models/peak-models-methods.qmd`.

The folder is laid out like the models in [reichlab/operational-models](https://github.com/reichlab/operational-models)
(e.g. `flu_flusion`) so it can be moved there as is.

## Running locally

From this directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --today_date=2026-10-07 --short_run
```

`main.py` fits the model and writes `output/model-output/UMass-peakGB/<reference_date>-UMass-peakGB.csv` (the reference
date is the Saturday on or after `--today_date`; `--short_run` uses 3 bags and 20 revision draws, for testing). It then
runs `Rscript plot.R <reference_date>`, which writes `output/plots/<reference_date>-UMass-peakGB.pdf`: one panel per
location, 6 per page, showing the NHSN data released on the Wednesday before the reference date, earlier seasons, the
peak-size forecast and the peak-week pmf (the figure function used in `analysis/peak-models/peak-models-slides.qmd`).
`output/` is gitignored in this repository.

A full run takes about 7 minutes on a laptop. On macOS, LightGBM needs libomp (`brew install libomp`, or set
`DYLD_FALLBACK_LIBRARY_PATH` to a directory that contains `libomp.dylib`).

## Requirements

- Python: `requirements.txt`.
- R (for `plot.R`): dplyr, tidyr, readr, ggplot2, patchwork, scales. In operational-models these go in a `renv.lock`
  generated as described in that repository's README, e.g.

```bash
Rscript -e "renv::install(c('dplyr', 'tidyr', 'readr', 'ggplot2', 'patchwork', 'scales'))"
```

`plot.R` downloads the NHSN release from the Reich Lab S3 bucket and location names from FluSight's
`auxiliary-data/locations.csv`.
