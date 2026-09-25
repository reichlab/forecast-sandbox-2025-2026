# peak_kcde

Kernel conditional density estimation / analog model (`idmodels.peak.PeakKCDEModel`).

Forecasts the seasonal targets `peak week inc flu hosp` (pmf over the 34 window Saturdays) and
`peak inc flu hosp` (23 quantiles). Training data are NHSN, ILINet x WHO-NREVSS percent positive and
FluSurv-NET seasons before the forecast season; revisions to recent NHSN data are simulated from NHSN data
vintages. Methods are described in `analysis/peak-models/peak-models-methods.qmd`.

## Running locally

From this directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --today_date=2026-01-07 --short_run
```

This writes `../../model-output/UMass-peak_kcde/<reference_date>-UMass-peak_kcde.csv`. On macOS, LightGBM needs
OpenMP (`brew install libomp`).

Hindcasts for a whole season (all peak models, sharing data loads) can be run from the repository root with
`python analysis/peak-models/run_hindcasts.py --season 2025/26`.
