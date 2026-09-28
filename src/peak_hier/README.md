# peak_hier

Bayesian hierarchical direct model (`idmodels.peak.PeakHierModel`, NumPyro): a discrete-time hazard for peak timing and a
truncated-normal model for the square root of the remaining log growth given the timing, with location and season
effects. The current season's effect is updated from right-censored observations of the season so far.

Forecasts the seasonal targets `peak week inc flu hosp` (pmf over the 34 window Saturdays) and
`peak inc flu hosp` (23 quantiles). Training data are NHSN, ILINet x WHO-NREVSS percent positive and
FluSurv-NET seasons before the forecast season; revisions to recent NHSN data are simulated from NHSN data
vintages. Methods are described in `analysis/peak-models/peak-models-methods.qmd`.

Development status: evaluated on ILINet development seasons only (`analysis/peak-models/validate_ilinet.py`). In
development a weaker likelihood tempering (`likelihood_weight` 0.02-0.05) scored better than the default 0.15.

## Running locally

From this directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --today_date=2026-01-07 --short_run
```

This writes `../../model-output/UMass-peak_hier/<reference_date>-UMass-peak_hier.csv`. A full fit (2 chains x 1000
iterations) takes roughly 15-20 minutes per season on a laptop.
