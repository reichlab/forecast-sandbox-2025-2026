# peak_hybrid

Hybrid of the gradient-boosted and hierarchical direct models (`idmodels.peak.PeakHybridModel`): the hierarchical model
(`peak_hier`) with `peak_gbqr`'s predictions as offsets for the already-peaked probability, the timing hazard and the
peak size. The hierarchical part adds season and location effects, including the update of the current season's effect
from the season observed so far. GBQR's offsets for the training rows are its out-of-bag predictions, so they are as
accurate as its real forecasts.

Forecasts the seasonal targets `peak week inc flu hosp` (pmf over the 34 window Saturdays) and
`peak inc flu hosp` (23 quantiles). Methods are described in `analysis/peak-models/peak-models-methods.qmd` (section
"Hybrid: GBQR offsets").

Development status: **paused** (idmodels 1dd468d). Scored on the ILINet development season 2018/19 only; the 2019/20
runs and the `current_update_weight` variants (`hybrid__cu01`, `hybrid__cu003`) are incomplete and on hold while the
GBQR model is locked in.

## Running locally

From this directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --today_date=2026-01-07 --short_run
```

This writes `../../model-output/UMass-peak_hybrid/<reference_date>-UMass-peak_hybrid.csv`. A full fit (GBQR bags plus
2 chains x 1000 MCMC iterations) takes roughly 30-40 minutes per season on a laptop.
