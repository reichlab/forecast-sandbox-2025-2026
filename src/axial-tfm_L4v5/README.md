# axial-tfm_L4v5

A simulation-pretrained axial transformer for flu hospital admissions. A small looped transformer
(81k parameters) that attends along each location's weeks and then across locations within each week is
pretrained on simulated two-strain SEIR seasons. It is then fine-tuned once per season on ILI+, FluSurv-NET
and the NHSN seasons completed before that season. The production model, `ens_fresh_hwr`, averages 15
networks: 5 seeds without NSSP and 10 with NSSP ED-visit inputs. The NSSP networks are weighted 0.75, 0.5,
0.5 and 0.25 at horizons 0 to 3 when NSSP is as current as NHSN; otherwise the 5 seeds without NSSP are used
alone. Newest NHSN weeks are corrected by a chain-ladder nowcast. The US median is set to the sum of the
state medians.

The model code, weights and full experiment log live in the axial repository
([nikhillebaka/axial](https://github.com/nikhillebaka/axial), currently private), as `src/axial`. This folder
is a thin wrapper, the way the gbqr models wrap `idmodels`. It was run at axial commit
`e0ea8651d29482c1e3be8dc9d3115838ebbdec36`.

# The name

`axial-tfm_L4v5`: the axial transformer (`tfm`), looped (`L4`: one axial block applied 4 times with shared
weights), pretrained on simulator v5 (`v5`: two-strain SEIR with summer circulation, wider drift and an NSSP-like
ED stream). In axial itself this model is `ens_fresh_hwr`: 15 networks fine-tuned at 70% real data, NSSP
networks weighted by horizon, nowcast inputs and a coherent national forecast.

# Targets in `model-output/axial-tfm_L4v5`

| target | reference dates | how |
|---|---|---|
| `wk inc flu hosp` | all 80 (2023/24 to 2025/26) | `ens_fresh_hwr` quantiles |
| `wk inc flu prop ed visits` | 53, from late April 2024 (point-in-time NSSP) | the ED head of the NSSP-input networks; the US is rescaled to CDC's national series |
| `peak inc flu hosp`, `peak week inc flu hosp` | 2025/26 (the hub's peak-week bins) | 500 season paths, each built by feeding the base networks' one-week-ahead draws back in as data |

Every forecast is a point-in-time backtest. NHSN and NSSP are as released before each reference date's
deadline (Wednesday 11:59 PM ET), and each season's networks are fine-tuned only on seasons that ended before
it. axial's leakage audit passes for each of these targets: on sampled dates, a rerun on data truncated at the deadline reproduces the output exactly. These are not live forecasts. 2023/24 to 2025/26 were
also used for design choices in axial, so scores on them are optimistic. The single held-out look at 2025/26
was taken with an earlier version (round 5, axial tag `round5-frozen`).

# To run locally

```bash
git clone https://github.com/nikhillebaka/axial ../../../axial
git -C ../../../axial checkout e0ea8651d29482c1e3be8dc9d3115838ebbdec36
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt   # axial's pinned environment (Python 3.13)
```

Then build axial's point-in-time data stores (gitignored there; see axial's `src/axial/README.md`):

```bash
cd ../../../axial
PYTHONPATH=src python -m axial.data.fetch
PYTHONPATH=src python -m axial.data.vintages
PYTHONPATH=src python -m axial.data.history
PYTHONPATH=src python -m axial.data.ilinet
PYTHONPATH=src python -m axial.data.fetch nssp delphi_hist
PYTHONPATH=src python -m axial.data.nssp
```

From this directory:

```bash
PYTHONPATH=../../../axial/src python main.py --today_date=2025-12-10
```

This writes `../../model-output/axial-tfm_L4v5/2025-12-13-axial-tfm_L4v5.csv`. `run-all-forecasts.sh` runs every
reference date in `hub-config/tasks.json`. Admissions and ED take a few seconds per date. The 2025/26 dates
also run the peak rollouts, about 4 minutes each.

# requirements.txt

`requirements.txt` is a copy of axial's `src/axial/requirements.lock` at the commit above.
