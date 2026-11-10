# axial-tfm_L4v5_r2

A revision of [`axial-tfm_L4v5`](../axial-tfm_L4v5/README.md): the same looped axial transformer (81k parameters)
and the same 15-network ensemble, retrained with three changes:
- missing and late newest weeks and dropouts are simulated in every pretraining and fine-tuning input, with a
  channel for how old each location's newest value is;
- week positions use rotary embeddings in place of the report-age channels;
- fine-tuning weights each location's loss by its size to the power 0.25.

The ensemble averages 5 seeds each of an admissions-only network and two NSSP-input networks. The NSSP networks are
weighted 0.75, 0.5, 0.5 and 0.25 at horizons 0 to 3 when NSSP is as current as NHSN; otherwise the admissions-only
networks are used alone. Newest NHSN weeks are corrected by a chain-ladder nowcast. The US median is set to the sum
of the state medians.

The model code, weights and full experiment log live in the axial repository
([nikhillebaka/axial](https://github.com/nikhillebaka/axial), currently private), as `src/axial`, where this model
is `ens_r2`. This folder is a thin wrapper, the way the gbqr models wrap `idmodels`. It was run at axial commit
`cd4b407a48efe85f4c571e84b9e961a04920f6de`.

# The name

`axial-tfm_L4v5_r2`: `axial-tfm_L4v5` (the axial transformer, looped 4 times, pretrained on simulator v5),
revision 2.

# Targets in `model-output/axial-tfm_L4v5_r2`

| target | reference dates | how |
|---|---|---|
| `wk inc flu hosp` | all 85 (2023/24 to 2025/26, and 2026/27 through 2026-11-07) | `ens_r2` quantiles |
| `wk inc flu prop ed visits` | 33: 2025-11-22 to 2026-05-30, and 2026/27 (the dates with a FluSight-baseline ED forecast here) | the ED head of the NSSP-input networks with an ED output; the US is rescaled to CDC's national series. Locations with no ED value in the current season are left out (Iowa in 2026/27) |
| `peak inc flu hosp` | 2025/26 and 2026/27 (over each season's peak-week bins in this hub) | 500 season paths, each built by feeding the admissions-only networks' one-week-ahead draws back in as data |
| `peak week inc flu hosp` | 2026/27 | the same paths |

2025/26 has no `peak week inc flu hosp`: the hub's peak-week bins are now the 2026/27 weeks, so 2025/26 peak-week
forecasts no longer validate. 2026/27 dates are added once FluSight-baseline has a forecast for them here.

Every forecast is a point-in-time backtest. NHSN and NSSP are as released before each reference date's
deadline (Wednesday 11:59 PM ET), and each season's networks are fine-tuned only on seasons that ended before
it. axial's leakage audit passes for each of these targets: on sampled dates, a rerun on data truncated at the
deadline reproduces the output exactly. These are not live forecasts. 2023/24 to 2025/26 were also used for design
choices in axial, so scores on them are optimistic. The 2026/27 forecasts were made after their deadlines (on
2026-11-09) with that season's networks, fine-tuned on data public by 2026-10-10. No design choice used 2026/27.

# To run locally

```bash
git clone https://github.com/nikhillebaka/axial ../../../axial
git -C ../../../axial checkout cd4b407a48efe85f4c571e84b9e961a04920f6de
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

This writes `../../model-output/axial-tfm_L4v5_r2/2025-12-13-axial-tfm_L4v5_r2.csv`. `run-all-forecasts.sh` runs
every reference date above. Admissions and ED take a few seconds per date. The 2025/26 and 2026/27 dates also run
the peak rollouts, a few minutes each.

# requirements.txt

`requirements.txt` is a copy of axial's `src/axial/requirements.lock` at the commit above.
