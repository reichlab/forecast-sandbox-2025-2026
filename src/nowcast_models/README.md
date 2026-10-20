# Nowcast comparison models

Retrospective 2025/26 forecasts that test whether nowcasting the most recent NHSN week before fitting improves GBQR
and SARIX forecasts (idmodels issue #34, branch `nowcast-preprocessing`).

Four base models, with the settings of the sandbox models of the same names, each run three ways on the same idmodels
commit, so that the nowcast is the only difference:

| Suffix | Nowcast |
|---|---|
| `_nc_none` | none (the control) |
| `_nc_scale` | the newest week divided by its completeness: hospitals reporting that week / the usual number (median of the 2-4 weeks before) |
| `_nc_state` | scaled as above, then corrected by a factor for the location, learned from this season's settled weeks |

Base models: `gbqr_1src` (NHSN), `gbqr_4src` (NHSN, NSSP, ILINet, FluSurv-NET), `AR6_pooled` and `AR2` (SARIX).
GBQR bags are seeded by the reference date, so the three variants of a GBQR model use the same bags.

The nowcasts read the number of hospitals reporting, and for `_nc_state` the history of admissions revisions, from
the Delphi Epidata archive, as published by the Wednesday release the forecast's data come from.

## Running

```bash
python -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python fetch_archives.py archives        # save the Delphi archives once
python main.py --model gbqr_1src_nc_state --today_date 2026-01-14 --archive_dir archives
sbatch submit-unity.sh                              # all 12 models x 28 dates on Unity
```

On a Mac without libomp, LightGBM needs
`DYLD_FALLBACK_LIBRARY_PATH=$PWD/.venv/lib/python3.12/site-packages/sklearn/.dylibs`.
