"""
Export the season-replay training rows used by relative-size-eda.qmd (the R document makes all figures and tables).

Rows are built with idmodels.peak.series.build_season_arrays / build_replay_rows, exactly as in the peak models, from
final ILINet (x percent positive) and FluSurv-NET data (the cached parquet written by validate_ilinet.py). NHSN is
deliberately not loaded: its seasons are the held-out test data. The ILINet series for Puerto Rico (72) and the US
Virgin Islands (78) are dropped before building rows (they are identically zero). Writes to analysis/peak-models/eda/:
  replay_rows.parquet     one row per (source, agg_level, location, season, season week t) with the model features
                          and the targets z and k (plus M = running max, eps, log_peak, at_zero, group)
  season_arrays.parquet   the season-aligned weekly values y (after the models' gap filling) and the running maximum
                          M_t for every week t, for every series

Usage (from the repository root, with an environment that has idmodels):
    DYLD_FALLBACK_LIBRARY_PATH=<venv>/lib/python3.12/site-packages/sklearn/.dylibs \
        python analysis/peak-models/relative_size_eda.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

from idmodels.peak.series import LOG_EPS, build_replay_rows, build_season_arrays, running_max

HERE = Path(__file__).resolve().parent
OUT = HERE / "eda"
ILI_PATH = HERE / "ilinet-validation" / "data-ilinet-flusurvnet.parquet"
DROP_ILINET_LOCATIONS = ["72", "78"]  # Puerto Rico, US Virgin Islands: all-zero series

W0, W1, REPLAY_START, MIN_OBS = 10, 43, 5, 25
ZERO_TOL = 1e-9
GROUPS = ["ILINet states", "ILINet national + regions", "FluSurv-NET"]
# columns that newer idmodels versions add to the replay rows; this export keeps the original model features only
EXTRA_MODEL_COLUMNS = ["sync_med_rel_max", "sync_med_g3", "sync_frac_past2", "sync_frac_half", "cum_vs_hist_total",
                       "cum_vs_hist_same_week"]


def group_of(source: pd.Series, agg_level: pd.Series) -> pd.Series:
    g = np.select(
        [(source == "ilinet") & (agg_level == "state"), source == "ilinet", source == "flusurvnet"],
        ["ILINet states", "ILINet national + regions", "FluSurv-NET"], "other")
    return pd.Series(g, index=source.index)


def load_rows(keep_extra=False):
    data = pd.read_parquet(ILI_PATH)
    data = data.loc[~((data["source"] == "ilinet") & data["location"].isin(DROP_ILINET_LOCATIONS))]
    arrays = build_season_arrays(data)
    rows = build_replay_rows(arrays, W0, W1, REPLAY_START, MIN_OBS)
    if not keep_extra:
        rows = rows.drop(columns=[c for c in EXTRA_MODEL_COLUMNS if c in rows.columns])
    rows["group"] = group_of(rows["source"], rows["agg_level"])
    rows["eps"] = rows["source"].map(LOG_EPS)
    rows["M"] = np.exp(rows["lm"]) - rows["eps"]
    rows["at_zero"] = rows["z"].abs() < ZERO_TOL
    rows["log_peak"] = rows["z"] + rows["lm"]
    return arrays, rows


def season_arrays_long(arrays) -> pd.DataFrame:
    """Long table of y and the running max M_t (as defined by idmodels.peak.series.running_max) for every week."""
    M = np.column_stack([running_max(arrays.y, t, W0)[0] for t in range(1, W1 + 1)])
    frames = []
    for j in range(W1):
        f = arrays.keys.copy()
        f["season_week"] = j + 1
        f["y"] = arrays.y[:, j]
        f["M"] = M[:, j]
        frames.append(f)
    out = pd.concat(frames, ignore_index=True)
    n_obs = np.sum(~np.isnan(arrays.y[:, W0 - 1:W1]), axis=1)
    complete = arrays.keys.assign(complete=n_obs >= MIN_OBS)
    return out.merge(complete, on=list(arrays.keys.columns))


def main():
    OUT.mkdir(exist_ok=True)
    arrays, rows = load_rows()
    rows.to_parquet(OUT / "replay_rows.parquet")
    season_arrays_long(arrays).to_parquet(OUT / "season_arrays.parquet")
    print(f"{len(rows)} replay rows, {rows.groupby(['source', 'agg_level', 'location', 'season']).ngroups} series")


if __name__ == "__main__":
    main()
