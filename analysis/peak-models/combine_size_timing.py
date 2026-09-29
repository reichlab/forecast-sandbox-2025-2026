"""
Combine the peak-size forecasts of one validation run with the peak-week forecasts of another (GBQR's size and timing
models are fit separately, so the combination equals a run with those size and timing feature sets).
Usage: python analysis/peak-models/combine_size_timing.py <name> <size run> <timing run>
"""
import sys
from pathlib import Path

import pandas as pd

F = Path(__file__).parent / "ilinet-validation" / "forecasts"
name, size_run, timing_run = sys.argv[1:4]
for season in ["2018-19", "2019-20"]:
    q = pd.read_parquet(F / f"{size_run}_{season}.parquet")
    p = pd.read_parquet(F / f"{timing_run}_{season}.parquet")
    out = pd.concat([q[q["output_type"] == "quantile"], p[p["output_type"] == "pmf"]], ignore_index=True)
    out.assign(model=name).to_parquet(F / f"{name}_{season}.parquet")
    print(name, season, len(out))
