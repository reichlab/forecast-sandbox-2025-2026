"""
Package analysis/peak-models/peak-scores.csv as scores.json for the peak score explorer
(analysis/peak-models/viz/scores.html).

Usage (from the repository root, after running score_peak_forecasts.py):
    python analysis/peak-models/viz/build_score_data.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
ROOT = HERE.parents[2]
LABELS = {"UMass-peak_gbqr": "GBQR", "UMass-peak_gbqr_offset": "GBQR (offset)", "UMass-peak_kcde": "KCDE",
          "UMass-peak_baseline": "Baseline"}
ORDER = ["UMass-peak_gbqr", "UMass-peak_gbqr_offset", "UMass-peak_kcde", "UMass-peak_baseline", "FluSight-ensemble",
         "NAU-vulPES", "UGA_flucast-Copycat", "PSI-PROF", "CU-ensemble", "FluSight-base_seasonal"]
DEFAULT_ON = {"UMass-peak_gbqr", "UMass-peak_gbqr_offset", "UMass-peak_baseline", "FluSight-ensemble"}
MIN_ROWS = 500  # drop models with only a handful of forecasts
METRICS = ["log_score", "rps", "prob_pm1", "wis_log", "wis", "cov50", "cov95"]
LOG_FLOOR = np.log(1e-4)


def main():
    s = pd.read_csv(HERE.parent / "peak-scores.csv", dtype={"location": str})
    counts = s["model"].value_counts()
    s = s.loc[s["model"].isin(counts[counts >= MIN_ROWS].index)].copy()
    # floor log scores at log(1e-4) so that zero-probability forecasts do not dominate means
    s["log_score"] = s["log_score"].clip(lower=LOG_FLOOR)
    models = [m for m in ORDER if m in set(s["model"])] + sorted(set(s["model"]) - set(ORDER))
    locs = pd.read_csv(ROOT / "auxiliary-data" / "locations.csv", dtype=str)
    rows = []
    for r in s.itertuples():
        vals = [None if pd.isna(getattr(r, m)) else round(float(getattr(r, m)), 5) for m in METRICS]
        rows.append([models.index(r.model), r.season, r.location, int(r.weeks_from_peak)] + vals)
    out = {
        "models": [{"id": m, "label": LABELS.get(m, m), "sandbox": m.startswith("UMass-"), "on": m in DEFAULT_ON}
                   for m in models],
        "metrics": METRICS,
        "columns": ["model", "season", "location", "wfp"] + METRICS,
        "rows": rows,
        "locations": [{"id": r.location, "name": r.location_name} for r in locs.itertuples()],
        "log_floor": LOG_FLOOR,
    }
    (HERE / "scores.json").write_text(json.dumps(out, separators=(",", ":"), allow_nan=False))
    print(len(rows), "rows,", (HERE / "scores.json").stat().st_size // 1024, "KB;", models)


if __name__ == "__main__":
    main()
