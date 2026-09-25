"""
Build the observed ("oracle") peak week and peak size for each season and location from the hub's weekly target data.

The peak is taken over the 34-week window of each season's peak-week pmf dates in hub-config/tasks.json (season
weeks 10-43). If several weeks share the maximum count the location's peak week is marked as tied; FluSight does not
score peak-week forecasts in that case.

Usage (from the repository root):
    python analysis/peak-models/build_peak_oracle.py
Writes analysis/peak-models/peak-oracle.csv with columns
    season, location, peak_week (date), peak_count, tied, n_weeks
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def season_windows(tasks_path: Path) -> dict[str, list[str]]:
    """Map each season (named by its first window date's year) to its list of pmf dates."""
    tasks = json.loads(tasks_path.read_text())
    windows = {}
    for mt in tasks["rounds"][0]["model_tasks"]:
        if mt["task_ids"]["target"]["required"] == ["peak week inc flu hosp"]:
            dates = mt["output_type"]["pmf"]["output_type_id"]["required"]
            y = int(dates[0][:4])
            windows[f"{y}/{str(y + 1)[2:]}"] = dates
    return windows


def main():
    oracle = pd.read_csv(ROOT / "target-data" / "oracle-output.csv", dtype={"location": str}, low_memory=False)
    wk = (oracle.loc[oracle["target"] == "wk inc flu hosp", ["target_end_date", "location", "oracle_value"]]
          .drop_duplicates(["target_end_date", "location"]))
    rows = []
    for season, dates in season_windows(ROOT / "hub-config" / "tasks.json").items():
        w = wk.loc[wk["target_end_date"].isin(dates)]
        for loc, g in w.groupby("location"):
            mx = g["oracle_value"].max()
            at_max = g.loc[g["oracle_value"] == mx, "target_end_date"].sort_values()
            rows.append({"season": season, "location": loc, "peak_week": at_max.iloc[0], "peak_count": mx,
                         "tied": len(at_max) > 1, "n_weeks": len(g)})
    out = pd.DataFrame(rows)
    # only seasons whose window is fully observed
    out = out.loc[out["n_weeks"] == 34]
    out.to_csv(Path(__file__).parent / "peak-oracle.csv", index=False)
    print(out.groupby("season").agg(n=("location", "size"), tied=("tied", "sum")))


if __name__ == "__main__":
    main()
