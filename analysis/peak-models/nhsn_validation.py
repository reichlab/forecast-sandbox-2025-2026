"""
NHSN test hindcasts for peak model configurations chosen on the ILINet development seasons (validate_ilinet.py).

For each reference date of a season's FluSight peak-week task (hub-config/tasks.json), the model is run exactly as in
real time: NHSN data and vintages as of the reference date, the NHSN revision simulation, training on seasons before
the forecast season. Influenza type counts (for the B-share features) are final WHO/NREVSS counts for training seasons
up to 2022/23 and Delphi clinical-lab counts after that; for the forecast season, the counts as published by the
reference date. Outputs go to analysis/peak-models/nhsn-hindcasts/UMass-peak_<name>/ (not model-output/); score with
    python analysis/peak-models/score_peak_forecasts.py --forecast_root analysis/peak-models/nhsn-hindcasts

Usage (from the repository root, with the idmodels environment):
    python analysis/peak-models/nhsn_validation.py --season 2023/24 --models gbqr__N1
"""
import argparse
import datetime
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
OUT = HERE / "nhsn-hindcasts"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "src"))

import epidata_ilinet as E  # noqa: E402
import validate_ilinet as V  # noqa: E402
from peak_common import make_run_config  # noqa: E402
from run_hindcasts import season_reference_dates  # noqa: E402


def _week_map(season: str) -> pd.DataFrame:
    from idmodels.peak.base import season_week_to_date

    return pd.DataFrame({"wk_end_date": pd.to_datetime([season_week_to_date(season, w) for w in range(1, 54)]),
                         "season_week": np.arange(1, 54)})


def final_strain(latest: pd.DataFrame) -> dict:
    """Final counts: WHO/NREVSS (iddata S3 file) through 2022/23, Delphi latest clinical-lab A/B for later seasons."""
    from idmodels.peak.base import season_of

    strain = dict(V.load_strain())
    latest = latest.dropna(subset=["location"]).copy()
    latest["season"] = [season_of(d.date()) for d in latest["wk_end_date"]]
    for season, g in latest.groupby("season"):
        if season <= "2022/23":
            continue
        g = g.merge(_week_map(season), on="wk_end_date")
        for loc, h in g.groupby("location"):
            arr = np.full((4, 53), np.nan)
            wk = h["season_week"].to_numpy().astype(int) - 1
            arr[0, wk], arr[1, wk] = h["total_a"].to_numpy(float), h["total_b"].to_numpy(float)
            strain[(loc, season)] = arr
    return strain


def main():
    from idmodels.peak.base import season_of

    parser = argparse.ArgumentParser()
    parser.add_argument("--season", required=True)
    parser.add_argument("--models", nargs="+", required=True)
    parser.add_argument("--skip_existing", action="store_true")
    args = parser.parse_args()

    lag_table, latest = E.nhsn_clinical_tables()
    strain = final_strain(latest)
    models = {name: V.load_model(name) for name in args.models}
    for ref in season_reference_dates(args.season):
        ref_date = datetime.date.fromisoformat(ref)
        todo = {n: m for n, m in models.items()
                if not (args.skip_existing and (OUT / f"UMass-peak_{n}" / f"{ref}-UMass-peak_{n}.csv").exists())}
        if not todo:
            continue
        t0 = time.time()
        inputs = next(iter(todo.values())).load_inputs(ref_date)
        season = season_of(ref_date)
        inputs.strain = E.strain_as_of(lag_table, ref_date, season, _week_map(season), strain)
        run_config = make_run_config(ref_date, OUT)
        for name, model in todo.items():
            df = model.forecast(inputs, run_config)
            out_dir = OUT / f"UMass-peak_{name}"
            out_dir.mkdir(parents=True, exist_ok=True)
            df.to_csv(out_dir / f"{ref}-UMass-peak_{name}.csv", index=False, na_rep="NA")
        print(f"{ref} {list(todo)}: {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
