import argparse
import datetime
import json
from pathlib import Path

import pandas as pd
from dateutil import relativedelta

from axial.eval.backtest import make_context
from axial.eval.splits import season_of
from axial.models import peaks, targets
from axial.models.transformer import FreshNsspHorizonWeightedR

HUB = Path(__file__).resolve().parents[2]
MODEL_ID = "axial-tfm_L4v5"
INC = "wk inc flu hosp"
# axial's NSSP store is point-in-time from Delphi's 2024 epiweek 16 issue on; before that it is a pseudo-vintage
# built from later (revised) values, so no ED forecast is written for earlier reference dates.
ED_FROM = pd.Timestamp("2024-04-21")
COLUMNS = ["location", "reference_date", "horizon", "target_end_date", "target", "output_type", "output_type_id", "value"]


def peak_bins() -> pd.DatetimeIndex:
    """The peak-week bins (Saturdays) in this hub's tasks.json."""
    tasks = json.loads((HUB / "hub-config" / "tasks.json").read_text())["rounds"][0]["model_tasks"]
    ids = next(m for m in tasks if m["task_ids"]["target"]["required"] == [peaks.PEAK_WEEK])["output_type"]["pmf"]["output_type_id"]["required"]
    return pd.DatetimeIndex(sorted(pd.to_datetime(ids)))


def weekly_rows(ctx, ref: pd.Timestamp) -> pd.DataFrame:
    """Admissions (always) and ED-visit proportion (point-in-time NSSP, all four horizons), horizons 0-3."""
    fc = FreshNsspHorizonWeightedR().forecast(ctx)
    parts = [targets.to_hub(fc[fc["horizon"].between(0, 3)], ref, INC)]
    ed = targets.ed_forecast(ctx) if ref >= ED_FROM else pd.DataFrame()
    if len(ed) and (ed.groupby("location")["horizon"].nunique() == 4).all():
        parts.append(targets.to_hub(ed, ref, targets.ED))
    return pd.concat(parts, ignore_index=True)


def peak_rows(ctx, ref: pd.Timestamp) -> pd.DataFrame | None:
    """Peak size and peak week over this hub's bins, for reference dates in the bins' season."""
    bins = peak_bins()
    season = season_of(ref)
    if season != season_of(bins[0]):
        return None
    peaks.WINDOWS[season] = (str(bins[0].date()), str(bins[-1].date()))  # the hub's bins, not axial's FluSight window
    pw, pk = peaks.peak_forecast(ctx)
    return peaks.to_hub(pw, pk, ref)


def write(rows: pd.DataFrame, ref: pd.Timestamp) -> Path:
    out = HUB / "model-output" / MODEL_ID / f"{ref.date()}-{MODEL_ID}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = rows.assign(horizon=rows["horizon"].astype("Int64"))[COLUMNS]
    rows.to_csv(out, index=False, na_rep="NA")
    return out


def main():
    """Generate axial-tfm_L4v5 forecasts for the reference date on or after --today_date."""
    p = argparse.ArgumentParser(description=main.__doc__)
    p.add_argument("--today_date", type=str, help="Date to use as effective model run date (YYYY-MM-DD)")
    try:
        today_date = datetime.date.fromisoformat(p.parse_args().today_date)
    except (TypeError, ValueError):  # if today_date is None or a bad format
        today_date = datetime.date.today()
    ref = pd.Timestamp(today_date + relativedelta.relativedelta(weekday=5))
    ctx = make_context(ref)
    rows = [weekly_rows(ctx, ref), peak_rows(ctx, ref)]
    print(write(pd.concat([r for r in rows if r is not None], ignore_index=True), ref))


if __name__ == "__main__":
    main()
