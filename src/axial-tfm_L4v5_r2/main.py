import argparse
import datetime
import json
from pathlib import Path

import pandas as pd
from dateutil import relativedelta

from axial.eval.backtest import make_context
from axial.eval.splits import season_of
from axial.models import peaks, targets
from axial.models.transformer import R2

HUB = Path(__file__).resolve().parents[2]
MODEL_ID = "axial-tfm_L4v5_r2"
INC = "wk inc flu hosp"
# The hub's 2025/26 peak-week bins (its tasks.json before the 2026/27 season was added); peak sizes are taken over them.
PEAK_SEASON, PEAK_BINS = "2025/26", ("2025-11-22", "2026-05-30")
COLUMNS = ["location", "reference_date", "horizon", "target_end_date", "target", "output_type", "output_type_id", "value"]


def ed_dates() -> set:
    """Reference dates for which FluSight-baseline has ED forecasts here; ED rows are written only for these."""
    files = sorted((HUB / "model-output" / "FluSight-baseline").glob("*.csv"))
    return {f.name[:10] for f in files if targets.ED in f.read_text()}


def weekly_rows(ctx, ref: pd.Timestamp) -> pd.DataFrame:
    """Admissions (always) and ED-visit proportion (dates with a FluSight-baseline ED forecast), horizons 0-3."""
    fc = R2().forecast(ctx)
    parts = [targets.to_hub(fc[fc["horizon"].between(0, 3)], ref, INC)]
    ed = targets.ed_forecast(ctx) if str(ref.date()) in ed_dates() else pd.DataFrame()
    if len(ed) and (ed.groupby("location")["horizon"].nunique() == 4).all():
        parts.append(targets.to_hub(ed, ref, targets.ED))
    return pd.concat(parts, ignore_index=True)


def peak_bins() -> pd.DatetimeIndex:
    """The peak-week bins (Saturdays) in this hub's tasks.json."""
    tasks = json.loads((HUB / "hub-config" / "tasks.json").read_text())["rounds"][0]["model_tasks"]
    ids = next(m for m in tasks if m["task_ids"]["target"]["required"] == [peaks.PEAK_WEEK])["output_type"]["pmf"]["output_type_id"]["required"]
    return pd.DatetimeIndex(sorted(pd.to_datetime(ids)))


def peak_rows(ctx, ref: pd.Timestamp) -> pd.DataFrame | None:
    """Both peak targets over the hub's bins, for reference dates in the bins' season (2026/27). For 2025/26, peak
    size only, over that season's bins: the hub's peak-week output_type_ids are now the 2026/27 weeks."""
    bins, season = peak_bins(), season_of(ref)
    if season == season_of(bins[0]):
        peaks.WINDOWS[season] = (str(bins[0].date()), str(bins[-1].date()))
        return peaks.to_hub(*peaks.peak_forecast(ctx), ref)
    if season != PEAK_SEASON:
        return None
    peaks.WINDOWS[PEAK_SEASON] = PEAK_BINS
    rows = peaks.to_hub(*peaks.peak_forecast(ctx), ref)
    return rows[rows["target"] == peaks.PEAK_INC]


def write(rows: pd.DataFrame, ref: pd.Timestamp) -> Path:
    out = HUB / "model-output" / MODEL_ID / f"{ref.date()}-{MODEL_ID}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = rows.assign(horizon=rows["horizon"].astype("Int64"))[COLUMNS]
    rows.to_csv(out, index=False, na_rep="NA")
    return out


def main():
    """Generate axial-tfm_L4v5_r2 forecasts for the reference date on or after --today_date."""
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
