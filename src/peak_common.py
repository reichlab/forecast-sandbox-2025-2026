"""Run configuration shared by the peak_* models (seasonal peak week and peak size of flu hospital admissions)."""
import datetime
from pathlib import Path

from dateutil import relativedelta
from iddata.enums import Disease
from idmodels.config import RunConfig

LOCATIONS = ["US", "01", "02", "04", "05", "06", "08", "09", "10", "11", "12", "13", "15", "16", "17", "18", "19", "20",
             "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38",
             "39", "40", "41", "42", "44", "45", "46", "47", "48", "49", "50", "51", "53", "54", "55", "56", "72"]
Q_LEVELS = [0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
            0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.975, 0.99]
Q_LABELS = ["0.01", "0.025", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5",
            "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95", "0.975", "0.99"]


def reference_date_from(today_date: str | None) -> datetime.date:
    """The Saturday on or after today_date (YYYY-MM-DD; defaults to today)."""
    try:
        today = datetime.date.fromisoformat(today_date)
    except (TypeError, ValueError):  # today_date is None or badly formatted
        today = datetime.date.today()
    return today + relativedelta.relativedelta(weekday=5)


def make_run_config(ref_date: datetime.date, output_root: Path) -> RunConfig:
    return RunConfig(
        disease=Disease.FLU,
        ref_date=ref_date,
        output_root=output_root,
        artifact_store_root=None,
        max_horizon=4,  # unused by the peak models
        states=LOCATIONS,
        hsas=[],
        q_levels=Q_LEVELS,
        q_labels=Q_LABELS)
