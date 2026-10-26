"""Reproduce the reference-date forecast plots used in
reading-literature/pymc-tutorial/AR2_SEASONAL_RESULTS.md and
AR6_FOURIER_POOLED_RESULTS.md, but reading directly from this model's own
production CSV output (model-output/UMass-WISAR6_fourthroot_adaptive_t/<date>-UMass-WISAR6_fourthroot_adaptive_t.csv)
rather than re-running the model.

Observed data is drawn from three sources so point-in-time correctness is
visible (see ../../WISAR6_REFERENCE_DATE_FIX.md):
  - "as of reference date": the archive snapshot main.py actually trains on
    (main._find_archive_snapshot), i.e. data published before the reference date;
  - "leaked snapshot": the archive file named for the reference date itself,
    published ~4 days later, which the pre-fix code trained on (--show-leaked);
  - "current (revised)": the hub's live target data, for evaluating the forecast.

Usage:

    uv run python plot_reference_date_forecast.py --reference-date=2025-04-05 \
        --locations=33,US --show-leaked --previous-ref=2d3651e
"""

from __future__ import annotations

import datetime
from pathlib import Path

import io
import subprocess

import click
import matplotlib.pyplot as plt
import pandas as pd

import main as wisar_main

MODEL_ABBR = "UMass-WISAR6_fourthroot_adaptive_t"
OUTPUT_ROOT = Path(__file__).resolve().parents[2] / "model-output" / MODEL_ABBR
HUB_ROOT = Path(__file__).resolve().parents[3] / "FluSight-forecast-hub"
TARGET_DATA_PATH = HUB_ROOT / "target-data" / "target-hospital-admissions.csv"
ARCHIVE_DIR = HUB_ROOT / "auxiliary-data" / "target-data-archive"
LOCATIONS_PATH = HUB_ROOT / "auxiliary-data" / "locations.csv"
FIGURES_DIR = Path(__file__).resolve().parent / "figures"


def load_forecast(reference_date: str) -> pd.DataFrame:
    path = OUTPUT_ROOT / f"{reference_date}-{MODEL_ABBR}.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"No forecast output found at {path}. Run "
            f"'uv run python main.py --today_date={reference_date}' first."
        )
    df = pd.read_csv(path)
    df["target_end_date"] = pd.to_datetime(df["target_end_date"])
    df["reference_date"] = pd.to_datetime(df["reference_date"])
    df["location"] = df["location"].astype(str)
    return df


def _clean(df: pd.DataFrame) -> pd.DataFrame:
    df["target_end_date"] = pd.to_datetime(df["target_end_date"])
    df["location"] = df["location"].astype(str)
    return df


def load_previous_forecast(reference_date: str, git_ref: str) -> pd.DataFrame | None:
    """The forecast CSV as committed at git_ref (e.g. the pre-fix version)."""
    rel = f"model-output/{MODEL_ABBR}/{reference_date}-{MODEL_ABBR}.csv"
    res = subprocess.run(["git", "show", f"{git_ref}:{rel}"], capture_output=True, text=True,
                         cwd=OUTPUT_ROOT.parents[1])
    return _clean(pd.read_csv(io.StringIO(res.stdout))) if res.returncode == 0 else None


def load_observed_context(reference_date: pd.Timestamp, loc: str, weeks_before: int = 10,
                           weeks_after: int = 6, path: Path = TARGET_DATA_PATH) -> pd.DataFrame:
    raw = pd.read_csv(path)
    raw["date"] = pd.to_datetime(raw["date"])
    raw["location"] = raw["location"].astype(str)
    sub = raw[raw["location"] == loc].sort_values("date")
    window_start = reference_date - pd.Timedelta(weeks=weeks_before)
    window_end = reference_date + pd.Timedelta(weeks=weeks_after)
    return sub[(sub["date"] >= window_start) & (sub["date"] <= window_end)]


def load_location_names() -> dict:
    loc_info = pd.read_csv(LOCATIONS_PATH)
    return dict(zip(loc_info["location"].astype(str), loc_info["location_name"]))


def plot_reference_date_forecast(forecast_df: pd.DataFrame, reference_date: pd.Timestamp,
                                  loc: str, loc_name: str, ax, show_leaked: bool = False,
                                  previous_df: pd.DataFrame | None = None) -> None:
    fdf = forecast_df[forecast_df["location"] == loc].copy()

    current = load_observed_context(reference_date, loc)
    ax.plot(current["date"], current["value"], color="0.6", lw=1, marker="o", ms=3,
            label="current (revised)")
    asof_path = wisar_main._find_archive_snapshot(reference_date.date())
    asof = load_observed_context(reference_date, loc, weeks_after=0, path=asof_path)
    asof = asof[asof["date"] < reference_date]
    ax.plot(asof["date"], asof["value"], color="black", lw=1.5, marker="o", ms=4,
            label=f"as of ref. date ({asof_path.stem.split('_')[-1]} file, trained on)")
    if show_leaked:
        leaked_path = ARCHIVE_DIR / f"target-hospital-admissions_{reference_date.date()}.csv"
        if leaked_path.exists():
            leaked = load_observed_context(reference_date, loc, weeks_after=0, path=leaked_path)
            leaked = leaked[leaked["date"] < reference_date]
            ax.plot(leaked["date"], leaked["value"], color="tab:red", lw=0, marker="x", ms=7,
                    mew=1.5, label=f"leaked snapshot ({reference_date.date()} file, pre-fix)")
    if previous_df is not None:
        prev = previous_df[(previous_df["location"] == loc)
                           & (previous_df["output_type_id"] == 0.5)].sort_values("target_end_date")
        ax.plot(prev["target_end_date"], prev["value"], color="tab:red", lw=1.5, ls="--",
                label="pre-fix median")
    ax.axvline(reference_date, color="gray", ls=":", lw=1, label="forecast origin")

    med = fdf[fdf["output_type_id"] == 0.5].sort_values("target_end_date")
    lo90 = fdf[fdf["output_type_id"] == 0.05].sort_values("target_end_date")
    hi90 = fdf[fdf["output_type_id"] == 0.95].sort_values("target_end_date")
    lo50 = fdf[fdf["output_type_id"] == 0.25].sort_values("target_end_date")
    hi50 = fdf[fdf["output_type_id"] == 0.75].sort_values("target_end_date")

    ax.plot(med["target_end_date"], med["value"], color="tab:orange", lw=2,
            marker="s", label="median (fixed)")
    ax.fill_between(med["target_end_date"], lo90["value"].values, hi90["value"].values,
                     color="tab:orange", alpha=0.15, label="90% interval")
    ax.fill_between(med["target_end_date"], lo50["value"].values, hi50["value"].values,
                     color="tab:orange", alpha=0.3, label="50% interval")

    ax.set_title(f"{loc}: {loc_name}")
    ax.set_ylabel("weekly hospitalizations")


@click.command()
@click.option("--reference-date", "reference_date_str", required=True,
              help="Reference date matching an existing model-output CSV (YYYY-MM-DD).")
@click.option("--locations", default="US,06,25,50",
              help="Comma-separated location codes to plot (default: US, California, "
                   "Massachusetts, Vermont -- matching the 4-location pymc-tutorial study).")
@click.option("--out", "out_name", default=None,
              help="Output filename (default: reference_date_forecast_<date>.png, "
                   "written under ./figures/).")
@click.option("--show-leaked", is_flag=True,
              help="Also plot the archive snapshot named for the reference date itself "
                   "(what the pre-fix code trained on).")
@click.option("--previous-ref", default=None,
              help="Git ref whose committed forecast CSV to overlay as the pre-fix median.")
def main(reference_date_str: str, locations: str, out_name: str | None, show_leaked: bool,
         previous_ref: str | None):
    reference_date = pd.Timestamp(reference_date_str)
    forecast_df = load_forecast(reference_date_str)
    previous_df = load_previous_forecast(reference_date_str, previous_ref) if previous_ref else None
    loc_names = load_location_names()
    loc_list = locations.split(",")

    n = len(loc_list)
    ncols = 2 if n > 1 else 1
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.5 * ncols, 4.5 * nrows), squeeze=False)
    axes_flat = axes.flatten()

    for ax, loc in zip(axes_flat, loc_list):
        plot_reference_date_forecast(forecast_df, reference_date, loc,
                                      loc_names.get(loc, loc), ax, show_leaked, previous_df)
    for ax in axes_flat[len(loc_list):]:
        ax.axis("off")

    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8, loc="lower center", ncol=4, frameon=False)
    fig.suptitle(f"WISAR6_fourthroot_adaptive_t forecast from reference date "
                 f"{reference_date.date()}")
    fig.autofmt_xdate()
    fig.tight_layout(rect=(0, 0.1, 1, 1))

    FIGURES_DIR.mkdir(exist_ok=True)
    out_name = out_name or f"reference_date_forecast_{reference_date_str}.png"
    out_path = FIGURES_DIR / out_name
    fig.savefig(out_path, dpi=150)
    click.echo(f"Saved {out_path}")


if __name__ == "__main__":
    main()
