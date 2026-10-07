"""
UMass-peakGB: gradient-boosted direct forecasts of the FluSight seasonal targets "peak week inc flu hosp" (pmf) and
"peak inc flu hosp" (quantiles).

Laid out like the models in reichlab/operational-models (e.g. flu_flusion) so the folder can be moved there as is: run
from this directory; the forecast is written to output/model-output/UMass-peakGB/ and plot.R then writes
output/plots/<reference_date>-UMass-peakGB.pdf.
"""

import datetime
import subprocess
from pathlib import Path

import click
from dateutil import relativedelta
from iddata.enums import Disease
from idmodels.config import PeakGBQRModelConfig, RunConfig
from idmodels.peak import PeakGBQRModel

MODEL_NAME = "peakGB"

LOCATIONS = ["US", "01", "02", "04", "05", "06", "08", "09", "10", "11", "12", "13", "15", "16", "17", "18", "19", "20",
             "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38",
             "39", "40", "41", "42", "44", "45", "46", "47", "48", "49", "50", "51", "53", "54", "55", "56", "72"]
Q_LEVELS = [0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
            0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.975, 0.99]
Q_LABELS = ["0.01", "0.025", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5",
            "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95", "0.975", "0.99"]


def build_model_config() -> PeakGBQRModelConfig:
    # peak size: core + SB (synchrony + burden) features; peak timing: core + holiday-week features.
    # Everything else (training sources ILINet + FluSurv-NET, 25 bags, revision simulation, pmf floor) is the default.
    return PeakGBQRModelConfig(
        model_name=MODEL_NAME,
        size_feature_groups=["core", "sb"],
        timing_feature_groups=["core", "holiday"],
    )


def build_run_config(ref_date: datetime.date, output_root: Path) -> RunConfig:
    return RunConfig(
        disease=Disease.FLU,
        ref_date=ref_date,
        output_root=output_root,
        artifact_store_root=None,
        max_horizon=4,  # unused by the peak models
        states=LOCATIONS,
        hsas=[],
        q_levels=Q_LEVELS,
        q_labels=Q_LABELS,
    )


@click.command()
@click.option("--today_date", type=str, required=False, help="Date to use as effective model run date (YYYY-MM-DD)")
@click.option("--short_run", is_flag=True, help="Run with reduced parameters for faster testing")
def main(today_date: str | None = None, short_run: bool = False):
    """Generate peak week and peak size forecasts from the peakGB model and plot them."""
    try:
        today_date = datetime.date.fromisoformat(today_date)
    except (TypeError, ValueError):  # if today_date is None or a bad format
        today_date = datetime.date.today()
    reference_date = today_date + relativedelta.relativedelta(weekday=5)

    model_config = build_model_config()
    if short_run:
        model_config.num_bags = 3
        model_config.num_revision_draws = 20
    run_config = build_run_config(reference_date, Path("output/model-output"))
    PeakGBQRModel(model_config).run(run_config)

    subprocess.run(["Rscript", "plot.R", str(reference_date)])


if __name__ == "__main__":
    main()
