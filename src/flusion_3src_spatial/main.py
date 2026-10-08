import datetime
import shutil
import subprocess
import sys
from pathlib import Path

import click
from dateutil import relativedelta

STATE_MODELS = ["gbqr_3src_spatial", "AR6_pooled"]
US_MODELS = ["gbqr_3src", "AR6_pooled"]


@click.command()
@click.option(
    "--today_date",
    type=str,
    required=False,
    help="Date to use as effective model run date (YYYY-MM-DD)",
)
def main(today_date: str | None = None):
    """Combine existing component forecasts into the flusion_3src_spatial ensemble.

    State-level forecasts blend gbqr_3src_spatial and AR6_pooled.
    US-level forecasts blend gbqr_3src (gbqr_3src_spatial has no national-level output,
    since directional wave spatial features only support a single aggregation level)
    and AR6_pooled. All components must already have published forecasts for the
    reference date in ../../model-output/ -- this model does not run any component
    itself.
    """
    try:
        today_date = datetime.date.fromisoformat(today_date)
    except (TypeError, ValueError):  # if today_date is None or a bad format
        today_date = datetime.date.today()
    reference_date = today_date + relativedelta.relativedelta(weekday=5)

    model_output_root = Path("../../model-output")
    missing = []
    for model_abbr in sorted(set(STATE_MODELS + US_MODELS)):
        expected = model_output_root / f"UMass-{model_abbr}" / f"{reference_date}-UMass-{model_abbr}.csv"
        if not expected.exists():
            missing.append(str(expected))
    if missing:
        print("ERROR: missing required component forecast(s) for reference date "
              f"{reference_date}:", file=sys.stderr)
        for m in missing:
            print(f"  {m}", file=sys.stderr)
        sys.exit(1)

    rscript_path = shutil.which("Rscript")
    if rscript_path is None:
        print("ERROR: Rscript not found in PATH", file=sys.stderr)
        print("Please ensure R is installed and in your PATH", file=sys.stderr)
        sys.exit(1)

    subprocess.run([rscript_path, "flusion_ensemble.R", str(reference_date)], check=True)


if __name__ == "__main__":
    main()
