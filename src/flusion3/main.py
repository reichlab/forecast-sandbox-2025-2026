"""
Main script to run the flusion3 ensemble model.

This model combines three components:
1. AR(6) pooled model (without Fourier)
2. AR(6) pooled model with Fourier seasonality
3. Gradient Boosted Quantile Regression (GBQR)

The three components are combined using equal-weighted quantile averaging.
"""

import click
import datetime
from dateutil import relativedelta
import subprocess
import sys
from typing import Optional


@click.command()
@click.option(
    "--today_date",
    type=str,
    required=False,
    help="Date to use as effective model run date (YYYY-MM-DD)",
)
@click.option(
    "--short_run",
    is_flag=True,
    help="Perform a short run for testing.",
)
def main(today_date: Optional[str] = None, short_run: bool = False):
    """Generate flu predictions from flusion3 ensemble (AR6 + AR6-Fourier + GBQR)."""

    try:
        today_date = datetime.date.fromisoformat(today_date)
    except (TypeError, ValueError):
        today_date = datetime.date.today()

    reference_date = today_date + relativedelta.relativedelta(weekday=5)

    if short_run:
        short_run_flag = ["--short_run"]
    else:
        short_run_flag = []

    print(f"Running flusion3 ensemble for reference date: {reference_date}")
    print("=" * 60)

    # Run AR6 without Fourier
    print("\n1/4: Running AR6 (no Fourier)...")
    subprocess.run(
        [sys.executable, "0_ar6_pooled.py", "--reference_date", str(reference_date)]
        + short_run_flag,
        check=True,
    )

    # Run AR6 with Fourier
    print("\n2/4: Running AR6 with Fourier seasonality...")
    subprocess.run(
        [sys.executable, "1_ar6_pooled_fourier.py", "--reference_date", str(reference_date)]
        + short_run_flag,
        check=True,
    )

    # Run GBQR
    print("\n3/4: Running GBQR...")
    subprocess.run(
        [sys.executable, "2_gbqr.py", "--reference_date", str(reference_date)]
        + short_run_flag,
        check=True,
    )

    # Run ensemble
    print("\n4/4: Creating ensemble...")
    subprocess.run(
        ["Rscript", "3_flusion3_ensemble.R", str(reference_date)],
        check=True,
    )

    print("\n" + "=" * 60)
    print("Flusion3 ensemble complete!")
    print(f"Output saved to: ../../model-output/UMass-flusion3/{reference_date}-UMass-flusion3.csv")


if __name__ == "__main__":
    main()
