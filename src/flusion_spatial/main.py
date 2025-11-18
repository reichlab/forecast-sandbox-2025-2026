import click
import datetime
from dateutil import relativedelta
import subprocess
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
    help="Perform a short run."
)
def main(today_date: Optional[str] = None, short_run: bool = False):
    """Generate flu predictions from flusion_spatial ensemble model."""
    try:
        today_date = datetime.date.fromisoformat(today_date)
    except (TypeError, ValueError):  # if today_date is None or a bad format
        today_date = datetime.date.today()
    reference_date = today_date + relativedelta.relativedelta(weekday=5)

    if short_run:
        short_run_flag = ["--short_run"]
    else:
        short_run_flag = []

    print(f"Running flusion_spatial ensemble for reference date: {reference_date}")
    print("=" * 60)

    print("\n[1/3] Running AR(6) pooled model...")
    subprocess.run(["python", "0_ar6_pooled.py",
                    "--reference_date", str(reference_date)] + short_run_flag,
                   check=True)

    print("\n[2/3] Running GBQR 3-source spatial model...")
    subprocess.run(["python", "1_gbqr_3src_spatial.py",
                    "--reference_date", str(reference_date)] + short_run_flag,
                   check=True)

    print("\n[3/3] Creating ensemble forecast...")
    subprocess.run(["Rscript", "2_flusion_spatial_ensemble.R", str(reference_date)],
                   check=True)

    print("\n" + "=" * 60)
    print(f"Ensemble forecast complete for {reference_date}")
    print(f"Output saved to: ../../model-output/UMass-flusion_spatial/")


if __name__ == "__main__":
    main()
