import click
import datetime
from dateutil import relativedelta
import subprocess
import shutil
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
    help="Perform a short run."
)
def main(today_date: Optional[str] = None, short_run: bool = False):
    """Generate flu predictions from flusion_spatial_mem ensemble model with MEM-based adaptive weighting."""
    try:
        today_date = datetime.date.fromisoformat(today_date)
    except (TypeError, ValueError):  # if today_date is None or a bad format
        today_date = datetime.date.today()
    reference_date = today_date + relativedelta.relativedelta(weekday=5)

    if short_run:
        short_run_flag = ["--short_run"]
    else:
        short_run_flag = []

    print(f"Running flusion_spatial_mem ensemble (MEM-based adaptive weighting) for reference date: {reference_date}")
    print("=" * 60)

    print("\n[1/4] Running AR(6) pooled model...")
    subprocess.run([sys.executable, "0_ar6_pooled.py",
                    "--reference_date", str(reference_date)] + short_run_flag,
                   check=True)

    print("\n[2/4] Running GBQR 3-source spatial model...")
    subprocess.run([sys.executable, "1_gbqr_3src_spatial2.py",
                    "--reference_date", str(reference_date)] + short_run_flag,
                   check=True)

    print("\n[3/4] Running GBQR 3-source spatial model...")
    subprocess.run([sys.executable, "2_gbqr_3src.py",
                    "--reference_date", str(reference_date)] + short_run_flag,
                   check=True)

    print("\n[4/4] Creating ensemble forecast...")
    # Find Rscript in PATH
    rscript_path = shutil.which("Rscript")
    if rscript_path is None:
        print("ERROR: Rscript not found in PATH", file=sys.stderr)
        print("Please ensure R is installed and in your PATH", file=sys.stderr)
        sys.exit(1)

    print(f"Using Rscript: {rscript_path}")
    subprocess.run([rscript_path, "3_flusion_mem_ensemble.R", str(reference_date)],
                   check=True)

    print("\n" + "=" * 60)
    print(f"Ensemble forecast complete for {reference_date}")
    print(f"Output saved to: ../../model-output/UMass-flusion_spatial_mem/")


if __name__ == "__main__":
    main()
