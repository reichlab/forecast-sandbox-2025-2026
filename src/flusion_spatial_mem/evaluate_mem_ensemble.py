#!/usr/bin/env python3
"""
Evaluation script for flusion_spatial_mem model.

This script generates forecasts for all Saturdays in a specified date range,
typically covering the last two flu seasons. This allows for retrospective
evaluation and comparison with other models like flusion_spatial2_prod.

Usage:
    python evaluate_mem_ensemble.py --start_date=2023-10-07 --end_date=2025-05-31
    python evaluate_mem_ensemble.py --start_date=2023-10-07 --end_date=2025-05-31 --short_run
"""

import click
import datetime
from dateutil import relativedelta
import subprocess
import sys
from pathlib import Path
from typing import List


def get_saturdays_in_range(start_date: datetime.date, end_date: datetime.date) -> List[datetime.date]:
    """
    Get all Saturdays within the specified date range.

    Args:
        start_date: Start of date range
        end_date: End of date range

    Returns:
        List of Saturday dates
    """
    # Find the first Saturday on or after start_date
    current = start_date
    # Move to the next Saturday (weekday 5)
    days_until_saturday = (5 - current.weekday()) % 7
    if days_until_saturday > 0:
        current = current + datetime.timedelta(days=days_until_saturday)

    saturdays = []
    while current <= end_date:
        saturdays.append(current)
        current = current + datetime.timedelta(days=7)

    return saturdays


@click.command()
@click.option(
    "--start_date",
    type=str,
    required=True,
    help="Start date for evaluation period (YYYY-MM-DD)",
)
@click.option(
    "--end_date",
    type=str,
    required=True,
    help="End date for evaluation period (YYYY-MM-DD)",
)
@click.option(
    "--short_run",
    is_flag=True,
    help="Perform short runs (fewer quantiles, less MCMC warmup)."
)
@click.option(
    "--dry_run",
    is_flag=True,
    help="Print dates that would be run without actually running the model."
)
def main(start_date: str, end_date: str, short_run: bool, dry_run: bool):
    """
    Generate forecasts for all Saturdays in the specified date range.

    This script is designed for retrospective evaluation of the flusion_spatial_mem
    model. It generates forecasts for each Saturday in the range, allowing comparison
    with actual outcomes and other models.

    Example date ranges for flu seasons:
    - 2023-2024 season: 2023-10-07 to 2024-05-18
    - 2024-2025 season: 2024-10-05 to 2025-05-17
    - Both seasons: 2023-10-07 to 2025-05-31
    """
    try:
        start = datetime.date.fromisoformat(start_date)
        end = datetime.date.fromisoformat(end_date)
    except ValueError as e:
        print(f"ERROR: Invalid date format: {e}", file=sys.stderr)
        print("Dates must be in YYYY-MM-DD format", file=sys.stderr)
        sys.exit(1)

    if start > end:
        print("ERROR: start_date must be before end_date", file=sys.stderr)
        sys.exit(1)

    # Get all Saturdays in the range
    saturdays = get_saturdays_in_range(start, end)

    print("=" * 70)
    print(f"Flusion Spatial MEM Ensemble Evaluation")
    print("=" * 70)
    print(f"Date range: {start_date} to {end_date}")
    print(f"Total Saturdays to process: {len(saturdays)}")
    print(f"Short run mode: {short_run}")
    print(f"Dry run mode: {dry_run}")
    print("=" * 70)

    if dry_run:
        print("\nSaturdays that would be processed:")
        for i, saturday in enumerate(saturdays, 1):
            print(f"  {i:2d}. {saturday}")
        print(f"\nTotal: {len(saturdays)} forecast dates")
        return

    # Track successes and failures
    successful = []
    failed = []

    # Run forecasts for each Saturday
    for i, saturday in enumerate(saturdays, 1):
        print(f"\n{'=' * 70}")
        print(f"Processing {i}/{len(saturdays)}: {saturday}")
        print(f"{'=' * 70}")

        # Build command
        cmd = [sys.executable, "main.py", f"--today_date={saturday}"]
        if short_run:
            cmd.append("--short_run")

        # Run the forecast
        try:
            result = subprocess.run(cmd, check=True, capture_output=False)
            successful.append(saturday)
            print(f"✓ Successfully completed forecast for {saturday}")
        except subprocess.CalledProcessError as e:
            failed.append(saturday)
            print(f"✗ FAILED to complete forecast for {saturday}", file=sys.stderr)
            print(f"  Error: {e}", file=sys.stderr)
            # Continue with next date even if this one failed

    # Print summary
    print("\n" + "=" * 70)
    print("EVALUATION SUMMARY")
    print("=" * 70)
    print(f"Total Saturdays processed: {len(saturdays)}")
    print(f"Successful: {len(successful)}")
    print(f"Failed: {len(failed)}")

    if failed:
        print("\nFailed dates:")
        for date in failed:
            print(f"  - {date}")

    print("\nOutput directory: ../../model-output/UMass-flusion_spatial_mem/")

    # Check if output directory exists and count files
    output_dir = Path("../../model-output/UMass-flusion_spatial_mem")
    if output_dir.exists():
        csv_files = list(output_dir.glob("*.csv"))
        print(f"CSV files in output directory: {len(csv_files)}")
    else:
        print("Output directory does not exist yet")

    print("=" * 70)

    # Exit with error if any forecasts failed
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
