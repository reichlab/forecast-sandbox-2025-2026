#!/usr/bin/env python3
"""
Run UMass-flusion3_weighted on all dates to match the full date range of other models.
"""

import subprocess
import sys
from datetime import datetime

# All dates from flusion model
ALL_DATES = [
    "2023-10-21", "2023-10-28", "2023-11-04", "2023-11-11", "2023-11-18", "2023-11-25",
    "2023-12-02", "2023-12-09", "2023-12-16", "2023-12-23", "2023-12-30", "2024-01-06",
    "2024-01-13", "2024-01-20", "2024-01-27", "2024-02-03", "2024-02-10", "2024-02-17",
    "2024-02-24", "2024-03-02", "2024-03-09", "2024-03-16", "2024-03-23", "2024-03-30",
    "2024-04-06", "2024-04-13", "2024-04-20", "2024-04-27", "2024-05-04", "2024-05-11",
    "2024-11-30", "2024-12-07", "2024-12-14", "2024-12-21", "2024-12-28",
    "2025-01-04", "2025-01-11", "2025-01-18", "2025-01-25", "2025-02-01", "2025-02-08",
    "2025-02-15", "2025-02-22", "2025-03-01", "2025-03-08", "2025-03-15", "2025-03-22",
    "2025-03-29", "2025-04-05", "2025-04-12", "2025-04-19", "2025-04-26",
    "2025-05-03", "2025-05-10", "2025-05-17"
]

# Dates already completed
COMPLETED_DATES = [
    "2023-10-21", "2023-10-28", "2023-11-04", "2023-11-11", "2023-11-18", "2023-11-25",
    "2023-12-02", "2023-12-09", "2023-12-16", "2023-12-23", "2023-12-30", "2024-01-06"
]

def run_forecast(date_str):
    """Run the model for a specific date."""
    cmd = [sys.executable, "main.py", "--today_date", date_str]

    print(f"\n{'='*60}")
    print(f"Running forecast for {date_str}")
    print(f"{'='*60}")

    result = subprocess.run(cmd, check=False)
    return result.returncode == 0

def main():
    # Find dates to run
    dates_to_run = [d for d in ALL_DATES if d not in COMPLETED_DATES]

    print("=" * 60)
    print("Running UMass-flusion3_weighted on all remaining dates")
    print("=" * 60)
    print(f"Total dates: {len(ALL_DATES)}")
    print(f"Already completed: {len(COMPLETED_DATES)}")
    print(f"To run: {len(dates_to_run)}")
    print(f"Estimated time: {len(dates_to_run) * 10} - {len(dates_to_run) * 15} minutes")
    print("=" * 60)

    success_count = 0
    failed_dates = []

    start_time = datetime.now()

    for i, date in enumerate(dates_to_run, 1):
        print(f"\nProgress: {i}/{len(dates_to_run)}")
        if run_forecast(date):
            success_count += 1
        else:
            print(f"WARNING: Forecast for {date} failed")
            failed_dates.append(date)

    end_time = datetime.now()
    elapsed = end_time - start_time

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Successful: {success_count}/{len(dates_to_run)}")
    print(f"Failed: {len(failed_dates)}/{len(dates_to_run)}")
    if failed_dates:
        print(f"Failed dates: {', '.join(failed_dates)}")
    print(f"Total time: {elapsed}")
    print(f"Average time per date: {elapsed / len(dates_to_run) if dates_to_run else 0}")

if __name__ == "__main__":
    main()
