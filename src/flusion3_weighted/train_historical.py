#!/usr/bin/env python3
"""
Train the flusion3_weighted model on historical dates.

This script runs the weighted ensemble for multiple historical dates
to build up a history of forecasts that can be used to learn weights.
"""

import subprocess
import sys
from datetime import datetime, timedelta

# Historical dates to run (full training set with more bags)
# Using 8 training dates and 3 test dates for more robust weight learning
TRAINING_DATES = [
    "2023-10-21",
    "2023-10-28",
    "2023-11-04",
    "2023-11-11",
    "2023-11-18",
    "2023-11-25",
    "2023-12-02",
    "2023-12-09",
]

# Test dates (after training period, where weights will be learned)
TEST_DATES = [
    "2023-12-16",
    "2023-12-23",
    "2023-12-30",
]

def run_forecast(date_str, short_run=False):
    """Run the weighted ensemble model for a specific date."""
    cmd = [sys.executable, "main.py", "--today_date", date_str]
    if short_run:
        cmd.append("--short_run")

    print(f"\n{'='*60}")
    print(f"Running forecast for {date_str}")
    print(f"{'='*60}")

    result = subprocess.run(cmd, check=False)
    return result.returncode == 0

def main():
    print("Training flusion3_weighted on historical data...")
    print(f"Training dates: {len(TRAINING_DATES)}")
    print(f"Test dates: {len(TEST_DATES)}")

    # Run training dates
    print("\n" + "="*60)
    print("PHASE 1: Training period (building history)")
    print("="*60)

    training_success = 0
    for date in TRAINING_DATES:
        if run_forecast(date, short_run=False):
            training_success += 1
        else:
            print(f"WARNING: Forecast for {date} failed")

    print(f"\nTraining phase: {training_success}/{len(TRAINING_DATES)} succeeded")

    # Run test dates (where weights will be learned)
    print("\n" + "="*60)
    print("PHASE 2: Test period (learning weights from history)")
    print("="*60)

    test_success = 0
    for date in TEST_DATES:
        if run_forecast(date, short_run=False):
            test_success += 1
        else:
            print(f"WARNING: Forecast for {date} failed")

    print(f"\nTest phase: {test_success}/{len(TEST_DATES)} succeeded")

    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Training: {training_success}/{len(TRAINING_DATES)} succeeded")
    print(f"Test: {test_success}/{len(TEST_DATES)} succeeded")
    print(f"Total: {training_success + test_success}/{len(TRAINING_DATES) + len(TEST_DATES)} succeeded")

if __name__ == "__main__":
    main()
