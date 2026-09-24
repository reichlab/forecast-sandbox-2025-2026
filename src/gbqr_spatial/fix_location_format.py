#!/usr/bin/env python3
"""
Fix location column format in UMass-gbqr_spatial forecast files.

Converts integer location values (1, 2, 4, 5, 6, ...) to string format
with leading zeros ('01', '02', '04', '05', '06', ...) to match the
expected format for hub-dashboard-predtimechart.

Usage:
    python fix_location_format.py [--dry-run] [--backup]
"""

import argparse
import pandas as pd
from pathlib import Path
import shutil
from datetime import datetime
import csv


def fix_location_format(csv_file: Path, dry_run: bool = False, backup: bool = False) -> dict:
    """
    Fix location column format in a single CSV file.

    Args:
        csv_file: Path to the CSV file
        dry_run: If True, only report what would be changed without modifying files
        backup: If True, create a backup before modifying the file

    Returns:
        Dictionary with status information
    """
    result = {
        'file': csv_file.name,
        'modified': False,
        'location_type_before': None,
        'location_type_after': None,
        'num_rows': 0,
        'error': None
    }

    try:
        # Read the CSV file
        df = pd.read_csv(csv_file)
        result['num_rows'] = len(df)

        # Check current location type
        if 'location' not in df.columns:
            result['error'] = "No 'location' column found"
            return result

        result['location_type_before'] = str(df['location'].dtype)

        # Check if location is already a string type
        if df['location'].dtype == 'object':
            # Check if it needs zero-padding
            sample_locs = df['location'].head(5).tolist()
            if all(isinstance(loc, str) and len(loc) == 2 for loc in sample_locs if loc != 'US'):
                result['error'] = "Already in correct format (string with leading zeros)"
                return result

        # Create backup if requested
        if backup and not dry_run:
            backup_path = csv_file.with_suffix('.csv.bak')
            shutil.copy2(csv_file, backup_path)
            result['backup_path'] = str(backup_path)

        # Convert location to string with leading zeros
        # Handle both numeric locations and 'US' if present
        df['location'] = df['location'].apply(
            lambda x: str(x).zfill(2) if str(x) != 'US' and str(x).isdigit() else str(x)
        )

        result['location_type_after'] = str(df['location'].dtype)
        result['modified'] = True

        # Save the modified file (unless dry-run)
        # Write with minimal quoting, then manually add quotes to location column only
        if not dry_run:
            # First save with minimal quoting
            df.to_csv(csv_file, index=False, quoting=csv.QUOTE_MINIMAL)

            # Read the file back and manually quote only the location column
            with open(csv_file, 'r') as f:
                lines = f.readlines()

            # Process each line to quote only the location column (first column)
            modified_lines = []
            for i, line in enumerate(lines):
                if i == 0:  # Header line
                    modified_lines.append(line)
                else:
                    parts = line.split(',', 1)  # Split only on first comma
                    if len(parts) == 2:
                        # Quote the location value if not already quoted
                        location = parts[0].strip()
                        if not location.startswith('"'):
                            location = f'"{location}"'
                        modified_lines.append(f'{location},{parts[1]}')
                    else:
                        modified_lines.append(line)

            # Write back
            with open(csv_file, 'w') as f:
                f.writelines(modified_lines)

        # Get sample of changed values
        result['sample_locations'] = df['location'].unique()[:10].tolist()

    except Exception as e:
        result['error'] = str(e)

    return result


def main():
    parser = argparse.ArgumentParser(
        description='Fix location column format in UMass-gbqr_spatial forecast files'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be changed without modifying files'
    )
    parser.add_argument(
        '--backup',
        action='store_true',
        help='Create backup files (.csv.bak) before modifying'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='../../model-output/UMass-gbqr_spatial',
        help='Directory containing the CSV files to fix (default: ../../model-output/UMass-gbqr_spatial)'
    )

    args = parser.parse_args()

    # Get the output directory
    output_dir = Path(args.output_dir)
    if not output_dir.exists():
        print(f"Error: Directory not found: {output_dir}")
        return 1

    # Find all CSV files
    csv_files = sorted(output_dir.glob('*.csv'))

    if not csv_files:
        print(f"No CSV files found in {output_dir}")
        return 1

    print(f"Found {len(csv_files)} CSV files in {output_dir}")
    if args.dry_run:
        print("DRY RUN MODE - No files will be modified")
    if args.backup:
        print("Backup mode enabled - .bak files will be created")
    print()

    # Process each file
    results = []
    for csv_file in csv_files:
        result = fix_location_format(csv_file, dry_run=args.dry_run, backup=args.backup)
        results.append(result)

    # Print summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)

    modified_count = sum(1 for r in results if r['modified'])
    error_count = sum(1 for r in results if r['error'])

    print(f"\nTotal files processed: {len(results)}")
    print(f"Files modified: {modified_count}")
    print(f"Files with errors: {error_count}")

    if error_count > 0:
        print("\nFiles with errors:")
        for result in results:
            if result['error']:
                print(f"  - {result['file']}: {result['error']}")

    if modified_count > 0:
        print(f"\nSuccessfully {'would modify' if args.dry_run else 'modified'} {modified_count} files")
        print("\nSample of changed files:")
        for result in [r for r in results if r['modified']][:5]:
            print(f"  - {result['file']}")
            print(f"    Type: {result['location_type_before']} -> {result['location_type_after']}")
            if 'sample_locations' in result:
                print(f"    Sample locations: {result['sample_locations'][:5]}")

    if args.dry_run and modified_count > 0:
        print(f"\nTo apply these changes, run without --dry-run flag")

    print()

    return 0


if __name__ == '__main__':
    exit(main())
