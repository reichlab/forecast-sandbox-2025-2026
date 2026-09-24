"""
Compare forecasts between gbqr_3src and gbqr_3src_spatial models.
Find dates where the models differ most substantially.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime

def load_forecast(model_name, reference_date):
    """Load forecast file for a given model and date."""
    file_path = Path(f"model-output/UMass-{model_name}/{reference_date}-UMass-{model_name}.csv")
    if file_path.exists():
        df = pd.read_csv(file_path, dtype={'location': str, 'output_type_id': str})
        return df
    return None

def compare_forecasts(df1, df2, reference_date):
    """Compare two forecast dataframes and compute difference metrics."""
    # Merge on common keys
    merged = df1.merge(
        df2,
        on=['reference_date', 'target', 'horizon', 'location', 'target_end_date', 'output_type', 'output_type_id'],
        suffixes=('_base', '_spatial')
    )

    # Only look at quantile forecasts
    merged = merged[merged['output_type'] == 'quantile'].copy()

    # Compute absolute differences
    merged['abs_diff'] = np.abs(merged['value_spatial'] - merged['value_base'])
    merged['rel_diff'] = merged['abs_diff'] / (merged['value_base'] + 1e-6)  # avoid division by zero

    # Aggregate metrics
    metrics = {
        'reference_date': reference_date,
        'mean_abs_diff': merged['abs_diff'].mean(),
        'median_abs_diff': merged['abs_diff'].median(),
        'max_abs_diff': merged['abs_diff'].max(),
        'mean_rel_diff': merged['rel_diff'].mean(),
        'std_abs_diff': merged['abs_diff'].std(),
        'q95_abs_diff': merged['abs_diff'].quantile(0.95),
        'n_forecasts': len(merged)
    }

    # Compute differences by horizon
    horizon_diffs = merged.groupby('horizon')['abs_diff'].mean().to_dict()
    for h, diff in horizon_diffs.items():
        metrics[f'mean_abs_diff_h{h}'] = diff

    # Compute median forecast differences (quantile 0.5)
    median_forecasts = merged[merged['output_type_id'] == '0.5'].copy()
    if len(median_forecasts) > 0:
        metrics['mean_abs_diff_median'] = median_forecasts['abs_diff'].mean()
        metrics['max_abs_diff_median'] = median_forecasts['abs_diff'].max()

    # Geographic variance - do some locations differ more?
    loc_diffs = merged.groupby('location')['abs_diff'].mean()
    metrics['geographic_variance'] = loc_diffs.var()
    metrics['max_location_diff'] = loc_diffs.max()

    return metrics, merged

def main():
    # Get all dates from one model
    model_dir = Path("model-output/UMass-gbqr_3src/")
    forecast_files = sorted(model_dir.glob("*.csv"))
    dates = [f.stem.split('-UMass-')[0] for f in forecast_files]

    print(f"Found {len(dates)} forecast dates to compare")
    print(f"Date range: {dates[0]} to {dates[-1]}\n")

    # Compare forecasts for each date
    results = []
    for date in dates:
        df_base = load_forecast("gbqr_3src", date)
        df_spatial = load_forecast("gbqr_3src_spatial", date)

        if df_base is not None and df_spatial is not None:
            metrics, _ = compare_forecasts(df_base, df_spatial, date)
            results.append(metrics)
            print(f"Processed {date}")

    # Create results dataframe
    results_df = pd.DataFrame(results)

    # Save full results
    results_df.to_csv("forecast_comparison_results.csv", index=False)
    print(f"\nSaved full results to forecast_comparison_results.csv")

    # Print summary statistics
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    print(f"\nOverall mean absolute difference: {results_df['mean_abs_diff'].mean():.2f}")
    print(f"Overall std of absolute difference: {results_df['mean_abs_diff'].std():.2f}")

    # Find dates with largest differences
    print("\n" + "="*80)
    print("TOP 10 DATES BY MEAN ABSOLUTE DIFFERENCE")
    print("="*80)
    top_dates = results_df.nlargest(10, 'mean_abs_diff')[
        ['reference_date', 'mean_abs_diff', 'median_abs_diff', 'max_abs_diff',
         'mean_abs_diff_median', 'geographic_variance']
    ]
    print(top_dates.to_string(index=False))

    # Find dates with largest geographic variance
    print("\n" + "="*80)
    print("TOP 10 DATES BY GEOGRAPHIC VARIANCE")
    print("="*80)
    top_geo = results_df.nlargest(10, 'geographic_variance')[
        ['reference_date', 'geographic_variance', 'mean_abs_diff',
         'max_location_diff']
    ]
    print(top_geo.to_string(index=False))

    # Find dates with largest differences in median forecasts
    print("\n" + "="*80)
    print("TOP 10 DATES BY MEDIAN FORECAST DIFFERENCE")
    print("="*80)
    top_median = results_df.nlargest(10, 'mean_abs_diff_median')[
        ['reference_date', 'mean_abs_diff_median', 'max_abs_diff_median']
    ]
    print(top_median.to_string(index=False))

    # Analyze by season
    results_df['year_month'] = pd.to_datetime(results_df['reference_date']).dt.to_period('M')
    monthly_avg = results_df.groupby('year_month')['mean_abs_diff'].mean().reset_index()
    print("\n" + "="*80)
    print("AVERAGE DIFFERENCE BY MONTH")
    print("="*80)
    print(monthly_avg.to_string(index=False))

    # Recommend dates for detailed analysis
    print("\n" + "="*80)
    print("RECOMMENDED DATES FOR DETAILED ANALYSIS")
    print("="*80)

    # Pick one early season date with high difference
    early_season = results_df[pd.to_datetime(results_df['reference_date']).dt.month.isin([10, 11, 12])]
    if len(early_season) > 0:
        top_early = early_season.nlargest(1, 'mean_abs_diff').iloc[0]
        print(f"\n1. Early season (Oct-Dec): {top_early['reference_date']}")
        print(f"   Mean abs diff: {top_early['mean_abs_diff']:.2f}")
        print(f"   Geographic variance: {top_early['geographic_variance']:.2f}")

    # Pick one peak season date with high difference
    peak_season = results_df[pd.to_datetime(results_df['reference_date']).dt.month.isin([1, 2])]
    if len(peak_season) > 0:
        top_peak = peak_season.nlargest(1, 'mean_abs_diff').iloc[0]
        print(f"\n2. Peak season (Jan-Feb): {top_peak['reference_date']}")
        print(f"   Mean abs diff: {top_peak['mean_abs_diff']:.2f}")
        print(f"   Geographic variance: {top_peak['geographic_variance']:.2f}")

    # Pick one late season date
    late_season = results_df[pd.to_datetime(results_df['reference_date']).dt.month.isin([3, 4, 5])]
    if len(late_season) > 0:
        top_late = late_season.nlargest(1, 'mean_abs_diff').iloc[0]
        print(f"\n3. Late season (Mar-May): {top_late['reference_date']}")
        print(f"   Mean abs diff: {top_late['mean_abs_diff']:.2f}")
        print(f"   Geographic variance: {top_late['geographic_variance']:.2f}")

if __name__ == "__main__":
    main()
