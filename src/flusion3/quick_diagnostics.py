"""
Quick diagnostic comparison of a single forecast date across models.
"""

import pandas as pd
import numpy as np
from pathlib import Path

MODEL_OUTPUT_DIR = Path("../../model-output")
TARGET_DATA_URL = "https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/target-data/target-hospital-admissions.csv"

def load_target_data():
    """Load observed hospitalization data."""
    df = pd.read_csv(TARGET_DATA_URL)
    df['date'] = pd.to_datetime(df['date'])
    return df

def load_single_forecast(model_name, ref_date):
    """Load a single forecast file."""
    file_path = MODEL_OUTPUT_DIR / model_name / f"{ref_date}-{model_name}.csv"
    if file_path.exists():
        df = pd.read_csv(file_path)
        df['reference_date'] = pd.to_datetime(ref_date)
        df['target_end_date'] = pd.to_datetime(df['target_end_date'])
        return df
    return None

def compute_wis_single_date(forecasts, targets):
    """Compute WIS for a single forecast date."""
    quantile_forecasts = forecasts[forecasts['output_type'] == 'quantile'].copy()
    quantile_forecasts['output_type_id'] = quantile_forecasts['output_type_id'].astype(float)

    merged = quantile_forecasts.merge(
        targets[['location', 'date', 'value']],
        left_on=['location', 'target_end_date'],
        right_on=['location', 'date'],
        how='inner',
        suffixes=('_forecast', '_observed')
    )

    quantiles = sorted(merged['output_type_id'].unique())
    results = []

    for (loc, target_date, horizon), group in merged.groupby(['location', 'target_end_date', 'horizon']):
        obs = group['value_observed'].iloc[0]

        # Median AE
        median_forecast = group[np.isclose(group['output_type_id'], 0.5)]['value_forecast'].values
        if len(median_forecast) > 0:
            median_ae = abs(obs - median_forecast[0])
        else:
            continue

        # Interval scores
        interval_scores = []
        K = len(quantiles) // 2

        for i in range(K):
            alpha = quantiles[i]
            lower_q = alpha
            upper_q = quantiles[-(i+1)]

            lower_val = group[np.isclose(group['output_type_id'], lower_q)]['value_forecast'].values
            upper_val = group[np.isclose(group['output_type_id'], upper_q)]['value_forecast'].values

            if len(lower_val) > 0 and len(upper_val) > 0:
                lower_val = lower_val[0]
                upper_val = upper_val[0]
                width = upper_val - lower_val

                if obs < lower_val:
                    penalty = (2 / alpha) * (lower_val - obs)
                elif obs > upper_val:
                    penalty = (2 / alpha) * (obs - upper_val)
                else:
                    penalty = 0

                interval_scores.append(width + penalty)

        wis = 0.5 * median_ae + 0.5 * np.mean(interval_scores) if interval_scores else median_ae

        results.append({
            'location': loc,
            'target_end_date': target_date,
            'horizon': horizon,
            'wis': wis,
            'mae': median_ae
        })

    return pd.DataFrame(results)

def main():
    ref_date = "2024-01-06"

    print("="*60)
    print(f"QUICK DIAGNOSTIC COMPARISON - {ref_date}")
    print("="*60)

    # Load target data
    targets = load_target_data()

    # Models to compare
    models = {
        'flusion': 'UMass-flusion',
        'flusion_fourier': 'UMass-flusion_fourier',
        'flusion3': 'UMass-flusion3'
    }

    results = {}
    for model_key, model_name in models.items():
        print(f"\nLoading {model_name}...")
        forecast = load_single_forecast(model_name, ref_date)
        if forecast is not None:
            wis_df = compute_wis_single_date(forecast, targets)
            results[model_key] = wis_df

            avg_wis = wis_df['wis'].mean()
            avg_mae = wis_df['mae'].mean()
            print(f"  Average WIS: {avg_wis:.2f}")
            print(f"  Average MAE: {avg_mae:.2f}")

    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print("\nOverall Performance (lower is better):")
    print(f"{'Model':<20} {'Mean WIS':>12} {'Mean MAE':>12}")
    print("-" * 45)
    for model_key in models.keys():
        if model_key in results:
            wis = results[model_key]['wis'].mean()
            mae = results[model_key]['mae'].mean()
            print(f"{model_key:<20} {wis:>12.2f} {mae:>12.2f}")

    print("\nBy Horizon:")
    for h in sorted(results['flusion']['horizon'].unique()):
        print(f"\n  Horizon {h}:")
        print(f"  {'Model':<20} {'Mean WIS':>12} {'Mean MAE':>12}")
        print("  " + "-" * 45)
        for model_key in models.keys():
            if model_key in results:
                h_df = results[model_key][results[model_key]['horizon'] == h]
                wis = h_df['wis'].mean()
                mae = h_df['mae'].mean()
                print(f"  {model_key:<20} {wis:>12.2f} {mae:>12.2f}")

if __name__ == "__main__":
    main()
