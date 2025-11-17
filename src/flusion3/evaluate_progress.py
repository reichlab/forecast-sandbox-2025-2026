"""
Evaluate completed flusion3 forecasts against baseline models.
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

def load_forecasts(model_name, ref_dates):
    """Load forecast data for specific reference dates."""
    dfs = []
    for ref_date in ref_dates:
        file_path = MODEL_OUTPUT_DIR / model_name / f"{ref_date}-{model_name}.csv"
        if file_path.exists():
            df = pd.read_csv(file_path)
            df['reference_date'] = pd.to_datetime(ref_date)
            df['target_end_date'] = pd.to_datetime(df['target_end_date'])
            dfs.append(df)

    if dfs:
        return pd.concat(dfs, ignore_index=True)
    return None

def compute_wis(forecasts, targets):
    """Compute WIS for forecasts."""
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

    for (loc, ref_date, target_date, horizon), group in merged.groupby(['location', 'reference_date', 'target_end_date', 'horizon']):
        obs = group['value_observed'].iloc[0]

        # Median AE
        median_forecast = group[np.isclose(group['output_type_id'], 0.5)]['value_forecast'].values
        if len(median_forecast) == 0:
            continue
        median_ae = abs(obs - median_forecast[0])

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
                width = upper_val[0] - lower_val[0]
                if obs < lower_val[0]:
                    penalty = (2 / alpha) * (lower_val[0] - obs)
                elif obs > upper_val[0]:
                    penalty = (2 / alpha) * (obs - upper_val[0])
                else:
                    penalty = 0
                interval_scores.append(width + penalty)

        wis = 0.5 * median_ae + 0.5 * np.mean(interval_scores) if interval_scores else median_ae

        results.append({
            'location': loc,
            'reference_date': ref_date,
            'target_end_date': target_date,
            'horizon': horizon,
            'wis': wis,
            'mae': median_ae
        })

    return pd.DataFrame(results)

def compute_coverage(forecasts, targets, alpha=0.95):
    """Compute prediction interval coverage."""
    lower_q = (1 - alpha) / 2
    upper_q = 1 - lower_q

    quantile_forecasts = forecasts[forecasts['output_type'] == 'quantile'].copy()
    quantile_forecasts['output_type_id'] = quantile_forecasts['output_type_id'].astype(float)

    lower = quantile_forecasts[np.isclose(quantile_forecasts['output_type_id'], lower_q)][
        ['location', 'reference_date', 'target_end_date', 'horizon', 'value']
    ].rename(columns={'value': 'lower'})

    upper = quantile_forecasts[np.isclose(quantile_forecasts['output_type_id'], upper_q)][
        ['location', 'reference_date', 'target_end_date', 'horizon', 'value']
    ].rename(columns={'value': 'upper'})

    intervals = lower.merge(upper, on=['location', 'reference_date', 'target_end_date', 'horizon'])

    coverage = intervals.merge(
        targets[['location', 'date', 'value']],
        left_on=['location', 'target_end_date'],
        right_on=['location', 'date'],
        how='inner'
    )

    coverage['covered'] = (coverage['value'] >= coverage['lower']) & (coverage['value'] <= coverage['upper'])
    coverage['interval_width'] = coverage['upper'] - coverage['lower']

    return coverage

def main():
    # Get completed flusion3 forecasts
    flusion3_files = sorted((MODEL_OUTPUT_DIR / "UMass-flusion3").glob("*.csv"))
    ref_dates = [f.stem.split("-UMass-flusion3")[0] for f in flusion3_files]

    print("="*70)
    print("FLUSION3 PROGRESS EVALUATION")
    print("="*70)
    print(f"\nEvaluating {len(ref_dates)} completed forecasts:")
    for date in ref_dates:
        print(f"  - {date}")

    # Load target data
    print("\nLoading target data...")
    targets = load_target_data()

    # Models to compare
    models = {
        'flusion': 'UMass-flusion',
        'flusion_fourier': 'UMass-flusion_fourier',
        'flusion3': 'UMass-flusion3'
    }

    # Load forecasts
    forecasts = {}
    wis_results = {}
    coverage_results = {}

    for model_key, model_name in models.items():
        print(f"\nLoading {model_name}...")
        fc = load_forecasts(model_name, ref_dates)
        if fc is not None:
            forecasts[model_key] = fc
            wis_results[model_key] = compute_wis(fc, targets)
            coverage_results[model_key] = compute_coverage(fc, targets)

    print("\n" + "="*70)
    print("OVERALL PERFORMANCE (lower is better)")
    print("="*70)
    print(f"\n{'Model':<20} {'Mean WIS':>12} {'Median WIS':>12} {'Mean MAE':>12} {'Coverage':>10}")
    print("-" * 70)

    for model_key in ['flusion', 'flusion_fourier', 'flusion3']:
        if model_key in wis_results:
            mean_wis = wis_results[model_key]['wis'].mean()
            median_wis = wis_results[model_key]['wis'].median()
            mean_mae = wis_results[model_key]['mae'].mean()
            coverage = coverage_results[model_key]['covered'].mean()
            print(f"{model_key:<20} {mean_wis:>12.2f} {median_wis:>12.2f} {mean_mae:>12.2f} {coverage:>10.3f}")

    print("\n" + "="*70)
    print("PERFORMANCE BY HORIZON")
    print("="*70)

    for h in sorted(wis_results['flusion']['horizon'].unique()):
        print(f"\nHorizon {h}:")
        print(f"{'Model':<20} {'Mean WIS':>12} {'Mean MAE':>12} {'Coverage':>10}")
        print("-" * 55)
        for model_key in ['flusion', 'flusion_fourier', 'flusion3']:
            if model_key in wis_results:
                h_wis = wis_results[model_key][wis_results[model_key]['horizon'] == h]
                h_cov = coverage_results[model_key][coverage_results[model_key]['horizon'] == h]
                if len(h_wis) > 0:
                    mean_wis = h_wis['wis'].mean()
                    mean_mae = h_wis['mae'].mean()
                    coverage = h_cov['covered'].mean() if len(h_cov) > 0 else 0
                    print(f"{model_key:<20} {mean_wis:>12.2f} {mean_mae:>12.2f} {coverage:>10.3f}")

    print("\n" + "="*70)
    print("KEY FINDINGS")
    print("="*70)

    # Calculate relative performance
    if 'flusion' in wis_results and 'flusion3' in wis_results:
        flusion_wis = wis_results['flusion']['wis'].mean()
        flusion3_wis = wis_results['flusion3']['wis'].mean()
        diff_pct = ((flusion3_wis - flusion_wis) / flusion_wis) * 100

        if diff_pct < 0:
            print(f"\n✓ flusion3 is {abs(diff_pct):.1f}% BETTER than flusion (lower WIS)")
        else:
            print(f"\n✗ flusion3 is {diff_pct:.1f}% WORSE than flusion (higher WIS)")

        print(f"  flusion WIS:  {flusion_wis:.2f}")
        print(f"  flusion3 WIS: {flusion3_wis:.2f}")

if __name__ == "__main__":
    main()
