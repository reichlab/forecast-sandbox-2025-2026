"""
Analyze geographic patterns in forecast differences between models.
Identifies which states benefit most from spatial features.

Usage:
    python analyze_geographic_differences.py --date 2025-02-15
"""

import click
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")

# State FIPS to name mapping
STATE_NAMES = {
    "01": "Alabama", "02": "Alaska", "04": "Arizona", "05": "Arkansas",
    "06": "California", "08": "Colorado", "09": "Connecticut", "10": "Delaware",
    "11": "DC", "12": "Florida", "13": "Georgia", "15": "Hawaii",
    "16": "Idaho", "17": "Illinois", "18": "Indiana", "19": "Iowa",
    "20": "Kansas", "21": "Kentucky", "22": "Louisiana", "23": "Maine",
    "24": "Maryland", "25": "Massachusetts", "26": "Michigan", "27": "Minnesota",
    "28": "Mississippi", "29": "Missouri", "30": "Montana", "31": "Nebraska",
    "32": "Nevada", "33": "New Hampshire", "34": "New Jersey", "35": "New Mexico",
    "36": "New York", "37": "North Carolina", "38": "North Dakota", "39": "Ohio",
    "40": "Oklahoma", "41": "Oregon", "42": "Pennsylvania", "44": "Rhode Island",
    "45": "South Carolina", "46": "South Dakota", "47": "Tennessee", "48": "Texas",
    "49": "Utah", "50": "Vermont", "51": "Virginia", "53": "Washington",
    "54": "West Virginia", "55": "Wisconsin", "56": "Wyoming", "72": "Puerto Rico"
}

# Regional groupings
REGIONS = {
    "Northeast": ["09", "23", "25", "33", "44", "50", "34", "36", "42"],
    "Midwest": ["17", "18", "26", "39", "55", "19", "20", "27", "29", "31", "38", "46"],
    "South": ["10", "11", "12", "13", "24", "37", "45", "51", "54", "01", "21", "28",
              "47", "05", "22", "40", "48", "72"],
    "West": ["04", "08", "16", "30", "32", "35", "49", "56", "02", "06", "15", "41",
             "53"]
}

def get_region(state_fips):
    """Get region for a state FIPS code."""
    for region, states in REGIONS.items():
        if state_fips in states:
            return region
    return "Other"

def load_forecast(model_name, reference_date):
    """Load forecast file."""
    file_path = Path(f"../model-output/UMass-{model_name}/{reference_date}-UMass-{model_name}.csv")
    if file_path.exists():
        df = pd.read_csv(file_path, dtype={'location': str, 'output_type_id': str})
        return df
    return None

def compare_forecasts(df_base, df_spatial, focus_quantile='0.5'):
    """Compare forecasts between models."""
    # Focus on median forecasts
    df_base_med = df_base[df_base['output_type_id'] == focus_quantile].copy()
    df_spatial_med = df_spatial[df_spatial['output_type_id'] == focus_quantile].copy()

    # Merge
    merged = df_base_med.merge(
        df_spatial_med,
        on=['location', 'horizon', 'target', 'target_end_date'],
        suffixes=('_base', '_spatial')
    )

    # Compute differences
    merged['abs_diff'] = np.abs(merged['value_spatial'] - merged['value_base'])
    merged['rel_diff'] = merged['abs_diff'] / (merged['value_base'] + 1)
    merged['diff'] = merged['value_spatial'] - merged['value_base']

    # Add state names and regions
    merged['state_name'] = merged['location'].map(STATE_NAMES)
    merged['region'] = merged['location'].apply(get_region)

    return merged

@click.command()
@click.option('--date', required=True, help='Reference date (YYYY-MM-DD)')
@click.option('--output_dir', default='outputs', help='Output directory')
def main(date, output_dir):
    """Analyze geographic patterns in forecast differences."""

    output_path = Path(output_dir) / date
    output_path.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*80}")
    print(f"GEOGRAPHIC ANALYSIS FOR {date}")
    print(f"{'='*80}\n")

    # Load forecasts
    print("Loading forecasts...")
    df_base = load_forecast("gbqr_3src", date)
    df_spatial = load_forecast("gbqr_3src_spatial", date)

    if df_base is None or df_spatial is None:
        print(f"ERROR: Could not load forecasts for {date}")
        return

    print(f"Base model locations: {df_base['location'].nunique()}")
    print(f"Spatial model locations: {df_spatial['location'].nunique()}\n")

    # Compare median forecasts
    print("Comparing median forecasts...")
    comparison = compare_forecasts(df_base, df_spatial)

    # Overall statistics
    print("\n" + "="*80)
    print("OVERALL STATISTICS")
    print("="*80)
    print(f"Mean absolute difference: {comparison['abs_diff'].mean():.2f}")
    print(f"Median absolute difference: {comparison['abs_diff'].median():.2f}")
    print(f"Max absolute difference: {comparison['abs_diff'].max():.2f}")
    print(f"Mean relative difference: {comparison['rel_diff'].mean():.2%}")

    # By horizon
    print("\n" + "="*80)
    print("DIFFERENCES BY HORIZON")
    print("="*80)
    horizon_stats = comparison.groupby('horizon').agg({
        'abs_diff': ['mean', 'median', 'max'],
        'rel_diff': 'mean'
    }).round(2)
    print(horizon_stats)

    # By region
    print("\n" + "="*80)
    print("DIFFERENCES BY REGION")
    print("="*80)
    region_stats = comparison.groupby('region').agg({
        'abs_diff': ['mean', 'median', 'max'],
        'rel_diff': 'mean',
        'location': 'count'
    }).round(2)
    print(region_stats)

    # Top states with largest differences
    state_diffs = comparison.groupby(['location', 'state_name', 'region'])['abs_diff'].mean().reset_index()
    state_diffs = state_diffs.sort_values('abs_diff', ascending=False)

    print("\n" + "="*80)
    print("TOP 15 STATES BY MEAN ABSOLUTE DIFFERENCE")
    print("="*80)
    print(state_diffs.head(15).to_string(index=False))

    print("\n" + "="*80)
    print("BOTTOM 15 STATES BY MEAN ABSOLUTE DIFFERENCE")
    print("="*80)
    print(state_diffs.tail(15).to_string(index=False))

    # Save detailed results
    comparison.to_csv(output_path / "forecast_comparison_detailed.csv", index=False)
    state_diffs.to_csv(output_path / "state_differences.csv", index=False)
    horizon_stats.to_csv(output_path / "horizon_differences.csv")
    region_stats.to_csv(output_path / "region_differences.csv")

    # Create visualizations
    print("\nCreating visualizations...")

    # 1. State-level differences bar chart
    fig, ax = plt.subplots(figsize=(14, 10))
    top_states = state_diffs.head(30)
    colors = [sns.color_palette()[i] for i in [REGIONS[region] for region in top_states['region']].index(top_states['region'].iloc[0])]

    # Color by region
    region_colors = {'Northeast': 0, 'Midwest': 1, 'South': 2, 'West': 3}
    colors = [sns.color_palette()[region_colors[r]] for r in top_states['region']]

    ax.barh(range(len(top_states)), top_states['abs_diff'], color=colors)
    ax.set_yticks(range(len(top_states)))
    ax.set_yticklabels([f"{row['state_name']} ({row['location']})"
                        for _, row in top_states.iterrows()], fontsize=9)
    ax.set_xlabel('Mean Absolute Difference in Median Forecast', fontsize=12)
    ax.set_title(f'Top 30 States by Forecast Difference - {date}', fontsize=14, fontweight='bold')
    ax.invert_yaxis()

    # Add legend for regions
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=sns.color_palette()[i], label=region)
                      for region, i in region_colors.items()]
    ax.legend(handles=legend_elements, loc='lower right')

    plt.tight_layout()
    plt.savefig(output_path / "state_differences.png", dpi=300)
    print(f"Saved: {output_path / 'state_differences.png'}")
    plt.close()

    # 2. Differences by horizon
    fig, ax = plt.subplots(figsize=(10, 6))
    horizon_means = comparison.groupby('horizon')['abs_diff'].mean()
    ax.bar(horizon_means.index, horizon_means.values, color=sns.color_palette()[0])
    ax.set_xlabel('Forecast Horizon', fontsize=12)
    ax.set_ylabel('Mean Absolute Difference', fontsize=12)
    ax.set_title(f'Forecast Differences by Horizon - {date}', fontsize=14, fontweight='bold')
    ax.set_xticks(horizon_means.index)
    plt.tight_layout()
    plt.savefig(output_path / "horizon_differences.png", dpi=300)
    print(f"Saved: {output_path / 'horizon_differences.png'}")
    plt.close()

    # 3. Regional comparison
    fig, ax = plt.subplots(figsize=(10, 6))
    region_means = comparison.groupby('region')['abs_diff'].mean().sort_values(ascending=False)
    ax.bar(range(len(region_means)), region_means.values,
           color=[sns.color_palette()[region_colors[r]] for r in region_means.index])
    ax.set_xticks(range(len(region_means)))
    ax.set_xticklabels(region_means.index, fontsize=12)
    ax.set_ylabel('Mean Absolute Difference', fontsize=12)
    ax.set_title(f'Forecast Differences by Region - {date}', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path / "region_differences.png", dpi=300)
    print(f"Saved: {output_path / 'region_differences.png'}")
    plt.close()

    # 4. Scatter plot: base vs spatial predictions for top differing states
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes = axes.flatten()

    top_6_states = state_diffs.head(6)
    for idx, (_, state_row) in enumerate(top_6_states.iterrows()):
        ax = axes[idx]
        state_data = comparison[comparison['location'] == state_row['location']]

        ax.scatter(state_data['value_base'], state_data['value_spatial'],
                  alpha=0.6, s=50, c=state_data['horizon'],
                  cmap='viridis')

        # Add diagonal line
        min_val = min(state_data['value_base'].min(), state_data['value_spatial'].min())
        max_val = max(state_data['value_base'].max(), state_data['value_spatial'].max())
        ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.3, linewidth=1)

        ax.set_xlabel('Base Model Forecast', fontsize=10)
        ax.set_ylabel('Spatial Model Forecast', fontsize=10)
        ax.set_title(f"{state_row['state_name']}\n(Mean diff: {state_row['abs_diff']:.1f})",
                    fontsize=11, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_path / "state_scatter_comparison.png", dpi=300)
    print(f"Saved: {output_path / 'state_scatter_comparison.png'}")
    plt.close()

    # 5. Heatmap of differences by state and horizon
    heatmap_data = comparison.pivot_table(
        values='abs_diff',
        index='state_name',
        columns='horizon',
        aggfunc='mean'
    )

    # Sort by mean difference
    heatmap_data['mean'] = heatmap_data.mean(axis=1)
    heatmap_data = heatmap_data.sort_values('mean', ascending=False).drop('mean', axis=1)
    heatmap_data = heatmap_data.head(30)  # Top 30 states

    fig, ax = plt.subplots(figsize=(10, 14))
    sns.heatmap(heatmap_data, annot=True, fmt='.1f', cmap='YlOrRd',
                cbar_kws={'label': 'Absolute Difference'}, ax=ax)
    ax.set_title(f'Forecast Differences by State and Horizon - {date}',
                fontsize=14, fontweight='bold')
    ax.set_xlabel('Forecast Horizon', fontsize=12)
    ax.set_ylabel('State', fontsize=12)
    plt.tight_layout()
    plt.savefig(output_path / "state_horizon_heatmap.png", dpi=300)
    print(f"Saved: {output_path / 'state_horizon_heatmap.png'}")
    plt.close()

    print(f"\nAll outputs saved to: {output_path}")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()
