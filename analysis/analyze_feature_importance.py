"""
Analyze and compare feature importance between gbqr_3src and gbqr_3src_spatial models.
Run this after the diagnostic jobs complete.

Usage:
    python analyze_feature_importance.py --date 2025-02-15
"""

import click
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import json

sns.set_style("whitegrid")

def categorize_feature(feature_name):
    """Categorize features by type."""
    if 'wave' in feature_name.lower():
        if 'velocity' in feature_name.lower():
            return 'spatial_velocity'
        elif 'aggregate' in feature_name.lower():
            return 'spatial_aggregate'
        else:
            # Extract direction if present
            directions = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
            for direction in directions:
                if f'_{direction}_' in feature_name or f'_{direction}' in feature_name:
                    return f'spatial_wave_{direction}'
            return 'spatial_wave_other'
    elif 'lag' in feature_name.lower() or '_t-' in feature_name:
        return 'temporal_lag'
    elif 'level' in feature_name.lower() or 'y_' in feature_name:
        return 'level'
    elif any(src in feature_name.lower() for src in ['flusurvnet', 'nhsn', 'ilinet']):
        return 'data_source'
    else:
        return 'other'

def load_feature_importance(model_dir, date):
    """Load feature importance from saved artifacts."""
    artifact_path = Path(model_dir) / "artifacts" / date

    if not artifact_path.exists():
        raise FileNotFoundError(f"Artifact directory not found: {artifact_path}")

    # Look for feature importance files (format may vary based on idmodels implementation)
    # Common patterns: feature_importance.csv, feat_importance.csv, or similar
    importance_files = list(artifact_path.glob("*importance*.csv")) + \
                      list(artifact_path.glob("*importance*.parquet"))

    if not importance_files:
        print(f"Warning: No feature importance files found in {artifact_path}")
        print(f"Contents: {list(artifact_path.glob('*'))}")
        return None

    # Load the first matching file
    importance_file = importance_files[0]
    print(f"Loading: {importance_file}")

    if importance_file.suffix == '.csv':
        df = pd.read_csv(importance_file)
    elif importance_file.suffix == '.parquet':
        df = pd.read_parquet(importance_file)
    else:
        raise ValueError(f"Unsupported file format: {importance_file.suffix}")

    return df

def aggregate_importance_by_category(df):
    """Aggregate feature importance by category."""
    df['category'] = df['feature'].apply(categorize_feature)
    category_importance = df.groupby('category')['importance'].sum().sort_values(ascending=False)
    return category_importance

def compare_top_features(df_base, df_spatial, top_n=30):
    """Compare top features between models."""
    # Normalize importance scores
    df_base['importance_norm'] = df_base['importance'] / df_base['importance'].sum()
    df_spatial['importance_norm'] = df_spatial['importance'] / df_spatial['importance'].sum()

    # Get top features from each model
    top_base = df_base.nlargest(top_n, 'importance_norm')
    top_spatial = df_spatial.nlargest(top_n, 'importance_norm')

    # Merge for comparison
    comparison = pd.merge(
        df_base[['feature', 'importance_norm']],
        df_spatial[['feature', 'importance_norm']],
        on='feature',
        how='outer',
        suffixes=('_base', '_spatial')
    ).fillna(0)

    comparison['diff'] = comparison['importance_norm_spatial'] - comparison['importance_norm_base']
    comparison['category'] = comparison['feature'].apply(categorize_feature)

    return comparison, top_base, top_spatial

@click.command()
@click.option('--date', required=True, help='Reference date (YYYY-MM-DD)')
@click.option('--base_dir', default='../src/gbqr_3src', help='Base model directory')
@click.option('--spatial_dir', default='../src/gbqr_3src_spatial', help='Spatial model directory')
@click.option('--output_dir', default='outputs', help='Output directory for plots and reports')
def main(date, base_dir, spatial_dir, output_dir):
    """Analyze and compare feature importance between models."""

    output_path = Path(output_dir) / date
    output_path.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*80}")
    print(f"FEATURE IMPORTANCE ANALYSIS FOR {date}")
    print(f"{'='*80}\n")

    # Load feature importance
    print("Loading feature importance data...")
    df_base = load_feature_importance(base_dir, date)
    df_spatial = load_feature_importance(spatial_dir, date)

    if df_base is None or df_spatial is None:
        print("ERROR: Could not load feature importance data.")
        print("\nExpected artifact structure:")
        print(f"  {base_dir}/artifacts/{date}/*importance*.csv")
        print(f"  {spatial_dir}/artifacts/{date}/*importance*.csv")
        return

    # Standardize column names
    importance_col = 'importance' if 'importance' in df_base.columns else 'gain'
    if importance_col not in df_base.columns:
        importance_col = df_base.columns[1]  # Assume second column is importance

    df_base = df_base.rename(columns={importance_col: 'importance'})
    df_spatial = df_spatial.rename(columns={importance_col: 'importance'})

    # Ensure 'feature' column exists
    if 'feature' not in df_base.columns:
        df_base['feature'] = df_base.iloc[:, 0]
    if 'feature' not in df_spatial.columns:
        df_spatial['feature'] = df_spatial.iloc[:, 0]

    print(f"Base model: {len(df_base)} features")
    print(f"Spatial model: {len(df_spatial)} features\n")

    # Analyze by category
    print("Aggregating importance by feature category...")
    cat_base = aggregate_importance_by_category(df_base)
    cat_spatial = aggregate_importance_by_category(df_spatial)

    # Create category comparison
    cat_comparison = pd.DataFrame({
        'base': cat_base,
        'spatial': cat_spatial
    }).fillna(0)
    cat_comparison['diff'] = cat_comparison['spatial'] - cat_comparison['base']

    print("\n" + "="*80)
    print("FEATURE IMPORTANCE BY CATEGORY")
    print("="*80)
    print(cat_comparison.to_string())

    # Save category comparison
    cat_comparison.to_csv(output_path / "category_importance.csv")

    # Compare top features
    print("\nComparing top features...")
    comparison, top_base, top_spatial = compare_top_features(df_base, df_spatial, top_n=30)

    # Print top features unique to spatial model
    spatial_only = comparison[comparison['importance_norm_base'] == 0].nlargest(20, 'importance_norm_spatial')
    print("\n" + "="*80)
    print("TOP 20 FEATURES UNIQUE TO SPATIAL MODEL")
    print("="*80)
    print(spatial_only[['feature', 'importance_norm_spatial', 'category']].to_string(index=False))

    # Save detailed comparison
    comparison.to_csv(output_path / "feature_comparison.csv", index=False)

    # Create visualizations
    print("\nCreating visualizations...")

    # 1. Category comparison bar chart
    fig, ax = plt.subplots(figsize=(12, 6))
    cat_comparison[['base', 'spatial']].plot(kind='bar', ax=ax, width=0.8)
    ax.set_title(f'Feature Importance by Category - {date}', fontsize=14, fontweight='bold')
    ax.set_xlabel('Feature Category', fontsize=12)
    ax.set_ylabel('Total Importance', fontsize=12)
    ax.legend(['Base Model', 'Spatial Model'])
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_path / "category_comparison.png", dpi=300)
    print(f"Saved: {output_path / 'category_comparison.png'}")
    plt.close()

    # 2. Top features comparison
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    # Base model top features
    top_base_plot = top_base.nlargest(20, 'importance_norm').copy()
    top_base_plot['category'] = top_base_plot['feature'].apply(categorize_feature)
    colors_base = [sns.color_palette()[0] if cat != 'spatial' else sns.color_palette()[1]
                   for cat in top_base_plot['category'].str.startswith('spatial')]

    ax1.barh(range(len(top_base_plot)), top_base_plot['importance_norm'], color=colors_base)
    ax1.set_yticks(range(len(top_base_plot)))
    ax1.set_yticklabels(top_base_plot['feature'], fontsize=8)
    ax1.set_xlabel('Normalized Importance', fontsize=11)
    ax1.set_title('Top 20 Features - Base Model', fontsize=12, fontweight='bold')
    ax1.invert_yaxis()

    # Spatial model top features
    top_spatial_plot = top_spatial.nlargest(20, 'importance_norm').copy()
    top_spatial_plot['category'] = top_spatial_plot['feature'].apply(categorize_feature)
    colors_spatial = [sns.color_palette()[1] if 'spatial' in cat else sns.color_palette()[0]
                      for cat in top_spatial_plot['category']]

    ax2.barh(range(len(top_spatial_plot)), top_spatial_plot['importance_norm'], color=colors_spatial)
    ax2.set_yticks(range(len(top_spatial_plot)))
    ax2.set_yticklabels(top_spatial_plot['feature'], fontsize=8)
    ax2.set_xlabel('Normalized Importance', fontsize=11)
    ax2.set_title('Top 20 Features - Spatial Model', fontsize=12, fontweight='bold')
    ax2.invert_yaxis()

    plt.tight_layout()
    plt.savefig(output_path / "top_features_comparison.png", dpi=300)
    print(f"Saved: {output_path / 'top_features_comparison.png'}")
    plt.close()

    # 3. Spatial features breakdown (for spatial model)
    spatial_features = df_spatial[df_spatial['feature'].apply(
        lambda x: categorize_feature(x).startswith('spatial')
    )].copy()
    spatial_features['category'] = spatial_features['feature'].apply(categorize_feature)
    spatial_cat = spatial_features.groupby('category')['importance'].sum().sort_values(ascending=False)

    if len(spatial_cat) > 0:
        fig, ax = plt.subplots(figsize=(10, 6))
        spatial_cat.plot(kind='bar', ax=ax, color=sns.color_palette()[1])
        ax.set_title(f'Spatial Features Breakdown - {date}', fontsize=14, fontweight='bold')
        ax.set_xlabel('Spatial Feature Type', fontsize=12)
        ax.set_ylabel('Total Importance', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(output_path / "spatial_features_breakdown.png", dpi=300)
        print(f"Saved: {output_path / 'spatial_features_breakdown.png'}")
        plt.close()

    # Generate summary report
    summary = {
        'date': date,
        'base_model_features': len(df_base),
        'spatial_model_features': len(df_spatial),
        'spatial_only_features': len(spatial_features),
        'top_category_base': cat_base.index[0],
        'top_category_spatial': cat_spatial.index[0],
        'spatial_feature_importance_pct': (spatial_features['importance'].sum() /
                                           df_spatial['importance'].sum() * 100),
    }

    with open(output_path / "summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Spatial features account for {summary['spatial_feature_importance_pct']:.1f}% "
          f"of total importance in spatial model")
    print(f"\nAll outputs saved to: {output_path}")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()
