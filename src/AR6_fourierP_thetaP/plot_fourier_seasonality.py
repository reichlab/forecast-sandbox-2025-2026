"""
Plot Fourier seasonality terms from AR6_fourierP_thetaP model.

This script runs the model and visualizes the seasonal pattern captured by
the pooled Fourier regression terms (shared across all locations).
"""

import datetime
from pathlib import Path

from dateutil import relativedelta
from iddata.ancillary.population import PopulationData
from iddata.enums import Disease, SourceType
from iddata.loader import DiseaseDataLoader
from idmodels.config import (
    PoolingStrategy,
    PowerTransform,
    RunConfig,
    SARIXFourierModelConfig,
)
from idmodels.sarix import SARIXFourierModel
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sarix import sarix


def main():
    """Run model and plot pooled Fourier seasonality curve."""

    # Use the date specified in the README
    today_date = datetime.date(2024, 1, 6)
    reference_date = today_date + relativedelta.relativedelta(weekday=5)

    print("Setting up model configuration...")
    model_config = SARIXFourierModelConfig(
        model_name="AR6_fourierP_thetaP",
        main_source=SourceType.NHSN,
        fit_locations_separately=False,
        p=6,
        P=0,
        d=0,
        D=0,
        season_period=1,
        power_transform=PowerTransform.FOURTH_ROOT,
        theta_pooling=PoolingStrategy.SHARED,
        sigma_pooling=PoolingStrategy.NONE,
        fourier_pooling=PoolingStrategy.SHARED,  # Pooled Fourier coefficients
        fourier_K=2,  # Number of Fourier harmonic pairs
        x=[],
        num_warmup=500,   # Reduced for faster testing
        num_samples=500,  # Reduced for faster testing
        num_chains=1,
    )

    run_config = RunConfig(
        disease=Disease.FLU,
        ref_date=reference_date,
        output_root=Path("../../model-output/"),
        artifact_store_root=None,
        max_horizon=4,
        states=[
            "US", "01", "02", "04", "05", "06", "08", "09", "10", "11",
            "12", "13", "15", "16", "17", "18", "19", "20", "21", "22",
            "23", "24", "25", "26", "27", "28", "29", "30", "31", "32",
            "33", "34", "35", "36", "37", "38", "39", "40", "41", "42",
            "44", "45", "46", "47", "48", "49", "50", "51", "53", "54",
            "55", "56", "72",
        ],
        hsas=[],
        q_levels=[
            0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
            0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70,
            0.75, 0.80, 0.85, 0.90, 0.95, 0.975, 0.99,
        ],
        q_labels=[
            "0.01", "0.025", "0.05", "0.1", "0.15", "0.2",
            "0.25", "0.3", "0.35", "0.4", "0.45", "0.5",
            "0.55", "0.6", "0.65", "0.7", "0.75", "0.8",
            "0.85", "0.9", "0.95", "0.975", "0.99",
        ],
    )

    print("Initializing model and loading data...")
    print(f"  Reference date: {reference_date}")
    print(f"  Number of locations: {len(run_config.states)}")
    print(f"  Fourier harmonics (K): {model_config.fourier_K}")
    print(f"  Fourier pooling: {model_config.fourier_pooling}")

    model = SARIXFourierModel(model_config)
    sources = model._build_sources(run_config)
    df = DiseaseDataLoader().load(sources=sources, as_of=run_config.ref_date, ancillary=[PopulationData()])
    df = model._filter_locations(df, run_config)
    df["unique_id"] = df["agg_level"] + df["location"]

    transform = model._build_transform()
    df = transform.apply(df)

    pipeline = model._build_feature_pipeline(run_config)
    df, feat_names = pipeline.apply(df)

    xy_colnames = model_config.x + ["inc_trans_cs"]
    df = df.query("wk_end_date >= '2022-10-01'").interpolate()
    batched_xy = df[xy_colnames].values.reshape(
        len(df["unique_id"].unique()), -1, len(xy_colnames)
    )

    extra_params = model._get_extra_sarix_params(df)

    print("Fitting SARIX model with pooled Fourier terms...")
    sarix_fit = sarix.SARIX(
        xy=batched_xy,
        p=model_config.p,
        d=model_config.d,
        P=model_config.P,
        D=model_config.D,
        season_period=model_config.season_period,
        transform="none",
        theta_pooling=model_config.theta_pooling,
        sigma_pooling=model_config.sigma_pooling,
        forecast_horizon=run_config.max_horizon,
        num_warmup=model_config.num_warmup,
        num_samples=model_config.num_samples,
        num_chains=model_config.num_chains,
        **extra_params,
    )

    print("\nExtracting and plotting pooled Fourier coefficients...")
    fourier_beta = sarix_fit.samples['fourier_beta']
    print(f"Fourier beta shape: {fourier_beta.shape}")

    K = model_config.fourier_K

    # Create day-of-year grid for plotting (full year)
    day_grid = np.arange(0, 365, 1)
    t_normalized = day_grid / 365.25

    # Calculate Fourier features for full year
    fourier_features_grid = []
    for k in range(1, K + 1):
        fourier_features_grid.append(np.sin(2 * np.pi * k * t_normalized))
        fourier_features_grid.append(np.cos(2 * np.pi * k * t_normalized))
    fourier_features_grid = np.stack(fourier_features_grid, axis=-1)  # (365, 2*K)

    # Calculate pooled seasonality curve
    fourier_beta_y = fourier_beta[:, -1, :]  # (num_samples, 2*K)

    # Compute smooth median curve from median coefficients
    beta_median = np.median(fourier_beta_y, axis=0)  # (2*K,)
    seasonality_median = fourier_features_grid @ beta_median  # (365,)

    # For credible intervals: compute all curves, then take pointwise percentiles
    seasonality_samples = fourier_features_grid @ fourier_beta_y.T  # (365, num_samples)
    seasonality_lower = np.percentile(seasonality_samples, 2.5, axis=1)  # (365,)
    seasonality_upper = np.percentile(seasonality_samples, 97.5, axis=1)  # (365,)

    # Convert day-of-year to month labels for x-axis
    dates = pd.date_range('2024-01-01', periods=365, freq='D')
    month_starts = [i for i, d in enumerate(dates) if d.day == 1]
    month_labels = [dates[i].strftime('%b') for i in month_starts]

    # Create plot
    print("\nCreating plot for pooled seasonality...")
    fig, ax = plt.subplots(figsize=(12, 6))

    ax.fill_between(day_grid, seasonality_lower, seasonality_upper,
                    alpha=0.3, color='blue', label='95% CI')
    ax.plot(day_grid, seasonality_median, 'b-', linewidth=2, label='Median')
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)

    ax.set_xlabel('Month', fontsize=12)
    ax.set_ylabel('Seasonal effect (pooled across locations)', fontsize=12)
    ax.set_title(
        f'Pooled Fourier Seasonal Pattern (K={K} harmonics)\n'
        f'Model: AR6_fourierP_thetaP | Reference Date: {reference_date}',
        fontsize=14,
        fontweight='bold',
    )
    ax.set_xticks(month_starts)
    ax.set_xticklabels(month_labels, fontsize=10)
    ax.tick_params(axis='y', labelsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, loc='best')

    output_path = Path('fourier_seasonality_pooled.png')
    print(f"\nSaving plot to: {output_path.absolute()}")
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print("✓ Plot saved successfully")

    print("\n" + "=" * 60)
    print("Summary Statistics:")
    print(f"  Median seasonal effect range: [{seasonality_median.min():.4f}, {seasonality_median.max():.4f}]")
    print(f"  Peak seasonality (day of year): {day_grid[np.argmax(seasonality_median)]}")
    print(f"  Trough seasonality (day of year): {day_grid[np.argmin(seasonality_median)]}")
    print("=" * 60)

    return seasonality_median, seasonality_lower, seasonality_upper


if __name__ == "__main__":
    seasonality_median, seasonality_lower, seasonality_upper = main()
