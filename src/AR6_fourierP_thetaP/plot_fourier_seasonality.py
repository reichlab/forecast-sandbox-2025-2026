"""
Plot Fourier seasonality terms from AR6_fourierP_thetaP model.

This script runs the model and visualizes the seasonal pattern captured by
the pooled Fourier regression terms (shared across all locations).
"""

import sys
from pathlib import Path
# Add sarix package to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent / 'sarix' / 'src'))

import datetime
from dateutil import relativedelta
from types import SimpleNamespace
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from idmodels.sarix import SARIXFourierModel


def main():
    """Run model and plot pooled Fourier seasonality curve."""

    # Use the date specified in the README
    today_date = datetime.date(2024, 1, 6)
    reference_date = today_date + relativedelta.relativedelta(weekday=5)

    print("Setting up model configuration...")
    model_config = SimpleNamespace(
        model_class = "sarix",
        model_name = "AR6_fourierP_thetaP",
        sources = ["nhsn"],
        fit_locations_separately = False,
        p = 6,
        P = 0,
        d = 0,
        D = 0,
        season_period = 1,
        power_transform = "4rt",
        theta_pooling="shared",
        sigma_pooling="none",
        fourier_pooling="shared",  # Pooled Fourier coefficients
        fourier_K = 2,  # Number of Fourier harmonic pairs
        x = []
    )

    run_config = SimpleNamespace(
        disease="flu",
        ref_date=reference_date,
        output_root=Path("../../model-output/"),
        artifact_store_root=None,
        max_horizon=4,
        locations=["US", "01", "02", "04", "05", "06", "08", "09", "10", "11",
                   "12", "13", "15", "16", "17", "18", "19", "20", "21", "22",
                   "23", "24", "25", "26", "27", "28", "29", "30", "31", "32",
                   "33", "34", "35", "36", "37", "38", "39", "40", "41", "42",
                   "44", "45", "46", "47", "48", "49", "50", "51", "53", "54",
                   "55", "56", "72"],
        q_levels = [0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
                    0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70,
                    0.75, 0.80, 0.85, 0.90, 0.95, 0.975, 0.99],
        q_labels = ['0.01', '0.025', '0.05', '0.1', '0.15', '0.2',
                    '0.25', '0.3', '0.35', '0.4', '0.45', '0.5',
                    '0.55', '0.6', '0.65', '0.7', '0.75', '0.8',
                    '0.85', '0.9', '0.95', '0.975', '0.99'],
        num_warmup = 500,   # Reduced for faster testing
        num_samples = 500,  # Reduced for faster testing
        num_chains = 1
    )

    print("Initializing and running model...")
    print(f"  Reference date: {reference_date}")
    print(f"  Number of locations: {len(run_config.locations)}")
    print(f"  Fourier harmonics (K): {model_config.fourier_K}")
    print(f"  Fourier pooling: {model_config.fourier_pooling}")

    # Load data and fit model to extract coefficients
    from iddata.loader import DiseaseDataLoader
    from iddata.utils import get_holidays
    from sarix import sarix

    # Load data (same as in SARIXFourierModel.run())
    print("\nLoading data...")
    fdl = DiseaseDataLoader()
    df = fdl.load_data(
        nhsn_kwargs={"as_of": run_config.ref_date, "disease": run_config.disease},
        sources=model_config.sources,
        power_transform=model_config.power_transform
    )
    df = df.loc[df["location"].isin(run_config.locations)]

    # Process data
    print("Processing data...")
    df = df.merge(
        get_holidays() \
            .query("holiday == 'Christmas Day'") \
            .drop(columns=["holiday", "date"]) \
            .rename(columns={"season_week": "xmas_week"}),
        how="left",
        on="season") \
    .assign(delta_xmas = lambda x: x["season_week"] - x["xmas_week"])
    df["xmas_spike"] = np.maximum(3 - np.abs(df["delta_xmas"]), 0)

    xy_colnames = model_config.x + ["inc_trans_cs"]
    df = df.query("wk_end_date >= '2022-10-01'").interpolate()
    batched_xy = df[xy_colnames].values.reshape(
        len(df["location"].unique()), -1, len(xy_colnames)
    )

    # Extract day-of-year
    day_of_year = df.groupby("location")["wk_end_date"].apply(
        lambda x: x.dt.dayofyear.values
    ).iloc[0]

    # Fit model
    print("Fitting SARIX model with pooled Fourier terms...")
    sarix_fit = sarix.SARIX(
        xy = batched_xy,
        p = model_config.p,
        d = model_config.d,
        P = model_config.P,
        D = model_config.D,
        season_period = model_config.season_period,
        transform="none",
        theta_pooling=model_config.theta_pooling,
        sigma_pooling=model_config.sigma_pooling,
        fourier_pooling=model_config.fourier_pooling,
        forecast_horizon = run_config.max_horizon,
        num_warmup = run_config.num_warmup,
        num_samples = run_config.num_samples,
        num_chains = run_config.num_chains,
        day_of_year = day_of_year,
        fourier_K = model_config.fourier_K
    )

    print("\nExtracting and plotting pooled Fourier coefficients...")

    # Extract Fourier coefficients
    # For pooled model: Shape should be (num_samples, n_vars, 2*K)
    # Since coefficients are shared across locations
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
    # fourier_beta: (num_samples, n_vars, 2*K) for pooled model
    # We want the coefficients for the response variable (inc_trans_cs), which is the last variable
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
    print(f"\nCreating plot for pooled seasonality...")

    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot credible interval
    ax.fill_between(day_grid, seasonality_lower, seasonality_upper,
                    alpha=0.3, color='blue', label='95% CI')

    # Plot median
    ax.plot(day_grid, seasonality_median, 'b-', linewidth=2, label='Median')

    # Add zero line
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)

    # Format axes
    ax.set_xlabel('Month', fontsize=12)
    ax.set_ylabel('Seasonal effect (pooled across locations)', fontsize=12)
    ax.set_title(f'Pooled Fourier Seasonal Pattern (K={K} harmonics)\nModel: AR6_fourierP_thetaP | Reference Date: {reference_date}',
                 fontsize=14, fontweight='bold')
    ax.set_xticks(month_starts)
    ax.set_xticklabels(month_labels, fontsize=10)
    ax.tick_params(axis='y', labelsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, loc='best')

    # Save plot
    output_path = Path('fourier_seasonality_pooled.png')
    print(f"\nSaving plot to: {output_path.absolute()}")
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"✓ Plot saved successfully")

    # Print summary statistics
    print("\n" + "="*60)
    print("Summary Statistics:")
    print(f"  Median seasonal effect range: [{seasonality_median.min():.4f}, {seasonality_median.max():.4f}]")
    print(f"  Peak seasonality (day of year): {day_grid[np.argmax(seasonality_median)]}")
    print(f"  Trough seasonality (day of year): {day_grid[np.argmin(seasonality_median)]}")
    print("="*60)

    print("\nDone! Generated plot:")
    print(f"  {output_path.absolute()}")
    print("="*60)

    return seasonality_median, seasonality_lower, seasonality_upper


if __name__ == "__main__":
    seasonality_median, seasonality_lower, seasonality_upper = main()
