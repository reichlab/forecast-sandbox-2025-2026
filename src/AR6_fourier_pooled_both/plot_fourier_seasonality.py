"""
Plot Fourier seasonality terms from AR6_fourier_pooled_both model.

This script runs the model and visualizes the seasonal patterns captured by
the Fourier regression terms for each state.
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
    """Run model and plot Fourier seasonality curves."""

    # Use a fixed date for reproducibility
    today_date = datetime.date(2025, 5, 9)
    reference_date = today_date + relativedelta.relativedelta(weekday=5)

    print("Setting up model configuration...")
    model_config = SimpleNamespace(
        model_class = "sarix",
        model_name = "AR6_fourier_pooled_both",
        sources = ["nhsn"],
        fit_locations_separately = False,
        p = 6,
        P = 0,
        d = 0,
        D = 0,
        season_period = 1,
        power_transform = "4rt",
        theta_pooling="shared",
        sigma_pooling="shared",
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

    # Skip the model.run() call - we refit the model below to extract coefficients
    # model = SARIXFourierModel(model_config)
    # model.run(run_config)

    print("\nExtracting Fourier coefficients from fitted model...")
    # Access the underlying SARIX object that was created during model.run()
    # We need to recreate it to get access to the samples

    # Actually, we need to modify the approach - let me create a custom run
    # that gives us access to the fitted model object
    from iddata.loader import DiseaseDataLoader
    from iddata.utils import get_holidays
    from sarixfourier import sarix_fourier

    # Load data (same as in SARIXFourierModel.run())
    fdl = DiseaseDataLoader()
    df = fdl.load_data(
        nhsn_kwargs={"as_of": run_config.ref_date, "disease": run_config.disease},
        sources=model_config.sources,
        power_transform=model_config.power_transform
    )
    df = df.loc[df["location"].isin(run_config.locations)]

    # Process data
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
    print("Fitting SARIX model with Fourier terms...")
    sarix_fit = sarix_fourier.SARIX(
        xy = batched_xy,
        p = model_config.p,
        d = model_config.d,
        P = model_config.P,
        D = model_config.D,
        season_period = model_config.season_period,
        transform="none",
        theta_pooling=model_config.theta_pooling,
        sigma_pooling=model_config.sigma_pooling,
        forecast_horizon = run_config.max_horizon,
        num_warmup = run_config.num_warmup,
        num_samples = run_config.num_samples,
        num_chains = run_config.num_chains,
        day_of_year = day_of_year,
        fourier_K = model_config.fourier_K
    )

    print("\nExtracting and plotting Fourier coefficients...")

    # Extract Fourier coefficients
    # Shape: (num_samples, n_locations, n_vars, 2*K)
    fourier_beta = sarix_fit.samples['fourier_beta']
    print(f"Fourier beta shape: {fourier_beta.shape}")

    # Get location names
    locations = run_config.locations
    n_locations = len(locations)
    K = model_config.fourier_K

    # Load location names from auxiliary data
    locations_df = pd.read_csv('/Users/nick/Documents/research-versioned/FluSight-forecast-hub/auxiliary-data/locations.csv')
    location_name_map = dict(zip(locations_df['location'], locations_df['location_name']))
    location_names = [location_name_map.get(loc, loc) for loc in locations]

    # Create day-of-year grid for plotting (full year)
    day_grid = np.arange(0, 365, 1)
    t_normalized = day_grid / 365.25

    # Calculate Fourier features for full year
    fourier_features_grid = []
    for k in range(1, K + 1):
        fourier_features_grid.append(np.sin(2 * np.pi * k * t_normalized))
        fourier_features_grid.append(np.cos(2 * np.pi * k * t_normalized))
    fourier_features_grid = np.stack(fourier_features_grid, axis=-1)  # (365, 2*K)

    # Calculate seasonality curves for each location
    # fourier_beta: (num_samples, n_locations, 1, 2*K)
    # Strategy:
    # - Median: use median coefficients for smooth curve
    # - Credible intervals: use all samples for correct uncertainty

    # Squeeze out the variable dimension first
    fourier_beta_reshaped = fourier_beta.squeeze(axis=2)  # (num_samples, n_locations, 2*K)

    # Compute smooth median curve from median coefficients
    beta_median = np.median(fourier_beta_reshaped, axis=0)  # (n_locations, 2*K)
    seasonality_median = np.einsum('tk,lk->lt', fourier_features_grid, beta_median)  # (n_locations, 365)

    # For credible intervals: compute all curves, then take pointwise percentiles
    # This gives correct uncertainty representation
    # Shape: (num_samples, n_locations, 365)
    seasonality_samples = np.einsum('tk,slk->slt', fourier_features_grid, fourier_beta_reshaped)
    seasonality_lower = np.percentile(seasonality_samples, 2.5, axis=0)  # (n_locations, 365)
    seasonality_upper = np.percentile(seasonality_samples, 97.5, axis=0)  # (n_locations, 365)

    # Convert day-of-year to month labels for x-axis
    dates = pd.date_range('2024-01-01', periods=365, freq='D')
    month_starts = [i for i, d in enumerate(dates) if d.day == 1]
    month_labels = [dates[i].strftime('%b') for i in month_starts]

    # Create plots
    print(f"\nCreating plots for {n_locations} locations...")

    # Create a geographically realistic compact layout
    # Format: location_code -> (row, col)
    geo_layout = {
        # Row 0: US centered at top
        'US': (0, 4),
        # Row 1: Alaska + Northern New England
        '02': (1, 0),  # AK
        '23': (1, 8),  # ME
        # Row 2: Pacific NW + Northern tier
        '53': (2, 0),  # WA
        '16': (2, 1),  # ID
        '30': (2, 2),  # MT
        '38': (2, 3),  # ND
        '27': (2, 4),  # MN
        '55': (2, 5),  # WI
        '26': (2, 6),  # MI
        '50': (2, 7),  # VT
        '33': (2, 8),  # NH
        # Row 3: West + Upper Midwest + Northeast
        '41': (3, 0),  # OR
        '56': (3, 1),  # WY
        '46': (3, 2),  # SD
        '19': (3, 3),  # IA
        '17': (3, 4),  # IL
        '18': (3, 5),  # IN
        '39': (3, 6),  # OH
        '36': (3, 7),  # NY
        '25': (3, 8),  # MA
        # Row 4: California + Mountain West + Midwest + Mid-Atlantic
        '06': (4, 0),  # CA
        '32': (4, 1),  # NV
        '49': (4, 2),  # UT
        '31': (4, 3),  # NE
        '29': (4, 4),  # MO
        '21': (4, 5),  # KY
        '54': (4, 6),  # WV
        '42': (4, 7),  # PA
        '44': (4, 8),  # RI
        # Row 5: Southwest + South Central + Mid-Atlantic
        '04': (5, 1),  # AZ
        '08': (5, 2),  # CO
        '20': (5, 3),  # KS
        '05': (5, 4),  # AR
        '47': (5, 5),  # TN
        '51': (5, 6),  # VA
        '34': (5, 7),  # NJ
        '09': (5, 8),  # CT
        # Row 6: Deep Southwest + Deep South + Southeast
        '35': (6, 1),  # NM
        '40': (6, 2),  # OK
        '22': (6, 3),  # LA
        '28': (6, 4),  # MS
        '01': (6, 5),  # AL
        '13': (6, 6),  # GA
        '37': (6, 7),  # NC
        '24': (6, 8),  # MD
        # Row 7: Texas + Southeast Coast + DC/DE
        '48': (7, 2),  # TX
        '45': (7, 5),  # SC
        '12': (7, 6),  # FL
        '11': (7, 7),  # DC
        '10': (7, 8),  # DE
        # Row 8: Hawaii + Puerto Rico
        '15': (8, 0),  # HI
        '72': (8, 6),  # PR
    }

    # Create figure with compact layout
    n_rows = 9
    n_cols = 9
    fig = plt.figure(figsize=(22, 18))
    gs = gridspec.GridSpec(n_rows, n_cols, figure=fig, hspace=0.4, wspace=0.3)

    for loc, loc_name, loc_median, loc_lower, loc_upper in zip(
        locations, location_names, seasonality_median, seasonality_lower, seasonality_upper
    ):
        # Get grid position for this location
        if loc not in geo_layout:
            print(f"Warning: No geo layout for location {loc}")
            continue

        row, col = geo_layout[loc]
        ax = fig.add_subplot(gs[row, col])

        # Plot credible interval
        ax.fill_between(day_grid, loc_lower, loc_upper, alpha=0.3, color='blue',
                        label='95% CI')

        # Plot median
        ax.plot(day_grid, loc_median, 'b-', linewidth=1.5, label='Median')

        # Add zero line
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

        # Format axes
        ax.set_title(f'{loc_name}', fontsize=10, fontweight='bold')
        ax.set_xlabel('Month', fontsize=8)
        ax.set_ylabel('Seasonal effect', fontsize=8)
        ax.set_xticks(month_starts)
        ax.set_xticklabels(month_labels, fontsize=7)
        ax.tick_params(axis='y', labelsize=7)
        ax.grid(True, alpha=0.3)

        # Add legend to US plot
        if loc == 'US':
            ax.legend(fontsize=7, loc='upper right')

    plt.suptitle(
        f'Fourier Seasonal Patterns by Location (K={K} harmonics)',
        fontsize=14, fontweight='bold', y=0.995
    )

    # Save plot
    output_path = Path('fourier_seasonality_by_location.png')
    print(f"\nSaving main plot to: {output_path.absolute()}")
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"✓ Plot saved successfully")

    # Also create a summary plot showing all locations on one plot
    print("\nCreating summary plot...")
    fig2, ax2 = plt.subplots(figsize=(12, 6))

    for idx, (loc, loc_median) in enumerate(zip(locations, seasonality_median)):
        alpha = 0.3 if loc != 'US' else 1.0
        linewidth = 0.5 if loc != 'US' else 2.0
        color = 'red' if loc == 'US' else 'blue'
        label = loc if loc == 'US' else None
        ax2.plot(day_grid, loc_median, alpha=alpha, linewidth=linewidth,
                color=color, label=label)

    ax2.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax2.set_xlabel('Month', fontsize=12)
    ax2.set_ylabel('Seasonal effect', fontsize=12)
    ax2.set_title(f'All Locations: Fourier Seasonal Patterns (K={K} harmonics)',
                 fontsize=14, fontweight='bold')
    ax2.set_xticks(month_starts)
    ax2.set_xticklabels(month_labels)
    ax2.grid(True, alpha=0.3)
    if any(loc == 'US' for loc in locations):
        ax2.legend()

    output_path2 = Path('fourier_seasonality_all_locations.png')
    print(f"Saving summary plot to: {output_path2.absolute()}")
    plt.savefig(output_path2, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"✓ Summary plot saved successfully")

    print("\n" + "="*60)
    print("Done! Generated plots:")
    print(f"  1. {output_path.absolute()}")
    print(f"  2. {output_path2.absolute()}")
    print("="*60)

    return seasonality_median, locations


if __name__ == "__main__":
    seasonality_median, locations = main()
