import click
import datetime
from dateutil import relativedelta
from pathlib import Path
from typing import Optional
import pandas as pd
import numpy as np
import torch
from chronos import ChronosPipeline

from iddata.loader import DiseaseDataLoader


@click.command()
@click.option(
    "--today_date",
    type=str,
    required=False,
    help="Date to use as effective model run date (YYYY-MM-DD)",
)
@click.option(
    "--model_size",
    type=click.Choice(["tiny", "mini", "small", "base", "large"], case_sensitive=False),
    default="small",
    help="Size of Chronos model to use (default: small)",
)
@click.option(
    "--context_length",
    type=int,
    default=104,
    help="Number of weeks of historical data to use (default: 104 = 2 years)",
)
@click.option(
    "--num_samples",
    type=int,
    default=200,
    help="Number of samples to generate for quantile estimation (default: 200)",
)
def main(
    today_date: Optional[str] = None,
    model_size: str = "small",
    context_length: int = 104,
    num_samples: int = 200,
):
    """Generate flu predictions using Chronos time series foundation model."""
    try:
        today_date = datetime.date.fromisoformat(today_date)
    except (TypeError, ValueError):  # if today_date is None or a bad format
        today_date = datetime.date.today()
    reference_date = today_date + relativedelta.relativedelta(weekday=5)

    # Configuration
    disease = "flu"
    output_root = Path("../../model-output/")
    team_abbr = "UMass"
    model_abbr = f"chronos_{model_size}"

    # All locations (US + states/territories, excluding AK, CT, HI)
    locations = [
        "US", "01", "02", "04", "05", "06", "08", "09", "10", "11",
        "12", "13", "15", "16", "17", "18", "19", "20", "21", "22",
        "23", "24", "25", "26", "27", "28", "29", "30", "31", "32",
        "33", "34", "35", "36", "37", "38", "39", "40", "41", "42",
        "44", "45", "46", "47", "48", "49", "50", "51", "53", "54",
        "55", "56", "72"
    ]

    # Quantile levels required by hub
    q_levels = [
        0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
        0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70,
        0.75, 0.80, 0.85, 0.90, 0.95, 0.975, 0.99
    ]
    q_labels = [
        '0.01', '0.025', '0.05', '0.1', '0.15', '0.2',
        '0.25', '0.3', '0.35', '0.4', '0.45', '0.5',
        '0.55', '0.6', '0.65', '0.7', '0.75', '0.8',
        '0.85', '0.9', '0.95', '0.975', '0.99'
    ]

    # Forecast horizons: -1, 0, 1, 2, 3
    horizons = [-1, 0, 1, 2, 3]
    prediction_length = 5  # Number of weeks to forecast (all 5 horizons)

    print(f"Loading Chronos-T5-{model_size} model...")
    model_name = f"amazon/chronos-t5-{model_size}"
    pipeline = ChronosPipeline.from_pretrained(
        model_name,
        device_map="cpu",  # Use "cuda" if GPU available
        dtype=torch.bfloat16,
    )

    print(f"Fetching historical flu hospitalization data...")
    # Fetch historical data up to reference_date
    loader = DiseaseDataLoader()
    target_data = loader.load_nhsn(
        disease=disease,
        as_of=reference_date,
        rates=False  # We want raw counts, not rates
    )

    # Convert date column and sort
    target_data['wk_end_date'] = pd.to_datetime(target_data['wk_end_date'])
    target_data = target_data.sort_values(['location', 'wk_end_date'])

    print(f"Generating forecasts for {len(locations)} locations...")
    all_forecasts = []

    for location in locations:
        # Get historical data for this location
        loc_data = target_data[target_data['location'] == location].copy()

        if len(loc_data) < context_length:
            print(f"Warning: Location {location} has only {len(loc_data)} weeks of data")

        # Use last context_length weeks (or all available data if less)
        context_data = loc_data.tail(context_length)

        # Extract hospitalization counts as time series
        # Use 'inc' column which contains incident hospitalizations
        # Handle NaN values by forward filling, then backward filling if needed
        inc_values = context_data['inc'].ffill().bfill().fillna(0)
        context_series = torch.tensor(inc_values.values, dtype=torch.float32)

        # Generate forecast samples
        forecast_samples = pipeline.predict(
            context=context_series,
            prediction_length=prediction_length,
            num_samples=num_samples,
        )

        # forecast_samples shape: (num_samples, prediction_length)
        # Convert to numpy for easier manipulation
        forecast_samples_np = forecast_samples.squeeze().numpy()

        # Ensure non-negative forecasts (counts can't be negative)
        forecast_samples_np = np.maximum(forecast_samples_np, 0)

        # Calculate quantiles for each horizon
        for horizon_idx, horizon in enumerate(horizons):
            # Get samples for this horizon
            horizon_samples = forecast_samples_np[:, horizon_idx]

            # Calculate quantiles
            quantile_values = np.quantile(horizon_samples, q_levels)

            # Calculate target_end_date (Saturday at the given horizon)
            if horizon >= 0:
                target_end_date = reference_date + datetime.timedelta(weeks=horizon)
            else:
                # Horizon -1 is the previous week
                target_end_date = reference_date + datetime.timedelta(weeks=horizon)

            # Create forecast rows for this location, horizon
            for q_level, q_label, q_value in zip(q_levels, q_labels, quantile_values):
                all_forecasts.append({
                    'reference_date': reference_date.strftime('%Y-%m-%d'),
                    'target': 'wk inc flu hosp',
                    'horizon': horizon,
                    'location': location,
                    'target_end_date': target_end_date.strftime('%Y-%m-%d'),
                    'output_type': 'quantile',
                    'output_type_id': q_label,
                    'value': float(q_value)
                })

    # Create DataFrame and save
    forecast_df = pd.DataFrame(all_forecasts)

    # Create output directory if needed
    output_dir = output_root / f"{team_abbr}-{model_abbr}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save forecast file
    output_file = output_dir / f"{reference_date.strftime('%Y-%m-%d')}-{team_abbr}-{model_abbr}.csv"
    forecast_df.to_csv(output_file, index=False)

    print(f"Forecast saved to: {output_file}")
    print(f"Generated {len(forecast_df)} forecast rows")


if __name__ == "__main__":
    main()
