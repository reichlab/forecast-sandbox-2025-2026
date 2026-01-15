"""
Build a static historical weather dataset for all US states.

This script fetches historical weather data from Open-Meteo's Archive API
(which uses ERA5 reanalysis) and aggregates it to weekly resolution.

The resulting dataset can be used as a static baseline, supplemented with
recent data from the Historical Forecast API for operational runs.

Usage:
    python build_historical_weather.py --start_year=2010 --end_year=2023
"""

import argparse
import time
from datetime import date
from pathlib import Path

import pandas as pd
import requests

from weather_utils import STATE_CENTROIDS, aggregate_daily_to_weekly


def fetch_historical_weather_era5(
    lat: float,
    lon: float,
    start_date: date,
    end_date: date,
    variables: list,
    timeout: int = 120,
    max_retries: int = 3
) -> pd.DataFrame:
    """
    Fetch historical weather from Open-Meteo Archive API (ERA5 reanalysis).

    This API has data from 1940-present and is more complete than the
    Historical Forecast API which only has data from 2022.
    """
    url = "https://archive-api.open-meteo.com/v1/archive"
    params = {
        "latitude": lat,
        "longitude": lon,
        "start_date": str(start_date),
        "end_date": str(end_date),
        "daily": ",".join(variables),
        "timezone": "UTC"
    }

    last_error = None
    for attempt in range(max_retries):
        try:
            response = requests.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            data = response.json()

            df = pd.DataFrame(data["daily"])
            df["time"] = pd.to_datetime(df["time"])
            return df

        except requests.exceptions.Timeout as e:
            last_error = e
            wait_time = 2 ** attempt
            print(f"  Timeout (attempt {attempt + 1}/{max_retries}), retrying in {wait_time}s...")
            time.sleep(wait_time)
        except requests.exceptions.HTTPError as e:
            last_error = e
            if e.response.status_code == 429:
                wait_time = 10 * (attempt + 1)
                print(f"  Rate limited (429), waiting {wait_time}s...")
                time.sleep(wait_time)
            else:
                print(f"  HTTP error: {e}")
                break
        except requests.RequestException as e:
            last_error = e
            print(f"  Request error: {e}")
            break

    print(f"  Failed after {max_retries} attempts: {last_error}")
    return None


def build_historical_dataset(
    start_year: int = 2010,
    end_year: int = 2023,
    output_path: str = "data/historical_weather.parquet",
    request_delay: float = 0.3
):
    """
    Build a complete historical weather dataset for all US states.
    """
    variables = ["temperature_2m_mean", "relative_humidity_2m_mean"]
    start_date = date(start_year, 1, 1)
    end_date = date(end_year, 12, 31)

    print(f"Building historical weather dataset")
    print(f"Date range: {start_date} to {end_date}")
    print(f"Locations: {len(STATE_CENTROIDS)}")
    print(f"Variables: {variables}")
    print("=" * 60)

    all_data = []

    for i, (loc, (lat, lon)) in enumerate(STATE_CENTROIDS.items()):
        print(f"[{i+1}/{len(STATE_CENTROIDS)}] Fetching {loc} ({lat}, {lon})...")

        df = fetch_historical_weather_era5(lat, lon, start_date, end_date, variables)

        if df is not None:
            # Aggregate to weekly
            weekly = aggregate_daily_to_weekly(df)
            weekly["location"] = loc
            all_data.append(weekly)
            print(f"  Success: {len(weekly)} weeks")
        else:
            print(f"  FAILED")

        # Delay between requests
        time.sleep(request_delay)

    if not all_data:
        print("No data retrieved!")
        return None

    # Combine all locations
    combined = pd.concat(all_data, ignore_index=True)

    # Rename columns
    rename_map = {
        "temperature_2m_mean": "weather_temp",
        "relative_humidity_2m_mean": "weather_humidity",
    }
    combined = combined.rename(columns=rename_map)

    # Ensure output directory exists
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Save as CSV (parquet requires pyarrow)
    if str(output_path).endswith('.parquet'):
        output_path = Path(str(output_path).replace('.parquet', '.csv'))
    combined.to_csv(output_path, index=False)
    print("=" * 60)
    print(f"Saved to {output_path}")
    print(f"Total rows: {len(combined)}")
    print(f"Date range: {combined['wk_end_date'].min()} to {combined['wk_end_date'].max()}")
    print(f"File size: {output_path.stat().st_size / 1024:.1f} KB")

    return combined


def main():
    parser = argparse.ArgumentParser(description="Build historical weather dataset")
    parser.add_argument("--start_year", type=int, default=2010, help="Start year")
    parser.add_argument("--end_year", type=int, default=2023, help="End year")
    parser.add_argument("--output", type=str, default="data/historical_weather.parquet",
                        help="Output file path")
    parser.add_argument("--delay", type=float, default=0.3,
                        help="Delay between API requests (seconds)")

    args = parser.parse_args()

    build_historical_dataset(
        start_year=args.start_year,
        end_year=args.end_year,
        output_path=args.output,
        request_delay=args.delay
    )


if __name__ == "__main__":
    main()
