"""
Utilities for fetching weather data from Open-Meteo API.

This module provides functions to:
1. Load static historical weather data from pre-built CSV
2. Fetch recent weather data from API (only for dates not in static data)
3. Fetch current weather forecasts (for prediction)
4. Aggregate daily data to weekly resolution
5. Cache API responses to avoid redundant calls
"""

import hashlib
import json
import os
from datetime import date, timedelta
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import requests


# Path to static historical weather data (built from ERA5 reanalysis)
STATIC_WEATHER_PATH = Path(__file__).parent / "data" / "historical_weather.csv"

# Cached static weather data (loaded once)
_static_weather_cache: Optional[pd.DataFrame] = None


def load_static_weather_data(path: Path = None) -> Optional[pd.DataFrame]:
    """
    Load the static historical weather dataset.

    This dataset was pre-built from Open-Meteo Archive API (ERA5 reanalysis)
    and covers 2010-2024 for all US states.

    Parameters
    ----------
    path : Path, optional
        Path to the CSV file. Defaults to STATIC_WEATHER_PATH.

    Returns
    -------
    pd.DataFrame or None
        DataFrame with columns: wk_end_date, weather_temp, weather_humidity, location
        Returns None if file doesn't exist.
    """
    global _static_weather_cache

    if path is None:
        path = STATIC_WEATHER_PATH

    # Return cached data if available
    if _static_weather_cache is not None:
        return _static_weather_cache

    if not path.exists():
        return None

    df = pd.read_csv(path)
    # Handle mixed date formats (some may have timestamps)
    df["wk_end_date"] = pd.to_datetime(df["wk_end_date"], format="mixed")
    _static_weather_cache = df
    return df


def get_static_weather_cutoff() -> Optional[date]:
    """
    Get the latest date in the static weather dataset.

    Returns
    -------
    date or None
        The latest wk_end_date in the static dataset, or None if no data.
    """
    static_df = load_static_weather_data()
    if static_df is None or static_df.empty:
        return None
    return static_df["wk_end_date"].max().date()


# Population-weighted centroids for US states (approximate)
# Format: FIPS code -> (latitude, longitude)
STATE_CENTROIDS = {
    "US": (39.8, -98.6),    # Geographic center of contiguous US
    "01": (33.0, -86.8),    # Alabama (Birmingham)
    "02": (61.2, -149.9),   # Alaska (Anchorage)
    "04": (33.4, -112.1),   # Arizona (Phoenix)
    "05": (34.7, -92.3),    # Arkansas (Little Rock)
    "06": (34.0, -118.2),   # California (Los Angeles)
    "08": (39.7, -105.0),   # Colorado (Denver)
    "09": (41.3, -72.9),    # Connecticut (New Haven)
    "10": (39.7, -75.5),    # Delaware (Wilmington)
    "11": (38.9, -77.0),    # DC (Washington)
    "12": (25.8, -80.2),    # Florida (Miami)
    "13": (33.7, -84.4),    # Georgia (Atlanta)
    "15": (21.3, -157.9),   # Hawaii (Honolulu)
    "16": (43.6, -116.2),   # Idaho (Boise)
    "17": (41.9, -87.6),    # Illinois (Chicago)
    "18": (39.8, -86.2),    # Indiana (Indianapolis)
    "19": (41.6, -93.6),    # Iowa (Des Moines)
    "20": (39.0, -94.6),    # Kansas (Kansas City area)
    "21": (38.2, -85.8),    # Kentucky (Louisville)
    "22": (30.0, -90.1),    # Louisiana (New Orleans)
    "23": (43.7, -70.3),    # Maine (Portland)
    "24": (39.3, -76.6),    # Maryland (Baltimore)
    "25": (42.4, -71.1),    # Massachusetts (Boston)
    "26": (42.3, -83.0),    # Michigan (Detroit)
    "27": (44.9, -93.3),    # Minnesota (Minneapolis)
    "28": (32.3, -90.2),    # Mississippi (Jackson)
    "29": (38.6, -90.2),    # Missouri (St. Louis)
    "30": (46.9, -110.4),   # Montana (Helena)
    "31": (41.3, -96.0),    # Nebraska (Omaha)
    "32": (36.2, -115.1),   # Nevada (Las Vegas)
    "33": (43.0, -71.5),    # New Hampshire (Manchester)
    "34": (40.1, -74.7),    # New Jersey (Trenton area)
    "35": (35.1, -106.6),   # New Mexico (Albuquerque)
    "36": (40.7, -74.0),    # New York (NYC)
    "37": (35.8, -78.6),    # North Carolina (Raleigh)
    "38": (46.8, -100.8),   # North Dakota (Bismarck)
    "39": (39.9, -82.9),    # Ohio (Columbus)
    "40": (35.5, -97.5),    # Oklahoma (Oklahoma City)
    "41": (45.5, -122.7),   # Oregon (Portland)
    "42": (40.0, -75.1),    # Pennsylvania (Philadelphia)
    "44": (41.8, -71.4),    # Rhode Island (Providence)
    "45": (34.0, -81.0),    # South Carolina (Columbia)
    "46": (43.5, -96.7),    # South Dakota (Sioux Falls)
    "47": (36.2, -86.8),    # Tennessee (Nashville)
    "48": (29.8, -95.4),    # Texas (Houston)
    "49": (40.8, -111.9),   # Utah (Salt Lake City)
    "50": (44.3, -72.6),    # Vermont (Montpelier area)
    "51": (37.5, -77.4),    # Virginia (Richmond)
    "53": (47.6, -122.3),   # Washington (Seattle)
    "54": (38.3, -81.6),    # West Virginia (Charleston)
    "55": (43.0, -89.4),    # Wisconsin (Madison)
    "56": (41.1, -104.8),   # Wyoming (Cheyenne)
    "72": (18.4, -66.1),    # Puerto Rico (San Juan)
}

# Weather variables to fetch for historical data
DEFAULT_WEATHER_VARS = [
    "temperature_2m_mean",
    "relative_humidity_2m_mean",
    "precipitation_sum",
    "shortwave_radiation_sum",
]

# Weather variables for forecasts (subset most relevant for flu)
DEFAULT_FORECAST_VARS = [
    "temperature_2m_mean",
    "relative_humidity_2m_mean",
]


def get_cache_path(cache_dir: str = ".weather_cache") -> Path:
    """Get or create the cache directory path."""
    cache_path = Path(cache_dir)
    cache_path.mkdir(exist_ok=True)
    return cache_path


def get_cache_key(lat: float, lon: float, start_date: date, end_date: date,
                  variables: list, api_type: str) -> str:
    """Generate a cache key for an API request."""
    key_data = f"{lat}_{lon}_{start_date}_{end_date}_{sorted(variables)}_{api_type}"
    return hashlib.md5(key_data.encode()).hexdigest()


def fetch_historical_weather(
    lat: float,
    lon: float,
    start_date: date,
    end_date: date,
    variables: list = None,
    cache_dir: str = ".weather_cache",
    timeout: int = 120,
    max_retries: int = 3
) -> Optional[pd.DataFrame]:
    """
    Fetch historical weather forecast data from Open-Meteo.

    Parameters
    ----------
    lat : float
        Latitude of location
    lon : float
        Longitude of location
    start_date : date
        Start date for data retrieval
    end_date : date
        End date for data retrieval
    variables : list, optional
        Weather variables to fetch. Defaults to DEFAULT_WEATHER_VARS.
    cache_dir : str
        Directory for caching API responses
    timeout : int
        Request timeout in seconds (default 120 for long date ranges)
    max_retries : int
        Number of retry attempts on failure

    Returns
    -------
    pd.DataFrame or None
        DataFrame with columns: date, and one column per weather variable
        Returns None if API call fails
    """
    import time as time_module

    if variables is None:
        variables = DEFAULT_WEATHER_VARS

    # Check cache first
    cache_path = get_cache_path(cache_dir)
    cache_key = get_cache_key(lat, lon, start_date, end_date, variables, "historical")
    cache_file = cache_path / f"{cache_key}.json"

    if cache_file.exists():
        with open(cache_file, "r") as f:
            data = json.load(f)
        return pd.DataFrame(data["daily"])

    # Fetch from API with retry logic
    url = "https://historical-forecast-api.open-meteo.com/v1/forecast"
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

            # Cache the response
            with open(cache_file, "w") as f:
                json.dump(data, f)

            df = pd.DataFrame(data["daily"])
            df["time"] = pd.to_datetime(df["time"])
            return df

        except requests.exceptions.Timeout as e:
            last_error = e
            wait_time = 2 ** attempt  # Exponential backoff: 1, 2, 4 seconds
            print(f"Timeout fetching weather (attempt {attempt + 1}/{max_retries}), retrying in {wait_time}s...")
            time_module.sleep(wait_time)
        except requests.exceptions.HTTPError as e:
            last_error = e
            if e.response.status_code == 429:
                # Rate limited - use longer backoff
                wait_time = 10 * (attempt + 1)  # 10, 20, 30 seconds
                print(f"Rate limited (429), waiting {wait_time}s before retry {attempt + 1}/{max_retries}...")
                time_module.sleep(wait_time)
            else:
                print(f"HTTP error fetching historical weather: {e}")
                break
        except requests.RequestException as e:
            last_error = e
            print(f"Error fetching historical weather: {e}")
            break  # Don't retry on other errors

    print(f"Failed to fetch historical weather after {max_retries} attempts: {last_error}")
    return None


def fetch_weather_forecast(
    lat: float,
    lon: float,
    forecast_days: int = 16,
    variables: list = None,
    cache_dir: str = ".weather_cache",
    timeout: int = 60,
    max_retries: int = 3
) -> Optional[pd.DataFrame]:
    """
    Fetch current weather forecast from Open-Meteo.

    Parameters
    ----------
    lat : float
        Latitude of location
    lon : float
        Longitude of location
    forecast_days : int
        Number of days to forecast (max 16)
    variables : list, optional
        Weather variables to fetch. Defaults to DEFAULT_WEATHER_VARS.
    cache_dir : str
        Directory for caching API responses
    timeout : int
        Request timeout in seconds
    max_retries : int
        Number of retry attempts on failure

    Returns
    -------
    pd.DataFrame or None
        DataFrame with columns: date, and one column per weather variable
        Returns None if API call fails
    """
    import time as time_module

    if variables is None:
        variables = DEFAULT_WEATHER_VARS

    # For forecasts, cache with today's date as part of key (forecasts change daily)
    today = date.today()
    cache_path = get_cache_path(cache_dir)
    cache_key = get_cache_key(lat, lon, today, today, variables, f"forecast_{forecast_days}")
    cache_file = cache_path / f"{cache_key}.json"

    if cache_file.exists():
        with open(cache_file, "r") as f:
            data = json.load(f)
        df = pd.DataFrame(data["daily"])
        df["time"] = pd.to_datetime(df["time"])
        return df

    # Fetch from API with retry logic
    url = "https://api.open-meteo.com/v1/forecast"
    params = {
        "latitude": lat,
        "longitude": lon,
        "daily": ",".join(variables),
        "timezone": "UTC",
        "forecast_days": min(forecast_days, 16)
    }

    last_error = None
    for attempt in range(max_retries):
        try:
            response = requests.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            data = response.json()

            # Cache the response
            with open(cache_file, "w") as f:
                json.dump(data, f)

            df = pd.DataFrame(data["daily"])
            df["time"] = pd.to_datetime(df["time"])
            return df

        except requests.exceptions.Timeout as e:
            last_error = e
            wait_time = 2 ** attempt
            print(f"Timeout fetching forecast (attempt {attempt + 1}/{max_retries}), retrying in {wait_time}s...")
            time_module.sleep(wait_time)
        except requests.exceptions.HTTPError as e:
            last_error = e
            if e.response.status_code == 429:
                wait_time = 10 * (attempt + 1)
                print(f"Rate limited (429), waiting {wait_time}s before retry {attempt + 1}/{max_retries}...")
                time_module.sleep(wait_time)
            else:
                print(f"HTTP error fetching forecast: {e}")
                break
        except requests.RequestException as e:
            last_error = e
            print(f"Error fetching weather forecast: {e}")
            break

    print(f"Failed to fetch weather forecast after {max_retries} attempts: {last_error}")
    return None


def aggregate_daily_to_weekly(
    df: pd.DataFrame,
    date_col: str = "time",
    agg_method: dict = None
) -> pd.DataFrame:
    """
    Aggregate daily weather data to weekly (Saturday-ending weeks).

    Parameters
    ----------
    df : pd.DataFrame
        Daily weather data with a date column
    date_col : str
        Name of the date column
    agg_method : dict, optional
        Aggregation method per column. Defaults to mean for temperature/humidity,
        sum for precipitation/radiation.

    Returns
    -------
    pd.DataFrame
        Weekly aggregated data with wk_end_date column
    """
    if agg_method is None:
        agg_method = {
            "temperature_2m_mean": "mean",
            "relative_humidity_2m_mean": "mean",
            "precipitation_sum": "sum",
            "shortwave_radiation_sum": "sum",
        }

    df = df.copy()
    df[date_col] = pd.to_datetime(df[date_col])

    # Find the Saturday ending each week
    # weekday(): Monday=0, Saturday=5
    df["wk_end_date"] = df[date_col] + pd.to_timedelta(
        (5 - df[date_col].dt.weekday) % 7, unit="D"
    )

    # Aggregate
    agg_dict = {col: method for col, method in agg_method.items() if col in df.columns}
    weekly = df.groupby("wk_end_date").agg(agg_dict).reset_index()

    return weekly


def fetch_weather_for_locations(
    locations: list,
    start_date: date,
    end_date: date,
    variables: list = None,
    include_forecast: bool = True,
    cache_dir: str = ".weather_cache",
    request_delay: float = 0.2,
    use_static_data: bool = True,
) -> pd.DataFrame:
    """
    Fetch weather data for multiple locations and combine into a single DataFrame.

    This function uses a hybrid approach:
    1. Load pre-built static data (2010-2024) from CSV for dates in the static dataset
    2. Fetch only recent data from the API for dates after the static cutoff
    3. Optionally append forecast data

    This reduces API calls from ~53 to just a few (only for recent weeks not in static data),
    and eliminates the ~75% missing values that occurred when using only the Historical
    Forecast API (which only has data from 2022+).

    Parameters
    ----------
    locations : list
        List of FIPS location codes (e.g., ["US", "01", "06"])
    start_date : date
        Start date for historical data
    end_date : date
        End date for historical data
    variables : list, optional
        Weather variables to fetch (only used for API calls, not static data)
    include_forecast : bool
        Whether to also fetch future forecasts
    cache_dir : str
        Directory for caching API responses
    request_delay : float
        Delay in seconds between API requests to avoid rate limiting
    use_static_data : bool
        If True, use pre-built static dataset for historical data.
        Set to False to fetch everything from API (not recommended).

    Returns
    -------
    pd.DataFrame
        Combined weather data with columns: location, wk_end_date, and weather variables
    """
    import time as time_module

    if variables is None:
        variables = DEFAULT_WEATHER_VARS

    all_data = []
    api_calls_made = 0

    # Try to load static data first
    static_df = None
    static_cutoff = None
    if use_static_data:
        static_df = load_static_weather_data()
        if static_df is not None:
            static_cutoff = get_static_weather_cutoff()
            print(f"Using static weather data (cutoff: {static_cutoff})")

    for loc in locations:
        if loc not in STATE_CENTROIDS:
            print(f"Warning: No centroid for location {loc}, skipping")
            continue

        # Use static data for this location if available
        if static_df is not None:
            loc_static = static_df[
                (static_df["location"] == loc) &
                (static_df["wk_end_date"] >= pd.Timestamp(start_date)) &
                (static_df["wk_end_date"] <= pd.Timestamp(end_date))
            ].copy()

            if not loc_static.empty:
                all_data.append(loc_static)

        lat, lon = STATE_CENTROIDS[loc]

        # Determine if we need to fetch recent data from API
        # Only fetch if end_date is after static_cutoff
        api_start_date = None
        if static_df is not None and static_cutoff is not None:
            if end_date > static_cutoff:
                # Fetch data from day after static cutoff to end_date
                api_start_date = static_cutoff + timedelta(days=1)
        else:
            # No static data, fetch everything from API
            api_start_date = start_date

        if api_start_date is not None and api_start_date <= end_date:
            # Fetch recent data from API
            hist_df = fetch_historical_weather(
                lat, lon, api_start_date, end_date, variables, cache_dir
            )
            api_calls_made += 1

            if hist_df is not None:
                hist_weekly = aggregate_daily_to_weekly(hist_df)
                hist_weekly["location"] = loc

                # Rename columns to match static data format
                rename_map = {
                    "temperature_2m_mean": "weather_temp",
                    "relative_humidity_2m_mean": "weather_humidity",
                    "precipitation_sum": "weather_precip",
                    "shortwave_radiation_sum": "weather_radiation",
                }
                hist_weekly = hist_weekly.rename(
                    columns={k: v for k, v in rename_map.items() if k in hist_weekly.columns}
                )
                all_data.append(hist_weekly)

            # Add delay between requests to avoid rate limiting
            if request_delay > 0:
                time_module.sleep(request_delay)

        # Fetch forecast data if requested
        if include_forecast:
            fcst_df = fetch_weather_forecast(lat, lon, 16, variables, cache_dir)
            api_calls_made += 1
            if fcst_df is not None:
                fcst_weekly = aggregate_daily_to_weekly(fcst_df)
                fcst_weekly["location"] = loc
                # Only keep forecast dates beyond end_date
                fcst_weekly = fcst_weekly[fcst_weekly["wk_end_date"] > pd.Timestamp(end_date)]

                # Rename columns
                rename_map = {
                    "temperature_2m_mean": "weather_temp",
                    "relative_humidity_2m_mean": "weather_humidity",
                    "precipitation_sum": "weather_precip",
                    "shortwave_radiation_sum": "weather_radiation",
                }
                fcst_weekly = fcst_weekly.rename(
                    columns={k: v for k, v in rename_map.items() if k in fcst_weekly.columns}
                )
                if len(fcst_weekly) > 0:
                    all_data.append(fcst_weekly)

            # Add delay after forecast request too
            if request_delay > 0:
                time_module.sleep(request_delay)

    if api_calls_made > 0:
        print(f"Made {api_calls_made} API calls for recent/forecast data")

    if not all_data:
        return pd.DataFrame()

    combined = pd.concat(all_data, ignore_index=True)

    # Remove duplicates (in case static and API overlap)
    combined = combined.drop_duplicates(subset=["location", "wk_end_date"], keep="last")
    combined = combined.sort_values(["location", "wk_end_date"]).reset_index(drop=True)

    return combined


def fetch_forecast_features_for_locations(
    locations: list,
    ref_date: date,
    forecast_horizons: list = None,
    variables: list = None,
    cache_dir: str = ".weather_cache",
    request_delay: float = 0.2
) -> pd.DataFrame:
    """
    Fetch weather forecast data and create horizon-specific features.

    This creates features like 'weather_temp_fcst_1wk' representing the
    forecasted temperature for 1 week ahead from the reference date.

    Parameters
    ----------
    locations : list
        List of FIPS location codes (e.g., ["US", "01", "06"])
    ref_date : date
        Reference date for the forecast (typically the model run date)
    forecast_horizons : list, optional
        Forecast horizons in weeks (e.g., [1, 2, 3]). Defaults to [1, 2, 3].
    variables : list, optional
        Weather variables to fetch. Defaults to DEFAULT_FORECAST_VARS.
    cache_dir : str
        Directory for caching API responses
    request_delay : float
        Delay in seconds between API requests to avoid rate limiting

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: location, and forecast features like
        weather_temp_fcst_1wk, weather_humidity_fcst_2wk, etc.
        One row per location.
    """
    import time as time_module

    if forecast_horizons is None:
        forecast_horizons = [1, 2, 3]

    if variables is None:
        variables = DEFAULT_FORECAST_VARS

    all_data = []

    for loc in locations:
        if loc not in STATE_CENTROIDS:
            print(f"Warning: No centroid for location {loc}, skipping forecast")
            continue

        lat, lon = STATE_CENTROIDS[loc]

        # Fetch 16-day forecast
        fcst_df = fetch_weather_forecast(lat, lon, 16, variables, cache_dir)

        # Add delay between requests
        if request_delay > 0:
            time_module.sleep(request_delay)

        if fcst_df is None:
            continue

        # Aggregate to weekly
        fcst_weekly = aggregate_daily_to_weekly(fcst_df)

        # Calculate target weeks for each horizon
        # Reference date is a Saturday, horizon 1 = next Saturday, etc.
        ref_saturday = pd.Timestamp(ref_date)
        # Ensure ref_date is a Saturday (adjust if needed)
        days_to_saturday = (5 - ref_saturday.weekday()) % 7
        ref_saturday = ref_saturday + pd.Timedelta(days=days_to_saturday)

        loc_data = {"location": loc}

        for horizon in forecast_horizons:
            target_date = ref_saturday + pd.Timedelta(weeks=horizon)

            # Find the forecast for this target week
            fcst_row = fcst_weekly[fcst_weekly["wk_end_date"] == target_date]

            for var in variables:
                # Map API variable names to feature names
                var_short = var.replace("temperature_2m_mean", "temp").replace(
                    "relative_humidity_2m_mean", "humidity"
                ).replace("precipitation_sum", "precip").replace(
                    "shortwave_radiation_sum", "radiation"
                )
                feat_name = f"weather_{var_short}_fcst_{horizon}wk"

                if len(fcst_row) > 0 and var in fcst_row.columns:
                    loc_data[feat_name] = fcst_row[var].values[0]
                else:
                    loc_data[feat_name] = np.nan

        all_data.append(loc_data)

    if not all_data:
        return pd.DataFrame()

    return pd.DataFrame(all_data)


def add_weather_features_to_df(
    df: pd.DataFrame,
    weather_df: pd.DataFrame,
    weather_cols: list = None,
    lags: list = None
) -> tuple:
    """
    Merge weather data into a surveillance DataFrame and create lagged features.

    Parameters
    ----------
    df : pd.DataFrame
        Surveillance data with 'location' and 'wk_end_date' columns
    weather_df : pd.DataFrame
        Weather data from fetch_weather_for_locations()
    weather_cols : list, optional
        Weather columns to use. Defaults to all weather_* columns.
    lags : list, optional
        Lag values to create (in weeks). Defaults to [0, 1, 2].

    Returns
    -------
    tuple
        (augmented_df, feature_names) where feature_names is the list of
        new weather feature column names
    """
    if weather_cols is None:
        weather_cols = [c for c in weather_df.columns if c.startswith("weather_")]

    if lags is None:
        lags = [0, 1, 2]

    # Ensure date columns are datetime
    df = df.copy()
    df["wk_end_date"] = pd.to_datetime(df["wk_end_date"])
    weather_df = weather_df.copy()
    weather_df["wk_end_date"] = pd.to_datetime(weather_df["wk_end_date"])

    # Merge weather data (lag 0)
    merge_cols = ["location", "wk_end_date"] + weather_cols
    df = df.merge(
        weather_df[merge_cols],
        on=["location", "wk_end_date"],
        how="left"
    )

    feature_names = []

    # Create lagged features
    for lag in lags:
        if lag == 0:
            # Already merged above
            for col in weather_cols:
                feature_names.append(col)
        else:
            # Create lagged version
            for col in weather_cols:
                lag_col = f"{col}_lag{lag}"
                # Shift within each location group
                df[lag_col] = df.groupby("location")[col].shift(lag)
                feature_names.append(lag_col)

    return df, feature_names


def add_forecast_features_to_df(
    df: pd.DataFrame,
    forecast_df: pd.DataFrame,
) -> tuple:
    """
    Merge forecast features into a surveillance DataFrame.

    Forecast features are location-specific but not time-varying (they represent
    the current forecast at model run time), so they are merged only on location.

    Parameters
    ----------
    df : pd.DataFrame
        Surveillance data with 'location' column
    forecast_df : pd.DataFrame
        Forecast data from fetch_forecast_features_for_locations()

    Returns
    -------
    tuple
        (augmented_df, feature_names) where feature_names is the list of
        new forecast feature column names
    """
    df = df.copy()

    # Get forecast feature columns (everything except 'location')
    forecast_cols = [c for c in forecast_df.columns if c != "location"]

    # Merge on location only
    df = df.merge(forecast_df, on="location", how="left")

    return df, forecast_cols


# Default noise standard deviations for weather variables (calibrated to typical forecast errors)
DEFAULT_NOISE_STD = {
    "weather_temp": 2.0,       # ~2°C typical 1-2 week forecast error
    "weather_humidity": 10.0,  # ~10% typical forecast error
    "weather_precip": 5.0,     # ~5mm (highly variable)
    "weather_radiation": 2.0,  # ~2 MJ/m²
}


def add_weather_features_with_leads(
    df: pd.DataFrame,
    weather_df: pd.DataFrame,
    weather_cols: list = None,
    lags: list = None,
    leads: list = None,
    noise_std: dict = None,
    random_seed: int = None,
) -> tuple:
    """
    Merge weather data and create lagged AND leading features with noise injection.

    Lagged and leading features are created in the weather dataframe FIRST,
    then merged with the surveillance data. This ensures that even the first
    surveillance observation can have lagged weather features (from prior weeks
    that exist in the weather data but not in the surveillance data).

    Leading features have Gaussian noise added to simulate forecast uncertainty.

    Parameters
    ----------
    df : pd.DataFrame
        Surveillance data with 'location' and 'wk_end_date' columns
    weather_df : pd.DataFrame
        Weather data from fetch_weather_for_locations(). Should include data
        from (earliest_surveillance - max_lag weeks) through (ref_date + max_lead weeks).
    weather_cols : list, optional
        Weather columns to use. Defaults to all weather_* columns.
    lags : list, optional
        Lag values to create in weeks (positive integers). Defaults to [1, 2].
    leads : list, optional
        Lead values to create in weeks (positive integers). Defaults to [1, 2].
    noise_std : dict, optional
        Standard deviation of Gaussian noise to add to leading features,
        keyed by weather column name. Defaults to DEFAULT_NOISE_STD.
        Noise is only added to non-NaN values.
    random_seed : int, optional
        Seed for reproducible noise generation.

    Returns
    -------
    tuple
        (augmented_df, feature_names) where feature_names is the list of
        all weather feature column names (current, lagged, and leading)
    """
    if weather_cols is None:
        weather_cols = [c for c in weather_df.columns if c.startswith("weather_")]

    if lags is None:
        lags = [1, 2]

    if leads is None:
        leads = [1, 2]

    if noise_std is None:
        noise_std = DEFAULT_NOISE_STD

    # Set up random number generator
    rng = np.random.default_rng(seed=random_seed)

    # Ensure date columns are datetime
    df = df.copy()
    df["wk_end_date"] = pd.to_datetime(df["wk_end_date"])
    weather_df = weather_df.copy()
    weather_df["wk_end_date"] = pd.to_datetime(weather_df["wk_end_date"])

    feature_names = list(weather_cols)  # Current week features

    # Create LAGGED features in weather_df FIRST (before merging)
    # This ensures surveillance data can get lagged weather even for its first week
    for lag in lags:
        for col in weather_cols:
            lag_col = f"{col}_lag{lag}"
            weather_df[lag_col] = weather_df.groupby("location")[col].shift(lag)
            feature_names.append(lag_col)

    # Create LEADING features in weather_df (before merging)
    for lead in leads:
        for col in weather_cols:
            lead_col = f"{col}_lead{lead}"
            # Shift backwards to get future values
            weather_df[lead_col] = weather_df.groupby("location")[col].shift(-lead)
            feature_names.append(lead_col)

    # Now merge ALL weather features (current, lagged, leading) with surveillance data
    merge_cols = ["location", "wk_end_date"] + feature_names
    df = df.merge(
        weather_df[merge_cols],
        on=["location", "wk_end_date"],
        how="left"
    )

    # Add noise to LEADING features to simulate forecast uncertainty
    for lead in leads:
        for col in weather_cols:
            lead_col = f"{col}_lead{lead}"

            has_value_mask = ~df[lead_col].isna()
            col_noise_std = noise_std.get(col, 2.0)

            n_with_value = has_value_mask.sum()
            if n_with_value > 0:
                noise = rng.normal(0, col_noise_std, size=n_with_value)
                df.loc[has_value_mask, lead_col] = df.loc[has_value_mask, lead_col] + noise

    return df, feature_names


# Convenience function for quick testing
def test_api():
    """Test the Open-Meteo API with a simple request."""
    print("Testing Open-Meteo API...")

    # Test historical
    lat, lon = STATE_CENTROIDS["25"]  # Massachusetts
    end_date = date.today() - timedelta(days=1)
    start_date = end_date - timedelta(days=28)

    print(f"\n1. Historical data for MA ({lat}, {lon})")
    print(f"   Date range: {start_date} to {end_date}")

    hist = fetch_historical_weather(lat, lon, start_date, end_date)
    if hist is not None:
        print(f"   Success! Shape: {hist.shape}")
        weekly = aggregate_daily_to_weekly(hist)
        print(f"   Weekly aggregated: {len(weekly)} weeks")
        print(weekly.to_string())
    else:
        print("   Failed to fetch historical data")

    # Test forecast
    print(f"\n2. Forecast data for MA")
    fcst = fetch_weather_forecast(lat, lon, 16)
    if fcst is not None:
        print(f"   Success! Shape: {fcst.shape}")
        print(f"   Date range: {fcst['time'].min()} to {fcst['time'].max()}")
    else:
        print("   Failed to fetch forecast data")

    print("\nAPI test complete!")


if __name__ == "__main__":
    test_api()
