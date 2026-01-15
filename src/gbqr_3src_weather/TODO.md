# GBQR Weather Model - TODO

## Current Status (2024-12-12)

**COMPLETED:** The weather data pipeline is fully implemented with correct retrospective evaluation support.

### What was done:
1. Built static historical weather dataset (2010 - Dec 2025) using Open-Meteo Archive API (ERA5)
2. Implemented hybrid data loading: static data + API only when needed
3. Corrected the leading feature logic: always use observed weather + noise for retrospective runs

### How weather data is sourced:

**For any reference date:**
1. Calculate the required date range: (earliest flu data) to (ref_date + max_lead weeks)
2. Check if static data covers the entire range
3. If YES (retrospective): Use 100% static data, add noise to leading features
4. If NO (current/live): Use static for dates it covers, fetch recent data from API

| Scenario | Example Date | Data Sources | API Calls |
|----------|--------------|--------------|-----------|
| Retrospective | 2023-12-01 | 100% static | 0 |
| Near-current | 2025-11-01 | 100% static | 0 |
| Beyond static | 2025-12-15 | Static + API | 53 |

### Key design decisions:
- **Leading features always get noise added** - This simulates forecast uncertainty for training
- **No "real" forecasts for retrospective runs** - For dates where we have observations, we use observations + noise rather than calling the Forecast API
- **LightGBM handles NaN natively** - Missing values at boundaries are fine

## Next Steps

### 1. Run Retrospective Evaluation on Unity

Unity scripts are ready:
- `submit-unity-parallel.sh` - SLURM array job for all 55 dates
- `run-all-forecasts-2023-2025.sh` - Sequential local run

**To run on Unity:**
```bash
# Copy the data directory to Unity
scp -r data/ unity:~/forecast-sandbox-2025-2026/src/gbqr_weather/

# Set up the virtual environment on Unity
ssh unity
cd forecast-sandbox-2025-2026/src/gbqr_weather
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Submit the array job
sbatch submit-unity-parallel.sh
```

### 2. Compare with gbqr_3src Baseline

After running retrospective forecasts, compare:
- WIS scores between `gbqr_3src` and `gbqr_3src_weather_alpha`
- Performance by horizon, location, and season phase
- Feature importance to see if weather variables are being used

## Files Overview

| File | Purpose |
|------|---------|
| `main.py` | WeatherGBQRModel - main entry point |
| `weather_utils.py` | Static data loading, API functions, feature creation |
| `build_historical_weather.py` | One-time script to build/extend static dataset |
| `data/historical_weather.csv` | Static ERA5 weather data (2010 - Dec 2025) |
| `plan.md` | Detailed implementation plan and research notes |
| `README.md` | User documentation |
| `.weather_cache/` | Cached API responses for recent data |

## Data Details

**Static historical dataset (`data/historical_weather.csv`):**
- Date range: 1997-01-04 to 2025-12-13 (covers all surveillance data back to 1997)
- Locations: 53 (US + all states/territories)
- Total rows: 80,083
- Columns: wk_end_date, weather_temp, weather_humidity, location
- Source: Open-Meteo Archive API (ERA5 reanalysis)
- **Zero missing values** when running retrospective forecasts

## Configuration

Current model config in `main.py`:
- `sources`: ["flusurvnet", "nhsn", "ilinet"] (matches gbqr_3src)
- `weather_variables`: temperature, humidity
- `weather_lags`: [1, 2] weeks (past observations)
- `weather_leads`: [1, 2] weeks (future "forecasts" - actually observations + noise)
- `weather_noise_std`: 2.0C for temp, 10% for humidity

## Extending the Static Dataset

To extend the static weather data to cover more recent dates:

```bash
cd src/gbqr_weather
source .venv/bin/activate

# Fetch data for a new date range
python -c "
import time
import pandas as pd
from datetime import date
from build_historical_weather import fetch_historical_weather_era5
from weather_utils import STATE_CENTROIDS, aggregate_daily_to_weekly

# Adjust these dates as needed
start_date = date(2025, 12, 1)
end_date = date(2026, 1, 15)  # ERA5 has ~5 day lag

variables = ['temperature_2m_mean', 'relative_humidity_2m_mean']
all_data = []

for i, (loc, (lat, lon)) in enumerate(STATE_CENTROIDS.items()):
    print(f'[{i+1}/53] {loc}...', end=' ', flush=True)
    time.sleep(1.5)
    df = fetch_historical_weather_era5(lat, lon, start_date, end_date, variables)
    if df is not None:
        weekly = aggregate_daily_to_weekly(df)
        weekly['location'] = loc
        all_data.append(weekly)
        print(f'{len(weekly)} weeks')

if all_data:
    new_df = pd.concat(all_data, ignore_index=True)
    new_df = new_df.rename(columns={
        'temperature_2m_mean': 'weather_temp',
        'relative_humidity_2m_mean': 'weather_humidity',
    })
    new_df['wk_end_date'] = pd.to_datetime(new_df['wk_end_date']).dt.strftime('%Y-%m-%d')

    existing = pd.read_csv('data/historical_weather.csv')
    existing['wk_end_date'] = pd.to_datetime(existing['wk_end_date'], format='mixed').dt.strftime('%Y-%m-%d')

    existing_max = existing['wk_end_date'].max()
    new_df = new_df[new_df['wk_end_date'] > existing_max]

    if len(new_df) > 0:
        combined = pd.concat([existing, new_df], ignore_index=True)
        combined = combined.sort_values(['location', 'wk_end_date']).reset_index(drop=True)
        combined.to_csv('data/historical_weather.csv', index=False)
        print(f'Added {len(new_df)} rows. Total: {len(combined)}')
"
```

## API Notes

**Open-Meteo APIs used:**

1. **Archive API** (`archive-api.open-meteo.com`) - ERA5 reanalysis, 1940-present
   - Used to build static dataset
   - ~5 day lag from present
   - Rate limits: be conservative (1.5s delay between requests)

2. **Historical Forecast API** (`historical-forecast-api.open-meteo.com`) - 2022-present
   - Used for recent data beyond static dataset (when ref_date is current)
   - Archives of past weather model runs

3. **Forecast API** (`api.open-meteo.com`) - Current forecasts up to 16 days
   - NOT used in current implementation
   - Could be added for true live/production forecasting
