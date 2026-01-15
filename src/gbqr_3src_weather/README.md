# GBQR Weather Model

GBQR flu forecasting model augmented with weather features from Open-Meteo.

## Overview

This model extends the base GBQR model by incorporating weather covariates:
- Temperature (daily mean)
- Relative humidity (daily mean)

Weather data is fetched from the [Open-Meteo API](https://open-meteo.com/), which provides:
- Historical forecast data from 2022 onwards
- Current forecasts up to 16 days ahead

## Feature Types

The model creates three types of weather features:

### 1. Current Week Weather
- `weather_temp` - Temperature for the current week
- `weather_humidity` - Humidity for the current week

### 2. Lagging Features (past weather)
These represent observed weather from previous weeks:
- `weather_temp_lag1` - Temperature from 1 week prior
- `weather_temp_lag2` - Temperature from 2 weeks prior
- `weather_humidity_lag1`, `weather_humidity_lag2` - etc.

### 3. Leading Features (future weather)
These represent weather for future weeks:
- `weather_temp_lead1` - Weather for 1 week ahead
- `weather_temp_lead2` - Weather for 2 weeks ahead
- `weather_humidity_lead1`, `weather_humidity_lead2` - etc.

## Handling Training/Inference Mismatch

A key challenge is that during training we have perfect observed weather values, but at
inference time we only have weather forecasts (which are noisier). This creates a
**covariate shift** problem that can lead to overconfident predictions.

### Our Approach: Noise Injection

To make the model robust to forecast uncertainty:

1. **Lagging features**: Always use observed values (no noise needed - these are known at forecast time)

2. **Leading features during training**: Use observed "future" weather + Gaussian noise
   - Temperature noise: σ ≈ 2°C (typical 1-2 week forecast error)
   - Humidity noise: σ ≈ 10% (typical forecast error)

3. **Leading features at inference**: Use actual forecasts from Open-Meteo API

This teaches the model not to rely too heavily on precise weather values, since at
inference time it will only have noisy forecasts.

### Research Background

This approach is informed by Yang et al. (2017) "The use of ambient humidity conditions
to improve influenza forecast" ([PMC5708837](https://pmc.ncbi.nlm.nih.gov/articles/PMC5708837/)),
which found that climatological humidity outperformed observed humidity for forecasting,
likely because day-to-day fluctuations introduce noise. Adding noise during training
is a standard technique to improve model robustness
([reference](https://machinelearningmastery.com/train-neural-networks-with-noise-to-reduce-overfitting/)).

## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Usage

```bash
# Run for a specific date
python main.py --today_date=2024-12-07

# Short run (fewer bags/quantiles, for testing)
python main.py --today_date=2024-12-07 --short_run
```

## Configuration

Weather features are controlled via `model_config` in `main.py`:

```python
model_config = SimpleNamespace(
    # ... standard GBQR config ...

    # Weather feature config
    use_weather_features=True,
    weather_variables=[
        "temperature_2m_mean",
        "relative_humidity_2m_mean",
    ],
    weather_lags=[1, 2],      # 1 and 2 weeks prior
    weather_leads=[1, 2],     # 1 and 2 weeks ahead
    weather_cache_dir=".weather_cache",

    # Noise injection for leading features
    weather_noise_std={
        "weather_temp": 2.0,       # ~2°C typical forecast error
        "weather_humidity": 10.0,  # ~10% typical forecast error
    },

    # Variables to use for forecasts at inference time
    forecast_variables=[
        "temperature_2m_mean",
        "relative_humidity_2m_mean",
    ],
)
```

## Available Weather Variables

From Open-Meteo API:
- `temperature_2m_mean` - Mean daily temperature
- `temperature_2m_max` / `temperature_2m_min` - Max/min temperature
- `apparent_temperature_mean` - "Feels like" temperature
- `relative_humidity_2m_mean` - Mean relative humidity
- `dewpoint_2m_mean` - Dewpoint temperature
- `precipitation_sum` - Total precipitation (mm)
- `shortwave_radiation_sum` - Solar radiation (MJ/m²)
- `uv_index_max` - Max UV index

## Files

- `main.py` - Entry point and WeatherGBQRModel class
- `weather_utils.py` - Open-Meteo API functions and utilities
- `plan.md` - Implementation plan and integration roadmap
- `.weather_cache/` - Cached API responses (created automatically)

## How It Works

1. Load surveillance data (same as base GBQR)
2. Fetch historical weather data for each state (using population-weighted centroids)
3. Aggregate daily weather to weekly (Saturday-ending weeks)
4. Merge weather data with surveillance data by location and date
5. Create lagging features (past weather - no noise)
6. Create leading features (future weather + noise for training rows)
7. For recent rows near the reference date, use actual forecasts instead of noisy observations
8. Pass augmented features to standard GBQR training/prediction

## Notes

- Weather data is cached locally in `.weather_cache/` to avoid redundant API calls
- Historical forecast API only has data from ~2022; earlier dates use median imputation
- Forecast API provides up to 16 days ahead (~2 weeks of weekly forecasts)
- Rate limits: 10,000 API calls/day (sufficient for our needs)
- Free for non-commercial research use (CC-BY 4.0 license)
- Random seed for noise is based on reference date for reproducibility

## API Clarification

Open-Meteo provides two relevant APIs:

1. **Historical Forecast API** (`historical-forecast-api.open-meteo.com`): Provides archived
   weather data assembled from past forecasts. This gives you what the forecast model predicted
   for each historical day, but concatenates only the initial forecast hours - essentially
   providing a continuous time series of "day 0" predictions, not multi-day forecast horizons.

2. **Forecast API** (`api.open-meteo.com`): Provides current weather forecasts up to 16 days
   ahead. Only works for today/future dates.

Neither API provides access to "what was the N-day forecast on historical date X" - this would
require a service that archives full forecast horizons over time. Our noise injection approach
addresses this limitation.

## Future Work

See `plan.md` for the roadmap to integrate this into `idmodels` as a proper subclass.
