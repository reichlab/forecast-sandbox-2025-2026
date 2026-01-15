# GBQR Weather Features - Implementation Plan

## Overview

This document outlines a plan to integrate weather forecast data from Open-Meteo as features in the GBQR flu forecasting model. The hypothesis is that weather variables (temperature, humidity, precipitation, solar radiation) may improve flu forecast accuracy due to known associations between weather and influenza transmission.

## Data Source

**Open-Meteo API** (https://open-meteo.com/)
- Free for non-commercial/research use (CC-BY 4.0 license)
- Rate limits: 10,000 calls/day, 5,000/hour, 600/min (sufficient for our needs)
- Two relevant endpoints:
  1. **Historical Forecast API**: Archived weather data from 2022+ assembled from past forecasts.
     Provides "day 0" forecast values for each historical date (what the model predicted for
     that day at the time), but does NOT provide access to historical forecast horizons.
  2. **Forecast API**: Current predictions up to 16 days ahead (for operational inference only)

**Important limitation**: The Historical Forecast API does not archive multi-day forecast horizons.
You cannot query "what was the 7-day forecast on Dec 7, 2023." This means forecast features
(e.g., `weather_temp_fcst_1wk`) are only available for operational runs, not retrospective analysis.

### Weather Variables

| Variable | API Name | Relevance to Flu |
|----------|----------|------------------|
| Temperature | `temperature_2m_mean` | Cold weather increases indoor crowding, virus stability |
| Relative Humidity | `relative_humidity_2m_mean` | Strong link to aerosol transmission |
| Precipitation | `precipitation_sum` | Affects behavior patterns |
| Solar Radiation | `shortwave_radiation_sum` | UV inactivation, Vitamin D proxy |

## Phase 1: Local Prototype (This Directory)

### Architecture

Create a self-contained model in `src/gbqr_weather/` that:
1. Fetches weather data from Open-Meteo
2. Merges it with flu surveillance data by location and date
3. Passes augmented features to the standard GBQR model

### Files

```
src/gbqr_weather/
├── plan.md              # This file
├── main.py              # Entry point, WeatherGBQRModel class
├── weather_utils.py     # Open-Meteo API functions
├── requirements.txt     # Dependencies
└── README.md            # Usage instructions
```

### Key Design Decisions

1. **Location mapping**: Use population-weighted state centroids for weather queries
   - Weather varies within states, but a representative point is sufficient
   - Could later extend to use multiple points or HSA-level centroids

2. **Temporal alignment**: Aggregate daily weather to weekly (Saturday-ending weeks)
   - Match MMWR epidemiological week structure
   - Use mean for temperature/humidity, sum for precipitation/radiation

3. **Feature engineering**:
   - Current week weather (for horizons 0+)
   - Lagged weather (1-4 weeks prior, for training signal)
   - Weather forecasts (for prediction horizons)

4. **Caching**: Cache API responses locally to avoid redundant calls
   - Store in `.weather_cache/` directory
   - Key by location + date range

### Implementation Approach

Rather than modifying `GBQRModel.run()` directly, we:
1. Load data using `DiseaseDataLoader` (same as base model)
2. Add weather columns to the DataFrame
3. Pass weather column names in `init_feats` to `create_features_and_targets()`
4. Use standard GBQR training/prediction pipeline

This requires copying some code from `GBQRModel.run()` but keeps the core model logic unchanged.

## Phase 2: Integration into idmodels

Once the prototype is validated, integrate into `idmodels` library:

### Option A: Config-Based (Preferred)

Add weather support to existing `GBQRModel` via config flags:

```python
model_config = SimpleNamespace(
    model_class = "gbqr",
    # ... existing config ...

    # Weather features (new)
    use_weather_features = True,
    weather_variables = ["temperature_2m_mean", "relative_humidity_2m_mean"],
    weather_lags = [0, 1, 2],  # weeks
    weather_cache_dir = ".weather_cache"
)
```

Changes needed in `idmodels`:
1. Add `weather_utils.py` with API functions
2. Modify `gbqr.py` lines 66-89 to optionally add weather features (similar to wave features pattern)
3. Add `requests` to dependencies

### Option B: Subclass

Create `WeatherGBQRModel` as a subclass in `idmodels`:

```python
# idmodels/gbqr_weather.py
from idmodels.gbqr import GBQRModel

class WeatherGBQRModel(GBQRModel):
    def run(self, run_config):
        # Override to add weather features
        ...
```

This is cleaner if weather features require substantially different logic.

### Recommended Path

1. **Now**: Implement Phase 1 prototype in this directory
2. **Validate**: Run retrospective forecasts, compare skill with/without weather
3. **If beneficial**: Implement Option A in `idmodels` (config-based, follows wave features pattern)
4. **If complex**: Fall back to Option B (subclass)

## Phase 3: Evaluation

### Metrics

- WIS (Weighted Interval Score) improvement vs. base GBQR
- Coverage calibration
- Skill by horizon, location, season phase

### Ablation Studies

- Which weather variables contribute most?
- Optimal lag structure?
- Does forecast weather help beyond lagged observed weather?

## Technical Notes

### API Call Efficiency

For 53 locations × ~150 weeks of historical data:
- One API call per location covers full date range
- Total: ~53 calls for training data + 53 for current forecasts = ~106 calls
- Well within rate limits

### Handling Missing Data

- Some early dates (pre-2022) won't have historical forecast data
- Options: use reanalysis data, impute, or restrict training to 2022+
- For now: restrict to data from 2022 onwards

### Handling Training/Inference Mismatch for Leading Weather Features

A key challenge is that during training we have observed weather values, but at inference time
we only have weather forecasts (which are noisier). This creates a **covariate shift** problem
that can lead to overconfident predictions if not handled properly.

#### Research Background

Yang et al. (2017) "The use of ambient humidity conditions to improve influenza forecast"
(https://pmc.ncbi.nlm.nih.gov/articles/PMC5708837/) found that:
- Climatological humidity **outperformed** observed humidity for forecasting
- Day-to-day humidity fluctuations introduce noise that degrades predictions
- Smoothed climatological values capture the seasonal signal that matters for flu

However, since our GBQR model already includes `season_week` as a covariate, seasonality
is already captured. The goal of adding weather features is to capture **deviations from
typical seasonal patterns** that might improve predictions.

#### Our Approach: Noise Injection During Training

To make the model robust to forecast uncertainty, we inject Gaussian noise into the
observed weather values during training. This teaches the model not to rely too heavily
on precise weather values, since at inference time it will only have noisy forecasts.

**Feature types:**
1. **Lagging features** (`weather_temp_lag1`, `weather_temp_lag2`): Weather from previous weeks
   - Always populated with observed values (available at both training and inference)
   - No noise injection needed (these are always "known" at forecast time)

2. **Leading features** (`weather_temp_lead1`, `weather_temp_lead2`): Weather for future weeks
   - During training: Use observed "future" weather + Gaussian noise
   - During inference: Use actual weather forecasts from Open-Meteo API
   - Noise level calibrated to typical forecast errors (~2°C for temperature, ~10% for humidity)

**Training data structure:**

| wk_end_date | weather_temp | weather_temp_lag1 | weather_temp_lead1 | weather_temp_lead2 |
|-------------|--------------|-------------------|--------------------|--------------------|
| 2023-01-07  | 2.3          | 1.8               | 3.1 + noise        | 4.2 + noise        |
| 2023-01-14  | 3.1          | 2.3               | 4.2 + noise        | 5.0 + noise        |
| ...         | ...          | ...               | ...                | ...                |
| 2024-12-07  | 4.5          | 3.8               | [forecast]         | [forecast]         |

For the most recent weeks where we don't have "future" observations to add noise to,
we use the actual forecasts directly (no noise), since these are what we'd use operationally.

**Noise calibration:**
- Temperature: σ ≈ 2°C (typical 1-2 week forecast error)
- Relative humidity: σ ≈ 10% (typical forecast error)
- Can be tuned via cross-validation

This approach:
- Allows retrospective evaluation (we can populate leading features for all historical dates)
- Makes the model robust to forecast uncertainty
- Avoids overconfident predictions from training on perfect observations

#### References

- Yang et al. (2017): https://pmc.ncbi.nlm.nih.gov/articles/PMC5708837/
- Shaman & Karspeck (2012): https://pmc.ncbi.nlm.nih.gov/articles/PMC3528592/
- Training with noise for robustness: https://machinelearningmastery.com/train-neural-networks-with-noise-to-reduce-overfitting/

### State Centroids

Use approximate population-weighted centroids:

| State | Lat | Lon | Notes |
|-------|-----|-----|-------|
| US | 39.8 | -98.6 | Geographic center (or pop-weighted) |
| 01 (AL) | 33.5 | -86.8 | Birmingham area |
| ... | ... | ... | ... |

Full mapping in `weather_utils.py`.

## Timeline

This plan does not include time estimates. Implementation order:
1. weather_utils.py (API functions, caching, state centroids)
2. main.py (WeatherGBQRModel)
3. Testing with single forecast date
4. Retrospective evaluation
5. Decision on idmodels integration

## Alternative Data Sources for Future Consideration

The current implementation uses Open-Meteo with noise injection to simulate forecast uncertainty.
If true historical forecast archives become important, here are alternative approaches:

### Option 1: Open-Meteo Previous Runs API

Open-Meteo has a "Previous Runs API" that archives forecasts with different lead time offsets
(1, 2, 3+ days). This would allow querying actual historical forecasts rather than using
noise injection.

- **URL**: https://open-meteo.com/en/docs/historical-forecast-api (see "Previous Model Runs" section)
- **Limitation**: Data collection started in early 2024, so historical coverage is limited
- **Advantage**: Free, same API structure as current implementation
- **Implementation**: Modify `fetch_historical_weather()` to use `previous_day` parameter

### Option 2: GFS/HRRR Archives via Herbie

The Herbie Python package provides access to raw NOAA weather model archives (GFS, HRRR)
dating back to 2015-2021. This gives true historical forecasts at multiple lead times.

- **GitHub**: https://github.com/blaylockbk/Herbie
- **Data sources**: AWS, Google Cloud, NOMADS, University of Utah Pando Archive
- **Advantage**: True historical forecasts with full lead time information
- **Disadvantage**: Requires processing GRIB files, more complex implementation
- **Implementation steps**:
  1. Install Herbie: `pip install herbie-data`
  2. Query GFS/HRRR for specific forecast initialization times
  3. Extract point data for state centroids
  4. Aggregate to weekly resolution

Example Herbie usage:
```python
from herbie import Herbie

# Get GFS forecast initialized on 2023-12-01 for 7-day lead time
H = Herbie("2023-12-01", model="gfs", fxx=168)  # 168 hours = 7 days
ds = H.xarray("TMP:2 m above ground")  # 2m temperature
```

### Option 3: Visual Crossing (Paid)

Visual Crossing offers true historical forecast queries via their `forecastBasisDate` parameter,
allowing you to retrieve "what was the 7-day forecast on date X."

- **URL**: https://www.visualcrossing.com/resources/documentation/weather-data/how-to-query-weather-forecasts-from-the-past-historical-forecasts/
- **Advantage**: Simple API, exactly what we need
- **Disadvantage**: Requires paid subscription for historical forecast access

### Recommendation

For now, the noise injection approach with Open-Meteo is sufficient for initial evaluation.
If weather features prove valuable and we need more rigorous validation:
1. First try Open-Meteo Previous Runs API (free, easy)
2. If more historical depth needed, implement Herbie/GFS approach
3. Visual Crossing as fallback if budget available

## References

- Open-Meteo API Docs: https://open-meteo.com/en/docs
- Historical Forecast API: https://open-meteo.com/en/docs/historical-forecast-api
- Herbie Documentation: https://herbie.readthedocs.io/
- GFS on AWS: https://registry.opendata.aws/noaa-gfs-bdp-pds/
- Influenza-weather associations: Shaman & Kohn (2009), Lowen et al. (2007)
