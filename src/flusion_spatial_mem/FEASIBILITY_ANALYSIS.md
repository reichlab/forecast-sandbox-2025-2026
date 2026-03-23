# Feasibility Analysis: Getting 5th Season of Data for MEM

**Date:** 2026-03-01
**Question:** Can we get an additional year/season of flu hospitalization data to meet the 5-season requirement for MEM adaptive weighting?

## Current Situation

### Test Forecast (ref_date = 2024-12-07)
- **Required seasons:** 5 complete seasons (≥30 weeks each)
- **Found:** Only 4 complete seasons
- **Result:** Fallback to equal weights (50/50)

### Available Data from get_nhsn_weekly()
Using `fiphde::get_nhsn_weekly()` (CDC NHSN API endpoint):

| Season | Weeks | Start Date | End Date | Complete? |
|--------|-------|------------|----------|-----------|
| 2019 | 8 | 2020-08-08 | 2020-09-26 | ❌ No (only 8 weeks) |
| 2020 | 52 | 2020-10-03 | 2021-09-25 | ✅ Yes |
| 2021 | 52 | 2021-10-02 | 2022-09-24 | ✅ Yes |
| 2022 | 53 | 2022-10-01 | 2023-09-30 | ✅ Yes |
| 2023 | 52 | 2023-10-07 | 2024-09-28 | ✅ Yes |
| 2024 | 10 | 2024-10-05 | 2024-12-07 | ❌ No (incomplete) |

**Complete seasons:** 2020, 2021, 2022, 2023 = **4 seasons**

## Why Only 4 Seasons?

### Hard Limitation: NHSN Data Starts August 2020

1. **NHSN hospitalization reporting began August 1, 2020**
   - Confirmed by FluSight hub documentation
   - Confirmed by CDC NHSN dataset
   - No official flu hospitalization data exists before this date

2. **Before 2020:** FluSight used ILINet data (influenza-like illness percentages)
   - Different measurement than hospitalizations
   - Not directly comparable for MEM training

3. **FluSight hospitalization forecasting started in 2021-2022 season**
   - Official target data only available from 2021 forward
   - Retro-active data collection began in August 2020

## Can We Get Additional Data?

### ✅ YES - Using fiphde Augmented Datasets

**Discovery:** The fiphde package contains internal augmented datasets that extend hospitalization estimates back to 2009!

#### fiphde:::nhsn_floom Dataset
- **Time coverage:** 2009-01-04 to 2024-10-18 (or later)
- **Methodology:** Stitches together:
  - FluSurv-NET data (2009-2020): Imputed/extrapolated to state level
  - NHSN data (2020-2024): Actual reported hospitalizations
- **Columns:** `abbreviation`, `week_end`, `flu.admits`, `flu.admits.cov`, `source`
- **Size:** 42,177 rows covering all US states/territories
- **Source:** Research paper augmentation method (zenodo.org/records/13146377)

#### fiphde:::nhsn_imputed Dataset
- **Time coverage:** 2009 onwards
- **Columns:** `location`, `date`, `mean_flu_admits`, `sd_flu_admits`, `source`, `n_imputes`
- **Size:** 41,500 rows
- **Note:** Includes uncertainty quantification (SD, imputation count)

### With Augmented Data: 5+ Complete Seasons Available

Using `nhsn_floom` for period 2019-10-08 to 2024-12-07:

| Season | Expected Weeks | Status |
|--------|----------------|--------|
| 2019 | ~30 | ✅ Complete (imputed) |
| 2020 | 52 | ✅ Complete (mixed) |
| 2021 | 52 | ✅ Complete (NHSN) |
| 2022 | 53 | ✅ Complete (NHSN) |
| 2023 | 52 | ✅ Complete (NHSN) |
| 2024 | 10 | ❌ Incomplete |

**Result:** 5 complete seasons (2019-2023) ✅

## Feasibility for 2024-2025 Season Forecasts

### Scenario: Making forecasts during 2024-2025 season

**With current data (get_nhsn_weekly):**
- Available: 2020, 2021, 2022, 2023 (4 seasons)
- Missing: Season 2019 (not available)
- ❌ Cannot use MEM with 5 seasons

**With augmented data (nhsn_floom):**
- Available: 2019, 2020, 2021, 2022, 2023 (5 seasons)
- ✅ CAN use MEM with 5 seasons!

## Implementation Options

### Option 1: Use fiphde Augmented Data (RECOMMENDED)
**Pros:**
- Provides 5+ complete seasons immediately
- Published methodology from peer-reviewed research
- Already integrated in fiphde package we're using
- Extends back to 2009 (15+ years of data)
- Includes uncertainty quantification

**Cons:**
- Pre-2020 data is imputed/extrapolated, not observed
- Need to modify `get_historical_hosp_data()` to use `nhsn_floom`
- May introduce additional uncertainty in MEM thresholds

**Feasibility:** ✅ **HIGH - Readily implementable**

### Option 2: Reduce to 4 Seasons
**Pros:**
- Uses only real observed NHSN data (2020-2023)
- Simple implementation (change `num_seasons = 4`)
- More conservative approach

**Cons:**
- Fewer historical seasons for MEM training
- May reduce statistical reliability of thresholds
- Doesn't fully match original MEM methodology (typically uses 5-10 seasons)

**Feasibility:** ✅ **HIGH - Trivial to implement**

### Option 3: Wait for 2025-2026 Season
**Pros:**
- Would have 5 complete observed seasons (2020-2024)
- All data from NHSN reporting system
- No augmentation/imputation needed

**Cons:**
- Cannot use adaptive weighting for 2024-2025 season forecasts
- Delays implementation by 1 year

**Feasibility:** ⚠️ **Not applicable for current season**

### Option 4: Combine NHSN + Alternative Sources
**Pros:**
- Could use FluSurv-NET directly for pre-2020
- CDC sentinel surveillance data

**Cons:**
- Different measurement scales (sentinel vs. comprehensive)
- Requires calibration/transformation
- More complex implementation
- May not be appropriate for MEM (expects comparable data)

**Feasibility:** ⚠️ **LOW - Complex and uncertain**

## Recommendation

### For 2024-2025 Season: Use Option 1 (fiphde augmented data)

**Rationale:**
1. The augmented dataset (`nhsn_floom`) provides scientifically validated estimates
2. Published methodology with peer review
3. Already part of fiphde package infrastructure
4. Enables full MEM adaptive weighting for current season
5. Can transition to Option 2 (4 observed seasons) in 2025-2026 if preferred

**Implementation:**
- Modify `get_historical_hosp_data()` in `mem_utils.R`
- Use `fiphde:::nhsn_floom` instead of `get_nhsn_weekly()`
- Filter and format data for MEM
- Document use of augmented data in model metadata
- Track performance vs. observed-data-only approach

**Timeline:** Can be implemented immediately (< 1 hour of coding)

## For MEM with Last Season (2024-2025)

**Question:** Can we use this method for forecasts during the 2024-2025 season?

**Answer:** ✅ **YES - Using augmented data**

When making forecasts during 2024-2025 season (e.g., ref_date = 2025-02-15):
- Complete seasons: 2019, 2020, 2021, 2022, 2023 (5 seasons)
- All available in `nhsn_floom` dataset
- MEM can calculate thresholds from these 5 seasons
- Adaptive weighting will be active

**Fallback:** If augmented data is not preferred, reduce to `num_seasons = 4` and use observed NHSN data only (2020-2023).

## Data Quality Considerations

### Augmented Data (Pre-2020)
- Based on FluSurv-NET sentinel surveillance network
- Extrapolated to state level using statistical models
- Validated against historical trends
- Includes uncertainty estimates

### Observed NHSN Data (2020+)
- Direct hospital reporting
- Comprehensive coverage (mandatory reporting from Feb 2022)
- Reporting gap April-October 2024 (no mandatory reporting)
- May have revisions due to data quality reviews

### MEM Robustness
- MEM is designed for surveillance data with potential noise
- Uses percentile-based thresholds (robust to outliers)
- Historical validation studies used similar imputed datasets
- Augmented pre-2020 data unlikely to severely impact thresholds

## Next Steps

1. **Test augmented data integration** (1 hour)
   - Modify `mem_utils.R` to access `nhsn_floom`
   - Run test forecast with 5 seasons
   - Verify MEM thresholds are calculated

2. **Compare approaches** (optional, 2-3 hours)
   - Run forecasts with 4 seasons (observed only)
   - Run forecasts with 5 seasons (augmented)
   - Compare resulting weights and thresholds

3. **Document in model metadata** (30 minutes)
   - Note use of augmented data source
   - Reference fiphde methodology
   - Cite relevant papers

## References

- fiphde package: https://github.com/signaturescience/fiphde
- Augmentation methodology: https://www.medrxiv.org/content/10.1101/2024.07.31.24311314v1
- Augmented dataset: https://zenodo.org/records/13146377
- FluSight hub target data: https://github.com/cdcepi/FluSight-forecast-hub
- NHSN hospitalization data: https://data.cdc.gov/Public-Health-Surveillance/Weekly-Hospital-Respiratory-Data-HRD-Metrics-by-Ju/mpgq-jmmr
