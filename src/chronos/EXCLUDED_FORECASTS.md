# Chronos Model: Excluded Forecast Dates

This document tracks forecast dates that have been excluded from the Chronos model submissions due to data inconsistencies or external factors.

## Excluded Dates

### 2025-01-25

**Reason**: Data inconsistencies following Presidential Inauguration period

The week following the U.S. Presidential Inauguration (January 20, 2025) showed unusual data patterns and reporting irregularities. The forecast for reference date 2025-01-25 exhibited inconsistencies with baseline models that suggest lingering data quality issues from the inauguration period disruptions.

**Impact**:
- Data vintage and reporting timing issues may affect forecast accuracy
- Baseline models show varying horizon patterns for dates in this period
- Evaluation comparability concerns due to data inconsistencies

**Status**: Excluded from evaluation due to data quality concerns

---

### 2025-04-26

**Reason**: Inconsistent horizon coverage in baseline models

The baseline AR6 model and other models exhibit unusual behavior for the 2025-04-26 reference date, including:
- Presence of horizon -1 (nowcast) in baseline models
- Missing horizon -1 in some model submissions
- Inconsistent horizon patterns compared to standard weekly submissions

This date appears to have special circumstances that affect horizon availability across all models, not just Chronos.

**Impact**:
- Chronos generated standard horizons (-1, 0, 1, 2)
- Baseline models show non-standard horizon patterns
- Direct comparison between models is not meaningful for this date

**Status**: Excluded from evaluation pending investigation of horizon inconsistencies

---

## Summary

**Total Chronos forecasts**: 55 originally generated
**Excluded forecasts**: 2 (2025-01-25, 2025-04-26)
**Available for evaluation**: 53 forecasts

The excluded forecasts remain in the git history but have been removed from the model submission to ensure evaluation metrics are calculated on comparable data.

## Forecast Coverage

After exclusions, Chronos provides:
- **2023-24 Season**: 30 forecasts (2023-10-21 to 2024-05-11)
- **2024-25 Season**: 23 forecasts (2024-11-30 to 2025-05-17, excluding 2025-01-25 and 2025-04-26)

---

**Last updated**: 2025-11-24
**Model**: UMass-chronos_small
**Maintained by**: Thomas Robacker (trobacker@umass.edu)
