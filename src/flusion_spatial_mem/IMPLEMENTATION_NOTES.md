# Implementation Notes: Flusion Spatial MEM

This document describes the implementation of the MEM-based adaptive ensemble weighting approach for the flusion_spatial_mem model.

## Overview

This model implements **Option 1** from the MEM integration research: **MEM-Based Adaptive Ensemble Weights**.

The key innovation is using the Moving Epidemic Method (MEM) to detect the current epidemic phase and intensity, then adaptively weighting ensemble components based on which models perform better in different epidemic phases.

## Architecture

### Component Models (unchanged from flusion_spatial2_prod)
1. **AR(6) pooled** (`0_ar6_pooled.py`) - SARIX model with shared AR coefficients
2. **GBQR 3-source spatial2** (`1_gbqr_3src_spatial2.py`) - GBQR with directional wave features for states
3. **GBQR 3-source** (`2_gbqr_3src.py`) - GBQR without spatial features for US-level

### Ensemble Script (NEW)
**`3_flusion_mem_ensemble.R`** - Implements MEM-based adaptive weighting

## Key Implementation Details

### Current Status (v0.2) - **MEM ADAPTIVE WEIGHTING FULLY IMPLEMENTED** ✅

The ensemble script now includes FULL MEM-based adaptive weighting:

✅ **Complete MEM Integration**
- `mem_utils.R`: Comprehensive helper functions for all MEM operations
- `get_historical_hosp_data()`: Fetches NHSN data using fiphde package
- `prepare_mem_data()`: Formats data for MEM package (seasons x weeks matrix)
- `calculate_mem_thresholds()`: Computes epidemic and intensity thresholds
- `determine_epidemic_phase()`: Classifies current epidemic status
- `get_phase_weights()`: Returns phase-specific ensemble weights
- `weighted_linear_pool()`: Applies adaptive weights to component forecasts

✅ **Adaptive Weighting Now Active**
- **NO LONGER using simple median** - now uses MEM-based weighted ensemble
- Weights adjust based on epidemic phase (baseline → very_high)
- Diagnostic output file tracks phase and weights for each forecast
- Automatic fallback to equal weights if MEM data unavailable

### MEM Workflow (FULLY IMPLEMENTED AND ACTIVE) ✅

```r
# 1. Calculate MEM thresholds from historical data
mem_thresholds <- calculate_mem_thresholds("US", ref_date, num_seasons = 5)

# 2. Determine current epidemic phase
current_phase <- determine_epidemic_phase("US", ref_date, mem_thresholds)

# 3. Get phase-specific weights
weights <- get_phase_weights(current_phase)
# Example outputs:
# - baseline:   ar6=0.5, gbqr=0.5
# - low:        ar6=0.45, gbqr=0.55
# - medium:     ar6=0.4, gbqr=0.6
# - high:       ar6=0.35, gbqr=0.65
# - very_high:  ar6=0.3, gbqr=0.7

# 4. Create weighted ensemble (NOW ACTIVE!)
ens_state <- weighted_linear_pool(state_dat, weights)
ens_us <- weighted_linear_pool(us_dat, weights)
ensemble <- bind_rows(ens_state, ens_us)

# 5. Save diagnostic information
# Saves phase, weights, and thresholds to diagnostic CSV
```

### Phase-Specific Weighting Strategy

| Phase | AR6 Weight | GBQR Weight | Rationale |
|-------|-----------|-------------|-----------|
| Baseline | 0.5 | 0.5 | Equal weights for stable periods |
| Low | 0.45 | 0.55 | Slight preference for GBQR as activity increases |
| Medium | 0.4 | 0.6 | Balanced but favoring GBQR for growth patterns |
| High | 0.35 | 0.65 | More weight to GBQR for rapid changes |
| Very High | 0.3 | 0.7 | Highest GBQR weight during peak activity |

**Note**: These weights are initial estimates and should be refined through validation on historical data.

## Implementation Status

### ✅ COMPLETED

1. **Data Pipeline Integration**
   - ✅ Implemented `get_historical_hosp_data()` using fiphde package
   - ✅ Created proper data formatting for MEM package (`prepare_mem_data()`)
   - ✅ Handle edge cases (insufficient historical data, missing values)

2. **MEM Calculation**
   - ✅ Implemented `calculate_mem_thresholds()` using mem R package
   - ✅ Implemented `determine_epidemic_phase()` for phase classification
   - ✅ Implemented `get_phase_weights()` for phase-specific weights

3. **Activation**
   - ✅ Replaced `simple_ensemble()` with `weighted_linear_pool()`
   - ✅ Added logging of phase and weights to diagnostic output file
   - ✅ Diagnostic CSV tracks phase, weights, and thresholds for each forecast

### 🔄 NEXT STEPS (Testing & Optimization)

1. **Installation & Setup**
   - [ ] Install fiphde package: `remotes::install_github('signaturescience/fiphde')`
   - [ ] Run `renv::restore()` to install all R dependencies
   - [ ] Verify Python environment has idmodels installed

2. **Initial Testing**
   - [ ] Run single test forecast: `python main.py --today_date=2024-12-07 --short_run`
   - [ ] Verify MEM thresholds are calculated correctly
   - [ ] Check diagnostic output file for phase detection
   - [ ] Inspect weights applied to ensemble

3. **MEM Validation**
   - [ ] Test MEM threshold calculation on historical data
   - [ ] Verify phase classification matches known epidemic periods
   - [ ] Compare MEM thresholds with CDC/public health definitions
   - [ ] Validate fiphde data quality and completeness

4. **Weight Optimization**
   - [ ] Run validation study on past 2-3 flu seasons
   - [ ] Calculate performance metrics (WIS, coverage, MAE) by phase
   - [ ] Optimize weights based on empirical performance
   - [ ] Consider location-specific vs. shared weights
   - [ ] Test sensitivity to weight variations

5. **Full Evaluation**
   - [ ] Run `evaluate_mem_ensemble.py` for retrospective forecasts (2023-2025)
   - [ ] Compare against flusion_spatial2_prod baseline
   - [ ] Analyze performance by epidemic phase
   - [ ] Create diagnostic plots showing phase transitions over time
   - [ ] Generate performance summary by phase and horizon

## Testing the MEM Implementation

The MEM-based adaptive weighting is now fully implemented and active!

```bash
# 1. Set up environment
cd src/flusion_spatial_mem
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
Rscript -e "renv::restore()"

# 2. Test single forecast
python main.py --today_date=2024-12-07 --short_run

# 3. Run evaluation (dry run to see dates)
python evaluate_mem_ensemble.py \
  --start_date=2024-10-01 \
  --end_date=2024-12-31 \
  --dry_run

# 4. Run actual evaluation
python evaluate_mem_ensemble.py \
  --start_date=2024-10-01 \
  --end_date=2024-12-31 \
  --short_run
```

## Technical Considerations

### Why Start with Simple Median?

1. **Baseline Comparison**: Allows direct comparison with flusion_spatial2_prod
2. **Pipeline Testing**: Verifies all components work before adding complexity
3. **Incremental Development**: MEM integration can be refined separately
4. **Risk Management**: Easy to roll back if MEM implementation needs revision

### MEM Package Dependencies

The `mem` package requires:
- R >= 3.4.0
- Packages: sm, RColorBrewer, rootSolve, reshape
- Historical seasonal data in specific format (matrix: seasons × weeks)

### Alternative Ensemble Methods to Consider

If simple weighted average proves suboptimal:
- **Quantile regression averaging** (QRA)
- **Stacking** with phase as a feature
- **Dynamic model averaging** (DMA)
- **Copula-based combination**

## References

### Moving Epidemic Method
- Vega T, Lozano JE, et al. (2013). Influenza surveillance in Europe: establishing epidemic thresholds by the moving epidemic method. *Influenza Other Respir Viruses*, 7(4):546-558.
- MEM R Package: https://cran.r-project.org/package=mem

### Ensemble Forecasting
- Reich NG, et al. (2019). A collaborative multi-model ensemble for real-time influenza season forecasting in the U.S. *PLOS Computational Biology*.
- Ray EL, Reich NG (2018). Prediction of infectious disease epidemics via weighted density ensembles. *PLOS Computational Biology*.

## Contact

For questions or contributions to this experimental model:
- Thomas Robacker (trobacker@umass.edu)
- Evan Ray (elray@umass.edu)
- Nicholas Reich (nick@umass.edu)
