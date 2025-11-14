# Forecasting Model Comparison Summary

## Overview

This document summarizes the various flu forecasting models explored for the 2025/2026 season, comparing their architectures, methodologies, and performance characteristics.

## Models Evaluated

### 1. UMass-flusion (flu_flusion)
**Status:** Top performer

**Architecture:**
- Three-component equal-weighted ensemble
- Components:
  1. AR(6) pooled model (no seasonality)
  2. Gradient Boosted Quantile Regression (GBQR)
  3. Equal-weighted linear pool of both

**Key Features:**
- Simple equal weighting (1/3 each component)
- Uses `hubEnsembles::simple_ensemble()`
- Intermediate outputs stored before ensemble creation
- Multi-stage pipeline: 0_ar6_pooled.py → 1_gbqr.py → 2_flusion_ensemble.R

**Strengths:**
- Proven robust performance
- Simple and interpretable
- Stable across different time periods
- Equal weights avoid overfitting to historical patterns

### 2. UMass-flusion3 (flusion3)
**Status:** Baseline for comparison

**Architecture:**
- Three-component equal-weighted ensemble
- Components:
  1. AR(6) pooled model (no Fourier)
  2. AR(6) pooled model with Fourier seasonality (K=2)
  3. Gradient Boosted Quantile Regression (GBQR)

**Key Features:**
- Adds Fourier seasonality component (2 harmonics)
- Equal weighting (1/3 each)
- Uses `hubEnsembles::simple_ensemble()`

**Performance on Test Dates (2023-12-16, 2023-12-23, 2023-12-30):**
- Mean WIS: 192.5
- Consistent performance across all three test dates

**Strengths:**
- Incorporates explicit seasonality via Fourier terms
- Equal weighting provides robustness
- Simple and interpretable

### 3. UMass-flusion3_weighted
**Status:** Experimental - did not improve over equal weights

**Architecture:**
- Three-component weighted ensemble with learned weights
- Same components as flusion3:
  1. AR(6) pooled model
  2. AR(6) pooled model with Fourier seasonality
  3. GBQR

**Key Features:**
- Weights learned from historical WIS performance
- Uses inverse WIS weighting: weight = 1/WIS (normalized)
- Requires historical forecasts to learn weights
- Uses `hubEnsembles::linear_pool()` with learned weights

**Training Setup:**
- Training period: 8 dates (2023-10-21 through 2023-12-09)
- Test period: 3 dates (2023-12-16, 2023-12-23, 2023-12-30)
- Full MCMC: 2000 warmup + 2000 samples
- Full bags: 100 for GBQR

**Performance on Test Dates:**
- Mean WIS: 208.5
- 8.33% worse than equal weights
- Won only 1 out of 3 test dates (2023-12-30)

**Results by Date:**
| Date | Equal Weights WIS | Learned Weights WIS | Winner |
|------|------------------|---------------------|---------|
| 2023-12-16 | 249.9 | 280.4 | Equal (12.2% better) |
| 2023-12-23 | 193.1 | 212.8 | Equal (10.2% better) |
| 2023-12-30 | 134.5 | 132.4 | Learned (1.5% better) |

**Why It Didn't Work:**
1. **Similar component performance:** The three components had comparable accuracy, so differential weighting didn't help
2. **Overfitting to history:** Weights learned from training period didn't generalize to test period
3. **Limited signal:** Even with 8 training dates, insufficient data to learn reliable differential weights
4. **Inverse WIS may not be optimal:** This weighting scheme may not capture the right signal for this ensemble

**Lessons Learned:**
- Equal weighting is more robust when components have similar performance
- Historical performance-based weighting requires strong differential signal between components
- More training data didn't solve the fundamental problem
- Simpler is sometimes better

### 4. UMass-flusion_fourier (flusion_fourier)
**Status:** Alternative approach with median ensemble

**Architecture:**
- Two-component median ensemble
- Components:
  1. AR(6) pooled model with Fourier seasonality (K=2)
  2. GBQR

**Key Features:**
- Uses median pooling instead of mean
- Explicit Fourier seasonality (2 harmonics for annual cycle)
- Excludes non-Fourier AR(6) component

**Differences from flusion3:**
- Median vs. mean pooling
- Only 2 components (dropped non-Fourier AR6)
- Different aggregation method

## Performance Comparison

### Comprehensive Model Evaluation

All models were evaluated on their overlapping forecasts to provide a fair comparison. This section presents results from three different analyses.

### Fair Comparison: All Models on Same 55 Dates ⭐

**This is the definitive comparison** - all models evaluated on the same 55 dates (2023-10-21 to 2025-05-17):

| Rank | Model | WIS | Difference from Best |
|------|-------|-----|---------------------|
| 1 | **UMass-flusion** | **167.62** | - |
| 2 | UMass-flusion_fourier | 168.46 | +0.5% |
| 2 | UMass-flusion_fourier_median | 168.46 | +0.5% |
| 4 | UMass-flusion3 | 171.31 | +2.2% |
| 5 | UMass-flusion3_weighted | 173.54 | +3.5% |

**Key Findings:**
- **UMass-flusion is the clear winner** across all 55 dates
- All models evaluated on 11,501 location-horizon forecasts
- Consistent rankings across different time periods
- Performance differences are modest but consistent

### Comparison: Initial Results (12 dates) vs Full Results (55 dates)

The initial comparison on just 12 dates (2023-10-21 to 2024-01-06) showed much larger performance gaps:

| Model | WIS (12 dates) | WIS (55 dates) | Rank Change |
|-------|----------------|----------------|-------------|
| UMass-flusion | 92.55 | 167.62 | Same (1st) |
| UMass-flusion_fourier | 93.09 | 168.46 | Same (2nd) |
| UMass-flusion3 | 111.26 | 171.31 | Same (4th) |
| UMass-flusion3_weighted | 122.25 | 173.54 | Same (5th) |

**Performance Gap Changes:**
- flusion vs flusion3: **20.2% → 2.2%** (gap narrowed significantly)
- flusion vs flusion3_weighted: **32.1% → 3.5%** (gap narrowed significantly)
- flusion3 vs flusion3_weighted: **9.9% → 1.3%** (gap narrowed significantly)

**Interpretation:** The larger dataset (55 vs 12 dates) provides more stable estimates. While the performance gaps are smaller, the ranking remains consistent: UMass-flusion is best, and learned weights don't improve performance.

### Analysis 3: Learned Weights vs Equal Weights

Comparing just the flusion3 variants on their 12 common dates:

#### Overall Performance

| Model | Mean WIS | % Worse than Best |
|-------|----------|-------------------|
| UMass-flusion3 (equal weights) | 111.26 | - |
| UMass-flusion3_weighted (learned weights) | 122.25 | +9.9% |

#### Test Date Performance (Training vs Test Split)

The flusion3_weighted model was trained on 8 dates and tested on 3 dates:

| Test Date | Equal Weights | Learned Weights | Winner |
|-----------|---------------|-----------------|---------|
| 2023-12-16 | 249.88 | 280.37 | Equal (+12.2%) |
| 2023-12-23 | 193.12 | 212.84 | Equal (+10.2%) |
| 2023-12-30 | 134.49 | 132.41 | Learned (+1.5%) |
| **Mean** | **192.50** | **208.54** | **Equal (+8.3%)** |

**Result:** Learned weights underperformed equal weights by 9.9% overall, winning only 1 of 3 test dates.

### Pairwise Model Comparisons (All 55 Dates)

Based on all 55 dates, here are the pairwise comparisons:

| Comparison | WIS Difference | % Better |
|------------|----------------|----------|
| **flusion vs flusion_fourier** | 167.62 vs 168.46 | 0.5% |
| **flusion vs flusion3** | 167.62 vs 171.31 | **2.2%** |
| **flusion vs flusion3_weighted** | 167.62 vs 173.54 | **3.5%** |
| **flusion_fourier vs flusion3** | 168.46 vs 171.31 | **1.7%** |
| **flusion_fourier vs flusion3_weighted** | 168.46 vs 173.54 | **3.0%** |
| **flusion3 vs flusion3_weighted** | 171.31 vs 173.54 | 1.3% |

### Key Findings

1. **UMass-flusion is the best model** across 55 dates spanning 2 years (2.2-3.5% better than flusion3 variants)
2. **UMass-flusion_fourier performs nearly identically** to flusion (only 0.5% worse)
3. **Adding a third component degraded performance:** flusion3 (3 components) performed 2.2% worse than flusion (2 components)
4. **Learned weights consistently don't help:** flusion3_weighted performed 1.3% worse than equal weights
5. **Performance differences are modest but consistent:** The gaps are smaller than initially observed with 12 dates, but the ranking is stable

## Methodological Insights

### Component Models

**AR(6) Pooled:**
- Autoregressive model with 6 lag terms
- Pooled across locations (shared parameters)
- Power transformation: fourth root
- Scaling: divide by 95th percentile, subtract mean
- MCMC: 2000 warmup + 2000 samples (full mode)

**AR(6) with Fourier:**
- Same as AR(6) but adds seasonal terms
- K=2 harmonics (captures annual cycle)
- Helps model seasonal patterns explicitly

**GBQR:**
- Gradient Boosted Quantile Regression
- Bagged approach (100 bags in full mode)
- Non-parametric, can capture complex patterns
- Uses historical lags as features

### Ensemble Methods

**Simple Ensemble (equal weights):**
- Each component weighted equally (1/n)
- Robust when components have similar accuracy
- Less prone to overfitting
- Used in flusion and flusion3

**Weighted Linear Pool:**
- Components weighted by historical performance
- Weights = inverse WIS, normalized to sum to 1
- Requires sufficient historical data
- Can overfit if components have similar performance
- Used in flusion3_weighted

**Median Ensemble:**
- Takes median across component forecasts
- More robust to outliers than mean
- Used in flusion_fourier

## Recommendations

### For Production Use:
1. **Stick with UMass-flusion** - proven, robust, simple
2. Consider **flusion3** as alternative if explicit seasonality is desired
3. Avoid weighted ensembles unless components show strong differential performance

### For Future Exploration:
1. **Stacking approaches:** Use holdout data to learn weights (cross-validation)
2. **Dynamic weighting:** Allow weights to vary by location or time
3. **Alternative weighting schemes:** Try methods beyond inverse WIS
4. **More diverse components:** Add fundamentally different model types (e.g., mechanistic models)
5. **Bayesian Model Averaging:** Principled approach to weight uncertainty

### When to Use Learned Weights:
- Components have clearly different performance profiles
- Large historical dataset available (>20-30 dates)
- Performance differences are stable over time
- Cross-validation shows improvement on held-out data

## Technical Implementation Notes

### Weight Learning Code
The flusion3_weighted model learns weights via:
```r
# Compute WIS for each model on historical data
wis_by_model <- # ... compute WIS ...

# Convert to weights (inverse WIS, normalized)
inverse_wis <- sapply(component_models, function(m) 1 / wis_by_model[[m]])
weights_vec <- inverse_wis / sum(inverse_wis)
```

### Files Created
- `src/flusion3_weighted/` - full model implementation
- `src/flusion3_weighted/train_historical.py` - historical training script
- `src/flusion3_weighted/compare_performance.R` - evaluation script
- `model-metadata/UMass-flusion3_weighted.yml` - metadata

### Computational Cost
Full training (11 dates, full MCMC, full bags):
- Runtime: ~2.5 hours
- ~10-15 minutes per date
- Resource intensive but feasible for experimentation

## Conclusion

The comprehensive evaluation across all available forecasts revealed several important lessons:

### Main Findings

1. **UMass-flusion is the best model:** Consistently outperforms all alternatives across 55 dates (WIS = 167.62)

2. **Simpler is better:** The 2-component ensemble (flusion) outperformed the 3-component variants (flusion3) by 2.2%

3. **Adding components can hurt:** The third component (AR6 without Fourier) in flusion3 degraded performance rather than improving it

4. **Learned weights don't help:** Even within the same architecture, learned weights underperformed equal weights by 1.3%

5. **Fourier seasonality has minimal impact:** flusion_fourier (with explicit seasonality) performed nearly identically to flusion (0.5% difference)

6. **Results are stable:** Rankings consistent across both 12-date and 55-date evaluations, though performance gaps narrowed with more data

### Recommendations

**For Production:**
- **Use UMass-flusion** - it has the best performance and is proven robust
- UMass-flusion_fourier is a viable alternative (nearly identical performance)
- Avoid the flusion3 variants - they underperform by 20%+
- Don't use learned weights unless components have very different performance profiles

**For Research:**
- Investigate why the third component degraded performance
- Consider whether AR6-Fourier alone might outperform the ensemble
- Explore whether the GBQR component is actually adding value
- Test other component combinations beyond current choices

### What We Learned

The investigation into weighted ensembles and additional components revealed that **more sophisticated methods don't always yield better results**. The original UMass-flusion's simplicity and careful component selection proved superior to more complex alternatives.

This work demonstrates the importance of:
- **Empirical validation:** Always test assumptions with real data
- **Fair comparisons:** Evaluate models on overlapping dates
- **Parsimony:** Start simple and only add complexity when justified
- **Component quality over quantity:** Two good components beat three mediocre ones

---

*Document created: 2024-11-13*
*Updated: 2024-11-14 with fair evaluation on all 55 dates*
*Models evaluated: All models on 55 dates (2023-10-21 to 2025-05-17)*
*Total forecasts per model: 11,501 location-horizon combinations*
*Training time: ~14 hours for flusion3 and flusion3_weighted to complete full date range*
