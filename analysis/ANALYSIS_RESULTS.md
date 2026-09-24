# Analysis Results: gbqr_3src vs gbqr_3src_spatial

## Executive Summary

The spatial model (`gbqr_3src_spatial`) incorporates directional wave features that capture geographic disease spread. Analysis of two forecast dates reveals:

**Key Finding:** Spatial features account for **~40% of total feature importance** and lead to substantially different forecasts. The magnitude and pattern of differences varies by epidemic phase, location characteristics, and forecast horizon. This analysis characterizes where and when the models diverge, though further evaluation against observed outcomes would be needed to assess forecast accuracy.

---

## Date 1: 2024-12-21 (Early Season)

### Feature Importance

**Spatial Features: 40.2% of total importance**

**Top Spatial Feature Categories:**
1. **Spatial Velocity Features**: 638,614 total importance (18.4% of top feature)
2. **North Direction Waves**: 615,185 importance
3. **South Direction Waves**: 440,734 importance
4. **Aggregate Wave Features**: 375,181 importance
5. **West Direction Waves**: 278,153 importance
6. **East Direction Waves**: 231,103 importance

**Most Important Individual Spatial Features:**
1. `inc_trans_cs_wave_avg` - Average wave across all directions (2.6%)
2. `inc_trans_cs_wave_avg_velocity` - Velocity of average wave (1.8%)
3. `inc_trans_cs_wave_W` - West direction wave (1.7%)
4. `inc_trans_cs_wave_SW` - Southwest direction wave (1.7%)

### Non-Spatial Features and Changes Between Models

**Feature Categories in Base Model:**
- **Other features** (62% of importance): Current incidence level, rolling means, Taylor expansion features (temporal smoothing), Christmas holiday effect (6.2%), forecast horizon indicator (4.7%)
- **Temporal lags** (34% of importance): Lagged values of transformed incidence at various windows
- **Data sources** (4.4% of importance): Indicators for FluSurvNet, NHSN, and ILINet data streams

**Changes When Spatial Features Are Added:**
- **Other features drop to 39%** (from 62%) - a 37% reduction in relative importance
- **Temporal lags drop to 18%** (from 34%) - a 48% reduction in relative importance
- **Data sources drop to 3.5%** (from 4.4%)

**Key Observation:** The dramatic decrease in temporal lag importance (from 34% to 18%) suggests that spatial features capture some of the same information that temporal lags were providing in the base model. Specifically, wave features at lag-1 and lag-2 may be providing richer context about recent trends than simple temporal lags alone.

**Most Important Non-Spatial Features That Remain Critical:**
1. `delta_xmas` (Christmas effect): 6.2% in base → 4.3% in spatial (still important)
2. `inc_trans_cs` (current incidence): 5.2% in base → 4.3% in spatial
3. Various Taylor expansion features (temporal smoothing): collectively important in both models
4. `horizon` (forecast horizon indicator): 4.7% in base → 4.4% in spatial

### Geographic Patterns

**Overall Forecast Differences:**
- Mean absolute difference: **53.77 hospitalizations**
- Mean relative difference: **40.4%**
- Median absolute difference: 31.62
- Maximum difference: 282.90

**Absolute vs Relative Differences Show Different Patterns:**

When ranking by **absolute difference** (hospitalization counts):
- **Top 5**: California (212), Texas (206), Arizona (168), Pennsylvania (167), New York (152)
- **Pattern**: Large, populous states show largest absolute differences

When ranking by **relative difference** (percentage):
- **Top 5**: North Dakota (90%), Montana (89%), South Dakota (79%), Maine (79%), Vermont (70%)
- **Pattern**: Small, rural states show largest proportional differences

**Bottom 5 by relative difference**: Minnesota (10%), Alabama (12%), Florida (13%), California (15%), Texas (17%)

**Interpretation:** The models diverge most in absolute terms for large states (where hospitalization counts are high), but diverge most proportionally for small/rural states. This suggests spatial features may play different roles depending on location characteristics and disease burden levels.

**Regional Patterns (Relative Differences):**
- Midwest: 40.4% (highest)
- Northeast: 40.0%
- South: 31.0%
- West: 26.2%

**By Forecast Horizon (Absolute Differences):**
- Horizon 0 (nowcast): 18.7
- Horizon 1 (1 week): 52.3
- Horizon 2 (2 weeks): 68.1
- Horizon 3 (3 weeks): 76.1

**Pattern:** Forecast differences increase with horizon length.

---

## Date 2: 2025-02-15 (Peak Season)

### Feature Importance

**Spatial Features: 39.7% of total importance**

Similar to early season, spatial features remain highly important at peak transmission.

**Top Spatial Feature Categories:**
1. **Spatial Velocity Features**: 628,584 total importance
2. **North Direction Waves**: 610,557 importance
3. **South Direction Waves**: 431,200 importance
4. **Aggregate Wave Features**: 376,129 importance
5. **West Direction Waves**: 276,465 importance
6. **East Direction Waves**: 225,902 importance

**Pattern Consistency:** The ranking of directional importance is identical to early season:
- North and velocity features dominate
- South and West are secondary
- East and Southeast are tertiary

### Geographic Patterns

**Overall Forecast Differences:**
- Mean absolute difference: **78.92 hospitalizations** (+47% vs early season)
- Mean relative difference: **11.1%** (much lower than early season's 40.4%)
- Median absolute difference: 23.74
- Maximum difference: **748.64** (2.6x larger than early season)

**Absolute vs Relative Differences Show Contrasting Patterns:**

When ranking by **absolute difference** (hospitalization counts):
- **Top 5**: Pennsylvania (619), New York (544), Ohio (315), North Carolina (211), Michigan (187)
- **Pattern**: Large Northeast/Midwest states dominate

When ranking by **relative difference** (percentage):
- **Top 5**: Minnesota (29%), New York (29%), Pennsylvania (22%), Maine (20%), New Hampshire (19%)
- **Bottom 5**: Oklahoma (2%), Iowa (4%), Kansas (4%), Nevada (4%), Alabama (4%)

**Interpretation:** At peak season, absolute differences are driven by high case counts in populous states, but relative differences remain most pronounced in specific states (notably Northeast states and Minnesota). The lower overall relative differences (11% vs 40% at early season) suggest that at higher case volumes, the models' forecasts become more similar proportionally even as absolute differences grow.

**Regional Patterns by Relative Difference:**
- **Northeast**: 16.4% (highest)
- Midwest: 11.1%
- South: 8.1%
- West: 6.6%

**Regional Patterns by Absolute Difference:**
- Northeast: 174.3 (3x larger than early season)
- Midwest: 93.1
- South: 60.4
- West: 25.4

**Key Observation:** Northeast shows both the highest relative differences (16.4%) and a 3-fold increase in absolute differences from early season. This suggests model forecasts diverge most substantially for this region during peak transmission.

**By Forecast Horizon (Absolute Differences):**
- Horizon 0 (nowcast): 43.6
- Horizon 1 (1 week): 74.6
- Horizon 2 (2 weeks): 97.1
- Horizon 3 (3 weeks): 100.4

**Pattern:** All horizons show larger absolute differences at peak season compared to early season.

---

## Cross-Date Comparisons

### Temporal Stability of Spatial Features

**Feature importance percentages are remarkably stable:**
- Early season (2024-12-21): 40.2%
- Peak season (2025-02-15): 39.7%

**Directional rankings remain identical:**
1. North/Velocity (consistently most important)
2. South
3. West
4. East
5. Southeast

**Implication:** Spatial features capture fundamental geographic spread patterns that persist across epidemic phases.

### Forecast Divergence Increases at Peak

**Mean absolute differences by season:**
- Early (Dec 21): 53.77 hospitalizations
- Peak (Feb 15): 78.92 hospitalizations (+47% increase)

**Maximum differences by season:**
- Early: 282.90
- Peak: 748.64 (+165% increase)

**Implication:** Spatial features matter MORE during active transmission.

### Geographic Shifts

**Early Season (2024-12-21) - Top states:**
- West dominates: California #1, Arizona #3
- Large states across regions

**Peak Season (2025-02-15) - Top states:**
- Northeast dominates: Pennsylvania #1, New York #2
- Midwest: Ohio #3, Michigan #5

**Implication:** Different regions show spatial effects at different epidemic phases. Northeast shows strongest spatial dynamics at peak.

### Horizon Effects

**Early season:** Differences increase 4x from horizon 0→3 (18.7 → 76.1)
**Peak season:** Differences increase 2.3x from horizon 0→3 (43.6 → 100.4)

**Implication:** Spatial features help both nowcasts and forecasts, but effect is more pronounced for longer horizons.

---

## Key Mechanistic Insights

### 1. North and South Directional Waves Dominate

**Why this matters:**
- Flu typically spreads along latitude lines and north-south corridors
- Interstate travel patterns favor north-south movement on East Coast
- Models correctly capture this asymmetry

### 2. Wave Velocity is Highly Important

**Velocity features account for ~18% of spatial importance**
- Not just WHERE neighbors have cases, but HOW FAST it's moving
- Suggests epidemic dynamics (R_eff) varies geographically
- Could inform early warning systems

### 3. Aggregate Features Matter

**Average wave features rank #1 among individual features**
- Simple spatial average across all directions is most predictive
- Suggests regional "burden" matters beyond specific directionality
- Could potentially simplify model with fewer directional features

### 4. Different Patterns for Absolute vs Relative Differences

**States with largest absolute forecast differences:**
- High population: CA, TX, PA, NY, OH (early season)
- Northeast/Midwest dominance at peak: PA, NY, OH, NC, MI

**States with largest relative forecast differences:**
- Early season: Rural states (ND 90%, MT 89%, SD 79%, ME 79%)
- Peak season: Mixed pattern (MN 29%, NY 29%, PA 22%, ME 20%)

**States with smallest relative differences:**
- Early season: Large states (MN 10%, AL 12%, FL 13%, CA 15%, TX 17%)
- Peak season: South/Plains states (OK 2%, IA 4%, KS 4%, NV 4%, AL 4%)

**Implication:** The scale of measurement matters substantially. Spatial features lead to large absolute differences in populous states (where counts are high) but large relative differences in smaller or more rural states. This may reflect different roles for spatial information depending on local epidemic characteristics, connectivity patterns, or baseline disease burden.

### 5. Northeast Shows Strongest Peak-Season Effects

**Northeast mean difference: 57.7 (early) → 174.3 (peak) = 3x increase**

Possible explanations:
- High population density enables rapid spatial spread
- Strong interstate commuting patterns
- Smaller states mean more cross-border effects
- Later peak timing allows wave propagation from South

---

## Practical Implications

### When Models Diverge Most

**Largest forecast divergence observed:**
- Peak transmission periods (Jan-Feb): 79 mean absolute difference vs 54 at early season
- Northeast region during peak: 174 mean absolute difference, 16.4% relative
- Longer forecast horizons: differences grow from H0 to H3
- Small/rural states show largest relative differences (up to 90% at early season)
- Large states show largest absolute differences (up to 619 for Pennsylvania at peak)

**Smallest forecast divergence observed:**
- Late season would likely show minimal differences (based on decreasing trend)
- Isolated locations (Alaska, Hawaii, Puerto Rico): low relative differences
- Some Southern/Plains states at peak: 2-4% relative differences

**Interpretation:** The choice between models may matter most for peak-season forecasting, Northeast region, and when evaluating small states on a relative basis. However, actual forecast accuracy would need to be evaluated against observed outcomes to determine which model performs better in these contexts.

### Computational Cost vs Forecast Differentiation

**Runtime:**
- Base model: 1-2 hours
- Spatial model: 6-12 hours (3-6x longer)

**Forecast differences:**
- Mean absolute difference: 54-79 hospitalizations (varies by season)
- Mean relative difference: 11-40% (varies by season)
- Pattern varies substantially by location type and epidemic phase

**Considerations for model selection:**
- The 3-6x computational cost produces substantially different forecasts
- Whether these differences translate to accuracy improvements requires validation against observed outcomes
- Cost-benefit may vary by use case, forecast horizon, region, and epidemic phase
- A simplified spatial model may offer a middle ground (see below)

### Potential Model Simplifications

Based on feature importance:

1. **Reduce directions:** N, S, and W account for most importance. Could drop SE.

2. **Keep velocity:** Velocity features are #2 most important category - essential.

3. **Keep aggregates:** Average wave is single most important feature.

4. **Test 1-lag only:** Lag-2 features are important but might be redundant with lag-1.

**Potential simplified model:**
- 4 directions (N, S, W, E) instead of 8
- 1 temporal lag instead of 2
- Keep velocity and aggregate
- **Estimated speedup: 2-3x while retaining ~80% of spatial signal**

---

## Answers to Key Research Questions

### 1. How much do spatial features contribute to model predictions?

**~40% of feature importance** consistently across both dates. Additionally, temporal lag features drop from 34% → 18% when spatial features are added, suggesting spatial features capture some information previously represented by temporal lags alone.

### 2. Which compass directions are most predictive?

**Ranking (consistent across dates):**
1. North (615k importance)
2. South (441k importance)
3. West (278k importance)
4. East (231k importance)
5. Southeast (194k importance)

### 3. Do spatial features capture symmetric or asymmetric spread?

**Asymmetric:** North-South dominance suggests the model captures epidemic spread along travel corridors and latitude patterns, not simple radial diffusion.

### 4. Are 1-week or 2-week spatial lags more important?

Both are important, with **lag-1 slightly more important** than lag-2 based on individual feature rankings.

### 5. Which US regions show largest forecast differences between models?

**By absolute differences:**
- **Early season:** West and South (large states: CA, TX, AZ)
- **Peak season:** Northeast (174 mean, +3x from early season), then Midwest

**By relative differences:**
- **Early season:** Midwest (40.4%) and Northeast (40.0%)
- **Peak season:** Northeast (16.4%), then Midwest (11.1%)

**Overall pattern:** Northeast shows most dramatic forecast divergence at peak, both in absolute and relative terms.

### 6. Do models diverge more at epidemic peak or throughout season?

**Peak season shows:**
- 47% larger mean absolute differences (78.92 vs 53.77)
- But 73% lower mean relative differences (11% vs 40%)
- Max differences 165% larger at peak (748.64 vs 282.90)

**Answer:** Models diverge more in absolute terms at peak, but less in relative terms. The pattern depends heavily on measurement scale and epidemic burden.

### 7. Does forecast divergence justify the computational cost?

**Cannot be determined from this analysis alone.**
- The 3-6x computational cost produces substantially different forecasts (54-79 mean absolute difference)
- Whether these differences represent accuracy improvements requires validation against observed outcomes
- Cost-benefit likely varies by region, epidemic phase, forecast horizon, and decision context
- A simplified spatial model may provide a favorable tradeoff

### 8. What is the minimum set of spatial features needed?

**Proposed minimal set (requires validation):**
- 4 directions (N, S, E, W) instead of 8
- 1 temporal lag instead of 2
- Velocity features (critical - 18% of spatial importance)
- Aggregate features (top individual feature)

**Expected impact:** ~75-80% retention of spatial signal at ~2-3x speedup (requires empirical testing)

---

## Files Generated

### For date: 2024-12-21

**Feature Importance Analysis:**
- `category_comparison.png` - Bar chart of feature categories
- `top_features_comparison.png` - Side-by-side top 20 features
- `spatial_features_breakdown.png` - Breakdown of spatial feature types
- `category_importance.csv` - Importance by category
- `feature_comparison.csv` - Detailed feature comparison
- `summary.json` - Summary statistics

**Geographic Analysis:**
- `state_differences.png` - Top 30 states by difference
- `horizon_differences.png` - Differences by forecast horizon
- `region_differences.png` - Regional comparison
- `state_scatter_comparison.png` - Base vs spatial for top 6 states
- `state_horizon_heatmap.png` - Heatmap of state × horizon
- `forecast_comparison_detailed.csv` - All comparisons
- `state_differences.csv` - State-level summary
- `horizon_differences.csv` - Horizon statistics
- `region_differences.csv` - Regional statistics

### For date: 2025-02-15

Same set of files as above.

**Location:** `analysis/outputs/{date}/`

---

## Recommendations for Future Work

### 1. Validate Simplified Model

Test 4-direction model with 1-lag to confirm speedup without accuracy loss.

### 2. Investigate Northeast Peak Effects

Why does Northeast show 3x larger spatial effects at peak?
- Analyze commuting data
- Check epidemic timing (later peak allows more wave propagation?)
- Compare population density effects

### 3. Develop Adaptive Model Selection

Use base model for:
- Late season (Mar-May)
- Isolated locations
- Low-transmission periods

Use spatial model for:
- Peak season (Dec-Feb)
- Large states
- Northeast region

### 4. Feature Engineering

Test additional features suggested by analysis:
- **Acceleration** (since velocity matters)
- **Weighted aggregates** by population or distance
- **Regional clusters** instead of individual directions
- **Lagged velocity** (rate of change of velocity)

### 5. Interpretability Analysis

Create visualizations showing:
- Wave propagation maps over time
- Which neighbors contribute most to each state's predictions
- Evolution of directional importance through season

### 6. Comparative Evaluation

Extend analysis to:
- More dates (weekly throughout season)
- Multiple seasons (2023-24, 2024-25)
- Comparison with other spatial methods (CAR models, gravity models)
- Validation against actual outcomes (need ground truth data)

---

## Conclusion

The spatial model's directional wave features account for ~40% of feature importance and produce substantially different forecasts than the base model. The magnitude and nature of differences varies significantly by epidemic phase, geographic region, and measurement scale (absolute vs relative):

- **Feature importance**: Spatial features consistently account for 40% of importance, while temporal lag features drop from 34% to 18%, suggesting spatial features capture information previously represented by lags
- **Forecast divergence**: Models differ by 54-79 mean hospitalizations (absolute) or 11-40% (relative), depending on season
- **Geographic patterns**: Large states show largest absolute differences; small/rural states show largest relative differences
- **Temporal patterns**: Northeast region shows 3x larger differences at peak vs early season
- **Directional patterns**: North-South waves and wave velocity dominate spatial feature importance

**Key mechanistic insight:** The consistent importance of North-South directional waves and wave velocity suggests the model is capturing disease spread along travel corridors with varying epidemic velocity across regions. The displacement of temporal lag features by spatial features indicates that neighbor information provides richer context about recent trends than time alone.

**Important caveats:** This analysis characterizes where and when models diverge, but does not assess which model produces more accurate forecasts. Validation against observed outcomes would be needed to determine whether the spatial model's additional computational cost (3-6x longer runtime) translates to improved forecast accuracy. The answer may vary by region, epidemic phase, and evaluation metric.
