# Analysis Plan: Understanding gbqr_3src_spatial Model Performance

## Overview

This analysis plan investigates why `gbqr_3src_spatial` outperforms `gbqr_3src` by examining:
1. Feature importance differences between models
2. Geographic patterns in forecast improvements
3. When spatial features matter most (by season and horizon)

## Key Differences Between Models

**gbqr_3src (Base Model):**
- 3 data sources: FluSurvNet, NHSN, ILINet
- Level features + temporal lags
- Pooled across locations

**gbqr_3src_spatial:**
- Same as base PLUS directional wave features:
  - 8 compass directions (N, NE, E, SE, S, SW, W, NW)
  - Temporal lags (1-2 weeks)
  - Maximum distance: 1500km
  - Includes wave velocity and aggregate features

## Selected Dates for Analysis

Based on forecast comparison analysis of 55 dates from 2023-2025:

### Date 1: **2024-12-21** (Early Season)
- Mean absolute difference: **43.67** hospitalizations
- Geographic variance: **1,731**
- Represents early epidemic spread phase

### Date 2: **2025-02-15** (Peak Season)
- Mean absolute difference: **67.77** hospitalizations
- Geographic variance: **9,287** (highest observed)
- Represents maximum spatial feature impact

**Why these dates?**
- Peak differences occur during active transmission (Jan-Feb 2025)
- Late season shows minimal differences (spatial features less important)
- These two dates capture early spread vs. peak transmission dynamics

---

## Step-by-Step Execution Plan

### Phase 1: Run Diagnostic Models on Unity Cluster

Both models need to be run with artifact saving enabled to capture feature importance.

#### 1.1 Prepare Virtual Environments

For each model directory:

```bash
# gbqr_3src
cd src/gbqr_3src
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# gbqr_3src_spatial
cd ../gbqr_3src_spatial
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

#### 1.2 Create Logs Directories

```bash
cd src/gbqr_3src
mkdir -p logs

cd ../gbqr_3src_spatial
mkdir -p logs
```

#### 1.3 Submit Jobs to Unity

**Base Model:**
```bash
cd src/gbqr_3src
sbatch submit-diagnostic-unity.sh
```

Expected runtime: ~1-2 hours (2 dates × 1 hour each)

**Spatial Model:**
```bash
cd src/gbqr_3src_spatial
sbatch submit-diagnostic-unity.sh
```

Expected runtime: ~6-12 hours (spatial features are computationally expensive)

#### 1.4 Monitor Job Progress

```bash
# Check job status
squeue -u $USER

# Monitor logs (replace JOBID with actual job ID)
tail -f logs/diagnostic-slurm-JOBID_0.out  # For 2024-12-21
tail -f logs/diagnostic-slurm-JOBID_1.out  # For 2025-02-15
```

#### 1.5 Verify Outputs

After jobs complete, verify artifacts were created:

```bash
# Check base model artifacts
ls -lh src/gbqr_3src/artifacts/2024-12-21/
ls -lh src/gbqr_3src/artifacts/2025-02-15/

# Check spatial model artifacts
ls -lh src/gbqr_3src_spatial/artifacts/2024-12-21/
ls -lh src/gbqr_3src_spatial/artifacts/2025-02-15/
```

Expected artifacts:
- Feature importance files (CSV or parquet)
- Model objects (pickled models)
- Potentially other diagnostic outputs

---

### Phase 2: Feature Importance Analysis

Analyzes which features drive predictions and how spatial features contribute.

#### 2.1 Setup Analysis Environment

```bash
cd analysis
python -m venv .venv
source .venv/bin/activate
pip install pandas numpy matplotlib seaborn click
```

#### 2.2 Run Feature Importance Analysis

For each date:

```bash
# Early season date
python analyze_feature_importance.py --date 2024-12-21

# Peak season date
python analyze_feature_importance.py --date 2025-02-15
```

#### 2.3 Outputs Generated

For each date, creates `analysis/outputs/{date}/`:
- `category_importance.csv` - Importance by feature category
- `feature_comparison.csv` - Detailed feature-level comparison
- `category_comparison.png` - Bar chart of category importance
- `top_features_comparison.png` - Side-by-side top 20 features
- `spatial_features_breakdown.png` - Breakdown of spatial feature types
- `summary.json` - Summary statistics

#### 2.4 Key Questions to Answer

1. **What percentage of importance do spatial features account for?**
   - Check `summary.json` → `spatial_feature_importance_pct`

2. **Which spatial directions are most important?**
   - Review `spatial_features_breakdown.png`
   - Look for directional patterns (N, S, E, W, etc.)

3. **Do spatial features matter more at peak vs. early season?**
   - Compare importance percentages between dates

4. **Which non-spatial features are still critical?**
   - Check `top_features_comparison.png` for overlap

---

### Phase 3: Geographic Analysis

Identifies which states benefit most from spatial features.

#### 3.1 Run Geographic Analysis

```bash
cd analysis

# Early season
python analyze_geographic_differences.py --date 2024-12-21

# Peak season
python analyze_geographic_differences.py --date 2025-02-15
```

#### 3.2 Outputs Generated

For each date, creates `analysis/outputs/{date}/`:
- `forecast_comparison_detailed.csv` - All location/horizon comparisons
- `state_differences.csv` - Mean difference by state
- `horizon_differences.csv` - Statistics by forecast horizon
- `region_differences.csv` - Statistics by US region
- `state_differences.png` - Top 30 states bar chart
- `horizon_differences.png` - Differences by horizon
- `region_differences.png` - Regional comparison
- `state_scatter_comparison.png` - Base vs. spatial for top 6 states
- `state_horizon_heatmap.png` - Heatmap of state × horizon differences

#### 3.3 Key Questions to Answer

1. **Which states show the largest improvements?**
   - Check `state_differences.csv` (top rows)
   - Review `state_differences.png`

2. **Are there regional patterns?**
   - Review `region_differences.png`
   - Check if certain regions benefit more (e.g., South vs. Northeast)

3. **Do improvements vary by forecast horizon?**
   - Check `horizon_differences.png`
   - See if spatial features help more for near-term vs. long-term forecasts

4. **Are high-transmission states benefiting most?**
   - Cross-reference top states with hospitalization levels
   - Check `state_scatter_comparison.png` for magnitude patterns

5. **Geographic spread patterns:**
   - Look for contiguous state clusters in top performers
   - Suggests spatial wave features capturing real geographic spread

---

### Phase 4: Synthesis and Interpretation

#### 4.1 Compare Across Dates

Create summary comparison:

```bash
cd analysis

# Create combined comparison (manual or script)
# Compare metrics between 2024-12-21 and 2025-02-15
```

Key comparisons:
- Spatial feature importance: early vs. peak
- Geographic variance: more pronounced at peak?
- Regional patterns: do they shift with epidemic phase?

#### 4.2 Hypothesis Testing

Based on results, evaluate hypotheses:

**H1: Spatial features matter most during active transmission**
- Evidence: Higher differences at peak season
- Compare mean absolute differences between dates

**H2: Border states benefit more from spatial features**
- Evidence: Check if states with more neighbors show larger improvements
- Review state difference rankings

**H3: Directional patterns reflect real disease spread**
- Evidence: Check if southward waves prominent (flu typically spreads south)
- Review directional feature importance in peak season

**H4: Spatial features improve nowcasts more than forecasts**
- Evidence: Compare differences by horizon (h=-1 vs. h=3)
- Review horizon statistics

#### 4.3 Create Summary Report

Document findings covering:

1. **Executive Summary**
   - Main finding: spatial features improve forecasts by X% during peak season
   - Key mechanism: directional wave features capture geographic spread

2. **Feature Importance Findings**
   - % importance from spatial features
   - Which directions/lags matter most
   - Temporal patterns (early vs. peak)

3. **Geographic Findings**
   - Top states benefiting from spatial features
   - Regional patterns
   - Relationship to epidemic intensity

4. **Practical Implications**
   - When to use spatial vs. base model
   - Computational cost vs. accuracy trade-off
   - Recommendations for future model development

---

## Expected Computational Requirements

### Unity Cluster Resources

**Base Model (gbqr_3src):**
- Memory: 32GB per job
- Cores: 8
- Time: 1-2 hours per date
- Total: ~2-4 hours for both dates

**Spatial Model (gbqr_3src_spatial):**
- Memory: 32GB per job
- Cores: 8
- Time: 3-6 hours per date (computing spatial features is expensive)
- Total: ~6-12 hours for both dates

**Analysis Scripts:**
- Can run locally (not computationally intensive)
- <5 minutes per script per date

### Disk Space

- Artifacts per model per date: ~100MB - 1GB (depends on idmodels implementation)
- Analysis outputs per date: ~5-10MB
- Total: ~2-5GB for complete analysis

---

## Troubleshooting

### Issue: No feature importance files found

**Symptoms:** Analysis script reports "No feature importance files found"

**Solutions:**
1. Check if `save_feat_importance=True` in diagnostic scripts
2. Verify `artifact_store_root` is set correctly
3. Check idmodels implementation - feature importance saving may be model-specific
4. Look for alternative artifact names (gain, importance, feat_imp, etc.)

### Issue: Job timeout on Unity

**Symptoms:** SLURM job killed due to time limit

**Solutions:**
1. Increase time limit in SBATCH directive (currently 2hr base, 6hr spatial)
2. Use `--short_run` flag to reduce quantiles and bags for testing
3. Request more cores to parallelize computation

### Issue: Out of memory errors

**Symptoms:** Job fails with memory allocation error

**Solutions:**
1. Increase `--mem` in SBATCH directive (currently 32GB)
2. Reduce `num_bags` in model config
3. Check if spatial features are memory-intensive (many neighbor calculations)

### Issue: Missing forecast files for selected dates

**Symptoms:** Geographic analysis can't find forecast CSVs

**Solutions:**
1. Check that reference dates match Saturday (forecasts use Saturday reference)
2. Verify model output directories: `model-output/UMass-{model}/`
3. Run models first if forecasts don't exist

---

## Next Steps After Analysis

1. **Validate Findings**
   - Extend to more dates if patterns unclear
   - Test on different seasons/years

2. **Model Improvements**
   - If certain directions dominate, simplify to fewer directions
   - If velocity matters, investigate acceleration features
   - If aggregate important, explore other spatial aggregations

3. **Computational Optimization**
   - Profile where spatial model spends time
   - Optimize spatial feature computation
   - Consider caching neighbor distances

4. **Documentation**
   - Write up findings for team
   - Update model metadata with insights
   - Create visualization dashboard for ongoing monitoring

---

## File Structure Summary

```
forecast-sandbox-2025-2026/
├── src/
│   ├── gbqr_3src/
│   │   ├── main.py                        # Standard run
│   │   ├── main_diagnostic.py             # With artifact saving
│   │   ├── submit-diagnostic-unity.sh     # Cluster submission
│   │   ├── artifacts/                     # Model outputs (created by job)
│   │   │   ├── 2024-12-21/
│   │   │   └── 2025-02-15/
│   │   └── logs/                          # SLURM logs
│   │
│   └── gbqr_3src_spatial/
│       ├── main.py                        # Standard run
│       ├── main_diagnostic.py             # With artifact saving
│       ├── submit-diagnostic-unity.sh     # Cluster submission
│       ├── artifacts/                     # Model outputs (created by job)
│       │   ├── 2024-12-21/
│       │   └── 2025-02-15/
│       └── logs/                          # SLURM logs
│
├── analysis/
│   ├── ANALYSIS_PLAN.md                   # This document
│   ├── analyze_feature_importance.py      # Feature analysis script
│   ├── analyze_geographic_differences.py  # Geographic analysis script
│   └── outputs/                           # Analysis outputs (created by scripts)
│       ├── 2024-12-21/
│       │   ├── *.png                      # Visualizations
│       │   ├── *.csv                      # Data tables
│       │   └── summary.json               # Summary stats
│       └── 2025-02-15/
│           └── (same structure)
│
├── model-output/                          # Existing forecasts (already generated)
│   ├── UMass-gbqr_3src/
│   └── UMass-gbqr_3src_spatial/
│
└── compare_model_forecasts.py             # Initial date selection script
```

---

## Timeline Estimate

**Phase 1: Cluster Jobs** (6-12 hours elapsed time)
- Submit jobs: 10 minutes
- Wait for completion: 6-12 hours (can run overnight)

**Phase 2: Feature Analysis** (30 minutes)
- Setup: 5 minutes
- Run both dates: 10 minutes
- Review outputs: 15 minutes

**Phase 3: Geographic Analysis** (30 minutes)
- Run both dates: 10 minutes
- Review outputs: 20 minutes

**Phase 4: Synthesis** (2-4 hours)
- Cross-date comparison: 1 hour
- Hypothesis testing: 1 hour
- Write summary report: 1-2 hours

**Total Active Time:** ~3-5 hours (excluding cluster wait time)
**Total Elapsed Time:** ~1-2 days (with overnight cluster jobs)

---

## Questions to Address in Final Report

1. What percentage of prediction improvement comes from spatial features?
2. Which compass directions are most predictive?
3. Do spatial features capture symmetric or asymmetric spread patterns?
4. Are 1-week or 2-week spatial lags more important?
5. Which US regions benefit most from spatial modeling?
6. Do spatial features help more at epidemic peak or throughout season?
7. Is the computational cost (3-6x longer runtime) justified by accuracy gains?
8. What is the minimum set of spatial features needed to capture most benefit?

---

## Contact

For questions about this analysis plan:
- Check idmodels documentation for artifact structure
- Review Unity cluster docs for job submission
- Contact Nick Reich or team members for domain expertise
