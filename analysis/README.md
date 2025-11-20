# Model Analysis: gbqr_3src vs gbqr_3src_spatial

This directory contains analysis tools to understand why the spatial model (`gbqr_3src_spatial`) outperforms the base model (`gbqr_3src`).

## Quick Start

### 1. Run Models with Diagnostics (on Unity cluster)

```bash
# Base model
cd ../src/gbqr_3src
sbatch submit-diagnostic-unity.sh

# Spatial model
cd ../src/gbqr_3src_spatial
sbatch submit-diagnostic-unity.sh
```

Wait for jobs to complete (~6-12 hours for spatial model).

### 2. Run Analysis (locally or on cluster)

```bash
cd analysis
./run_full_analysis.sh
```

This will analyze both dates (2024-12-21 and 2025-02-15) and generate all outputs.

**Or run individually:**

```bash
# Feature importance analysis
python analyze_feature_importance.py --date 2025-02-15

# Geographic analysis
python analyze_geographic_differences.py --date 2025-02-15
```

## Files

### Analysis Scripts
- `analyze_feature_importance.py` - Compare feature importance between models
- `analyze_geographic_differences.py` - Analyze geographic patterns in forecast differences
- `run_full_analysis.sh` - Run all analyses for both dates

### Documentation
- `ANALYSIS_PLAN.md` - **Detailed analysis plan and execution guide**
- `README.md` - This file (quick reference)

### Model Files
- `../src/gbqr_3src/main_diagnostic.py` - Base model with artifact saving
- `../src/gbqr_3src_spatial/main_diagnostic.py` - Spatial model with artifact saving
- `../src/gbqr_3src/submit-diagnostic-unity.sh` - Cluster submission script (base)
- `../src/gbqr_3src_spatial/submit-diagnostic-unity.sh` - Cluster submission script (spatial)

## Output Structure

After running analyses, outputs will be in `outputs/{date}/`:

```
outputs/
├── 2024-12-21/
│   ├── Feature Importance Analysis:
│   │   ├── category_importance.csv
│   │   ├── feature_comparison.csv
│   │   ├── category_comparison.png
│   │   ├── top_features_comparison.png
│   │   ├── spatial_features_breakdown.png
│   │   └── summary.json
│   │
│   └── Geographic Analysis:
│       ├── forecast_comparison_detailed.csv
│       ├── state_differences.csv
│       ├── horizon_differences.csv
│       ├── region_differences.csv
│       ├── state_differences.png
│       ├── horizon_differences.png
│       ├── region_differences.png
│       ├── state_scatter_comparison.png
│       └── state_horizon_heatmap.png
│
└── 2025-02-15/
    └── (same structure)
```

## Key Questions Answered

### Feature Importance Analysis
1. What % of importance comes from spatial features?
2. Which compass directions are most predictive?
3. Are wave velocity features important?
4. How does feature importance change from early to peak season?

### Geographic Analysis
1. Which states benefit most from spatial features?
2. Are there regional patterns?
3. Do improvements vary by forecast horizon?
4. Are high-transmission states benefiting more?

## Analysis Dates

**2024-12-21** (Early Season)
- Mean forecast difference: 43.67 hospitalizations
- Geographic variance: 1,731
- Captures early epidemic spread

**2025-02-15** (Peak Season)
- Mean forecast difference: 67.77 hospitalizations
- Geographic variance: 9,287
- Maximum spatial feature impact

These dates were selected by analyzing all 55 available forecast dates (see `../compare_model_forecasts.py`).

## Dependencies

```bash
pip install pandas numpy matplotlib seaborn click
```

## Troubleshooting

**"No feature importance files found"**
- Ensure diagnostic models ran successfully
- Check artifact directories exist:
  - `../src/gbqr_3src/artifacts/{date}/`
  - `../src/gbqr_3src_spatial/artifacts/{date}/`
- Look for `*importance*.csv` files in artifact directories

**"Could not load forecasts"**
- Ensure forecast CSV files exist in:
  - `../model-output/UMass-gbqr_3src/`
  - `../model-output/UMass-gbqr_3src_spatial/`

**Memory issues**
- Analysis scripts are lightweight, should run on any machine
- Model runs need 32GB RAM (handled by Unity cluster)

## For More Details

See **ANALYSIS_PLAN.md** for:
- Complete execution guide
- Detailed methodology
- Interpretation guidelines
- Troubleshooting
- Expected results
