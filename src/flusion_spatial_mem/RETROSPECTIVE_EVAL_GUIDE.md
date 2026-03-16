# Retrospective Evaluation Guide for flusion_spatial_mem

## Quick Start

The model is fully set up and ready for retrospective evaluation!

### Run Full Evaluation (85 forecasts)

```bash
cd src/flusion_spatial_mem
./run_retrospective_evaluation.sh
```

### Test with Short Run (faster, for development)

```bash
cd src/flusion_spatial_mem
./run_retrospective_evaluation.sh --short-run
```

## What Will Happen

The evaluation script will:

1. **Generate forecasts for 85 Saturdays** (Oct 2023 - May 2025)
   - 2023-2024 flu season: ~33 weeks
   - 2024-2025 flu season: ~32 weeks
   - Summer 2024: ~20 weeks

2. **For each forecast date:**
   - Fetch NHSN hospitalization data available at that time
   - Calculate MEM thresholds using 3-4 complete past seasons
   - Determine current epidemic phase (baseline/low/medium/high/very_high)
   - Apply adaptive ensemble weights based on phase
   - Generate component forecasts (AR6, GBQR models)
   - Create weighted ensemble forecast
   - Save outputs and diagnostics

3. **Outputs saved to:**
   - `../../model-output/UMass-flusion_spatial_mem/YYYY-MM-DD-UMass-flusion_spatial_mem.csv`
   - `../../model-output/UMass-flusion_spatial_mem/YYYY-MM-DD-diagnostic-info.csv`

## MEM Behavior During Retrospective Evaluation

**Key Feature: Truly Retrospective - No Look-Ahead Bias**

- MEM thresholds are recalculated for EACH forecast date
- Only uses data available UP TO that forecast date
- Example timeline:
  - **2023-10-07**: Uses seasons 2020, 2021, 2022 (3 complete seasons)
  - **2024-01-06**: Uses seasons 2020, 2021, 2022, 2023 (4 complete seasons)
  - **2024-10-05**: Uses seasons 2020, 2021, 2022, 2023 (4 complete seasons)
  - **2025-05-17**: Uses seasons 2020, 2021, 2022, 2023, 2024 (can use parts of 2024)

**Automatic Adaptation:**
- System tries 5 seasons first, falls back to 4 if insufficient data
- Early forecasts may use 3 seasons (minimum for MEM)
- Falls back to equal weights (50/50) if MEM fails

**No Pre-Training Required:**
- MEM doesn't need "warming up"
- Each forecast is independent
- Thresholds adapt naturally over time as more data becomes available

## Expected Runtime

**Full Evaluation (85 forecasts):**
- AR6 model: ~20 seconds per forecast
- GBQR models: ~2-3 minutes per forecast
- MEM calculation: ~10-30 seconds per forecast
- **Total: ~3-4 hours for all 85 forecasts**

**Short Run (85 forecasts):**
- AR6: ~10 seconds per forecast
- GBQR: ~1-2 minutes per forecast
- MEM: ~10-30 seconds per forecast
- **Total: ~2-2.5 hours**

## Verification Steps

After completion, check:

```bash
# Count output files (should be 85)
ls -1 ../../model-output/UMass-flusion_spatial_mem/*.csv | grep -v diagnostic | wc -l

# Count diagnostic files (should be 85)
ls -1 ../../model-output/UMass-flusion_spatial_mem/*diagnostic*.csv | wc -l

# View sample diagnostic info to see MEM weights evolution
head -5 ../../model-output/UMass-flusion_spatial_mem/*-diagnostic-info.csv
```

## Analysis After Evaluation

Once forecasts are complete, you can:

1. **Compare with flusion_spatial2_prod:**
   - Both models use same components (AR6, GBQR)
   - flusion_spatial2_prod uses simple median
   - flusion_spatial_mem uses MEM-adaptive weights

2. **Analyze weight evolution:**
   - Extract diagnostic files showing phase and weights over time
   - See how weights changed during epidemic peaks
   - Correlate phases with actual outbreak intensity

3. **Performance metrics:**
   - WIS (Weighted Interval Score)
   - Coverage at different quantile levels
   - Stratify by epidemic phase for targeted analysis

## Troubleshooting

**If forecasts fail:**
1. Check Python environment is activated: `source .venv/bin/activate`
2. Verify R packages installed: `Rscript -e "renv::status()"`
3. Check logs in terminal output for specific errors
4. Try single forecast first: `python main.py --today_date=2023-10-07`

**If MEM shows as unavailable:**
- Check diagnostic files for `mem_available = FALSE`
- This is expected for very early dates with <3 complete seasons
- Model will use equal weights (50/50) as fallback

## Next Steps

After retrospective evaluation:
1. Compare scores against flusion_spatial2_prod
2. Analyze if adaptive weighting improves performance
3. Consider parameter tuning (weight schedules, lookback windows)
4. Decide if model is competitive for production use
