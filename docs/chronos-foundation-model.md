# Chronos Time Series Foundation Model

## Overview

This document describes the implementation and use of Amazon's Chronos time series foundation model for flu forecasting in the 2025-2026 Reich Lab forecast hub.

## What are Time Series Foundation Models?

Time series foundation models represent a paradigm shift in forecasting, analogous to how large language models (LLMs) transformed natural language processing. These models are:

- **Pre-trained on massive datasets**: Trained on thousands of diverse time series from different domains
- **Zero-shot capable**: Can forecast new time series without requiring training or fine-tuning
- **Transferable**: Leverage patterns learned from broad time series data to new forecasting tasks
- **Probabilistic**: Generate uncertainty estimates naturally through sampling

## Chronos: Learning the Language of Time Series

### Architecture

Chronos is built on the T5 (Text-to-Text Transfer Transformer) encoder-decoder architecture, originally designed for natural language tasks. The key innovation is treating time series values as "tokens" similar to words in text:

1. **Tokenization**: Continuous time series values are discretized into a finite vocabulary of bins
2. **Sequence Modeling**: The binned values are treated as token sequences
3. **Autoregressive Generation**: The model generates future tokens (forecasts) autoregressively
4. **De-tokenization**: Predicted tokens are converted back to continuous values

### Model Variants

Chronos comes in five sizes, offering different accuracy-speed trade-offs:

| Variant | Parameters | Typical Use Case |
|---------|-----------|------------------|
| chronos-t5-tiny | 8M | Rapid experimentation, embedded systems |
| chronos-t5-mini | 20M | Fast inference with reasonable accuracy |
| chronos-t5-small | 46M | **Default: Best balance** |
| chronos-t5-base | 200M | Higher accuracy applications |
| chronos-t5-large | 710M | Maximum accuracy, research settings |

Our implementation defaults to `small` (46M parameters) for a good balance of accuracy and inference speed.

## Implementation Details

### Data Pipeline

```
NHSN Hospitalization Data (iddata)
    ↓
Weekly aggregated counts by location
    ↓
Last 104 weeks (2 years) used as context
    ↓
Chronos model generates samples
    ↓
Quantiles computed from samples
    ↓
Hubverse format output
```

### Key Characteristics

1. **No Data Transformations**: Unlike traditional models (AR, SARIX, GBQR) that apply power transforms, scaling, and centering, Chronos receives raw hospitalization counts. The model learns appropriate normalization internally.

2. **Location-Level Forecasting**: Forecasts are generated independently for each location. There is no pooling or sharing of information across locations, unlike our pooled AR models.

3. **Sampling-Based Quantiles**: The model generates 200 trajectory samples, from which the 23 required quantile levels are computed empirically.

4. **Context Window**: By default, 104 weeks (2 years) of historical data provide context. This can be adjusted via the `--context_length` parameter.

### Computational Requirements

- **First Run**: Downloads pre-trained weights (~200MB for small model)
- **Inference Time**: 5-15 minutes on CPU for small model, all locations
- **Memory**: ~2-4GB RAM for small model
- **GPU Support**: Can use CUDA for faster inference (requires code modification)

## Comparison with Traditional Models

### Advantages

1. **No Training Required**: Zero-shot capability means no parameter estimation phase
2. **Robust to Data Quality**: Pre-training on diverse data makes it robust to missing values, outliers
3. **Automatic Feature Learning**: No need for manual feature engineering or covariate selection
4. **Captures Complex Patterns**: Transformer architecture can model long-range dependencies
5. **Quick to Deploy**: Can be applied to new forecasting problems rapidly

### Potential Limitations

1. **Domain Specificity**: May not capture flu-specific patterns as well as trained models
2. **Black Box Nature**: Less interpretable than statistical models like AR or SARIX
3. **Computational Overhead**: Requires more compute than simple statistical models
4. **No Covariate Integration**: Cannot easily incorporate external signals (weather, search trends, etc.)
5. **Fixed Architecture**: Cannot customize model structure for problem-specific needs

## Usage Examples

### Basic Usage (Default Settings)

```bash
cd src/chronos
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

python main.py --today_date=2025-01-04
```

### Advanced Usage

```bash
# Use larger model for better accuracy
python main.py --today_date=2025-01-04 --model_size=base

# Increase context window to 3 years
python main.py --today_date=2025-01-04 --context_length=156

# More samples for smoother quantile estimates
python main.py --today_date=2025-01-04 --num_samples=500
```

## Performance Expectations

Based on published Chronos benchmarks and our understanding of flu forecasting:

### Expected Strengths
- Short-term forecasts (horizons 0-1) where recent trends dominate
- Handling data irregularities and reporting artifacts
- Generalizing across different flu seasons
- Producing well-calibrated prediction intervals

### Potential Weaknesses
- Long-range forecasts (horizon 3) where domain knowledge is valuable
- Capturing flu-specific seasonality patterns
- Leveraging location-specific characteristics
- Responding quickly to regime changes (new variants, behavioral shifts)

## Ensemble Opportunities

Chronos can be combined with traditional models to leverage complementary strengths:

1. **Chronos + AR Models**: Foundation model flexibility + statistical model interpretability
2. **Chronos + GBQR**: Zero-shot capability + covariate integration
3. **Multi-Model Ensemble**: Include Chronos in existing ensemble frameworks

## Research and References

### Primary Reference

Ansari, A. F., Stella, L., Turkmen, C., Zhang, X., Mercado, P., Shen, H., ... & Wang, Y. (2024).
**Chronos: Learning the Language of Time Series.**
*arXiv preprint arXiv:2403.07815.*

### Key Findings from Paper

- Outperforms specialized deep learning models on many benchmarks
- Strong zero-shot performance competitive with trained models
- Scales well: larger models generally perform better
- Robust across different time series characteristics (frequency, length, domain)

### Other Time Series Foundation Models

For context, other notable models in this space include:

- **TimesFM** (Google Research): 200M parameter decoder-only model
- **Lag-Llama**: Built on LLaMA architecture, probabilistic forecasting focus
- **Moirai** (Salesforce): Universal time series transformer
- **TimeGPT** (Nixtla): Commercial API-based model

Chronos was chosen for its:
- Strong open-source support and documentation
- Active development and community
- Balance of performance and accessibility
- Clear benchmarking results

## Running on Unity HPC Cluster

The Chronos model includes SLURM scripts for efficient execution on the UMass Unity cluster.

### Setup on Unity

```bash
# SSH to Unity
ssh your_username@unity.rc.umass.edu

# Navigate to chronos directory
cd ~/forecast-sandbox-2025-2026/src/chronos

# Load Python and create environment
module load python/3.11
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Create logs directory
mkdir -p logs
```

### Running Jobs

**Single Forecast**:
```bash
# Default settings (today, small model)
sbatch run_chronos.slurm

# Specific configuration
sbatch --export=TODAY_DATE=2025-01-04,MODEL_SIZE=base run_chronos.slurm
```

**Multiple Dates (Array Job)**:
```bash
# Edit DATES array in run_chronos_array.slurm, then:
sbatch run_chronos_array.slurm
```

### Resource Requirements by Model Size

| Model | CPUs | Memory | Time Limit | Typical Runtime |
|-------|------|--------|------------|-----------------|
| tiny  | 4    | 8GB    | 10 min     | ~2 min          |
| mini  | 8    | 10GB   | 10 min     | ~3 min          |
| small | 8    | 12GB   | 15 min     | ~5 min          |
| base  | 8    | 16GB   | 20 min     | ~8 min          |
| large | 16   | 32GB   | 45 min     | ~25 min         |

### Measured Performance

Based on testing with Chronos-T5-Base (200M parameters) on Unity CPU nodes:
- Average runtime: **8:10 minutes** (std: 5 seconds)
- Test dates: 2024-12-28, 2025-01-04, 2025-01-11
- Configuration: 8 CPUs, 16GB RAM
- Output: 6,095 forecast rows per run

### GPU Acceleration (Optional)

For faster inference, GPU support is available:

```bash
# Modify SLURM script:
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
```

Also update `main.py` line 89 to use CUDA:
```python
device_map="cuda",  # Change from "cpu"
```

Expected GPU speedup: 3-5x faster than CPU.

### Monitoring and Output

```bash
# Check job status
squeue -u $USER

# View logs
tail -f logs/chronos_<job_id>.out

# Output location
ls ../../model-output/UMass-chronos_<model_size>/
```

See `SLURM_README.md` in the `src/chronos/` directory for complete documentation, troubleshooting, and advanced usage.

## Future Directions

Potential enhancements to explore:

1. **Fine-Tuning**: Adapt Chronos to flu-specific patterns (requires additional development)
2. **Hybrid Models**: Combine Chronos with traditional models in novel ways
3. **Multi-Horizon Optimization**: Tune context length and sampling for each horizon separately
4. **Ensemble Weighting**: Determine optimal weight for Chronos in ensemble combinations
5. **GPU Acceleration**: Deploy on GPU for real-time forecasting workflows

## Technical Support

For questions or issues:

- **Model Implementation**: Thomas Robacker (trobacker@umass.edu)
- **Chronos Library**: https://github.com/amazon-science/chronos-forecasting
- **Hub Structure**: Reich Lab team

## Appendix: Command-Line Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--today_date` | String | Today | Model run date (YYYY-MM-DD) |
| `--model_size` | Choice | small | Model variant (tiny/mini/small/base/large) |
| `--context_length` | Integer | 104 | Weeks of historical data (2 years default) |
| `--num_samples` | Integer | 200 | Samples for quantile estimation |

## Version History

- **v1.0** (2025-01-04): Initial implementation with chronos-t5-small as default
  - Zero-shot forecasting for weekly flu hospitalizations
  - Support for all model sizes
  - Configurable context window and sampling
  - Integration with iddata pipeline
