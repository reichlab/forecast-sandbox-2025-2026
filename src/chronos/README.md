# Chronos Time Series Foundation Model for Influenza Forecasting

This model uses Amazon's Chronos, a pre-trained time series foundation model, to generate probabilistic forecasts of weekly incident flu hospitalizations.

## Model Overview

Chronos is a family of pre-trained time series forecasting models based on language model architectures (T5). The models treat time series as sequences of tokens and leverage patterns learned from training on a large corpus of diverse time series data. This enables zero-shot forecasting without requiring model training on the specific flu hospitalization data.

Key features:
- **Zero-shot forecasting**: No training or fine-tuning required
- **Probabilistic predictions**: Generates quantile forecasts through sampling
- **Foundation model approach**: Leverages pre-trained knowledge from diverse time series

## Model Sizes Available

The model supports multiple sizes with different trade-offs between accuracy and speed:
- `tiny`: 8M parameters - fastest, lowest accuracy
- `mini`: 20M parameters
- `small`: 46M parameters (default) - good balance
- `base`: 200M parameters
- `large`: 710M parameters - highest accuracy, slowest

## Configuration Options

- `--today_date`: Effective model run date (YYYY-MM-DD). Defaults to today.
- `--model_size`: Which Chronos variant to use (tiny/mini/small/base/large). Default: small
- `--context_length`: Number of weeks of historical data to use. Default: 104 (2 years)
- `--num_samples`: Number of samples for quantile estimation. Default: 200

## Running Locally

To test this model locally, run the following from this directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

# Run with default settings (small model)
python main.py --today_date=2025-01-04

# Run with a larger model for better accuracy
python main.py --today_date=2025-01-04 --model_size=base

# Run with more historical context
python main.py --today_date=2025-01-04 --context_length=156  # 3 years
```

The first run will download the pre-trained model (~200MB for small), which is cached for subsequent runs.

## Model Methodology

1. **Data Preparation**: Historical weekly flu hospitalization counts from NHSN are fetched using the iddata library
2. **Context Window**: The most recent `context_length` weeks of data are provided to the model
3. **Forecasting**: Chronos generates `num_samples` trajectory samples for 4 weeks ahead (horizons -1, 0, 1, 2, 3)
4. **Quantile Estimation**: The 23 required quantile levels are computed from the forecast samples
5. **Output**: Forecasts are saved in Hubverse format to `model-output/UMass-chronos_{size}/`

## Data Sources

- **NHSN**: Weekly incident flu hospitalizations from the National Healthcare Safety Network

## Computational Requirements

- **CPU**: Runs on CPU by default. Typical runtime: 5-15 minutes for small model
- **GPU**: Can use CUDA-enabled GPU for faster inference (modify `device_map` in main.py)
- **Memory**: ~2-4GB RAM for small model, more for larger variants
- **Disk**: ~200MB for model weights (downloaded on first run)

## References

- Ansari, A. F., Stella, L., Turkmen, C., Zhang, X., Mercado, P., Shen, H., ... & Wang, Y. (2024). Chronos: Learning the Language of Time Series. arXiv preprint arXiv:2403.07815.
- GitHub: https://github.com/amazon-science/chronos-forecasting

## Author

Thomas Robacker (trobacker@umass.edu)
