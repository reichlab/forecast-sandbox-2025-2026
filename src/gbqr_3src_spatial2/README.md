# GBQR model for influenza with spatial features

This is a GBQR model with standard flusion options, HNSN, FluSurvNet, ILINet data, and with the directional wave features added. This differs from `gbqr_3src_spatial` model in that this model only uses the 4 cardinal directions (not 8) and 1 lag (not 2). Trying to see if that reduces complexity and computation of the model while keeping the value of the spatial features.

# To run locally without Docker

To test this out locally, run the following with this directory, e.g. `gbqr_3src_spatial2` as your working directory.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --today_date=2024-01-06 --short_run
```

This should result in a model output file and a pdf with a plot under the directory specified in `main.py`.

# requirements.txt 

`requirements.txt` was generated according to [README.md](../README.md). 
