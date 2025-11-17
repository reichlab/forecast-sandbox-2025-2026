# GBQR model for influenza

This is just a vanilla gbqr model to use as a baseline. HNSN, FluSurvNet and ILINet data, standard flusion options. 

# To run locally without Docker

To test this out locally, run the following with this directory, e.g. `gbqr` as your working directory.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --today_date=2024-01-06 --short_run
```

This should result in a model output file and a pdf with a plot under the directory specified in `main.py`.

# requirements.txt 

`requirements.txt` was generated according to [README.md](../README.md). 
