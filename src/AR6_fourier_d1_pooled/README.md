# AR p=6, d=1 model for influenza

# To run locally without Docker

To test this out locally, run the following with this directory, e.g. `AR6_fourier_d1_pooled` as your working directory.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python main.py --today_date=2024-01-06
```

This should result in a model output file under the directory specified in `main.py`.

# requirements.txt and renv.lock details

`requirements.txt` was generated according to [README.md](../README.md). 