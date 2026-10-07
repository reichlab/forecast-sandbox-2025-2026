# gbqr_3src

GBQR flu model with standard flusion options, fit jointly to all 53 locations on NHSN,
FluSurv-NET, and ILINet data (fourth-root transform). Same configuration as
`flusion_spatial2_prod/2_gbqr_3src.py`. It supplies the national (US) GBQR forecast for the
`gbqr_3src_spatial` x AR ensembles in `../gbqr_ar_ensembles`, since `gbqr_3src_spatial` is
state-level only.

```bash
python main.py --today_date=2024-01-06 --short_run   # quick local test
sbatch submit-unity-parallel.sh                      # all hub reference dates on Unity
```

Set up with `python -m venv .venv && .venv/bin/pip install -r requirements.txt`
(idmodels v2.1.0, the same pin as `gbqr_3src_spatial` and `flusion_spatial2_prod`).
