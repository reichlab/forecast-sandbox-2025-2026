# peak_gbqr_sb

GBQR direct peak model (`idmodels.peak.PeakGBQRModel`) with separate feature sets for its two parts:
peak size uses **SB (synchrony + burden)** features, peak timing uses **core** features
(`size_feature_groups`, `timing_feature_groups`; see `analysis/peak-models/peak-models-methods.qmd`, "Versions and names").
Otherwise identical to `peak_gbqr` (training data, NHSN revision simulation, output format).

Chosen from the ILINet development runs and the NHSN 2023/24–2024/25 test runs; the hindcasts in
`model-output/UMass-peak_gbqr_sb/` for 2023/24–2025/26 were produced with `analysis/peak-models/nhsn_validation.py`
(validation name `gbqr__N1`).

Run locally as for `peak_gbqr`: `python main.py --today_date=2026-01-07 --short_run`.
