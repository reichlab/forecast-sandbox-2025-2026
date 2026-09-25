"""
Score peak week (pmf) and peak size (quantile) forecasts against analysis/peak-models/peak-oracle.csv.

Scores, for each (model, reference_date, location):
  peak week:  log score log p(W*) (natural log; higher is better), ranked probability score (RPS, lower is better),
              and the probability placed within +/- 1 week of W*.
  peak size:  weighted interval score (WIS) of the quantiles on the count scale and on the log(count + 1) scale,
              and 50% / 95% interval coverage.
Locations with tied peaks are excluded from peak-week scores.

Usage (from the repository root):
    python analysis/peak-models/score_peak_forecasts.py [--hub_root ../FluSight-forecast-hub] [--models ...]
Writes analysis/peak-models/peak-scores.csv (one row per forecast) and prints summaries.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
WEEK_TARGET = "peak week inc flu hosp"
SIZE_TARGET = "peak inc flu hosp"


def season_of_ref(ref: str) -> str:
    y = int(ref[:4]) if int(ref[5:7]) >= 7 else int(ref[:4]) - 1
    return f"{y}/{str(y + 1)[2:]}"


def wis(q_levels: np.ndarray, q_values: np.ndarray, y: float) -> float:
    """WIS computed as the mean pinball loss over the quantile levels, times 2 (equivalent to Bracher et al. 2021
    when the levels are symmetric and include the median)."""
    return 2 * np.mean([(float(y <= v) - tau) * (v - y) for tau, v in zip(q_levels, q_values)])


def load_forecasts(model_dirs: list[Path]) -> pd.DataFrame:
    frames = []
    for d in model_dirs:
        for f in sorted(d.glob("*.csv")):
            df = pd.read_csv(f, dtype={"location": str, "output_type_id": str}, low_memory=False)
            df = df.loc[df["target"].isin([WEEK_TARGET, SIZE_TARGET])]
            if len(df):
                df["model"] = d.name
                frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def score(forecasts: pd.DataFrame, oracle: pd.DataFrame, windows: dict[str, list[str]]) -> pd.DataFrame:
    forecasts = forecasts.assign(season=forecasts["reference_date"].map(season_of_ref))
    rows = []
    for (model, ref, loc, season), g in forecasts.groupby(["model", "reference_date", "location", "season"]):
        truth = oracle.loc[(oracle["season"] == season) & (oracle["location"] == loc)]
        if len(truth) == 0:
            continue
        truth = truth.iloc[0]
        rec = {"model": model, "reference_date": ref, "location": loc, "season": season,
               "weeks_from_peak": (pd.Timestamp(ref) - pd.Timestamp(truth["peak_week"])).days // 7}

        pmf = g.loc[g["target"] == WEEK_TARGET].set_index("output_type_id")["value"]
        if len(pmf) and not truth["tied"]:
            dates = windows[season]
            p = pmf.reindex(dates).fillna(0.0).to_numpy()
            p = p / p.sum()
            i_true = dates.index(truth["peak_week"])
            rec["log_score"] = np.log(max(p[i_true], 1e-300))
            cdf = np.cumsum(p)
            obs_cdf = (np.arange(len(dates)) >= i_true).astype(float)
            rec["rps"] = np.sum((cdf - obs_cdf) ** 2)
            rec["prob_pm1"] = p[max(i_true - 1, 0):i_true + 2].sum()
            rec["pmf_dates_submitted"] = int(pmf.index.isin(dates).sum())

        q = g.loc[(g["target"] == SIZE_TARGET) & (g["output_type"] == "quantile")]
        if len(q):
            q = q.assign(tau=q["output_type_id"].astype(float)).sort_values("tau")
            taus, vals, y = q["tau"].to_numpy(), q["value"].to_numpy(dtype=float), float(truth["peak_count"])
            rec["wis"] = wis(taus, vals, y)
            rec["wis_log"] = wis(taus, np.log1p(np.maximum(vals, 0)), np.log1p(y))
            rec["rel_error_median"] = np.interp(0.5, taus, vals) / max(y, 1) - 1
            for lo, hi, name in [(0.25, 0.75, "cov50"), (0.025, 0.975, "cov95")]:
                if lo in taus and hi in taus:
                    rec[name] = float(np.interp(lo, taus, vals) <= y <= np.interp(hi, taus, vals))
        rows.append(rec)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--models", nargs="*", default=None,
                        help="model directory names under model-output/ (default: all UMass-peak_*)")
    parser.add_argument("--hub_root", default=None,
                        help="optional path to a FluSight-forecast-hub clone whose peak forecasts to also score")
    args = parser.parse_args()

    oracle = pd.read_csv(HERE / "peak-oracle.csv", dtype={"location": str})
    tasks = json.loads((ROOT / "hub-config" / "tasks.json").read_text())
    windows = {}
    for mt in tasks["rounds"][0]["model_tasks"]:
        if mt["task_ids"]["target"]["required"] == [WEEK_TARGET]:
            dates = mt["output_type"]["pmf"]["output_type_id"]["required"]
            windows[season_of_ref(dates[0])] = dates

    names = args.models or sorted(p.name for p in (ROOT / "model-output").glob("UMass-peak_*"))
    dirs = [ROOT / "model-output" / n for n in names]
    if args.hub_root:
        dirs += [d for d in sorted((Path(args.hub_root) / "model-output").iterdir()) if d.is_dir()]
    forecasts = load_forecasts(dirs)
    scores = score(forecasts, oracle, windows)
    scores.to_csv(HERE / "peak-scores.csv", index=False)

    ours = scores.loc[scores["model"].str.startswith("UMass-peak_")]
    cols = ["log_score", "rps", "prob_pm1", "wis", "wis_log", "cov50", "cov95"]
    print(ours.groupby(["season", "model"])[cols].mean().round(3).to_string())
    print(ours.groupby("model")[cols].mean().round(3).to_string())


if __name__ == "__main__":
    main()
