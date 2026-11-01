"""
Score run_eval.sh output against target-data/oracle-output.csv.

Reports quantile-based WIS (from the 7 quantile levels) relative to the
no-weather variant, overall / by horizon / by season, with a bootstrap over
reference dates, plus 50% and 95% interval coverage.

Usage: python score_eval.py OUTPUT_ROOT
"""

import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ORACLE = Path(__file__).parents[2] / "target-data" / "oracle-output.csv"
BASE = "none"


def load(root):
    frames = []
    for f in glob.glob(f"{root}/model-output/UMass-gbqr_3src_gefs_*/*.csv"):
        d = pd.read_csv(f, dtype={"location": str})
        d["variant"] = Path(f).parent.name.replace("UMass-gbqr_3src_gefs_", "")
        frames.append(d)
    d = pd.concat(frames)
    d = d[(d.output_type == "quantile") & (d.target == "wk inc flu hosp")]
    orc = pd.read_csv(ORACLE, dtype={"location": str})
    orc = orc[orc.target == "wk inc flu hosp"].drop_duplicates(["location", "target_end_date"])
    return d.merge(orc[["location", "target_end_date", "oracle_value"]], on=["location", "target_end_date"])


def score(d):
    q = d.output_type_id.astype(float)
    y, x = d.oracle_value, d.value
    d = d.assign(pl=2 * np.where(y < x, (1 - q) * (x - y), q * (y - x)))
    keys = ["variant", "reference_date", "location", "horizon"]
    s = d.groupby(keys).pl.mean().rename("wis").to_frame()
    wide = d.assign(q=q).pivot_table(index=keys, columns="q", values="value")
    yv = d.groupby(keys).oracle_value.first()
    s["cov50"] = (yv >= wide[0.25]) & (yv <= wide[0.75])
    s["cov95"] = (yv >= wide[0.025]) & (yv <= wide[0.975])
    return s.reset_index()


def rel_wis(s, by=None, n_boot=2000, seed=0):
    """Ratio of summed WIS (variant / none) on matched forecasts, bootstrap over reference dates."""
    w = s.pivot_table(index=["reference_date", "location", "horizon"], columns="variant", values="wis").dropna()
    groups = [("all", w)] if by is None else list(w.groupby(level=by))
    rng = np.random.default_rng(seed)
    rows = []
    for g, wg in groups:
        per_date = wg.groupby(level="reference_date").sum()
        dates = per_date.index.to_numpy()
        boot = rng.integers(len(dates), size=(n_boot, len(dates)))
        for v in [c for c in w.columns if c != BASE]:
            num, den = per_date[v].to_numpy(), per_date[BASE].to_numpy()
            b = num[boot].sum(1) / den[boot].sum(1)
            rows.append({"group": g, "variant": v, "rel_wis": num.sum() / den.sum(),
                         "lo95": np.quantile(b, 0.025), "hi95": np.quantile(b, 0.975),
                         "n_dates": len(dates), "n": len(wg)})
    return pd.DataFrame(rows)


def main(root):
    s = score(load(root))
    s = s[s.location != "US"]
    s["season"] = np.where(pd.to_datetime(s.reference_date) < "2024-08-01", "2023/24", "2024/25")
    pd.set_option("display.width", 200)
    print("Forecasts per variant:", s.groupby("variant").size().to_dict())
    print("\nMean WIS / coverage:")
    print(s.groupby("variant")[["wis", "cov50", "cov95"]].mean().round(3).to_string())
    print("\nRelative WIS vs no-weather (bootstrap 95% over reference dates):")
    print(rel_wis(s).round(3).to_string(index=False))
    print("\nBy horizon:")
    print(rel_wis(s, by="horizon").round(3).to_string(index=False))
    s2 = s.set_index(["reference_date", "location", "horizon"])
    for season in ["2023/24", "2024/25"]:
        print(f"\nSeason {season}:")
        print(rel_wis(s2[s2.season == season].reset_index()).round(3).to_string(index=False))
    s.to_csv(Path(root) / "scores.csv", index=False)


if __name__ == "__main__":
    main(sys.argv[1])
