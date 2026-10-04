"""
Data for the revision-simulation slides: the revision vectors that idmodels.peak.revision.RevisionModel learns from the
NHSN releases available on one forecast date, and the revised draws it produces for one location on that date.

Writes eda/revision_vectors.parquet (one row per vector x lag) and eda/revision_example.parquet (reported, final and
draw values for the example location). Run from the repository root with the idmodels environment:
    python analysis/peak-models/revision_example.py --ref_date 2025-01-11 --location US
"""
import argparse
import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from idmodels.peak.revision import RevisionModel, apply_revisions, load_nhsn_vintages

HERE = Path(__file__).parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref_date", default="2025-01-11")
    parser.add_argument("--location", default="US")
    parser.add_argument("--num_draws", type=int, default=200)
    args = parser.parse_args()
    ref = datetime.date.fromisoformat(args.ref_date)

    vint = load_nhsn_vintages(ref)
    rm = RevisionModel(max_lag=10).fit(vint)
    vec = pd.DataFrame(rm.vectors, columns=range(rm.max_lag))
    vec["vector"] = np.arange(len(vec))
    vec["stratum"] = rm.strata
    vec = vec.melt(id_vars=["vector", "stratum"], var_name="lag", value_name="rho")
    vec["ref_date"] = args.ref_date
    vec.to_parquet(HERE / "eda" / "revision_vectors.parquet", index=False)

    # the example location: the release in use on the forecast date, its revised draws, and the final values
    latest = vint.loc[(vint["as_of"] == vint["as_of"].max()) & (vint["location"] == args.location)]
    latest = latest.dropna(subset=["inc"]).sort_values("wk_end_date")
    counts = latest["inc"].to_numpy(float)[None, :]
    last_week = np.array([counts.shape[1] - 1])
    rho = rm.sample(counts[:, -1], args.num_draws, np.random.default_rng(42))
    revised = apply_revisions(counts, last_week, rho)[0]  # (num_draws, weeks)
    final = load_nhsn_vintages(datetime.date.today())
    final = final.loc[(final["as_of"] == final["as_of"].max()) & (final["location"] == args.location)]
    final = final.set_index("wk_end_date")["inc"]

    dates = latest["wk_end_date"].to_numpy()
    rows = [pd.DataFrame({"date": dates, "series": "reported", "draw": -1, "value": counts[0]}),
            pd.DataFrame({"date": dates, "series": "final", "draw": -1, "value": final.reindex(dates).to_numpy()})]
    for d in range(args.num_draws):
        rows.append(pd.DataFrame({"date": dates, "series": "draw", "draw": d, "value": revised[d]}))
    out = pd.concat(rows, ignore_index=True)
    out["ref_date"], out["location"], out["n_vectors"] = args.ref_date, args.location, len(rm.vectors)
    out["stratum"] = int(rm._stratum(counts[:, -1])[0])
    out.to_parquet(HERE / "eda" / "revision_example.parquet", index=False)
    print(f"{len(rm.vectors)} vectors; example stratum {out['stratum'].iloc[0]}")


if __name__ == "__main__":
    main()
