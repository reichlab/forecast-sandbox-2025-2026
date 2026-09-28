"""
Equal-weight ensembles of saved ILINet development forecasts (validate_ilinet.py): peak week pmfs are averaged (linear
pool) and peak size quantiles are averaged level by level (Vincentization), as planned for peak_ensemble. Writes
forecasts/<name>_<season>.parquet with model <name>, then rescore with validate_ilinet.py --score_only.

Usage: python analysis/peak-models/ensemble_validation.py ens_gbqr_hier gbqr@rtrev hier__w05@rtrev
"""
import sys
from pathlib import Path

import pandas as pd

F = Path(__file__).parent / "ilinet-validation" / "forecasts"


def main():
    name, members = sys.argv[1], sys.argv[2:]
    for season in ["2018-19", "2019-20"]:
        fc = pd.concat([pd.read_parquet(F / f"{m}_{season}.parquet") for m in members], ignore_index=True)
        keys = ["season", "reference_date", "last_obs_week", "location", "output_type", "output_type_id"]
        # only forecasts that every member made
        n = fc.groupby(keys)["model"].transform("nunique")
        ens = fc.loc[n == len(members)].groupby(keys, as_index=False)["value"].mean().assign(model=name)
        pmf = ens["output_type"] == "pmf"
        ens.loc[pmf, "value"] = ens.loc[pmf, "value"] / ens.loc[pmf].groupby(["reference_date", "location"])[
            "value"].transform("sum")
        ens.to_parquet(F / f"{name}_{season}.parquet")
        print(name, season, len(ens))


if __name__ == "__main__":
    main()
