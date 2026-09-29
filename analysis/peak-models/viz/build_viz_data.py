"""
Collect peak forecasts, observed weekly admissions and observed peaks into data.json for the peak forecast viewer
(analysis/peak-models/viz/index.html).

Forecasts come from the sandbox peak models (model-output/UMass-peak_*) and, optionally, from the real-time peak
forecasts submitted to the CDC FluSight hub (a local clone of cdcepi/FluSight-forecast-hub).

Usage (from the repository root):
    python analysis/peak-models/viz/build_viz_data.py [--hub_root ../FluSight-forecast-hub]
"""
import argparse
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).parent
# final candidate models; labels are "<size features> / <timing features>" for the GBQR models
SANDBOX_MODELS = [("UMass-peak_gbqr_sb", "GBQR SB / core"), ("UMass-peak_gbqr_core_hol", "GBQR core / core+hol"),
                  ("UMass-peak_baseline", "Baseline")]
HUB_FIRST = ["FluSight-ensemble", "FluSight-base_seasonal"]  # listed first among hub models
DEFAULT_ON = {"UMass-peak_gbqr_sb", "UMass-peak_gbqr_core_hol", "UMass-peak_baseline", "FluSight-ensemble"}
MIN_HUB_REFS = 8  # hub models need at least this many reference dates in a season to be shown
WEEK_TARGET, SIZE_TARGET = "peak week inc flu hosp", "peak inc flu hosp"


def season_of(date: str) -> str:
    y = int(date[:4]) if int(date[5:7]) >= 7 else int(date[:4]) - 1
    return f"{y}/{str(y + 1)[2:]}"


def read_peak_rows(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
    else:
        df = pd.read_csv(path, dtype={"location": str, "output_type_id": str}, low_memory=False)
    df = df.loc[df["target"].isin([WEEK_TARGET, SIZE_TARGET])]
    if len(df):
        df = df.assign(location=df["location"].astype(str), output_type_id=df["output_type_id"].astype(str),
                       reference_date=df["reference_date"].astype(str).str[:10])
    return df


def add_as_of_data(seasons: dict) -> None:
    """
    For each season and forecast date, the weekly admissions of every location as published in the latest NHSN (or,
    before 2024-11-15, HHS) data release on or before the forecast date, i.e. the data the models used. Stored as
    seasons[season]["asof"][ref][location] = counts per window date (null where not yet reported).
    """
    import datetime

    from idmodels.peak.revision import load_nhsn_vintages

    last = max((r for s in seasons.values() for r in s["refs"]), default=None)
    if last is None:
        return
    v = load_nhsn_vintages(datetime.date.fromisoformat(last))
    v = v.assign(wk=pd.to_datetime(v["wk_end_date"]).dt.strftime("%Y-%m-%d"))
    as_ofs = sorted(v["as_of"].unique())
    by_as_of = {a: g for a, g in v.groupby("as_of")}
    for s in seasons.values():
        s["asof"], s["asof_date"] = {}, {}
        for ref in s["refs"]:
            ref_ts = pd.Timestamp(ref)
            avail = [a for a in as_ofs if a <= ref_ts]
            if not avail:
                continue
            g = by_as_of[avail[-1]]
            g = g.loc[g["wk"].isin(s["dates"]) & (pd.to_datetime(g["wk"]) < ref_ts)]
            s["asof"][ref] = {loc: [None if pd.isna(x) else int(round(x)) for x in
                                    h.drop_duplicates("wk", keep="last").set_index("wk")["inc"].reindex(s["dates"])]
                              for loc, h in g.groupby("location")}
            s["asof_date"][ref] = str(pd.Timestamp(avail[-1]).date())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hub_root", default=str(ROOT.parent / "FluSight-forecast-hub"))
    parser.add_argument("--no_vintages", action="store_true",
                        help="skip the as-published weekly admissions (NHSN data vintages) for each forecast date")
    parser.add_argument("--extra", nargs=3, action="append", default=[], metavar=("MODEL_ID", "LABEL", "DIR"),
                        help="an additional sandbox model version to include, e.g. an earlier run kept elsewhere")
    args = parser.parse_args()

    tasks = json.loads((ROOT / "hub-config" / "tasks.json").read_text())
    windows = {}
    for mt in tasks["rounds"][0]["model_tasks"]:
        if mt["task_ids"]["target"]["required"] == [WEEK_TARGET]:
            dates = mt["output_type"]["pmf"]["output_type_id"]["required"]
            windows[season_of(dates[0])] = dates

    locs = pd.read_csv(ROOT / "auxiliary-data" / "locations.csv", dtype=str)
    oracle = pd.read_csv(ROOT / "target-data" / "oracle-output.csv", dtype={"location": str}, low_memory=False)
    wk = (oracle.loc[oracle["target"] == "wk inc flu hosp", ["target_end_date", "location", "oracle_value"]]
          .drop_duplicates(["target_end_date", "location"]))
    peaks = pd.read_csv(HERE.parent / "peak-oracle.csv", dtype={"location": str})

    seasons = {}
    for season, dates in windows.items():
        w = wk.loc[wk["target_end_date"].isin(dates)]
        obs = {loc: [None if pd.isna(v) else float(v)
                     for v in g.set_index("target_end_date")["oracle_value"].reindex(dates)]
               for loc, g in w.groupby("location")}
        truth = {r.location: {"week": dates.index(r.peak_week), "count": float(r.peak_count), "tied": bool(r.tied)}
                 for r in peaks.loc[peaks["season"] == season].itertuples()}
        seasons[season] = {"dates": dates, "obs": obs, "truth": truth, "forecasts": {}, "models": []}

    # model directories: sandbox models, then hub models with peak forecasts
    model_dirs = [(m, lbl, "sandbox", ROOT / "model-output" / m) for m, lbl in SANDBOX_MODELS]
    model_dirs += [(m, lbl, "sandbox", Path(d)) for m, lbl, d in args.extra]
    hub_out = Path(args.hub_root) / "model-output"
    if hub_out.exists():
        hub = sorted(d.name for d in hub_out.iterdir() if d.is_dir())
        hub = [m for m in HUB_FIRST if m in hub] + [m for m in hub if m not in HUB_FIRST]
        model_dirs += [(m, m, "hub", hub_out / m) for m in hub]

    models = []
    for model_id, label, group, d in model_dirs:
        refs_by_season: dict[str, int] = {}
        per_model = {}
        for f in sorted(list(d.glob("*.csv")) + list(d.glob("*.parquet"))):
            df = read_peak_rows(f)
            if len(df) == 0:
                continue
            ref = df["reference_date"].iloc[0]
            season = season_of(ref)
            if season not in seasons:
                continue
            dates = seasons[season]["dates"]
            by_loc = {}
            for loc, g in df.groupby("location"):
                entry = {}
                pmf = g.loc[(g["target"] == WEEK_TARGET) & (g["output_type"] == "pmf")]
                if len(pmf):
                    p = pmf.set_index("output_type_id")["value"].astype(float)
                    p = p[~p.index.duplicated()].reindex(dates).fillna(0.0)
                    entry["p"] = [int(round(v * 1e4)) for v in p]  # probability x 10,000
                q = g.loc[(g["target"] == SIZE_TARGET) & (g["output_type"] == "quantile")]
                if len(q):
                    q = q.assign(tau=q["output_type_id"].astype(float)).drop_duplicates("tau").sort_values("tau")
                    if len(q) == 23:
                        entry["q"] = [int(round(float(v))) for v in q["value"]]
                if entry:
                    by_loc[loc] = entry
            if by_loc:
                per_model.setdefault(season, {})[ref] = by_loc
                refs_by_season[season] = refs_by_season.get(season, 0) + 1
        kept = [s for s, n in refs_by_season.items() if group == "sandbox" or n >= MIN_HUB_REFS]
        if not kept:
            continue
        for season in kept:
            seasons[season]["models"].append(model_id)
            for ref, by_loc in per_model[season].items():
                seasons[season]["forecasts"].setdefault(ref, {})[model_id] = by_loc
        models.append({"id": model_id, "label": label, "group": group, "on": model_id in DEFAULT_ON})

    for s in seasons.values():
        s["refs"] = sorted(s["forecasts"])
    if not args.no_vintages:
        add_as_of_data(seasons)
    out = {
        "models": models,
        "p_scale": 1e4,
        "q_levels": [0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
                     0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99],
        "locations": [{"id": r.location, "abbr": r.abbreviation, "name": r.location_name} for r in locs.itertuples()],
        "seasons": seasons,
    }
    text = json.dumps(out, separators=(",", ":"), allow_nan=False)
    (HERE / "data.json").write_text(text)
    # the published page loads data-<content hash>.json, so browsers and caches never serve an older data file;
    # write that copy and point index.html at it (publish index.html with the hashed file)
    import hashlib
    import re

    name = f"data-{hashlib.sha256(text.encode()).hexdigest()[:10]}.json"
    for old in HERE.glob("data-*.json"):
        old.unlink()
    (HERE / name).write_text(text)
    page = HERE / "index.html"
    page.write_text(re.sub(r'fetch\("data[^"]*\.json"', f'fetch("{name}"', page.read_text()))
    print("page now loads", name)
    print({k: (len(v["refs"]), len(v["models"])) for k, v in seasons.items()},
          (HERE / "data.json").stat().st_size // 1024, "KB")
    print([m["id"] for m in models])


if __name__ == "__main__":
    main()
