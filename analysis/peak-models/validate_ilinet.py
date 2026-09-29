"""
Model development set: forecast the seasonal peak of state-level ILINet (ILI x percent positive) in "real time" for
the last two pre-COVID seasons (default 2018/19 and 2019/20), without using any NHSN data.

For each validation season S, a model is trained on every season before S of ILINet (national, HHS regions, states)
and FluSurv-NET, and then forecasts the peak week and peak size of each state's ILINet series at each reference date
r of season S, using only the weeks reported by r: ILINet for the week ending on Saturday w is published the following
Friday, so at reference date r the last reported week is the one ending r - 7. Forward chaining means the 2019/20 fit
also sees 2018/19, but nothing from S or later is used in training, and within S nothing after the last reported week.

The data are the current (final) versions: iddata has no ILINet vintages, so unlike the NHSN hindcasts there is no
revision uncertainty here (a single "draw" of the current season is used). Peak size is in ILINet units and is scored
on the log scale (see score()).

Usage (from the repository root, with the idmodels environment):
    python analysis/peak-models/validate_ilinet.py --models baseline gbqr kcde
    python analysis/peak-models/validate_ilinet.py --score_only
Forecasts go to analysis/peak-models/ilinet-validation/forecasts/<model>_<season>.parquet, scores to
analysis/peak-models/ilinet-validation/scores.csv.
"""
import argparse
import datetime
import importlib.util
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
OUT = HERE / "ilinet-validation"
sys.path.insert(0, str(ROOT / "src"))

from peak_common import Q_LEVELS  # noqa: E402

TARGET_SOURCE = "ilinet"
# offset for scoring peak size on the log scale, in ILINet units (state peaks are typically 0.5-5)
SCORE_LOG_EPS = 0.01
# reference dates: season weeks 14 (late October) through 42 (mid May), comparable to the FluSight rounds
REF_SEASON_WEEKS = range(14, 43)
# ILINet for Puerto Rico and the US Virgin Islands is unreliable (zero for whole seasons; sparse lab testing), so these
# are neither training series (PeakModelConfig.exclude_training_series) nor validation targets
EXCLUDED_LOCATIONS = ["72", "78"]
# ILINet revision model (--revisions): offset for log revision ratios, in ILINet units, and the lag-0 SD of the
# fallback perturbation used while fewer than 200 revision vectors are available (the first reference dates of
# 2018/19, since clinical-lab vintages start in October 2018)
ILINET_REVISION_OFFSET = 0.01
ILINET_FALLBACK_SD = 0.25

MODEL_CLASSES = {"baseline": "PeakBaselineModel", "gbqr": "PeakGBQRModel", "gbqr_offset": "PeakGBQRModel",
                 "kcde": "PeakKCDEModel", "hier": "PeakHierModel", "hybrid": "PeakHybridModel"}


def load_model(name: str):
    """
    A model by name. `name` is either a sandbox model directory (src/peak_<name>/main.py, whose build_model_config()
    gives the configuration) or <base>__<variant>, where the variant is defined in VARIANTS below.
    """
    import idmodels.peak as peak

    base, _, variant = name.partition("__")
    spec = importlib.util.spec_from_file_location(f"peak_{base}_main", ROOT / "src" / f"peak_{base}" / "main.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    cfg = mod.build_model_config()
    cfg.model_name = f"peak_{name}"
    if variant and name not in VARIANTS:
        raise ValueError(f"unknown variant {name}")
    for k, v in VARIANTS.get(name, {}).items():
        if not hasattr(cfg, k):
            raise ValueError(f"{name}: config has no field {k}")
        setattr(cfg, k, v)
    return getattr(peak, MODEL_CLASSES[base])(cfg)


# named configuration variants, <base>__<variant>: field overrides of the base model's configuration
VARIANTS: dict[str, dict] = {
    # synchrony and burden-to-date features (idmodels.peak.series.SYNC_BURDEN_FEATURES)
    "gbqr__sb": dict(sync_burden_features=True),
    "hier__sb": dict(sync_burden_features=True),
    # hierarchical-model ablations and settings
    "hier__tv": dict(time_varying_coefs=True),  # feature effects vary linearly with the current week
    "hier__nocur": dict(current_season_update=False),  # no censored current-season update of the season effect
    "hier__noseason": dict(season_effects=False),  # no season effects at all
    "hier__wide": dict(effect_scale=1.0),  # wider prior on the effect SDs
    "hier__w05": dict(likelihood_weight=0.05),  # stronger tempering (~2 effective observations per series)
    "hier__w40": dict(likelihood_weight=0.4),  # weaker tempering (~16 per series)
    "hier__w05c12": dict(likelihood_weight=0.05, current_update_components=[1, 2]),  # + no censored update of comp 0
    "hier__w05c12sb": dict(likelihood_weight=0.05, current_update_components=[1, 2], sync_burden_features=True),
    "hier__w02c12": dict(likelihood_weight=0.02, current_update_components=[1, 2]),
    # GBQR feature-set runs (reported-at-t synchrony; separate size / timing feature groups)
    "gbqr__R1": dict(sync_reported_only=True, size_feature_groups=["core", "sb"], timing_feature_groups=["core", "sb"]),
    "gbqr__R2": dict(sync_reported_only=True, size_offset=True, size_feature_groups=["core", "sb"],
                     timing_feature_groups=["core", "bshare"]),
    "gbqr__R3": dict(sync_reported_only=True, size_feature_groups=["core", "sb", "bshare"],
                     timing_feature_groups=["core", "latlon"]),
    "gbqr__R4": dict(sync_reported_only=True, size_feature_groups=["core", "sb", "trend"],
                     timing_feature_groups=["core", "recession"]),
    "gbqr__R5": dict(sync_reported_only=True, size_feature_groups=["core", "sb", "h3"],
                     timing_feature_groups=["core", "holiday"]),
    "gbqr__R6": dict(sync_reported_only=True, size_feature_groups=["core", "sb", "recession"],
                     timing_feature_groups=["core", "bshare", "latlon", "recession", "holiday"]),
    "gbqr__R7": dict(sync_reported_only=True, size_feature_groups=["core", "sb", "latlon"],
                     timing_feature_groups=["core", "sb", "bshare"]),
    "gbqr__R8": dict(sync_reported_only=True, size_offset=True, size_feature_groups=["core", "sb", "bshare"],
                     timing_feature_groups=["core", "sb", "latlon", "recession"]),
    "gbqr__R9": dict(sync_reported_only=True, size_feature_groups=["core", "sb", "trend", "bshare"],
                     timing_feature_groups=["core", "sb", "bshare", "latlon", "recession"]),
    "gbqr__R10": dict(sync_reported_only=True, size_offset=True, size_feature_groups=["core"],
                      timing_feature_groups=["core", "sb", "bshare", "latlon", "recession", "holiday"]),
    # GBQR configurations carried to the NHSN test seasons (nhsn_validation.py); carried-forward synchrony
    "gbqr__N1": dict(size_feature_groups=["core", "sb"], timing_feature_groups=["core"]),
    "gbqr__N2": dict(size_feature_groups=["core", "sb", "latlon"], timing_feature_groups=["core", "holiday"]),
    "gbqr__N4": dict(size_feature_groups=["core"], timing_feature_groups=["core", "holiday"]),
    "gbqr__N3": dict(size_feature_groups=["core", "sb", "bshare"],
                     timing_feature_groups=["core", "bshare", "latlon", "recession", "holiday"]),
    "hybrid__w02": dict(likelihood_weight=0.02),
    "hybrid__feat": dict(hybrid_keep_features=True),  # also keep the hierarchical model's own feature terms
    "hybrid__cu01": dict(current_update_weight=0.01),  # weaker censored current-season update
    "hybrid__cu003": dict(current_update_weight=0.003),
    "hier__w05d": dict(likelihood_weight=0.05, wsm_dummies=True),  # + indicators for weeks since max 0..3
    "hier__w05dtv": dict(likelihood_weight=0.05, wsm_dummies=True, time_varying_coefs=True),
    "hier__smoke": dict(num_warmup=30, num_samples=30, num_chains=1, origin_stride=3, num_posterior_draws=30,
                        progress_bar=True),
}


def load_strain() -> dict:
    """Final weekly A, B, A(H1), A(H3) positives by (geography, season) (WHO/NREVSS via the iddata S3 file; cached
    under eda/strain/ by the exploratory analysis)."""
    from relative_size_predictors import load_strain as _load

    return _load()


def load_data() -> pd.DataFrame:
    """ILINet and FluSurv-NET (final data), cached locally."""
    path = OUT / "data-ilinet-flusurvnet.parquet"
    if path.exists():
        return pd.read_parquet(path)
    from iddata.ancillary.population import PopulationData
    from iddata.loader import DiseaseDataLoader
    from iddata.sources.flusurvnet import FluSurvNetDataSource
    from iddata.sources.ilinet import ILINetDataSource

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df = DiseaseDataLoader().load(sources=[ILINetDataSource(), FluSurvNetDataSource()],
                                      as_of=datetime.date.today(), ancillary=[PopulationData()])
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path)
    return df


def season_targets(data: pd.DataFrame, season: str, cfg):
    """Final state and national ILINet series for `season` (SeasonArrays) and their peaks (the "oracle")."""
    from idmodels.peak.series import build_season_arrays, season_peaks

    cur = data.loc[(data["source"] == TARGET_SOURCE) & (data["season"] == season)
                   & data["agg_level"].isin(["state", "national"]) & ~data["location"].isin(EXCLUDED_LOCATIONS)]
    arrays = build_season_arrays(cur)
    w = arrays.y[:, cfg.window_start_week - 1:cfg.window_end_week]
    n_obs = np.sum(~np.isnan(w), axis=1)
    arrays = arrays.subset(n_obs >= cfg.min_window_obs)
    peak, peak_week = season_peaks(arrays.y, cfg.window_start_week, cfg.window_end_week)
    w = arrays.y[:, cfg.window_start_week - 1:cfg.window_end_week]
    tied = np.sum(np.isclose(w, peak[:, None]), axis=1) > 1
    oracle = arrays.keys[["location", "agg_level"]].assign(season=season, peak=peak, peak_week=peak_week, tied=tied)
    return arrays, oracle


def realtime_reported(data: pd.DataFrame, season: str, locations: list[str], ref_date, table: pd.DataFrame):
    """
    The season-aligned ILINet series of `locations` as published by ref_date (see epidata_ilinet.py), shape
    (n_loc, 53); weeks not yet published are missing.
    """
    import epidata_ilinet as E

    ili = data.loc[(data["source"] == TARGET_SOURCE) & data["agg_level"].isin(["state", "national"])]
    final = ili[["location", "wk_end_date", "inc"]].assign(wk_end_date=lambda d: pd.to_datetime(d["wk_end_date"]))
    asof = E.as_of_values(table, ref_date, final)
    weeks = ili.loc[ili["season"] == season, ["wk_end_date", "season_week"]].drop_duplicates()
    weeks = weeks.assign(wk_end_date=pd.to_datetime(weeks["wk_end_date"]))
    asof = asof.merge(weeks, on="wk_end_date")
    wide = asof.pivot_table(index="location", columns="season_week", values="inc", aggfunc="last")
    return wide.reindex(index=locations, columns=np.arange(1, 54, dtype=float)).to_numpy(dtype=float)


def run_model(name: str, data: pd.DataFrame, seasons: list[str], realtime: bool = False,
              revisions: bool = False) -> None:
    """
    Fit and forecast every reference date of each season. With realtime, the current season is the ILINet series as
    published by each reference date (Delphi Epidata vintages); otherwise the final series truncated at the last
    reported week. Results are labeled <name>@rt in the realtime case. With revisions (implies realtime), the model
    is also applied to simulated final versions of the current season from an ILINet revision model fit to the
    vintages known at each reference date (label <name>@rtrev).
    """
    from idmodels.peak.base import season_week_to_date

    model = load_model(name)
    cfg = model.model_config
    realtime = realtime or revisions
    label = f"{name}@rtrev" if revisions else (f"{name}@rt" if realtime else name)
    table = None
    if realtime:
        import epidata_ilinet as E

        table = E.lag_table()
    strain = load_strain()
    for season in seasons:
        out_path = OUT / "forecasts" / f"{label}_{season.replace('/', '-')}.parquet"
        t0 = time.time()
        # no NHSN rows exist in `data`; fit() further restricts training to seasons before `season`
        model.fit(data, season, Q_LEVELS, strain=strain)
        print(f"{name} {season}: fit {time.time() - t0:.0f}s", flush=True)
        if hasattr(model, "mcmc_stats_"):
            print(f"{name} {season}: MCMC {json.dumps(model.mcmc_stats_)}", flush=True)
        arrays, _ = season_targets(data, season, cfg)
        locations = arrays.keys["location"].tolist()
        weeks = np.arange(cfg.window_start_week, cfg.window_end_week + 1)
        week_dates = [str(season_week_to_date(season, int(w))) for w in weeks]
        frames = []
        for r_sw in REF_SEASON_WEEKS:
            ref_date = season_week_to_date(season, r_sw)
            last_obs = r_sw - 1  # the week ending ref_date - 7
            if realtime:
                reported = realtime_reported(data, season, locations, ref_date, table)
            else:
                reported = arrays.y.copy()
            reported[:, last_obs:] = np.nan  # column j holds season week j + 1
            ok = ~np.all(np.isnan(reported), axis=1)
            locs = [loc for loc, o in zip(locations, ok) if o]
            rev = None
            if revisions:
                from idmodels.peak.revision import RevisionModel

                rev = RevisionModel(max_lag=E.MAX_LAG, strata_bounds=(), offset=ILINET_REVISION_OFFSET,
                                    fallback_sd=ILINET_FALLBACK_SD).fit(E.revision_vintages(table, ref_date))
            cur_strain = strain
            if realtime:
                ili = data.loc[(data["source"] == TARGET_SOURCE) & (data["season"] == season)]
                week_map = ili[["wk_end_date", "season_week"]].drop_duplicates().assign(
                    wk_end_date=lambda d: pd.to_datetime(d["wk_end_date"]))
                cur_strain = E.strain_as_of(table, ref_date, season, week_map, strain)
            pmf, q = model.predict_series(reported[ok], locs, TARGET_SOURCE, nat_loc="US" if "US" in locs else None,
                                          revision_model=rev, rng=np.random.default_rng(r_sw), season=season,
                                          strain=cur_strain)
            pmf_df = pd.DataFrame(pmf, columns=week_dates).assign(location=locs).melt(
                id_vars="location", var_name="output_type_id", value_name="value").assign(output_type="pmf")
            q_df = pd.DataFrame(q, columns=[str(x) for x in Q_LEVELS]).assign(location=locs).melt(
                id_vars="location", var_name="output_type_id", value_name="value").assign(output_type="quantile")
            frames.append(pd.concat([pmf_df, q_df]).assign(reference_date=str(ref_date), last_obs_week=last_obs))
        out = pd.concat(frames, ignore_index=True).assign(model=label, season=season)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out.to_parquet(out_path)
        print(f"{label} {season}: {len(REF_SEASON_WEEKS)} reference dates, {time.time() - t0:.0f}s total", flush=True)


def wis(taus: np.ndarray, vals: np.ndarray, y: float) -> float:
    return 2 * np.mean((np.asarray(y <= vals, dtype=float) - taus) * (vals - y))


def score(data: pd.DataFrame) -> pd.DataFrame:
    """
    Scores per (model, season, reference_date, location), as in score_peak_forecasts.py: peak week log score (natural
    log), RPS and probability within +/- 1 week (tied peaks excluded); peak size WIS on the log(y + 0.01) scale, and
    50% / 95% interval coverage.
    """
    from idmodels.config import PeakBaselineModelConfig
    from idmodels.peak.base import season_week_to_date

    cfg = PeakBaselineModelConfig(model_name="scoring")
    files = sorted((OUT / "forecasts").glob("*.parquet"))
    fc = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    oracles = []
    for season in fc["season"].unique():
        _, o = season_targets(data, season, cfg)
        o["peak_date"] = [str(season_week_to_date(season, int(w))) for w in o["peak_week"]]
        oracles.append(o)
    oracle = pd.concat(oracles).set_index(["season", "location"])
    taus = np.asarray(Q_LEVELS)
    rows = []
    for (model, season, ref, loc), g in fc.groupby(["model", "season", "reference_date", "location"]):
        if (season, loc) not in oracle.index:
            continue
        tr = oracle.loc[(season, loc)]
        rec = {"model": model, "season": season, "reference_date": ref, "location": loc,
               "agg_level": tr["agg_level"], "last_obs_week": g["last_obs_week"].iloc[0],
               "weeks_from_peak": g["last_obs_week"].iloc[0] + 1 - int(tr["peak_week"])}
        pmf = g.loc[g["output_type"] == "pmf"].set_index("output_type_id")["value"]
        if not tr["tied"]:
            dates = sorted(pmf.index)
            p = pmf.reindex(dates).to_numpy()
            i = dates.index(tr["peak_date"])
            rec["log_score"] = np.log(max(p[i], 1e-300))
            rec["rps"] = np.sum((np.cumsum(p) - (np.arange(len(p)) >= i)) ** 2)
            rec["prob_pm1"] = p[max(i - 1, 0):i + 2].sum()
        q = g.loc[g["output_type"] == "quantile"].assign(tau=lambda d: d["output_type_id"].astype(float))
        q = q.sort_values("tau")
        vals, y = q["value"].to_numpy(), float(tr["peak"])
        rec["wis_log"] = wis(taus, np.log(np.maximum(vals, 0) + SCORE_LOG_EPS), np.log(y + SCORE_LOG_EPS))
        rec["cov50"] = float(np.interp(0.25, taus, vals) <= y <= np.interp(0.75, taus, vals))
        rec["cov95"] = float(np.interp(0.025, taus, vals) <= y <= np.interp(0.975, taus, vals))
        rec["pit_size"] = float(np.interp(y, vals, taus, left=0.0, right=1.0))
        rows.append(rec)
    return pd.DataFrame(rows)


def summarize(scores: pd.DataFrame) -> None:
    st = scores.loc[scores["agg_level"] == "state"]
    cols = ["log_score", "rps", "prob_pm1", "wis_log", "cov50", "cov95"]
    pd.set_option("display.width", 200)
    print("\nStates, by season:")
    print(st.groupby(["season", "model"])[cols].mean().round(3).to_string())
    print("\nStates, both seasons:")
    overall = st.groupby("model")[cols].mean()
    if "baseline" in overall.index:
        overall["rel_wis_log"] = overall["wis_log"] / overall.loc["baseline", "wis_log"]
        overall["rel_rps"] = overall["rps"] / overall.loc["baseline", "rps"]
    print(overall.sort_values("log_score", ascending=False).round(3).to_string())
    bins = pd.cut(st["weeks_from_peak"], [-99, -9, -5, -2, 1, 4, 99],
                  labels=["<=-9", "-8..-5", "-4..-2", "-1..1", "2..4", ">=5"])
    print("\nLog score by weeks from peak (current week - peak week):")
    print(st.groupby([bins, "model"], observed=True)["log_score"].mean().unstack().round(2).to_string())
    print("\nlog-WIS by weeks from peak:")
    print(st.groupby([bins, "model"], observed=True)["wis_log"].mean().unstack().round(3).to_string())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--models", nargs="*", default=[])
    parser.add_argument("--seasons", nargs="+", default=["2018/19", "2019/20"])
    parser.add_argument("--score_only", action="store_true")
    parser.add_argument("--skip_existing", action="store_true")
    parser.add_argument("--revisions", action="store_true",
                        help="realtime, plus revision simulation from ILINet vintages (label <name>@rtrev)")
    parser.add_argument("--realtime", action="store_true",
                        help="forecast from ILINet as published at each reference date (Delphi Epidata vintages)")
    args = parser.parse_args()

    data = load_data()
    assert (data["source"] != "nhsn").all()
    if not args.score_only:
        for name in args.models:
            label = f"{name}@rtrev" if args.revisions else (f"{name}@rt" if args.realtime else name)
            seasons = [s for s in args.seasons if not (args.skip_existing and (
                OUT / "forecasts" / f"{label}_{s.replace('/', '-')}.parquet").exists())]
            if seasons:
                run_model(name, data, seasons, realtime=args.realtime, revisions=args.revisions)
    scores = score(data)
    scores.to_csv(OUT / "scores.csv", index=False)
    summarize(scores)


if __name__ == "__main__":
    main()
