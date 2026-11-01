"""
GBQR (3 sources) with GEFS weather anomaly features.

Variants (see weather_features.FEATURE_SETS):
    none        no weather (baseline, same pipeline)
    obs         observed temp / precipitable-water anomalies for the last two weeks
    obs_fc      + GEFS forecast anomalies for the next three weeks
    obs_oracle  + observed (perfect-foresight) anomalies for the next three weeks; diagnostic only

Written against idmodels before the April 2026 refactor (idmodels.preprocess),
like src/gbqr_3src_weather.
"""

import datetime
from pathlib import Path
from types import SimpleNamespace

import click
import numpy as np
import pandas as pd
from dateutil import relativedelta
from iddata.loader import DiseaseDataLoader
from idmodels.gbqr import GBQRModel
from idmodels.preprocess import create_features_and_targets
from idmodels.utils import build_save_path

from weather_features import add_weather_features

STATES = [
    "US", "01", "02", "04", "05", "06", "08", "09", "10", "11",
    "12", "13", "15", "16", "17", "18", "19", "20", "21", "22",
    "23", "24", "25", "26", "27", "28", "29", "30", "31", "32",
    "33", "34", "35", "36", "37", "38", "39", "40", "41", "42",
    "44", "45", "46", "47", "48", "49", "50", "51", "53", "54",
    "55", "56", "72",
]
Q_FULL = [0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
          0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.975, 0.99]
Q_SHORT = [0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975]


class GEFSWeatherGBQRModel(GBQRModel):
    """GBQR with weather features merged in before feature/target construction."""

    def run(self, run_config):
        mc = self.model_config
        fdl = DiseaseDataLoader()
        df = fdl.load_data(
            nhsn_kwargs={"as_of": run_config.ref_date, "disease": run_config.disease},
            ilinet_kwargs={"scale_to_positive": False},
            flusurvnet_kwargs={"burden_adj": False},
            sources=mc.sources,
            power_transform=mc.power_transform,
        )
        # same location filter as gbqr_3src (drops ILINet HHS-region rows)
        df = df.loc[df["location"].isin(run_config.states) & (df["agg_level"] != "hsa")]
        df["wk_end_date"] = pd.to_datetime(df["wk_end_date"])

        df, weather_feats = add_weather_features(df, mc.weather_variant, ref_date=run_config.ref_date)
        if weather_feats:
            n_missing = df[weather_feats].isna().mean().round(3)
            print(f"Weather features ({mc.weather_variant}): missing fraction {n_missing.max()}")

        init_feats = ["inc_trans_cs", "season_week", "log_pop"] + weather_feats
        df, feat_names = create_features_and_targets(
            df=df,
            incl_level_feats=mc.incl_level_feats,
            max_horizon=run_config.max_horizon,
            curr_feat_names=init_feats,
        )
        df = df.query("season_week >= 5 and season_week <= 45")

        df_test = df.loc[df.wk_end_date == df.wk_end_date.max()].copy()
        df_train = df.loc[~df["delta_target"].isna().values]
        preds_df = self._train_gbq_and_predict(run_config, df_train, df_test, feat_names)

        save_path = build_save_path(root=run_config.output_root, run_config=run_config, model_config=mc)
        preds_df.to_csv(save_path, index=False)


@click.command()
@click.option("--today_date", type=str, required=False, help="Effective run date (YYYY-MM-DD)")
@click.option("--variant", type=click.Choice(["none", "obs", "obs_fc", "obs_oracle"]), default="obs_fc")
@click.option("--num_bags", type=int, default=100)
@click.option("--short_q", is_flag=True, help="Use 7 quantile levels instead of 23")
@click.option("--output_root", type=str, default="../../model-output")
@click.option("--artifact_root", type=str, default=None, help="Save feature importance here if given")
def main(today_date, variant, num_bags, short_q, output_root, artifact_root):
    try:
        today_date = datetime.date.fromisoformat(today_date)
    except (TypeError, ValueError):
        today_date = datetime.date.today()
    reference_date = today_date + relativedelta.relativedelta(weekday=5)

    model_config = SimpleNamespace(
        model_class="gbqr",
        model_name=f"gbqr_3src_gefs_{variant}",
        incl_level_feats=True,
        num_bags=num_bags,
        bag_frac_samples=0.7,
        reporting_adj=False,
        sources=["flusurvnet", "nhsn", "ilinet"],
        fit_locations_separately=False,
        power_transform="4rt",
        weather_variant=variant,
        save_feat_importance=artifact_root is not None,
    )
    q = Q_SHORT if short_q else Q_FULL
    run_config = SimpleNamespace(
        disease="flu",
        ref_date=reference_date,
        output_root=Path(output_root),
        artifact_store_root=Path(artifact_root) if artifact_root else None,
        max_horizon=4,
        states=STATES,
        hsas=[],
        q_levels=q,
        q_labels=[str(x) for x in q],
    )
    print(f"Running GEFSWeatherGBQRModel ({variant}) for reference date {reference_date}")
    GEFSWeatherGBQRModel(model_config).run(run_config)


if __name__ == "__main__":
    main()
