"""
GBQR model with weather features from Open-Meteo.

This module extends the base GBQRModel to include weather covariates
(temperature, humidity, precipitation, solar radiation) as features.
"""

import click
import datetime
from datetime import timedelta
from dateutil import relativedelta
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from iddata.loader import DiseaseDataLoader
from idmodels.gbqr import GBQRModel
from idmodels.preprocess import create_features_and_targets
from idmodels.utils import build_save_path

from weather_utils import (
    fetch_weather_for_locations,
    add_weather_features_with_leads,
    get_static_weather_cutoff,
    STATE_CENTROIDS,
    DEFAULT_NOISE_STD,
)


class WeatherGBQRModel(GBQRModel):
    """
    GBQR model augmented with weather features.

    This class overrides the run() method to inject weather features
    into the training and test data before model fitting.
    """

    def run(self, run_config):
        """
        Load flu data, add weather features, generate predictions, and save.

        This method is largely copied from GBQRModel.run() with the addition
        of weather feature integration.
        """
        # ===== Data Loading (same as base class) =====
        if self.model_config.reporting_adj:
            ilinet_kwargs = None
            flusurvnet_kwargs = None
        else:
            ilinet_kwargs = {"scale_to_positive": False}
            flusurvnet_kwargs = {"burden_adj": False}

        valid_sources = ["flusurvnet", "nhsn", "ilinet", "nssp"]
        if not np.isin(np.array(self.model_config.sources), valid_sources).all():
            raise ValueError(
                "For GBQR, the only supported data sources are 'nhsn', 'flusurvnet', 'ilinet', or 'nssp'."
            )

        if all(src in self.model_config.sources for src in ["nhsn", "nssp"]):
            raise ValueError("Only one of 'nhsn' or 'nssp' may be selected as a data source.")

        fdl = DiseaseDataLoader()
        if "nhsn" in self.model_config.sources:
            df = fdl.load_data(
                nhsn_kwargs={"as_of": run_config.ref_date, "disease": run_config.disease},
                ilinet_kwargs=ilinet_kwargs,
                flusurvnet_kwargs=flusurvnet_kwargs,
                sources=self.model_config.sources,
                power_transform=self.model_config.power_transform,
            )
        elif "nssp" in self.model_config.sources:
            df = fdl.load_data(
                nssp_kwargs={"as_of": run_config.ref_date, "disease": run_config.disease},
                ilinet_kwargs=ilinet_kwargs,
                flusurvnet_kwargs=flusurvnet_kwargs,
                sources=self.model_config.sources,
                power_transform=self.model_config.power_transform,
            )

        if (run_config.states == []) & (run_config.hsas == []):
            raise ValueError("User must request a non-empty set of locations to forecast for.")

        if (run_config.states != []) & (run_config.hsas != []):
            raise NotImplementedError(
                "Functionality for simultaneously forecasting state- and hsa-level locations is not yet implemented."
            )

        df_states = df.loc[(df["location"].isin(run_config.states)) & (df["agg_level"] != "hsa")]
        df_hsas = df.loc[(df["location"].isin(run_config.hsas)) & (df["agg_level"] == "hsa")]
        df = pd.concat([df_states, df_hsas], join="inner", axis=0)

        # ===== WEATHER FEATURE INTEGRATION (NEW) =====
        weather_feat_names = []
        if getattr(self.model_config, "use_weather_features", False):
            df, weather_feat_names = self._add_weather_features(df, run_config)
            print(f"Added {len(weather_feat_names)} weather features: {weather_feat_names}")

        # ===== Feature Creation =====
        if run_config.disease == "flu":
            init_feats = ["inc_trans_cs", "season_week", "log_pop"]
        elif run_config.disease == "covid":
            init_feats = ["inc_trans_cs", "log_pop"]

        # Add weather features to initial feature list
        init_feats = init_feats + weather_feat_names

        df, feat_names = create_features_and_targets(
            df=df,
            incl_level_feats=self.model_config.incl_level_feats,
            max_horizon=run_config.max_horizon,
            curr_feat_names=init_feats,
        )

        # Keep only in-season rows
        if run_config.disease == "flu":
            df = df.query("season_week >= 5 and season_week <= 45")

        # ===== Train/Test Split =====
        df_test = df.loc[df.wk_end_date == df.wk_end_date.max()].copy()
        df_train = df.loc[~df["delta_target"].isna().values]

        # ===== Model Training and Prediction =====
        if self.model_config.fit_locations_separately:
            locations = df_test["location"].unique()
            preds_df = [
                self._train_gbq_and_predict(run_config, df_train, df_test, feat_names, location)
                for location in locations
            ]
            preds_df = pd.concat(preds_df, axis=0)
        else:
            preds_df = self._train_gbq_and_predict(run_config, df_train, df_test, feat_names)

        # ===== Save Output =====
        save_path = build_save_path(
            root=run_config.output_root, run_config=run_config, model_config=self.model_config
        )
        preds_df.to_csv(save_path, index=False)

    def _add_weather_features(self, df, run_config):
        """
        Fetch weather data and merge it with the surveillance DataFrame.

        This adds three types of weather features:
        1. Current week weather: weather_temp, weather_humidity
        2. Lagged weather: weather_temp_lag1, weather_temp_lag2 (past weeks)
        3. Leading weather: weather_temp_lead1, weather_temp_lead2 (future weeks)

        Weather data sourcing logic:
        - If all required dates (through ref_date + max_lead weeks) are within
          the static historical dataset, use 100% static data (no API calls)
        - If some dates are beyond static data:
          - For past dates: fetch from Historical Forecast API
          - For future dates: leave as NaN (model handles missing values)

        Leading features always have Gaussian noise added to simulate forecast
        uncertainty, whether using observed data (retrospective) or actual
        forecasts (live/production).

        Parameters
        ----------
        df : pd.DataFrame
            Surveillance data from DiseaseDataLoader
        run_config : SimpleNamespace
            Run configuration with ref_date, etc.

        Returns
        -------
        tuple
            (augmented_df, weather_feature_names)
        """
        # Get weather config from model_config
        weather_vars = getattr(
            self.model_config,
            "weather_variables",
            ["temperature_2m_mean", "relative_humidity_2m_mean"],
        )
        weather_lags = getattr(self.model_config, "weather_lags", [1, 2])
        weather_leads = getattr(self.model_config, "weather_leads", [1, 2])
        cache_dir = getattr(self.model_config, "weather_cache_dir", ".weather_cache")
        noise_std = getattr(self.model_config, "weather_noise_std", DEFAULT_NOISE_STD)

        # Determine locations to fetch weather for
        locations = df["location"].unique()
        locations_with_centroids = [loc for loc in locations if loc in STATE_CENTROIDS]

        if len(locations_with_centroids) < len(locations):
            missing = set(locations) - set(locations_with_centroids)
            print(f"Warning: No weather centroids for locations: {missing}")

        # Determine date range for weather data
        min_date = df["wk_end_date"].min()
        max_lag = max(weather_lags) if weather_lags else 0
        max_lead = max(weather_leads) if weather_leads else 0

        # Start date: earliest flu data minus lag buffer
        # Static weather data covers 1997+ (ERA5 reanalysis)
        start_date = (min_date - timedelta(weeks=max_lag + 1)).date()

        # End date: ref_date + lead weeks (to populate leading features)
        ref_date = run_config.ref_date
        if isinstance(ref_date, datetime.datetime):
            ref_date = ref_date.date()
        end_date = ref_date + timedelta(weeks=max_lead)

        # Check if we're doing a retrospective run (all data available in static dataset)
        static_cutoff = get_static_weather_cutoff()
        today = datetime.date.today()
        is_retrospective = (static_cutoff is not None and end_date <= static_cutoff)

        if is_retrospective:
            print(f"Retrospective run: all weather data available in static dataset")
            print(f"Date range needed: {start_date} to {end_date}")
            print(f"Static data available through: {static_cutoff}")
        else:
            print(f"Current/live run: may need API calls for recent data")
            print(f"Date range needed: {start_date} to {end_date}")
            print(f"Static data available through: {static_cutoff}")
            print(f"Today's date: {today}")

        print(f"Fetching weather data for {len(locations_with_centroids)} locations")

        # Fetch weather data (hybrid: static + API as needed)
        weather_df = fetch_weather_for_locations(
            locations=locations_with_centroids,
            start_date=start_date,
            end_date=end_date,
            variables=weather_vars,
            include_forecast=False,  # We don't want Forecast API data mixed in
            cache_dir=cache_dir,
            use_static_data=True,
        )

        if weather_df.empty:
            print("Warning: No weather data retrieved")
            return df, []

        # Get weather column names (after renaming in fetch function)
        weather_cols = [c for c in weather_df.columns if c.startswith("weather_")]

        # Generate random seed based on reference date for reproducibility
        import calendar
        random_seed = int(calendar.timegm(run_config.ref_date.timetuple()))

        print(f"Creating weather features with lags={weather_lags}, leads={weather_leads}")
        print(f"Noise std for leading features: {noise_std}")

        # Add weather features with lagged and leading features + noise injection
        df, all_feat_names = add_weather_features_with_leads(
            df=df,
            weather_df=weather_df,
            weather_cols=weather_cols,
            lags=weather_lags,
            leads=weather_leads,
            noise_std=noise_std,
            random_seed=random_seed,
        )

        # Report missing weather data (LightGBM handles NaN natively)
        total_rows = len(df)
        for col in all_feat_names:
            if col in df.columns:
                n_missing = df[col].isna().sum()
                if n_missing > 0:
                    pct_missing = 100 * n_missing / total_rows
                    print(f"  {col}: {n_missing} missing ({pct_missing:.1f}%) - LightGBM will handle NaN")

        return df, all_feat_names


@click.command()
@click.option(
    "--today_date",
    type=str,
    required=False,
    help="Date to use as effective model run date (YYYY-MM-DD)",
)
@click.option(
    "--short_run",
    is_flag=True,
    help="Perform a short run with fewer bags and quantiles.",
)
def main(today_date: str, short_run: bool):
    """Generate flu predictions from GBQR model with weather features."""
    try:
        today_date = datetime.date.fromisoformat(today_date)
    except (TypeError, ValueError):
        today_date = datetime.date.today()
    reference_date = today_date + relativedelta.relativedelta(weekday=5)

    model_config = SimpleNamespace(
        model_class="gbqr",
        model_name="gbqr_3src_weather",
        # Standard GBQR config
        incl_level_feats=True,
        num_bags=100,
        bag_frac_samples=0.7,
        reporting_adj=False,
        sources=["flusurvnet", "nhsn", "ilinet"],
        fit_locations_separately=False,
        power_transform="4rt",
        # Weather feature config
        use_weather_features=True,
        weather_variables=[
            "temperature_2m_mean",
            "relative_humidity_2m_mean",
        ],
        weather_lags=[1, 2],      # 1 and 2 weeks prior (lagging features)
        weather_leads=[1, 2],     # 1 and 2 weeks ahead (leading features)
        weather_cache_dir=".weather_cache",
        # Noise injection for leading features (simulates forecast uncertainty)
        weather_noise_std={
            "weather_temp": 2.0,       # ~2°C typical 1-2 week forecast error
            "weather_humidity": 10.0,  # ~10% typical forecast error
        },
    )

    run_config = SimpleNamespace(
        disease="flu",
        ref_date=reference_date,
        output_root=Path("../../model-output"),
        artifact_store_root=None,
        save_feat_importance=False,
        max_horizon=4,
        states=[
            "US", "01", "02", "04", "05", "06", "08", "09", "10", "11",
            "12", "13", "15", "16", "17", "18", "19", "20", "21", "22",
            "23", "24", "25", "26", "27", "28", "29", "30", "31", "32",
            "33", "34", "35", "36", "37", "38", "39", "40", "41", "42",
            "44", "45", "46", "47", "48", "49", "50", "51", "53", "54",
            "55", "56", "72",
        ],
        hsas=[],
        q_levels=[
            0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
            0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70,
            0.75, 0.80, 0.85, 0.90, 0.95, 0.975, 0.99,
        ],
        q_labels=[
            "0.01", "0.025", "0.05", "0.1", "0.15", "0.2",
            "0.25", "0.3", "0.35", "0.4", "0.45", "0.5",
            "0.55", "0.6", "0.65", "0.7", "0.75", "0.8",
            "0.85", "0.9", "0.95", "0.975", "0.99",
        ],
    )

    if short_run:
        run_config.q_levels = [0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975]
        run_config.q_labels = ["0.025", "0.1", "0.25", "0.5", "0.75", "0.9", "0.975"]
        model_config.num_bags = 10

    print(f"Running WeatherGBQRModel for reference date: {reference_date}")
    model = WeatherGBQRModel(model_config)
    model.run(run_config)
    print("Done!")


if __name__ == "__main__":
    main()
