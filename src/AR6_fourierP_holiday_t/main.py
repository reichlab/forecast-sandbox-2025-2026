import datetime
from pathlib import Path

import click
from dateutil import relativedelta
from iddata.enums import Disease
from idmodels.config import (
    PoolingStrategy,
    PowerTransform,
    RunConfig,
    SARIXFourierModelConfig,
    SourceType,
)
from idmodels.sarix import SARIXFourierModel

# AR6_fourierP_thetaP_cal_xm1 (src/AR6_pooled_holiday, variant fp_xm1) with Student's t innovations. Everything but
# the innovation distribution is as in fp_xm1: Fourier seasonality (K = 2, shared across locations) and holiday
# effects for the week before Christmas week and the Christmas and New Year's weeks, shared across locations
HOLIDAYS_XM1 = ["christmas-1", "christmas", "new_year"]
CALENDAR_EFFECT_SD_PRIOR_SCALE = 0.25
CALENDAR_EFFECT_MU_PRIOR_SCALE = 0.5
FOURIER_K = 2
FOURIER_POOLING = PoolingStrategy.SHARED

# variant -> (model_name, innovation_df_pooling). The t innovations are standardized so that sigma is the innovation
# sd, with degrees of freedom above 2 and sarix's default prior on them (see sarix's SARIX docs)
VARIANTS = {
    "tdf_shared": ("AR6_fP_tP_cal_xm1_tdf_shared", PoolingStrategy.SHARED),
    "tdf_partial": ("AR6_fP_tP_cal_xm1_tdf_partial", PoolingStrategy.PARTIAL),
}

# last reference date main() will run: 2023/24 to 2025/26 (LAST_VISIBLE_WEEK in sarix-flu-dev's holdout.py is
# 2026-06-20). 2026/27 is held out; raise this only when it is unblinded
LAST_REFERENCE_DATE = datetime.date(2026, 5, 30)

STATES = ["US", "01", "02", "04", "05", "06", "08", "09", "10", "11", "12", "13", "15", "16", "17", "18", "19", "20",
          "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38",
          "39", "40", "41", "42", "44", "45", "46", "47", "48", "49", "50", "51", "53", "54", "55", "56", "72"]


def build_configs(variant: str, reference_date: datetime.date, short_run: bool = False,
                  output_root: Path = Path("../../model-output")) -> tuple[SARIXFourierModelConfig, RunConfig]:
    """Model and run configs for one variant and reference date."""
    if reference_date.weekday() != 5:
        raise ValueError(f"reference_date must be a Saturday, got {reference_date} ({reference_date:%A})")
    model_name, innovation_df_pooling = VARIANTS[variant]
    model_config = SARIXFourierModelConfig(
        model_name=model_name,
        main_source=SourceType.NHSN,
        fit_locations_separately=False,
        power_transform=PowerTransform.FOURTH_ROOT,
        p=6,
        P=0,
        d=0,
        D=0,
        season_period=1,
        theta_pooling=PoolingStrategy.SHARED,
        sigma_pooling=PoolingStrategy.NONE,
        x=[],
        num_warmup=2000,
        num_samples=2000,
        num_chains=1,
        fourier_K=FOURIER_K,
        fourier_pooling=FOURIER_POOLING,
        calendar_holidays=HOLIDAYS_XM1,
        calendar_pooling=PoolingStrategy.SHARED,
        calendar_effect_sd_prior_scale=CALENDAR_EFFECT_SD_PRIOR_SCALE,
        calendar_effect_mu_prior_scale=CALENDAR_EFFECT_MU_PRIOR_SCALE,
        innovation_dist="t",
        innovation_df_pooling=innovation_df_pooling)

    run_config = RunConfig(
        disease=Disease.FLU,
        ref_date=reference_date,
        output_root=output_root,
        artifact_store_root=None,
        max_horizon=4,
        states=STATES,
        hsas=[],
        q_levels=[0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
                  0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70,
                  0.75, 0.80, 0.85, 0.90, 0.95, 0.975, 0.99],
        q_labels=["0.01", "0.025", "0.05", "0.1", "0.15", "0.2",
                  "0.25", "0.3", "0.35", "0.4", "0.45", "0.5",
                  "0.55", "0.6", "0.65", "0.7", "0.75", "0.8",
                  "0.85", "0.9", "0.95", "0.975", "0.99"])

    if short_run:
        run_config.q_levels = [0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975]
        run_config.q_labels = ["0.025", "0.1", "0.25", "0.5", "0.75", "0.9", "0.975"]
        model_config.num_warmup = 200
        model_config.num_samples = 200

    return model_config, run_config


@click.command()
@click.option("--variant", type=click.Choice(list(VARIANTS)), required=True,
              help="Pooling of the t innovations' degrees of freedom across locations: tdf_shared or tdf_partial")
@click.option("--today_date", type=click.DateTime(["%Y-%m-%d"]), required=True,
              help="Date to use as effective model run date (YYYY-MM-DD)")
@click.option("--short_run", is_flag=True, help="Run with reduced parameters for faster testing")
def main(variant: str, today_date: datetime.datetime, short_run: bool = False):
    """Generate flu hospitalization forecasts from AR(6) SARIX models with Fourier seasonality, holiday effects and
    t innovations."""
    reference_date = today_date.date() + relativedelta.relativedelta(weekday=5)
    if reference_date > LAST_REFERENCE_DATE:
        raise click.BadParameter(f"reference date {reference_date} is after {LAST_REFERENCE_DATE}, in the held-out "
                                 "seasons", param_hint="--today_date")

    model_config, run_config = build_configs(variant, reference_date, short_run)
    model = SARIXFourierModel(model_config)
    model.run(run_config)
    if model.sarix_fit.num_divergences > 0:
        print(f"WARNING: {model.sarix_fit.num_divergences} divergent transitions "
              f"({model_config.model_name}, {reference_date})")


if __name__ == "__main__":
    main()
