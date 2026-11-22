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
    SARIXModelConfig,
    SourceType,
)
from idmodels.sarix import SARIXFourierModel, SARIXModel

# Christmas and New Year's weeks were elevated in every season in NHSN (2022/23, 2023/24), FluSurv-NET and ILI+
# (2010/11 - 2022/23); the week before Christmas week ("christmas-1") in most seasons but not 2022/23; the week after
# New Year's week ("new_year+1") modestly in NHSN and FluSurv-NET. Thanksgiving showed no consistent effect
# (see sarix/holiday-effects/01-eda-holiday-weeks.ipynb in the sarix-flu-dev repo)
HOLIDAYS = ["christmas", "new_year"]
HOLIDAYS_XM1 = ["christmas-1", "christmas", "new_year"]
HOLIDAYS_XM1_NY1 = ["christmas-1", "christmas", "new_year", "new_year+1"]

# priors on the holiday effects, which are log multipliers. mu: the common effect per holiday ('partial'),
# Normal(0, 0.5), allowing multipliers up to about e^1 = 2.7 at two standard deviations. sd: the spread of
# location effects around it ('partial') or of the effects around 0 ('unshared' / 'shared'), HalfCauchy(0.25);
# states' holiday multipliers had an interquartile range of about +-0.15 on the log scale
CALENDAR_EFFECT_SD_PRIOR_SCALE = 0.25
CALENDAR_EFFECT_MU_PRIOR_SCALE = 0.5

# variant -> (model_name, calendar_pooling, calendar_holidays, fourier). "none" is the no-holiday baseline: the same
# model as AR6_pooled, but run with the same idmodels/sarix versions as the holiday variants so that comparisons isolate
# the holiday effects. The fp_ variants add AR6_fourierP_thetaP's Fourier seasonality (K = 2 harmonics, shared across
# locations); fp_none is their no-holiday baseline on the same versions
VARIANTS = {
    "none": ("AR6_pooled_cal_none", None, [], False),
    "shared": ("AR6_pooled_cal_shared", PoolingStrategy.SHARED, HOLIDAYS, False),
    "partial": ("AR6_pooled_cal_partial", PoolingStrategy.PARTIAL, HOLIDAYS, False),
    "xm1": ("AR6_pooled_cal_xm1", PoolingStrategy.SHARED, HOLIDAYS_XM1, False),
    "xm1_ny1": ("AR6_pooled_cal_xm1_ny1", PoolingStrategy.SHARED, HOLIDAYS_XM1_NY1, False),
    "fp_none": ("AR6_fourierP_thetaP_cal_none", None, [], True),
    "fp_xm1": ("AR6_fourierP_thetaP_cal_xm1", PoolingStrategy.SHARED, HOLIDAYS_XM1, True),
}

# Fourier settings of AR6_fourierP_thetaP (src/AR6_fourierP_thetaP/main.py)
FOURIER_K = 2
FOURIER_POOLING = PoolingStrategy.SHARED

# last reference date main() will run: 2023/24 to 2025/26, whose forecast targets end by 2026-06-20
# (LAST_VISIBLE_WEEK in sarix-flu-dev's holdout.py; 2024/25 and 2025/26 were unblinded on 2026-11-20). 2026/27 is held
# out; raise this only when it is unblinded
LAST_REFERENCE_DATE = datetime.date(2026, 5, 30)

STATES = ["US", "01", "02", "04", "05", "06", "08", "09", "10", "11", "12", "13", "15", "16", "17", "18", "19", "20",
          "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38",
          "39", "40", "41", "42", "44", "45", "46", "47", "48", "49", "50", "51", "53", "54", "55", "56", "72"]


def build_configs(variant: str, reference_date: datetime.date, short_run: bool = False,
                  output_root: Path = Path("../../model-output")) -> tuple[SARIXModelConfig, RunConfig]:
    """Model and run configs for one variant and reference date (also used by the sarix-flu-dev notebooks)."""
    if reference_date.weekday() != 5:
        raise ValueError(f"reference_date must be a Saturday, got {reference_date} ({reference_date:%A})")
    model_name, calendar_pooling, calendar_holidays, fourier = VARIANTS[variant]
    config_class, fourier_args = SARIXModelConfig, {}
    if fourier:
        config_class, fourier_args = SARIXFourierModelConfig, {"fourier_K": FOURIER_K,
                                                               "fourier_pooling": FOURIER_POOLING}
    model_config = config_class(
        **fourier_args,
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
        calendar_holidays=calendar_holidays,
        calendar_pooling=calendar_pooling,
        calendar_effect_sd_prior_scale=CALENDAR_EFFECT_SD_PRIOR_SCALE,
        calendar_effect_mu_prior_scale=CALENDAR_EFFECT_MU_PRIOR_SCALE)

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


def make_model(model_config: SARIXModelConfig) -> SARIXModel:
    """The idmodels model for a config from build_configs: SARIXFourierModel for the fp_ variants."""
    if isinstance(model_config, SARIXFourierModelConfig):
        return SARIXFourierModel(model_config)
    return SARIXModel(model_config)


@click.command()
@click.option("--variant", type=click.Choice(list(VARIANTS)), required=True,
              help="Holiday effects and model: none (baseline), shared, partial, xm1, xm1_ny1, or with Fourier "
                   "seasonality fp_none or fp_xm1 (see VARIANTS)")
@click.option("--today_date", type=click.DateTime(["%Y-%m-%d"]), required=True,
              help="Date to use as effective model run date (YYYY-MM-DD)")
@click.option("--short_run", is_flag=True, help="Run with reduced parameters for faster testing")
def main(variant: str, today_date: datetime.datetime, short_run: bool = False):
    """Generate flu hospitalization forecasts from AR(6) pooled SARIX models with and without holiday effects."""
    reference_date = today_date.date() + relativedelta.relativedelta(weekday=5)
    if reference_date > LAST_REFERENCE_DATE:
        raise click.BadParameter(f"reference date {reference_date} is after {LAST_REFERENCE_DATE}, in the held-out "
                                 "seasons", param_hint="--today_date")

    model_config, run_config = build_configs(variant, reference_date, short_run)
    model = make_model(model_config)
    model.run(run_config)
    if model.sarix_fit.num_divergences > 0:
        print(f"WARNING: {model.sarix_fit.num_divergences} divergent transitions "
              f"({model_config.model_name}, {reference_date})")


if __name__ == "__main__":
    main()
