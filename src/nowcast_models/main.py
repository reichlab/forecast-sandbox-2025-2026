"""
GBQR and SARIX models with and without a nowcast of the most recent NHSN week (idmodels.nowcast, idmodels issue #34).

Each base model is run three ways, all on the same idmodels commit, so that the only difference is the nowcast:
  <base>_nc_none   no nowcast (the control)
  <base>_nc_scale  the newest week divided by its completeness (hospitals reporting / the usual number)
  <base>_nc_state  the scaled value corrected by a factor for its location, learned from this season's settled weeks

Usage, from this directory (LightGBM needs DYLD_FALLBACK_LIBRARY_PATH on a Mac without libomp; see README.md):
  python main.py --model gbqr_1src_nc_state --today_date 2026-01-14 [--archive_dir archives]
"""

import datetime
from pathlib import Path

import click
from dateutil import relativedelta
from iddata.enums import Disease
from idmodels.config import (
    GBQRModelConfig,
    PoolingStrategy,
    PowerTransform,
    RunConfig,
    SARIXModelConfig,
    ScaleByCompletenessConfig,
    ScaleByStateConfig,
    SourceType,
)
from idmodels.gbqr import GBQRModel
from idmodels.sarix import SARIXModel

# archives of NHSN admissions and hospitals reporting: files in archive_dir (written by fetch_archives.py), or the
# Delphi Epidata archive (fetched by each run) if archive_dir is None
NOWCASTS = {
    "nc_none": lambda archive_dir: None,
    "nc_scale": lambda archive_dir: ScaleByCompletenessConfig(
        reporting_path=archive_dir and archive_dir / "nhsn_hosp_reporting.parquet"
    ),
    "nc_state": lambda archive_dir: ScaleByStateConfig(
        reporting_path=archive_dir and archive_dir / "nhsn_hosp_reporting.parquet",
        admissions_path=archive_dir and archive_dir / "nhsn_admissions.parquet",
    ),
}


def gbqr_config(name, nowcast, supplementary_sources):
    return GBQRModelConfig(
        model_name=name,
        main_source=SourceType.NHSN,
        supplementary_sources=supplementary_sources,
        incl_level_feats=True,
        num_bags=100,
        bag_frac_samples=0.7,
        reporting_adj=False,
        fit_locations_separately=False,
        power_transform=PowerTransform.FOURTH_ROOT,
        nowcast=nowcast,
    )


def sarix_config(name, nowcast, p, theta_pooling):
    return SARIXModelConfig(
        model_name=name,
        main_source=SourceType.NHSN,
        fit_locations_separately=False,
        power_transform=PowerTransform.FOURTH_ROOT,
        p=p,
        P=0,
        d=0,
        D=0,
        season_period=1,
        theta_pooling=theta_pooling,
        sigma_pooling=PoolingStrategy.NONE,
        x=[],
        num_warmup=2000,
        num_samples=2000,
        num_chains=1,
        nowcast=nowcast,
    )


# base models: the settings of the sandbox models with these names
BASES = {
    "gbqr_1src": (GBQRModel, lambda name, nc: gbqr_config(name, nc, [])),
    "gbqr_4src": (
        GBQRModel,
        lambda name, nc: gbqr_config(name, nc, [SourceType.NSSP, SourceType.ILINET, SourceType.FLUSURVNET]),
    ),
    "AR6_pooled": (SARIXModel, lambda name, nc: sarix_config(name, nc, 6, PoolingStrategy.SHARED)),
    "AR2": (SARIXModel, lambda name, nc: sarix_config(name, nc, 2, PoolingStrategy.NONE)),
}
MODELS = [f"{base}_{nc}" for base in BASES for nc in NOWCASTS]

STATES = ["US", "01", "02", "04", "05", "06", "08", "09", "10", "11", "12", "13", "15", "16", "17", "18", "19", "20",
          "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38",
          "39", "40", "41", "42", "44", "45", "46", "47", "48", "49", "50", "51", "53", "54", "55", "56", "72"]
Q_LEVELS = [0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75,
            0.80, 0.85, 0.90, 0.95, 0.975, 0.99]
Q_LABELS = ["0.01", "0.025", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6",
            "0.65", "0.7", "0.75", "0.8", "0.85", "0.9", "0.95", "0.975", "0.99"]


@click.command()
@click.option("--model", type=click.Choice(MODELS), required=True, help="Model to run")
@click.option("--today_date", type=str, required=False, help="Date to use as effective model run date (YYYY-MM-DD)")
@click.option(
    "--archive_dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    required=False,
    help="Directory with the Delphi archives saved by fetch_archives.py (default: fetch them from Delphi)",
)
def main(model: str, today_date: str | None = None, archive_dir: Path | None = None):
    """Generate flu predictions from one of the nowcast comparison models."""
    try:
        today_date = datetime.date.fromisoformat(today_date)
    except (TypeError, ValueError):  # if today_date is None or a bad format
        today_date = datetime.date.today()
    reference_date = today_date + relativedelta.relativedelta(weekday=5)

    base, nc = model.rsplit("_nc_", 1)
    model_class, make_config = BASES[base]
    model_config = make_config(model, NOWCASTS[f"nc_{nc}"](archive_dir))
    run_config = RunConfig(
        disease=Disease.FLU,
        ref_date=reference_date,
        output_root=Path("../../model-output"),
        artifact_store_root=None,
        max_horizon=4,
        states=STATES,
        hsas=[],
        q_levels=Q_LEVELS,
        q_labels=Q_LABELS,
    )
    model_class(model_config).run(run_config)


if __name__ == "__main__":
    main()
