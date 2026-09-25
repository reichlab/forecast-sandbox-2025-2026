import sys
from pathlib import Path

import click
from idmodels.config import PeakKCDEModelConfig
from idmodels.peak import PeakKCDEModel

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from peak_common import make_run_config, reference_date_from  # noqa: E402


def build_model_config() -> PeakKCDEModelConfig:
    return PeakKCDEModelConfig(model_name="peak_kcde")


@click.command()
@click.option("--today_date", type=str, required=False, help="Date to use as effective model run date (YYYY-MM-DD)")
@click.option("--short_run", is_flag=True, help="Run with reduced parameters for faster testing")
def main(today_date: str | None = None, short_run: bool = False):
    """Generate peak week and peak size forecasts from the peak_kcde model."""
    model_config = build_model_config()
    if short_run:
        model_config.max_tuning_iter = 10
    run_config = make_run_config(reference_date_from(today_date), Path(__file__).resolve().parents[2] / "model-output")
    PeakKCDEModel(model_config).run(run_config)


if __name__ == "__main__":
    main()
