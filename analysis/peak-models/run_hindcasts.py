"""
Run the peak_* models for every reference date of one season in hub-config/tasks.json, writing hub-formatted files
to model-output/UMass-peak_<model>/.

Each model is fit once per season (training uses only earlier seasons, so the fit does not change within a season)
and the data loaded for a reference date are shared across models. NHSN data are loaded as of each reference date,
so the current-season inputs are exactly what was available in real time.

Usage (from the repository root, with an environment that has idmodels installed):
    python analysis/peak-models/run_hindcasts.py --season 2025/26 --models baseline gbqr kcde
"""
import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from peak_common import make_run_config  # noqa: E402

MODEL_CLASSES = {"baseline": "PeakBaselineModel", "gbqr": "PeakGBQRModel", "gbqr_offset": "PeakGBQRModel",
                 "kcde": "PeakKCDEModel"}


def load_model_config(name: str):
    spec = importlib.util.spec_from_file_location(f"peak_{name}_main", ROOT / "src" / f"peak_{name}" / "main.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.build_model_config()


def season_reference_dates(season: str) -> list[str]:
    tasks = json.loads((ROOT / "hub-config" / "tasks.json").read_text())
    for mt in tasks["rounds"][0]["model_tasks"]:
        if mt["task_ids"]["target"]["required"] != ["peak week inc flu hosp"]:
            continue
        dates = mt["output_type"]["pmf"]["output_type_id"]["required"]
        y = int(dates[0][:4])
        if season == f"{y}/{str(y + 1)[2:]}":
            return mt["task_ids"]["reference_date"]["optional"]
    raise ValueError(f"no peak week task for season {season}")


def main():
    import datetime

    import idmodels.peak as peak

    parser = argparse.ArgumentParser()
    parser.add_argument("--season", required=True)
    parser.add_argument("--models", nargs="+", default=["baseline", "gbqr", "kcde"], choices=list(MODEL_CLASSES))
    parser.add_argument("--output_root", default=str(ROOT / "model-output"))
    parser.add_argument("--skip_existing", action="store_true")
    args = parser.parse_args()

    models = {name: getattr(peak, MODEL_CLASSES[name])(load_model_config(name)) for name in args.models}
    for ref in season_reference_dates(args.season):
        ref_date = datetime.date.fromisoformat(ref)
        run_config = make_run_config(ref_date, Path(args.output_root))
        todo = {n: m for n, m in models.items()
                if not (args.skip_existing and (Path(args.output_root) / f"UMass-peak_{n}" /
                                                f"{ref}-UMass-peak_{n}.csv").exists())}
        if not todo:
            continue
        t0 = time.time()
        inputs = next(iter(todo.values())).load_inputs(ref_date)
        for name, model in todo.items():
            t1 = time.time()
            df = model.forecast(inputs, run_config)
            out_dir = Path(args.output_root) / f"UMass-peak_{name}"
            out_dir.mkdir(parents=True, exist_ok=True)
            df.to_csv(out_dir / f"{ref}-UMass-peak_{name}.csv", index=False, na_rep="NA")
            print(f"{ref} {name}: {time.time() - t1:.0f}s", flush=True)
        print(f"{ref} total {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
