"""Orchestrator for the gbqr_3src_spatial_wisar(_nssp) pipeline: runs all
component models and the R ensembling/rbind step with one command.

Components (gbqr, gbqr_us, wisar, nssp run concurrently - none of them
depend on each other's output, only the final ensemble step does):
  1. gbqr/    - GBQR 3-source spatial model (state-level forecasts; no US)
  2. gbqr_us/ - GBQR 3-source model, no spatial features (US-level forecast)
  3. wisar/   - WISAR6_fourthroot_adaptive_t (state + US)
  4. nssp/    - gbqr_4src_nssp (independent model, different target: "wk inc
     flu prop ed visits"; only valid from 2025/26 onward)
  5. 3_ensemble.R - simple_ensemble() of gbqr (state) + gbqr_us (US) with
     wisar, published alone to ../../model-output/UMass-gbqr_3src_spatial_wisar/,
     then row-bound with nssp's output (when present for this reference
     date) and published to
     ../../model-output/UMass-gbqr_3src_spatial_wisar_nssp/

Usage (from this directory):
    python main.py --today_date=2024-12-07
    python main.py --today_date=2024-12-07 --short_run

Per-component logs are written to logs/<component>-<reference_date>.log
(also gitignored scratch space) so concurrent output doesn't interleave on
the terminal; each is also tailed live with a "[<component>]" prefix.

See README.md for environment setup (three Python venvs: idmodels' pandas~=2.0
pin conflicts with WISAR6's pandas>=3.0/PyMC stack, so gbqr/gbqr_us share one
venv, wisar/ has its own, and nssp/ has its own since it pins a different
idmodels commit than gbqr/gbqr_us).
"""

import datetime
import os
import shutil
import subprocess
import sys
import threading
from pathlib import Path
from typing import Optional

import click
from dateutil import relativedelta

HERE = Path(__file__).resolve().parent
GBQR_VENV_PYTHON = HERE / ".venv-gbqr" / "bin" / "python"
WISAR_VENV_PYTHON = HERE / "wisar" / ".venv" / "bin" / "python"
NSSP_VENV_PYTHON = HERE / ".venv-nssp" / "bin" / "python"

# NSSP's reporting lag makes some holiday-adjacent reference dates produce
# invalid output; gbqr_4src_nssp's own run-all-forecasts.sh skips these.
# When the nssp step fails or is skipped for other dates, 3_ensemble.R falls
# back to publishing the hosp ensemble alone under the combined model name.
NSSP_KNOWN_INVALID_DATES = {"2025-11-29", "2025-12-27", "2026-01-03"}

# gbqr, gbqr_us, and nssp all fit LightGBM models, which by default use every
# available core; running three of them at once would oversubscribe the
# machine (the same reason Unity's submit scripts pin these explicitly). Cap
# each to a share of the available cores, leaving wisar's own MCMC chains
# (controlled separately by WISAR6_MCMC_CHAINS/CORES) out of this split.
_GBQR_STEP_COUNT = 3  # gbqr, gbqr_us, nssp
_CPU_COUNT = os.cpu_count() or 4
_THREADS_PER_GBQR_STEP = max(1, _CPU_COUNT // _GBQR_STEP_COUNT)


class ComponentRun:
    """One component's subprocess: launches immediately, streams its output
    live to both a per-component log file and the console (prefixed), and
    exposes wait()/returncode like subprocess.Popen."""

    def __init__(self, label: str, cmd: list[str], cwd: Path, log_path: Path, env: Optional[dict] = None):
        self.label = label
        self.log_path = log_path
        full_env = {**os.environ, "PYTHONUNBUFFERED": "1", **(env or {})}
        self.log_file = open(log_path, "w")
        self.proc = subprocess.Popen(
            cmd, cwd=cwd, env=full_env,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1,
        )
        self._thread = threading.Thread(target=self._pump, daemon=True)
        self._thread.start()

    def _pump(self):
        for line in self.proc.stdout:
            self.log_file.write(line)
            self.log_file.flush()
            print(f"[{self.label}] {line}", end="")

    def wait(self) -> int:
        returncode = self.proc.wait()
        self._thread.join()
        self.log_file.close()
        return returncode


def _gbqr_step_env():
    # idmodels' bagging uses joblib with n_jobs=-1, which ignores the BLAS-level
    # OMP/MKL/OPENBLAS env vars entirely - LOKY_MAX_CPU_COUNT is the one joblib's
    # default (loky) backend actually respects to cap how many workers -1 means.
    n = str(_THREADS_PER_GBQR_STEP)
    return {
        "OMP_NUM_THREADS": n, "MKL_NUM_THREADS": n, "OPENBLAS_NUM_THREADS": n,
        "LOKY_MAX_CPU_COUNT": n,
    }


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
    help="Perform a short run (gbqr/gbqr_us/nssp only; wisar/ has no short-run mode).",
)
def main(today_date: Optional[str] = None, short_run: bool = False):
    """Generate flu predictions from the gbqr_3src_spatial_wisar(_nssp) pipeline."""
    try:
        today_date = datetime.date.fromisoformat(today_date)
    except (TypeError, ValueError):  # if today_date is None or a bad format
        today_date = datetime.date.today()
    reference_date = today_date + relativedelta.relativedelta(weekday=5)

    short_run_flag = ["--short_run"] if short_run else []

    print(f"Running gbqr_3src_spatial_wisar(_nssp) pipeline for reference date: {reference_date}")
    print(f"({_CPU_COUNT} CPUs detected; capping each of gbqr/gbqr_us/nssp to {_THREADS_PER_GBQR_STEP} threads)")
    print("=" * 60)

    for python_bin, label in [
        (GBQR_VENV_PYTHON, "gbqr/gbqr_us"),
        (WISAR_VENV_PYTHON, "wisar"),
        (NSSP_VENV_PYTHON, "nssp"),
    ]:
        if not python_bin.exists():
            print(f"ERROR: {label} venv not found at {python_bin}", file=sys.stderr)
            print("See README.md for environment setup.", file=sys.stderr)
            sys.exit(1)

    logs_dir = HERE / "logs"
    logs_dir.mkdir(exist_ok=True)
    skip_nssp = str(reference_date) in NSSP_KNOWN_INVALID_DATES

    print("\nLaunching gbqr, gbqr_us, wisar"
          + ("" if skip_nssp else ", nssp") + " concurrently...")
    if skip_nssp:
        print(f"  (nssp skipped: {reference_date} is a known-invalid NSSP reference date, "
              "holiday reporting lag)")

    runs = [
        ComponentRun(
            "gbqr", [str(GBQR_VENV_PYTHON), "main.py", "--today_date", str(reference_date)] + short_run_flag,
            cwd=HERE / "gbqr", log_path=logs_dir / f"gbqr-{reference_date}.log", env=_gbqr_step_env(),
        ),
        ComponentRun(
            "gbqr_us", [str(GBQR_VENV_PYTHON), "main.py", "--today_date", str(reference_date)] + short_run_flag,
            cwd=HERE / "gbqr_us", log_path=logs_dir / f"gbqr_us-{reference_date}.log", env=_gbqr_step_env(),
        ),
        ComponentRun(
            "wisar", [str(WISAR_VENV_PYTHON), "main.py", "--today_date", str(reference_date)],
            cwd=HERE / "wisar", log_path=logs_dir / f"wisar-{reference_date}.log",
        ),
    ]
    if not skip_nssp:
        runs.append(ComponentRun(
            "nssp", [str(NSSP_VENV_PYTHON), "main.py", "--today_date", str(reference_date)] + short_run_flag,
            cwd=HERE / "nssp", log_path=logs_dir / f"nssp-{reference_date}.log", env=_gbqr_step_env(),
        ))

    required = {"gbqr", "gbqr_us", "wisar"}  # nssp is optional; its absence just means hosp-only output
    failed_required = []
    for run in runs:
        returncode = run.wait()
        if returncode == 0:
            print(f"\n[{run.label}] finished OK")
        elif run.label in required:
            print(f"\n[{run.label}] FAILED (exit {returncode}) - see {run.log_path}", file=sys.stderr)
            failed_required.append(run.label)
        else:
            print(f"\n[{run.label}] WARNING: failed (exit {returncode}) - see {run.log_path}; "
                  "the combined submission will publish the hosp ensemble alone for this date.")

    if failed_required:
        print(f"\nERROR: required component(s) failed: {', '.join(failed_required)}", file=sys.stderr)
        sys.exit(1)

    print("\nCreating ensemble forecast...")
    rscript_path = shutil.which("Rscript")
    if rscript_path is None:
        print("ERROR: Rscript not found in PATH", file=sys.stderr)
        print("Please ensure R is installed and in your PATH", file=sys.stderr)
        sys.exit(1)

    print(f"Using Rscript: {rscript_path}")
    subprocess.run(
        [rscript_path, "3_ensemble.R", str(reference_date)],
        check=True, cwd=HERE,
    )

    print("\n" + "=" * 60)
    print(f"Pipeline complete for {reference_date}")
    print("Hosp-only output:  ../../model-output/UMass-gbqr_3src_spatial_wisar/")
    print("Combined submission: ../../model-output/UMass-gbqr_3src_spatial_wisar_nssp/")


if __name__ == "__main__":
    main()
