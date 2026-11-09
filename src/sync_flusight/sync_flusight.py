#!/usr/bin/env python3
"""
Sync the data this sandbox borrows from cdcepi/FluSight-forecast-hub:

  1. target-data/oracle-output.csv  <- FluSight's target-data/oracle-output.csv (copied as-is)
  2. model-output/FluSight-baseline <- FluSight's FluSight-baseline forecasts, filtered to what this hub accepts:
                                       quantile output, horizons 0-3, the two weekly targets, and the hub's locations

Run weekly during the season (after FluSight's Wednesday target-data update / Thursday baseline) from the hub root:

    python3 src/sync_flusight/sync_flusight.py                      # season 2026-2027 (default)
    python3 src/sync_flusight/sync_flusight.py --since 2025-11-22   # any other window of reference dates

Baseline files that already exist locally are left alone unless --overwrite is given. Only reference dates listed in
hub-config/tasks.json and not in the future are fetched. Standard library only.
"""
import argparse
import csv
import datetime as dt
import io
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

RAW = "https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main"
BASELINE = "FluSight-baseline"
WEEKLY_TARGETS = ("wk inc flu hosp", "wk inc flu prop ed visits")
HORIZONS = {"0", "1", "2", "3"}


def fetch(url: str) -> bytes | None:
    try:
        with urllib.request.urlopen(url, timeout=120) as resp:
            return resp.read()
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        raise


def hub_reference_dates_and_locations(hub: Path) -> tuple[list[str], set[str]]:
    tasks = json.loads((hub / "hub-config" / "tasks.json").read_text())
    ref_dates, locations = set(), set()
    for mt in tasks["rounds"][0]["model_tasks"]:
        tid = mt["task_ids"]
        if (tid["target"]["required"] or [None])[0] not in WEEKLY_TARGETS:
            continue
        for key, acc in (("reference_date", ref_dates), ("location", locations)):
            acc.update((tid[key]["required"] or []) + (tid[key]["optional"] or []))
    return sorted(ref_dates), locations


def sync_oracle(hub: Path) -> None:
    data = fetch(f"{RAW}/target-data/oracle-output.csv")
    if data is None:
        sys.exit("FluSight oracle-output.csv not found")
    header = data.split(b"\n", 1)[0].decode()
    expected = '"as_of","target","target_end_date","location","horizon","output_type","output_type_id","oracle_value"'
    if header.strip() != expected:
        sys.exit(f"unexpected oracle-output.csv header: {header!r}")
    out = hub / "target-data" / "oracle-output.csv"
    out.write_bytes(data)
    as_of = data.split(b"\n", 2)[1].split(b",", 1)[0].decode()
    n_rows = data.count(b"\n") - 1
    print(f"oracle-output.csv: updated (as_of {as_of}, {n_rows} rows)")


def sync_baseline(hub: Path, since: str, overwrite: bool) -> None:
    ref_dates, locations = hub_reference_dates_and_locations(hub)
    today = dt.date.today().isoformat()
    out_dir = hub / "model-output" / BASELINE
    for ref_date in (d for d in ref_dates if since <= d <= today):
        out = out_dir / f"{ref_date}-{BASELINE}.csv"
        if out.exists() and not overwrite:
            continue
        data = fetch(f"{RAW}/model-output/{BASELINE}/{ref_date}-{BASELINE}.csv")
        if data is None:
            print(f"{out.name}: not yet available upstream")
            continue
        reader = csv.DictReader(io.StringIO(data.decode()))
        rows = [r for r in reader
                if r["target"] in WEEKLY_TARGETS and r["output_type"] == "quantile"
                and r["horizon"] in HORIZONS and r["location"] in locations]
        if not rows:
            print(f"{out.name}: upstream file has no usable rows, skipped")
            continue
        with open(out, "w", newline="") as fp:
            writer = csv.DictWriter(fp, fieldnames=reader.fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        print(f"{out.name}: wrote {len(rows)} rows")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hub", default=".", help="hub root (default: current directory)")
    parser.add_argument("--since", default="2026-10-10", help="first reference date to sync baseline for")
    parser.add_argument("--overwrite", action="store_true", help="re-download baseline files that already exist")
    parser.add_argument("--skip-oracle", action="store_true")
    parser.add_argument("--skip-baseline", action="store_true")
    args = parser.parse_args()
    hub = Path(args.hub).resolve()
    if not (hub / "hub-config" / "tasks.json").exists():
        sys.exit(f"{hub} is not a hub root")
    if not args.skip_oracle:
        sync_oracle(hub)
    if not args.skip_baseline:
        sync_baseline(hub, args.since, args.overwrite)


if __name__ == "__main__":
    main()
