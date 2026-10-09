#!/usr/bin/env python3
"""
Summarize a Slurm account's past jobs to answer: "what's the max total number of CPU cores, and memory per CPU core,
you anticipate?"

Reads job history with sacct (standard library only; no extra modules needed) and reports:
  1. Concurrent CPU cores in use by the account: the all-time peak, plus how much of the time usage was at or below
     various levels (a short burst can set the peak).
  2. Cores per job.
  3. Memory per core, both requested (allocated memory / allocated cores) and actually used (peak RSS / cores).

Usage (on a Unity login node):
  python3 slurm_usage_summary.py                       # your PI account(s) found with sacctmgr, last 365 days
  python3 slurm_usage_summary.py --account pi_nick_umass_edu --start 2025-07-01
  python3 slurm_usage_summary.py --account pi_x --exclude-gpu --csv jobs.csv

To see all users' jobs in an account you must be the account's PI/coordinator (otherwise sacct shows only your own).
"""

import argparse
import csv
import datetime as dt
import getpass
import re
import subprocess
import sys
from collections import defaultdict

FIELDS = ["JobID", "User", "Account", "Partition", "State", "Start", "End", "ElapsedRaw", "AllocCPUS", "AllocTRES",
          "MaxRSS"]
UNITS = {"K": 1 / 1024**2, "M": 1 / 1024, "G": 1.0, "T": 1024.0, "P": 1024.0**2}


def to_gb(s: str) -> float:
    """'16G', '512000K', '1.5T' -> GB; '' -> 0."""
    m = re.match(r"^([\d.]+)([KMGTP]?)", s or "")
    if not m:
        return 0.0
    return float(m.group(1)) * UNITS.get(m.group(2) or "M", 1 / 1024)  # sacct's default unit is MB


def tres(s: str) -> dict:
    out = {}
    for part in (s or "").split(","):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def parse_time(s: str):
    try:
        return dt.datetime.fromisoformat(s)
    except ValueError:
        return None  # Unknown / None


def pct(values, p):
    """p-th percentile (0-100) of a list, nearest-rank."""
    if not values:
        return float("nan")
    v = sorted(values)
    return v[min(len(v) - 1, max(0, int(round(p / 100 * len(v) + 0.5)) - 1))]


def default_accounts():
    user = getpass.getuser()
    res = subprocess.run(["sacctmgr", "-nP", "show", "assoc", f"user={user}", "format=account"],
                         capture_output=True, text=True)
    accts = sorted({a.strip() for a in res.stdout.split() if a.strip()})
    pi = [a for a in accts if a.startswith("pi_")]
    return pi or accts


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--account", help="Slurm account(s), comma-separated (default: your pi_* accounts)")
    ap.add_argument("--start", default=(dt.date.today() - dt.timedelta(days=365)).isoformat(),
                    help="earliest job start date, YYYY-MM-DD (default: 365 days ago)")
    ap.add_argument("--end", default="now", help="latest date (default: now)")
    ap.add_argument("--exclude-gpu", action="store_true", help="ignore jobs that were allocated GPUs")
    ap.add_argument("--csv", help="also write one row per job to this csv file")
    args = ap.parse_args()

    accounts = args.account.split(",") if args.account else default_accounts()
    if not accounts:
        sys.exit("No Slurm account found; pass --account")
    cmd = ["sacct", "--allusers", f"--accounts={','.join(accounts)}", f"--starttime={args.start}",
           f"--endtime={args.end}", "--parsable2", "--noheader", f"--format={','.join(FIELDS)}"]
    print("Running:", " ".join(cmd), file=sys.stderr)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit(res.stderr)

    # one allocation line per job, followed by its steps (JobID.batch, JobID.0, ...), which carry MaxRSS
    jobs, step_rss = {}, defaultdict(float)
    for line in res.stdout.splitlines():
        r = dict(zip(FIELDS, line.split("|")))
        jid = r["JobID"]
        if "." in jid:
            step_rss[jid.split(".")[0]] = max(step_rss[jid.split(".")[0]], to_gb(r["MaxRSS"]))
            continue
        start = parse_time(r["Start"])
        if start is None or int(r["ElapsedRaw"] or 0) == 0:
            continue  # never ran
        end = parse_time(r["End"]) or dt.datetime.now()
        t = tres(r["AllocTRES"])
        cpus = int(r["AllocCPUS"] or 0)
        gpus = int(t.get("gres/gpu", 0) or 0)
        if cpus == 0 or (args.exclude_gpu and gpus > 0):
            continue
        jobs[jid] = dict(user=r["User"], account=r["Account"], partition=r["Partition"], state=r["State"].split()[0],
                         start=start, end=end, cpus=cpus, gpus=gpus, mem_gb=to_gb(t.get("mem", "")),
                         hours=int(r["ElapsedRaw"]) / 3600)
    for jid, j in jobs.items():
        j["maxrss_gb"] = step_rss.get(jid, 0.0)
        j["mem_per_cpu_req"] = j["mem_gb"] / j["cpus"]
        j["mem_per_cpu_used"] = j["maxrss_gb"] / j["cpus"] if j["maxrss_gb"] else None

    if not jobs:
        sys.exit("No jobs found for these accounts and dates.")
    js = list(jobs.values())
    print(f"\nAccounts: {', '.join(accounts)}   period: {args.start} to {args.end}")
    print(f"Jobs that ran: {len(js)}   users: {len({j['user'] for j in js})}   "
          f"core-hours: {sum(j['cpus'] * j['hours'] for j in js):,.0f}"
          + ("   (GPU jobs excluded)" if args.exclude_gpu else ""))

    # 1. concurrent cores: sweep over start/end events, time-weighted
    events = sorted([(j["start"], j["cpus"]) for j in js] + [(j["end"], -j["cpus"]) for j in js])
    level, peak, peak_at, last, dur = 0, 0, None, None, defaultdict(float)
    for t, d in events:
        if last is not None and t > last:
            dur[level] += (t - last).total_seconds()
        level += d
        last = t
        if level > peak:
            peak, peak_at = level, t
    total = sum(dur.values())
    busy = {k: v for k, v in dur.items() if k > 0}
    busy_total = sum(busy.values())

    def level_at_fraction(frac, d, tot):
        acc = 0.0
        for k in sorted(d):
            acc += d[k]
            if acc >= frac * tot:
                return k
        return max(d)

    print("\n1. Total CPU cores in use at the same time (all jobs in the account)")
    print(f"   peak: {peak} cores (at {peak_at:%Y-%m-%d %H:%M})")
    print("   share of time at or below:  " + "  ".join(
        f"{int(f * 100)}%: {level_at_fraction(f, dur, total)}" for f in (0.5, 0.9, 0.95, 0.99)))
    if busy_total:
        print("   same, counting only time when something was running:  " + "  ".join(
            f"{int(f * 100)}%: {level_at_fraction(f, busy, busy_total)}" for f in (0.5, 0.9, 0.95, 0.99)))
        print(f"   average while running: {sum(k * v for k, v in busy.items()) / busy_total:.0f} cores; "
              f"something was running {100 * busy_total / total:.0f}% of the period")

    # 2. cores per job
    cpus = [j["cpus"] for j in js]
    print("\n2. Cores per job")
    print(f"   median {pct(cpus, 50)}, 95th percentile {pct(cpus, 95)}, max {max(cpus)}")

    # 3. memory per core
    req = [j["mem_per_cpu_req"] for j in js if j["mem_gb"] > 0]
    used = [j["mem_per_cpu_used"] for j in js if j["mem_per_cpu_used"] is not None]
    print("\n3. Memory per core (GB)")
    if req:
        print(f"   requested (allocated mem / cores): median {pct(req, 50):.1f}, 95th pct {pct(req, 95):.1f}, "
              f"max {max(req):.1f}")
    if used:
        print(f"   actually used (peak RSS / cores):  median {pct(used, 50):.2f}, 95th pct {pct(used, 95):.2f}, "
              f"max {max(used):.2f}   ({len(used)} jobs with memory records)")
    peak_mem_jobs = sorted(js, key=lambda j: -j["maxrss_gb"])[:3]
    print("   largest jobs by peak memory: " + "; ".join(
        f"{j['maxrss_gb']:.0f} GB on {j['cpus']} cores ({j['user']}, {j['partition']})" for j in peak_mem_jobs))

    # by user, to see who drives the peak
    print("\nBy user: jobs, core-hours, max cores in one job, max GB/core used")
    by_user = defaultdict(list)
    for j in js:
        by_user[j["user"]].append(j)
    for u, uj in sorted(by_user.items(), key=lambda kv: -sum(j["cpus"] * j["hours"] for j in kv[1])):
        mu = [j["mem_per_cpu_used"] for j in uj if j["mem_per_cpu_used"] is not None]
        print(f"   {u:<20} {len(uj):>6} {sum(j['cpus'] * j['hours'] for j in uj):>12,.0f} "
              f"{max(j['cpus'] for j in uj):>6} {max(mu) if mu else float('nan'):>8.2f}")

    if args.csv:
        with open(args.csv, "w", newline="") as f:
            w = csv.writer(f)
            cols = ["user", "account", "partition", "state", "start", "end", "hours", "cpus", "gpus", "mem_gb",
                    "maxrss_gb", "mem_per_cpu_req", "mem_per_cpu_used"]
            w.writerow(["jobid"] + cols)
            for jid, j in jobs.items():
                w.writerow([jid] + [j[c] for c in cols])
        print(f"\nWrote {args.csv}")


if __name__ == "__main__":
    main()
