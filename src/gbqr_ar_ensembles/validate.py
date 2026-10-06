"""Validate GBQR x AR ensemble model-output files.

Checks every file in model-output/<ensemble>/: same keys as its components,
value = mean of the two components (state GBQR + AR for states, US GBQR + AR
for US), no NA/negative values, monotone quantiles, and reference dates /
locations / quantile levels allowed by hub-config/tasks.json.

Usage (from this directory):
  python3 validate.py "<ensemble_id> <gbqr_id> <ar_id> [us_gbqr_id]" ...
"""
import csv
import json
import os
import sys
from collections import defaultdict

root = "../.."
mo = f"{root}/model-output"
tasks = json.load(open(f"{root}/hub-config/tasks.json"))
allowed = defaultdict(set)
req_q = set()
for rnd in tasks["rounds"]:
    for mt in rnd["model_tasks"]:
        for k in ("location", "reference_date"):
            v = mt["task_ids"][k]
            allowed[k] |= {str(x) for x in (v.get("required") or []) + (v.get("optional") or [])}
        if "quantile" not in mt["output_type"]:
            continue
        q = mt["output_type"]["quantile"]["output_type_id"]
        req_q |= {float(x) for x in (q.get("required") or []) + (q.get("optional") or [])}


def load(path):
    return {(r["location"], int(r["horizon"]), float(r["output_type_id"])): (float(r["value"]), r)
            for r in csv.DictReader(open(path))}


ok = True
for spec in sys.argv[1:]:
    ens, gbqr, ar, *rest = spec.split()
    us_gbqr = rest[0] if rest else gbqr
    files = sorted(os.listdir(f"{mo}/{ens}"))
    maxdiff, issues, nrows = 0.0, [], set()
    for f in files:
        d = f[:10]
        E = load(f"{mo}/{ens}/{f}")
        nrows.add(len(E))
        G, A = load(f"{mo}/{gbqr}/{d}-{gbqr}.csv"), load(f"{mo}/{ar}/{d}-{ar}.csv")
        U = load(f"{mo}/{us_gbqr}/{d}-{us_gbqr}.csv") if us_gbqr != gbqr else G
        expected_keys = {k for k in A if (k in U if k[0] == "US" else k in G)}
        if set(E) != expected_keys:
            issues.append(f"{d}: key mismatch")
        for k, (v, r) in E.items():
            g = U if k[0] == "US" else G
            if k in g and k in A:
                maxdiff = max(maxdiff, abs(v - (g[k][0] + A[k][0]) / 2))
            if v != v or v < 0:
                issues.append(f"{d} {k} bad value {v}")
            if r["reference_date"] != d or r["output_type"] != "quantile":
                issues.append(f"{d} {k} bad row")
            if r["location"] not in allowed["location"]:
                issues.append(f"{d} bad location {r['location']}")
        if d not in allowed["reference_date"]:
            issues.append(f"{d} not a hub reference_date")
        if {k[2] for k in E} != req_q:
            issues.append(f"{d}: quantile levels differ from tasks.json")
        by_loc_h = defaultdict(list)
        for (loc, h, q), (v, _) in E.items():
            by_loc_h[(loc, h)].append((q, v))
        for key, qv in by_loc_h.items():
            vals = [v for _, v in sorted(qv)]
            if any(b < a - 1e-9 for a, b in zip(vals, vals[1:])):
                issues.append(f"{d} {key} non-monotone")
    print(f"{ens}: {len(files)} files, rows/file {sorted(nrows)}, "
          f"max |ens - mean(components)| = {maxdiff:.2e}, issues = {len(issues)}", issues[:5])
    ok &= not issues
print("ALL OK" if ok else "PROBLEMS FOUND")
sys.exit(0 if ok else 1)
