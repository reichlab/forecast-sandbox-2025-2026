#!/bin/bash
# Submit the three leak-fix rerun arrays (rerun-leakfix-season-*.sh) as a
# chain: each season starts once the previous one has finished (afterany, so
# one failed date doesn't block the rest). With %10 per array, at most
# 10 tasks x 4 CPUs = 40 cores are in use at a time.
#
# Run from src/WISAR6_fourthroot_adaptive_t on a login node:
#   ./submit-rerun-leakfix.sh                  # all three seasons
#   ./submit-rerun-leakfix.sh 2024-25 2025-26  # just these seasons
#
# Smoke test first (one date, ~1 h):
#   sbatch --array=0 rerun-leakfix-season-2023-24.sh

set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs

seasons=("$@")
[ ${#seasons[@]} -eq 0 ] && seasons=(2023-24 2024-25 2025-26)

prev=""
for season in "${seasons[@]}"; do
    dep=()
    [ -n "$prev" ] && dep=(--dependency="afterany:${prev}")
    prev=$(sbatch --parsable "${dep[@]}" "rerun-leakfix-season-${season}.sh")
    echo "Submitted ${season}: job ${prev}${dep:+ (${dep[0]})}"
done

echo
echo "Monitor:  squeue -u \$USER"
echo "Afterwards: sacct -j <jobid> --format=JobID%20,State,Elapsed,TotalCPU,MaxRSS"
