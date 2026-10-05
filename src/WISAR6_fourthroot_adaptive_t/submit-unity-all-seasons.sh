#!/bin/bash
# Submit the three WISAR6_fourthroot_adaptive_t season arrays as a chain: each season starts once
# the previous one has finished (afterany, so one failed date doesn't block
# the rest). With %10 per array, at most 10 tasks x 4 CPUs = 40 cores are in
# use at a time. Dates that already have output are skipped by wisar6-task.sh.
#
# Run from src/WISAR6_fourthroot_adaptive_t on a login node:  ./submit-unity-all-seasons.sh
# Rough total wall time if 10 tasks run at once: unknown until the first
# Student-t tasks finish (WISAR6_fourthroot: ~1.5 h + 2 h + 3 h).

set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs

prev=""
for season in 2023-24 2024-25 2025-26; do
    dep=()
    [ -n "$prev" ] && dep=(--dependency="afterany:${prev}")
    prev=$(sbatch --parsable "${dep[@]}" "submit-unity-season-${season}.sh")
    echo "Submitted ${season}: job ${prev}${dep:+ (${dep[0]})}"
done

echo
echo "Monitor:  squeue -u \$USER"
echo "Afterwards: sacct -j <jobid> --format=JobID%20,State,Elapsed,TotalCPU,MaxRSS"
