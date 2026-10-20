#!/bin/bash
# Run the SARIX nowcast comparison models (AR6_pooled and AR2, 3 variants each) for the 28 dates of 2025/26 locally,
# N_PARALLEL at a time. Skips forecasts that already exist. Run from this directory.
N_PARALLEL=${N_PARALLEL:-4}
dates=$(sed -n '/^dates=(/,/^)/p' submit-unity.sh | grep -o '20[0-9-]*')
models="AR6_pooled_nc_none AR6_pooled_nc_scale AR6_pooled_nc_state AR2_nc_none AR2_nc_scale AR2_nc_state"
for m in $models; do for d in $dates; do echo "$m $d"; done; done |
  xargs -P "$N_PARALLEL" -L 1 bash -c '
    m=$0; d=$1
    ref=$(date -j -v+3d -f %F "$d" +%F)
    out="../../model-output/UMass-$m/$ref-UMass-$m.csv"
    [ -f "$out" ] && { echo "skip $m $d"; exit 0; }
    # SIP strips DYLD_* from protected binaries such as bash, so set it here, right before python
    export DYLD_FALLBACK_LIBRARY_PATH=$PWD/.venv/lib/python3.12/site-packages/sklearn/.dylibs
    export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
    export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=2"
    start=$(date +%s)
    .venv/bin/python main.py --model "$m" --today_date "$d" --archive_dir archives > "logs/local-$m-$d.log" 2>&1
    echo "$m $d exit $? $(( $(date +%s) - start ))s"
  '
echo "all done"
