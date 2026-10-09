#!/bin/bash
# Build all GBQR x AR pairwise ensembles from existing model-output files.
# No model fitting happens here; it runs in well under a minute.
# Run from this directory: bash run_all.sh
set -euo pipefail

module load r/4.4.0
# reuse the renv library from flusion_spatial2_prod (has hubEnsembles/hubUtils)
export R_LIBS_USER="$(cd ../flusion_spatial2_prod && pwd)/renv/library/linux-ubuntu-noble/R-4.4/x86_64-pc-linux-gnu"
echo "R: $(which Rscript) ($(Rscript --version 2>&1 | head -1))"

# output_model_id  gbqr_component  ar_component  [us_gbqr_component]
pairs=(
  "UMass-gbqr_4src_ar6p   UMass-gbqr_4src  UMass-AR6_pooled"
  "UMass-gbqr_4src_ar6fp  UMass-gbqr_4src  UMass-AR6_fourierP_thetaP"
  "UMass-gbqr_4src_wisar  UMass-gbqr_4src  UMass-WISAR6_fourthroot_adaptive_t"
  "UMass-gbqr_nih10_ar6p  UMass-gbqr_5src_nihflu_otid10  UMass-AR6_pooled"
  "UMass-gbqr_nih10_ar6fp UMass-gbqr_5src_nihflu_otid10  UMass-AR6_fourierP_thetaP"
  "UMass-gbqr_nih10_wisar UMass-gbqr_5src_nihflu_otid10  UMass-WISAR6_fourthroot_adaptive_t"
  "UMass-gbqr_3src_ar6p   UMass-gbqr_3src  UMass-AR6_pooled"
  "UMass-gbqr_3src_ar6fp  UMass-gbqr_3src  UMass-AR6_fourierP_thetaP"
  "UMass-gbqr_3src_wisar  UMass-gbqr_3src  UMass-WISAR6_fourthroot_adaptive_t"
)

for p in "${pairs[@]}"; do
  # shellcheck disable=SC2086
  Rscript ensemble_pair.R $p
done
