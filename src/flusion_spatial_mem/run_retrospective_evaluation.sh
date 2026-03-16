#!/bin/bash
#
# Retrospective Evaluation for flusion_spatial_mem Model
#
# This script generates forecasts for all Saturdays from Oct 2023 to May 2025,
# covering the last two complete flu seasons. This allows for performance
# comparison with flusion_spatial2_prod and other models.
#
# Usage:
#   ./run_retrospective_evaluation.sh              # Full evaluation
#   ./run_retrospective_evaluation.sh --short-run  # Short run for testing
#

set -e  # Exit on error

# Activate virtual environment
if [ ! -d ".venv" ]; then
    echo "ERROR: Virtual environment not found. Please run setup first."
    exit 1
fi

source .venv/bin/activate

# Date range for FluSight forecast season (not year-round)
# 2023-2024 season was: 2023-10-21 to 2024-05-04
# 2024-2025 season: 2024-11-30 to 2025-05-17
START_DATE="2024-11-30"  # Start of 2024-2025 season
END_DATE="2025-05-17"    # End of 2024-2025 season

# Parse command line arguments
SHORT_RUN_FLAG=""
if [ "$1" == "--short-run" ]; then
    SHORT_RUN_FLAG="--short_run"
    echo "Running in SHORT RUN mode (faster, for testing)"
fi

# Run the evaluation
echo "========================================"
echo "Flusion Spatial MEM Retrospective Evaluation"
echo "========================================"
echo "Date range: $START_DATE to $END_DATE"
echo "Output: ../../model-output/UMass-flusion_spatial_mem/"
echo ""
echo "Note: MEM thresholds will be calculated adaptively for each"
echo "forecast date using only data available at that time."
echo "========================================"
echo ""

python evaluate_mem_ensemble.py \
    --start_date="$START_DATE" \
    --end_date="$END_DATE" \
    $SHORT_RUN_FLAG

echo ""
echo "Retrospective evaluation complete!"
echo "Check output directory: ../../model-output/UMass-flusion_spatial_mem/"
echo ""
echo "Next steps:"
echo "  1. Review diagnostic files to see MEM phase evolution"
echo "  2. Compare performance against flusion_spatial2_prod"
echo "  3. Analyze adaptive weight behavior across epidemic phases"
