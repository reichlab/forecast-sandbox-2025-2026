#!/bin/bash
# Convenience script to run complete analysis for both dates

set -e  # Exit on error

echo "=========================================="
echo "RUNNING FULL ANALYSIS"
echo "=========================================="

# Check if virtual environment exists
if [ ! -d ".venv" ]; then
    echo "Creating virtual environment..."
    python -m venv .venv
    source .venv/bin/activate
    echo "Installing dependencies..."
    pip install pandas numpy matplotlib seaborn click
else
    source .venv/bin/activate
fi

# Dates to analyze
DATES=("2024-12-21" "2025-02-15")

for date in "${DATES[@]}"; do
    echo ""
    echo "=========================================="
    echo "ANALYZING DATE: $date"
    echo "=========================================="

    echo ""
    echo "1. Feature Importance Analysis..."
    python analyze_feature_importance.py --date "$date"

    echo ""
    echo "2. Geographic Analysis..."
    python analyze_geographic_differences.py --date "$date"

    echo ""
    echo "Completed analysis for $date"
    echo "Results saved to: outputs/$date/"
done

echo ""
echo "=========================================="
echo "ALL ANALYSES COMPLETE"
echo "=========================================="
echo ""
echo "Review results in:"
for date in "${DATES[@]}"; do
    echo "  - analysis/outputs/$date/"
done
echo ""
