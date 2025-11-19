#!/bin/bash
# Test script to verify Unity environment is set up correctly
# Run this before submitting the SLURM job

echo "=========================================="
echo "Testing Unity Environment Setup"
echo "=========================================="

# Check if we're in the right directory
if [ ! -f "main.py" ]; then
    echo "ERROR: Run this script from src/flusion_spatial directory"
    exit 1
fi

# Test 1: Check virtual environment
echo ""
echo "[1/5] Checking virtual environment..."
if [ -d ".venv" ]; then
    echo "✓ Virtual environment found"
    source .venv/bin/activate
    echo "✓ Virtual environment activated"
else
    echo "✗ Virtual environment not found"
    echo "  Run: python -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt"
    exit 1
fi

# Test 2: Check Python packages
echo ""
echo "[2/5] Checking Python packages..."
python -c "import click; import idmodels" 2>/dev/null
if [ $? -eq 0 ]; then
    echo "✓ Required Python packages installed"
else
    echo "✗ Missing Python packages"
    echo "  Run: pip install -r requirements.txt"
    exit 1
fi

# Test 3: Load R module and check
echo ""
echo "[3/5] Checking R availability..."
module load R/4.4.1 2>/dev/null
if [ $? -ne 0 ]; then
    echo "⚠ Could not load R module (might work differently on Unity)"
fi

which Rscript >/dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✓ Rscript found: $(which Rscript)"
    echo "  Version: $(Rscript --version 2>&1 | head -1)"
else
    echo "✗ Rscript not found in PATH"
    echo "  Try: module load R/4.4.1"
    exit 1
fi

# Test 4: Check R packages
echo ""
echo "[4/5] Checking R packages..."
Rscript -e "library(dplyr); library(hubData); library(hubEnsembles)" 2>/dev/null
if [ $? -eq 0 ]; then
    echo "✓ Required R packages installed"
else
    echo "✗ Missing R packages"
    echo "  Required packages: dplyr, hubData, hubEnsembles, readr, idforecastutils"
    exit 1
fi

# Test 5: Test Python can find Rscript
echo ""
echo "[5/5] Testing Python subprocess can find Rscript..."
python -c "import shutil; rscript = shutil.which('Rscript'); print(f'✓ Python can find Rscript: {rscript}') if rscript else exit(1)"
if [ $? -ne 0 ]; then
    echo "✗ Python cannot find Rscript"
    echo "  This is the issue! R module environment not visible to Python"
    exit 1
fi

# Test 6: Quick model test
echo ""
echo "=========================================="
echo "Environment OK! Ready to run jobs."
echo "=========================================="
echo ""
echo "To test a single forecast before submitting all jobs:"
echo "  python main.py --today_date=2024-12-28 --short_run"
echo ""
echo "To submit all jobs:"
echo "  sbatch submit-unity-parallel.sh"
echo ""
