#!/bin/bash
# Script to find available R modules on Unity

echo "Searching for R modules on Unity..."
echo ""
echo "Method 1: module spider R"
module spider R 2>&1 | head -50
echo ""
echo "=========================================="
echo ""
echo "Method 2: module avail R"
module avail R 2>&1
echo ""
echo "=========================================="
echo ""
echo "Method 3: Looking for common R module names"
for version in R/4.4 R/4.3 R/4.2 R/4.1 R/4 R; do
    echo -n "Trying $version... "
    module load $version 2>/dev/null && echo "✓ FOUND" && module unload $version || echo "not found"
done
