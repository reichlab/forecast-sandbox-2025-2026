#!/bin/bash

# Directory containing your model-output CSV files
CSV_DIR="../../model-output/UMass-flusion"

# Loop through all CSV files in the directory
for file in "$CSV_DIR"/*.csv; do
    # Remove all double quotes and overwrite the file
    sed -i '' 's/"//g' "$file"
    echo "Cleaned: $file"
done
