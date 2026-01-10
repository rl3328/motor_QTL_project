#!/bin/bash

# Script to find regions missing RDS files
# Usage: ./find_missing_rds.sh [output_folder] [bed_file] [output_file]

# Set default values
OUTPUT_FOLDER=${1:-"output"}
BED_FILE=${2:-"/mnt/lustre/lab/gwang/home/rl3328/gpQTL_project/TWAS/hg38_1362_blocks.bed"}
OUTPUT_FILE=${3:-"missing_regions.tsv"}

# Check if output folder exists
if [ ! -d "$OUTPUT_FOLDER" ]; then
    echo "Error: Output folder '$OUTPUT_FOLDER' does not exist"
    exit 1
fi

# Check if BED file exists
if [ ! -f "$BED_FILE" ]; then
    echo "Error: BED file '$BED_FILE' does not exist"
    exit 1
fi

# Create temporary files
TEMP_RDS_REGIONS=$(mktemp)
TEMP_ALL_REGIONS=$(mktemp)

# Extract region names from RDS files
echo "Extracting regions from RDS files in $OUTPUT_FOLDER..."
find "$OUTPUT_FOLDER" -name "*.rds" -type f | \
    grep -o 'chr[^_]*_[0-9]*_[0-9]*' | \
    sort -u > "$TEMP_RDS_REGIONS"

# Extract all region names from BED file (skip header if it starts with #)
echo "Extracting all regions from BED file..."
awk 'BEGIN{OFS="\t"} !/^#/ {print $4}' "$BED_FILE" | \
    sort -u > "$TEMP_ALL_REGIONS"

# Find regions that are in BED file but not in RDS files
echo "Finding missing regions..."
comm -23 "$TEMP_ALL_REGIONS" "$TEMP_RDS_REGIONS" > "$OUTPUT_FILE"

# Count results
RDS_COUNT=$(wc -l < "$TEMP_RDS_REGIONS")
TOTAL_COUNT=$(wc -l < "$TEMP_ALL_REGIONS")
MISSING_COUNT=$(wc -l < "$OUTPUT_FILE")

echo "Results:"
echo "  Total regions in BED file: $TOTAL_COUNT"
echo "  Regions with RDS files: $RDS_COUNT"
echo "  Missing RDS regions: $MISSING_COUNT"
echo "  Output written to: $OUTPUT_FILE"

# Clean up temporary files
rm "$TEMP_RDS_REGIONS" "$TEMP_ALL_REGIONS"

echo "Done!"