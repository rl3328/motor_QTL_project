#!/bin/bash

# Script to check for regions with stderr files but missing stdout or rds files
# Usage: ./check_missing_files.sh /path/to/folder

# Check if folder path is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <folder_path>"
    echo "Example: $0 /data/regions"
    exit 1
fi

FOLDER="$1"

# Check if folder exists
if [ ! -d "$FOLDER" ]; then
    echo "Error: Folder '$FOLDER' does not exist"
    exit 1
fi

# Output files
MISSING_STDOUT="regions_missing_stdout.txt"
MISSING_RDS="regions_missing_rds.txt"

# Clear output files if they exist
> "$MISSING_STDOUT"
> "$MISSING_RDS"

echo "Scanning folder: $FOLDER"
echo "================================"

# Find all stderr files
stderr_files=$(find "$FOLDER" -name "*.stderr" -type f)

if [ -z "$stderr_files" ]; then
    echo "No stderr files found in $FOLDER"
    exit 0
fi

# Counter variables
total_stderr=0
missing_stdout=0
missing_rds=0

# Process each stderr file
while IFS= read -r stderr_file; do
    ((total_stderr++))
    
    # Extract filename and directory
    filename=$(basename "$stderr_file")
    dir_name=$(dirname "$stderr_file")
    
    # Extract region (chr*_*_*) from filename
    region=$(echo "$filename" | grep -oP 'chr\w+_\d+_\d+')
    
    if [ -z "$region" ]; then
        echo "Warning: Could not extract region from $filename"
        continue
    fi
    
    # Extract prefix (everything before the region)
    prefix=$(echo "$filename" | sed "s/\.$region\..*//")
    
    # Construct expected filenames
    stdout_file="$dir_name/${prefix}.${region}.twas.tsv.stdout"
    rds_file="$dir_name/${prefix}.${region}.twas_data.rds"
    
    # Check for corresponding stdout file
    if [ ! -f "$stdout_file" ]; then
        echo "$region" >> "$MISSING_STDOUT"
        ((missing_stdout++))
    fi
    
    # Check for corresponding rds file
    if [ ! -f "$rds_file" ]; then
        echo "$region" >> "$MISSING_RDS"
        ((missing_rds++))
    fi
    
done <<< "$stderr_files"

# Print summary
echo ""
echo "Summary:"
echo "================================"
echo "Total stderr files found: $total_stderr"
echo "Regions with stderr but no stdout: $missing_stdout"
echo "Regions with stderr but no rds: $missing_rds"
echo ""
echo "Results saved to:"
echo "  - $MISSING_STDOUT"
echo "  - $MISSING_RDS"

# Show first few entries if available
if [ $missing_stdout -gt 0 ]; then
    echo ""
    echo "First 5 regions missing stdout:"
    head -5 "$MISSING_STDOUT"
fi

if [ $missing_rds -gt 0 ]; then
    echo ""
    echo "First 5 regions missing rds:"
    head -5 "$MISSING_RDS"
fi