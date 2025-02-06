#!/bin/bash

# Check if arguments are provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <file_pattern_or_files>"
    echo "Example 1: $0 *_sorted_filtered.bam"
    echo "Example 2: $0 file1_sorted_filtered.bam file2_sorted_filtered.bam"
    exit 1
fi

# Iterate over all arguments (can be patterns or specific filenames)
for sorted_bam_file in "$@"; do
    # Check if the file exists
    if [ -f "$sorted_bam_file" ]; then
        dedup_bam_file="${sorted_bam_file/_sorted_filtered.bam/_filtered_dedup.bam}"
        dedup_log="${sorted_bam_file/_sorted_filtered.bam/_filtered_dedup.log}"
        echo "Processing: $sorted_bam_file"
        umi_tools dedup --method unique --stdin="$sorted_bam_file" --log="$dedup_log" --paired > "$dedup_bam_file"
        echo "Deduplicated file created: $dedup_bam_file"
    else
        echo "File not found: $sorted_bam_file"
    fi
done
