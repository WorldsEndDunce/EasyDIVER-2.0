#!/bin/bash

outdir=$1
counts_dir="$outdir/counts.aa" # $outdir/counts or $outdir/counts.aa

# Get the maximum round
max_round=0
for file in "$counts_dir"/*-out_counts.aa.txt; do
    filename=$(basename "$file")
    round=${filename%-out_counts.aa.txt}

    if [ $round -gt $max_round ]; then
        max_round=$round
    fi
done

# Check if there are any files matching the format "*-in_counts.aa.txt"
if [ ! -n "$(find "$counts_dir" -name '*-in_counts.aa.txt' -print -quit)" ]; then
    # Cases 1A and 1B: Loop up to max_round - 1
    for ((i = 1; i < max_round; i++)); do
        # Run the modified_counts_bash.sh script with the appropriate arguments
        if [ ! -n "$(find "$counts_dir" -name '*-neg_counts.aa.txt' -print -quit)" ]; then
          python3 ./scripts_enrichments/modified_counts.py -out "$counts_dir/$i-out_counts.aa.txt" -res "./scripts_enrichments/$i-simple_res.txt"
        else
          python3 ./scripts_enrichments/modified_counts.py -out "$counts_dir/$i-out_counts.aa.txt" -neg "$counts_dir/$(($i + 1))-neg_counts.aa.txt" -res "./scripts_enrichments/$i-simple_res.txt"
        fi
    done
else
    # Case 2A and 2B: Loop up to max_round
    for ((i = 1; i <= max_round; i++)); do
        # Run the modified_counts_bash.sh script with the appropriate arguments
        if [ ! -n "$(find "$counts_dir" -name '*-neg_counts.aa.txt' -print -quit)" ]; then
          python3 ./scripts_enrichments/modified_counts.py -in "$counts_dir/$i-in_counts.aa.txt" -out "$counts_dir/$i-out_counts.aa.txt" -res "./scripts_enrichments/$i-simple_res.txt"
        else
          python3 ./scripts_enrichments/modified_counts.py -in "$counts_dir/$i-in_counts.aa.txt" -out "$counts_dir/$i-out_counts.aa.txt" -neg "$counts_dir/$i-neg_counts.aa.txt" -res "./scripts_enrichments/$i-simple_res.txt"
        fi
    done
fi

