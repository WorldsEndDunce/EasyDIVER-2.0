#!/bin/bash

python ./scripts_enrichments/modified_counts.py $1 $2 $3 $4 $5 $6 $7 $8 # .. when running easydiver

#head -n3 $4 > final.txt
#(tail -n+5 $4 | sort -k4 -n -r )  >> final.txt
#mv final.txt $4;

# Sort
#header=$(head -n3 "$4")
#echo "$header" > sorted.txt
#tail -n+4 "$4" | sort -k10 -nr -k1,1 >> sorted.txt
