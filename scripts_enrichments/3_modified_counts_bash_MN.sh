#!/bin/bash

python ../3_modified_counts_MN.py $1 $2 $3 $4 $5 $6

head -n3 $4 > final.txt
(tail -n+4 $4 | sort -k9 -n -r )  >> final.txt
mv final.txt $4;

