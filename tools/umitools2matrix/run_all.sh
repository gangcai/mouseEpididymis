#!/bin/bash
python step1_combine.py
bash step2_sort.sh
python step3_convert2matrix.py
rm gene_combined.tmp.tsv 
rm gene_combined_sorted.tmp.tsv 
echo "completed"
