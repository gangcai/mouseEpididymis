#!/bin/bash
prefix="gene"
sort -k 1,1 "$prefix"_combined.tmp.tsv  > "$prefix"_combined_sorted.tmp.tsv
