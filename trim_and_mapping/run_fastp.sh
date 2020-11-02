#!/bin/bash
sample=${arg1}
fastp -i ../"$sample"_R1.fq.gz -I ../"$sample"_R2.fq.gz -o "$sample"_Clean_R1.fastq.gz -O "$sample"_Clean_R2.fastq.gz --json "$sample"_fastp.json --html "$sample"_fastp.html --thread 6
echo "completed"
