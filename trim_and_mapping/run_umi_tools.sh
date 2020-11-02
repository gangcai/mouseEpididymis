#! /bin/env bash
# Step 1: Prepare data
sample_name=$1
fq_dir="../fq/"
fq1="$fq_dir"/"$sample_name"_Clean_R1.fastq.gz # Barcode (12) and UMI at 5p' (8)
fq2="$fq_dir"/"$sample_name"_Clean_R2.fastq.gz # cDNA
zcat $fq1 | head -2
pattern="CCCCCCCCCCCCNNNNNNNN"
genome_dir="/home/db/public/annotation/SoftwareIndex/STAR/Mouse/genomeDir/"
gtf="/home/db/public/annotation/Gencode/mouse/gencode.vM25.primary_assembly.annotation.gtf"

# Step 2: Identify correct cell barcodes
umi_tools whitelist --stdin $fq1 \
                    --bc-pattern=$pattern \
                    --set-cell-number=5000 \
		    --error-correct-threshold=2 \
                    --log2stderr > "$sample_name"_whitelist.txt;

# Step 3: Extract barcdoes and UMIs and add to read names
umi_tools extract --bc-pattern=$pattern \
                  --stdin $fq1 \
                  --stdout "$sample_name"_1_extracted.fastq.gz \
                  --read2-in $fq2 \
                  --read2-out="$sample_name"_2_extracted.fastq.gz \
                  --filter-cell-barcode \
		  --error-correct-cell \
                  --whitelist="$sample_name"_whitelist.txt; 
# Step 4: Map reads
STAR --runThreadN 15 \
     --genomeDir $genome_dir \
     --readFilesIn "$sample_name"_2_extracted.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 1 \
     --outSAMtype BAM SortedByCoordinate \
     --sjdbGTFfile $gtf \
     --outFileNamePrefix $sample_name;

     
# Step 5: Assign reads to genes
featureCounts -F GTF \
	      -t exon \
	      -g gene_name \
	      -a $gtf \
              -o "$sample_name"_gene_counts.txt \
              -R BAM "$sample_name"Aligned.sortedByCoord.out.bam \
              -T 15;            
samtools sort "$sample_name"Aligned.sortedByCoord.out.bam.featureCounts.bam -o "$sample_name"_assigned_sorted.bam;
samtools index "$sample_name"_assigned_sorted.bam;
              
# Step 6: Count UMIs per gene per cell
umi_tools count --per-gene --gene-tag=XT --edit-distance-threshold 1 --umi-separator "_" --per-cell -I "$sample_name"_assigned_sorted.bam -S "$sample_name"_counts.tsv.gz
