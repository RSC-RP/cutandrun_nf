#!/bin/bash

# Create toy dataset
BAMS=$(ls -1 results/bowtie2/*.bam )
REGION="chr17:39957839-40353938"

for bam in $(echo "$BAMS")
do
    bam_file=$(basename $bam)
    echo $bam
    samtools index $bam
    samtools view -bh $bam $REGION > test_data/bams/${bam_file%.bam}_chr17.bam
done

for bam in $(ls -1 test_data/bams/*.bam)
do  
    dir="test_data/fastqs"
    file_name=$(basename $bam)
    file_name="$dir/${file_name%.bam}"
    echo $file_name
    samtools collate -u -O $bam | samtools fastq -1 ${file_name}_R1.fastq -2 ${file_name}_R2.fastq -0 /dev/null -s /dev/null -n
done

gzip test_data/fastqs/*.fastq