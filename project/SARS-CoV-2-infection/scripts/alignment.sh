#! /bin/bash

# STAR Aligner (v. 2.5.2b): mapping raw sequencing data (fastq.gz) to reference genome
# Input: fastq.gz file
# Output: BAM file
# Script command line usage: sh alignment.sh SAMPLE

InputDir=#Replace with input directory
OutputDir=#Replace with output directory
GenomeDir=genomes/GRCh38/starIndex_gencode24

# Run STAR on single end fastq.gz RNA-seq samples
STAR --readFilesIn ${InputDir}/${1}.fastq.gz --readFilesCommand zcat --genomeDir ${GenomeDir} --outFilterMultimapNmax 1 --outFilterIntronMotifs RemoveNoncanonical --outFilterMismatchNmax 5 --alignSJDBoverhangMin 6 --alignSJoverhangMin 6 --outFilterType BySJout --alignIntronMin 25 --alignIntronMax 1000000 --outSAMstrandField intronMotif --outSAMunmapped Within --runThreadN 24 --outStd SAM --alignMatesGapMax 1000000 | samtools sort -@4 -O bam -o ${OutputDir}/${1}.sorted.bam
samtools index ${OutputDir}/${1}.sorted.bam ${OutputDir}/${1}.sorted.bai