#! /bin/bash -l

##########################################
# PART 1: Run for each sample separately #
##########################################

# StringTie v2.1.1: assembles RNA-seq aligned data into transcripts w/ gene abundance & isoform info
# Script command line usage: sbatch StringTie_Server.sh SAMPLE_1

InputDir=bam
OutputDir=counts
GFF3=genomes/GRCh38/starIndex_gencode24/gencode.v24.annotation.gff3

# Run StringTie to assemble RNA-seq sample transcript info w/ ref genome (gene ID, gene name, gene length, # exons, coverage, FPKM, TPM)
# -G= reference genome
# -A= create separate file w/ gene abundance data
mkdir ${OutputDir}/${1}
stringtie ${InputDir}/${1}.sorted.bam -o ${OutputDir}/${1}/StringTie_transcripts.gtf -G ${GFF3} -A ${OutputDir}/${1}/StringTie_gene_abund.tab

echo "End of Part 1!"