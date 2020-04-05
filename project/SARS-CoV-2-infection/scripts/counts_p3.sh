#! /bin/bash

###########################################################################
# PART 3: Re-running StringTie for each sample using the output of Part 2 #
###########################################################################

# StringTie v2.1.1: assembles RNA-seq aligned data into transcripts w/ gene abundance & isoform info
# Script command line usage: sbatch counts_p2.sh SAMPLE

OutputDir=counts
InputDir=bam

mkdir ${OutputDir}/final/${1}

# Rerun stringtie to also estimate transcript abundances (-eB) & generate read coverage tables
# Use merged GTF file (containing non-redundant transcripts across all samples) as ref genome (-G)
stringtie ${InputDir}/${1}.sorted.bam -o ${OutputDir}/final/${1}/StringTie_transcripts_filtered.gtf -eB -G ${OutputDir}/stringtie_merged_transcripts.gtf -A ${OutputDir}/final/${1}/StringTie_gene_abund_filtered.tab

echo "Part 3 is done!"