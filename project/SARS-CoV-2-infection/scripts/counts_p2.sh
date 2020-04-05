#! /bin/bash

###############################################################
# PART 2: Run after Part 1 has been completed for all samples #
###############################################################

# StringTie v2.1.1: assembles RNA-seq aligned data into transcripts w/ gene abundance & isoform info
# Script command line usage: sbatch counts_p2.sh

OutputDir=counts
GFF3=genomes/GRCh38/starIndex_gencode24/gencode.v24.annotation.gff3

# Run stringtie --merge to create a global, unified transcriptome (isoforms) across sample data
# I.e. generates non-redundant set of transcripts
# Create .txt file w/ directories to StringTie GTF results for each sample that you want to merge
cd ${OutputDir}
ls */StringTie_transcripts.gtf > assembly_GTF_list.txt
cat assembly_GTF_list.txt
stringtie --merge -o stringtie_merged_transcripts.gtf -G ${GFF3} assembly_GTF_list.txt

# Analyze accuracy of StringTie via GFFCompare (i.e. assessing quality of StringTie's transcript abundance predictions) 
# Compare transcripts generated from ref genome to ref genome annotations
gffcompare -r ${GFF3} -o gffcompare stringtie_merged_transcripts.gtf
cat gffcompare.stats

mkdir final

echo "Part 2 is done!"
