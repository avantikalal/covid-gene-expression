##############################################

### Modify this section for your local environment!!### 

# Set working directory
setwd("proj_2_rbp_data")

# List input data files

## NCBI # Downloaded from https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.
ref_file="reference/sequence.fasta.txt"
gff_file = "ncbi/ref/GCF_009858895.2_ASM985889v3_genomic.gff"

## ATtRACT # Downloaded from https://attract.cnic.es/download
rbp_file = "ATtRACT/ATtRACT_db.txt"
pwm_file = "ATtRACT/pwm.txt"

################################################

# Import requirements
library(data.table)
library(seqinr)
library(Biostrings)
library(TFBSTools)
library(foreach)
library(doParallel)
library(rtracklayer)
source("covid-gene-expression/project/human-RBP-analysis/rbp_functions.R")

registerDoParallel(detectCores())

################################################

# Read input data

# Reference genome
ref = readFasta(ref_file)

# RBPs
rbp = fread(rbp_file)

# PWMs
pwm = readPWMsFromFasta(pwm_file)

# GFF
gff = readGFF("ncbi/ref/GCF_009858895.2_ASM985889v3_genomic.gff")

##############################################

# Select only human RBPs
rbp = rbp[Organism == "Homo_sapiens",]

# See experiment types
rbp[, .N, by=Experiment_description][order(N, decreasing = T),]

# See source databases
rbp[, .N, by=Database][order(N, decreasing = T),]

# How many unique proteins are present?
rbp[, length(unique(Gene_name))]

# Select unique RBP-PWM mappings
rbp_pwm = unique(rbp[, .(Gene_name, Matrix_id)])

# How many proteins have multiple PWMs?
rbp_pwm[, .N, by=Gene_name][order(N, decreasing = T),]

##############################################

# Select PWMs that match these RBPs
pwm = pwm[names(pwm) %in% rbp_pwm[, Matrix_id]]

# Scan the genome with these PWMs
sites = ScanGenomeWithPWMs(refString, names(ref), pwm)

##############################################
# Match binding sites with protein name
sites = merge(sites, rbp_pwm, by="Matrix_id")

# Find duplicate binding sites for the same protein and choose the one with the highest score
sites = sites[sites[, .I[score == max(score)], by=.(start, end, strand, Gene_name)]$V1]
sites = sites[!duplicated(sites[, .(start, end, strand, Gene_name)]),]

# sort by position
sites = sites[order(start),]

# Add site length
sites[, len:=nchar(seq)]

# Save sites
save(sites, file="output/sites.RData")

##############################################

# Annotate sites

sites_ranges = GRanges(sites)
site_overlaps = findOverlaps(sites_ranges, gff, ignore.strand=F, type="any")

annotated_sites = sites_ranges[queryHits(hits)]
annotated_sites$type = gff$type[subjectHits(hits)]
annotated_sites$gene = gff$gene[subjectHits(hits)]

# Save sites
save(annotated_sites, file="output/annotated_sites.RData")

