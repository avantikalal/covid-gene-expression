##############################################

# List input data files

## NCBI # Downloaded from https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.
ref_file="reference/sequence.fasta.txt"

################################################

# Import requirements
library(data.table)
library(seqinr)
library(Biostrings)
library(TFBSTools)
library(foreach)
library(doParallel)
source("covid-gene-expression/project/human-RBP-analysis/rbp_functions.R")

registerDoParallel(detectCores())

################################################

# Read input data

# Reference genome
print("Reading reference genome")
ref = readFasta(ref_file)

# RBPs
print("Reading filtered RBPs")
load("output/filtered_rbp.RData")

# PWMs
print("Reading filtered PWMs")
load("output/filtered_pwm.RData")

##############################################

# Scan the genome with these PWMs
print("Scanning the reference genome for PWMs")
sites = ScanSeqWithPWMs(ref[[1]], pwm, names(ref))
print(paste0("Found ", nrow(sites), " sites."))

##############################################

# Select unique RBP-PWM mappings
print("Mapping RBPs to PWMs")
rbp_pwm = unique(rbp[, .(Gene_name, Matrix_id)])

# Match binding sites with protein name
print("Matching sites to RBPs")
sites = merge(sites, rbp_pwm, by="Matrix_id")

# Find duplicate binding sites for the same protein and choose the one with the highest score
print("Removing duplicate sites")
initial_n = nrow(sites)
sites = sites[sites[, .I[score == max(score)], by=.(start, end, strand, Gene_name)]$V1]
sites = sites[!duplicated(sites[, .(start, end, strand, Gene_name)]),]
print(paste0("Reduced number of sites from ", initial_n, " to ", nrow(sites)))

# sort by position
print("sorting sites")
sites = sites[order(start),]

# Add site length
print("adding site length")
sites[, len:=nchar(seq)]

# Save sites
print("saving sites")
save(sites, file="output/sites.RData")


