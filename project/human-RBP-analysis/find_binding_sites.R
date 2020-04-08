##############################################

### Modify this section for your local environment!!### 

# Set working directory
setwd("/covid-omics")

# List input data files

## NCBI # Downloaded from https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.
ref_file="reference/sequence.fasta.txt"

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
source("covid-gene-expression/project/human-RBP-analysis/rbp_functions.R")

registerDoParallel(detectCores())

################################################

# Read input data

# Reference genome
print("Reading reference genome")
ref = readFasta(ref_file)

# RBPs
print("Reading RBPs")
rbp = fread(rbp_file)

# PWMs
print("Reading PWMs")
pwm = readPWMsFromFasta(pwm_file)

##############################################

# Select only human RBPs
print("Filtering human RBPs")
initial_n = nrow(rbp)
rbp = rbp[Organism == "Homo_sapiens",]
print(paste0("Reduced number of RBPs from ", initial_n, " to ", nrow(rbp)))

# See experiment types
print("Types of experiments")
print(rbp[, .N, by=Experiment_description][order(N, decreasing = T),])

# See source databases
print("Source databases")
print(rbp[, .N, by=Database][order(N, decreasing = T),])

# How many unique proteins are present?
print("Number of unique human RBPs")
print(rbp[, length(unique(Gene_name))])

# Select unique RBP-PWM mappings
print("Mapping RBPs to PWMs")
print(rbp_pwm = unique(rbp[, .(Gene_name, Matrix_id)]))

# How many proteins have multiple PWMs?
print("Number of PWMs per protein")
print(rbp_pwm[, .N, by=Gene_name][order(N, decreasing = T),])

##############################################

# Select PWMs that match these RBPs
print("Selecting PWMs that match to human RBPs")
initial_n = nrow(pwm)
pwm = pwm[names(pwm) %in% rbp_pwm[, Matrix_id]]
print(paste0("Reduced number of PWMs from ", initial_n, " to ", nrow(pwm)))

# Scan the genome with these PWMs
print("Scanning the reference genome for PWMs")
sites = ScanGenomeWithPWMs(refString, names(ref), pwm)
print(paste0("Found ", nrow(sites), " sites."))

##############################################
# Match binding sites with protein name
print("Matching sites to RBPs")
sites = merge(sites, rbp_pwm, by="Matrix_id")

# Find duplicate binding sites for the same protein and choose the one with the highest score
print("Removing duplicate sites")
sites = sites[sites[, .I[score == max(score)], by=.(start, end, strand, Gene_name)]$V1]
sites = sites[!duplicated(sites[, .(start, end, strand, Gene_name)]),]

# sort by position
print("sorting sites")
sites = sites[order(start),]

# Add site length
print("adding site length")
sites[, len:=nchar(seq)]

# Save sites
print("saving sites")
save(sites, file="output/sites.RData")


