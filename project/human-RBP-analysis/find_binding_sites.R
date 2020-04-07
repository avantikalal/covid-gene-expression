##############################################

### Modify this section for your local environment!!### 

# Set working directory
setwd("proj_2_rbp_data")

# List input data files

## NCBI # This file is the reference genome for SARS-CoV-2. It was downloaded from https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.
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

################################################

# Define functions

# Function to read PWMs
readPWMsFromFasta = function (pwm_file) {
  # Read all lines from PWM file
  lines = readLines(pwm_file)
  # Find header lines, start and end of each PWM
  ind = which(substr(lines, 1L, 1L) == ">")
  nseq = length(ind)
  start = ind + 1
  end = ind - 1
  end = c(end[-1], length(lines))
  # Split PWMs
  pwms = lapply(seq_len(nseq), function(i) strsplit(lines[start[i]:end[i]], "\t"))
  # Format as numeric matrix
  pwms = lapply(pwms, function(x) matrix(as.numeric(unlist(x)), ncol=4, byrow=T))
  # Name with PWM ID
  names(pwms) = lapply(seq_len(nseq), function(i) {
    firstword <- strsplit(lines[ind[i]], "\t")[[1]][1]
    substr(firstword, 2, nchar(firstword))
  })
  return(pwms)
}

##############################################

# Read input data

# Reference genome
ref = read.fasta(ref_file, as.string = T, forceDNAtolower = F)
refString = DNAString(ref[[1]][[1]])
#refString = RNAString(gsub("T", "U", ref[[1]][[1]]))

# RBPs
rbp = fread(rbp_file)

# PWMs
pwm = readPWMsFromFasta(pwm_file)

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
sites = list()
for(i in 1:length(pwm)){
  # Read PWM ID
  id = as.character(names(pwm)[i])
  # Save PWM as PWMatrix class
  curr_pwm = PWMatrix(profileMatrix=matrix(c(pwm[[i]]), byrow=TRUE, nrow=4, dimnames=list(c("A", "C", "G", "T"))), ID=id, name="Unknown")
  # Scan genome
  curr_sites = searchSeq(curr_pwm, refString, min.score="95%", strand="*", seqname = names(ref))
  if(length(curr_sites) > 0){
    # Save in list
    curr_sites = as.data.table(writeGFF3(curr_sites))
    curr_sites[, seq:= curr_sites[, tstrsplit(attributes, split=";|=", perl=T)][, V6]]
    curr_sites[, Matrix_id:= id]
    sites[[id]] = curr_sites 
  }
}
sites = rbindlist(sites)

##############################################
# Match protein name
sites = merge(sites, rbp_pwm, by="Matrix_id")

# Add site length
sites[, len:=nchar(seq)]

# Save sites
# save(sites, file="output/sites.RData")

##############################################



