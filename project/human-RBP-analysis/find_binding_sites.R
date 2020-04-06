##############################################

### Modify this section for your local environment!!### 

# Set working directory
setwd("~/ngc-workspaces/mnt-covid-omics/")

# List input data files

## RBPDB
exp_table = "rbpdb/RBPDB_v1.3.1_experiments_human_2012-11-21.tdt"
protexp_table = "rbpdb/RBPDB_v1.3.1_protExp_human_2012-11-21.tdt"
prot_table = "rbpdb/RBPDB_v1.3.1_proteins_human_2012-11-21.tdt"
## The above three files were downloaded from http://rbpdb.ccbr.utoronto.ca/downloads/RBPDB_v1.3.1_human_2012-11-21_TDT.zip. Unzip the downloaded directory to get the tables.
pwm_dir = "rbpdb/matrices_human/PWMDir"
# The pwm_dir contains PWM files downloaded from http://rbpdb.ccbr.utoronto.ca/downloads/matrices_human.zip.

## NCBI
ref_file="ncbi/ref/sequence.fasta.txt"
# This file is the reference genome for SARS-CoV-2. It was downloaded from https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.

################################################

# Import requirements
library(data.table)
library(TFBSTools)
library(seqinr)
library(Biostrings)

# Define functions

# Funation to load PWMs
loadPWM = function(PWMfile, id, name){
  mat=as.matrix(read.table(PWMfile))
  pwm = PWMatrix(profileMatrix=matrix(c(mat), byrow=TRUE, nrow=4, dimnames=list(c("A", "C", "G", "T"))), ID=as.character(id), name=name)
  return(pwm)
}

# Function to scan sequence with PWM
scanGenome = function(ref_file, pwm, expID, protID, min.score="80%", strand="*"){
  # Read the genome sequence from FASTA file
  ref = read.fasta(ref_file, as.string = T, forceDNAtolower = F)
  refString=DNAString(ref[[1]][[1]])
  # Search the genome with the given PWM
  siteset = searchSeq(pwm, refString, min.score=min.score, strand=strand, seqname = names(ref))
  # Format the results as a table
  siteset = as.data.table(writeGFF3(siteset))
  # Add attributes to the table
  siteset_attributes = siteset[, tstrsplit(attributes, split=";|=", perl=T)]
  siteset[, tf:= siteset_attributes[, V2]]
  siteset[, seq:= siteset_attributes[, V6]]
  # Add experiment and protein IDs to the table
  siteset[, expID:= expID] 
  siteset[, protID:= protID] 
  return(siteset)
}

# Load tables from RBPDB
exp = fread(exp_table, header=F, col.names = c("id", "pmID", "exptype", "notes", "sequence_motif", "SELEX_file", 
                                               "aligned_SELEX_file", "aligned_motif_file", "PWM_file", "PFM_file", 
                                               "logo_file", "invivo_notes", "invivo_file", "secondary_structure", "flag"),
            na.strings = "\\N")
protexp = fread(protexp_table, header=F, col.names = c("protID", "expID", "homolog", "id"), na.strings = "\\N")
prot = fread(prot_table, header=F, col.names = c("id", "annotID", "createDate", "updateDate", "geneName", "geneDesc", 
                                                 "species", "taxID", "domains", "aliases", "flag", "flagNote", "PDBIDs"), 
             na.strings = "\\N")

# Select only the experiments for which PWM files are available
exp_with_pwm = exp[!is.na(PWM_file)]

# What types of experiments are these?
exp_with_pwm[, .N, by=exptype]

# Match these experiment to their protein IDs
exp_with_pwm = merge(exp_with_pwm, protexp, by.x="id", by.y="expID", all.x=T)

# Match these experiment to their protein names and other protein information
exp_with_pwm = merge(exp_with_pwm, prot, by.x="protID", by.y="id", all.x=T)

# For each experiment, load the matched PWM file
pwms = list()
for(i in 1:nrow(exp_with_pwm)){
  path_to_pwm_file = paste0(pwm_dir, "/", exp_with_pwm[i, PWM_file])
  expID = exp_with_pwm[i, id]
  protName = exp_with_pwm[i, geneName]
  pwms[[i]] = loadPWM(PWMfile = path_to_pwm_file, id = expID, name = protName)
}

# Scan the virus genome with each PWM to identify RBP-binding sites on the viral genome
sitesets = list()
for(i in 1:nrow(exp_with_pwm)){
  pwm = pwms[[i]]
  expID = exp_with_pwm[i, id]
  protID = exp_with_pwm[i, protID]
  sitesets[[i]] = scanGenome(ref_file = ref_file, pwm = pwm, expID=expID, protID=protID)
}
sitesets = rbindlist(sitesets)

# Save sitesets
save(sitesets, file="sitesets.RData")

