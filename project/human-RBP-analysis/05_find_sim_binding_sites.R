
# Import requirements
library(data.table)
library(Biostrings)
library(foreach)
library(doParallel)
library(rtracklayer)
source("covid-gene-expression/project/human-RBP-analysis/rbp_functions.R")

registerDoParallel(detectCores())

##############################################
# Load simulated data
print("Loading simulated sequences")
load("output/sim_genomes.RData")
load("output/sim_t_utrs.RData")
load("output/sim_f_utrs.RData")

# Load binding sites
print("Loading binding sites")
load("output/sites.RData")
load("output/annotated_sites.RData")

# Load fitered PWMs
print("Loading PWMs")
load("output/filtered_pwm.RData")

# Load filtered RBPs
print("Loading RBPs")
load("output/filtered_rbp.RData")
##############################################

# Scan simulated genomes

# List all PWMs for RBPs that were found in the genome
print("Filtering PWMs for RBPs that were found in the genome")
cand_pwm_ids = rbp[Gene_name %in% sites[, Gene_name], unique(Matrix_id)]
cand_pwm = pwm[names(pwm) %in% cand_pwm_ids]
print(paste0("Reduced number of PWMs from ", length(pwm), " to ", length(cand_pwm)))
 
# Scan the simulated genomes with these PWMs
print("Scanning simulated genomes with candidate PWMs")
sim_genome_sites = ScanSeqWithPWMs(sim_genomes, cand_pwm)
print(paste0("obtained ", nrow(sim_genome_sites), " binding sites from 1000 simulated genomes"))
 
# Save sites
print("saving binding sites on simulated genomes")
save(sim_genome_sites, file="output/sim_genome_sites_raw.RData") 

################################################

# Scan simulated 3' UTRs
print("Extracting binding sites in the 3'UTR")
t_utr_sites = annotated_sites[annotated_sites$type == "three_prime_UTR",]
print(paste0("Found ", nrow(t_utr_sites), " sites on the 3'UTR."))

# List candidate pwms - PWMs for RBPs whose binding sites were found
print("Filtering PWMs of RBPs that were found in the 3'UTR")
cand_pwm_ids = rbp[Gene_name %in% t_utr_sites$Gene_name, unique(Matrix_id)]
cand_pwm = pwm[names(pwm) %in% cand_pwm_ids]
print(paste0("Reduced number of PWMs from ", length(pwm), " to ", length(cand_pwm)))

# Scan the simulated genomes with these PWMs
print("Scanning simulated 3'UTRs with candidate PWMs")
sim_t_utr_sites = ScanSeqWithPWMs(sim_t_utrs, cand_pwm)
print(paste0("obtained ", nrow(sim_t_utr_sites), " binding sites from 1000 simulated 3'UTRs"))

# Save sites
print("saving binding sites on simulated 3'UTRs")
save(sim_t_utr_sites, file="output/sim_t_utr_sites_raw.RData") 

################################################

# Scan simulated 5' UTRs
print("Extracting binding sites in the 5'UTR")
f_utr_sites = annotated_sites[annotated_sites$type == "five_prime_UTR",]

# List candidate pwms - PWMs for RBPs whose binding sites were found
print("Filtering PWMs for RBPs that were found in the 5'UTR")
cand_pwm_ids = rbp[Gene_name %in% f_utr_sites$Gene_name, unique(Matrix_id)]
cand_pwm = pwm[names(pwm) %in% cand_pwm_ids]
print(paste0("Reduced number of PWMs from ", length(pwm), " to ", length(cand_pwm)))

# Scan the simulated genomes with these PWMs
print("Scanning simulated 5'UTRs with candidate PWMs")
sim_f_utr_sites = ScanSeqWithPWMs(sim_f_utrs, cand_pwm)
print(paste0("obtained ", nrow(sim_f_utr_sites), " binding sites from 1000 simulated 5'UTRs"))

# Save sites
print("Saving binding sites on simulated 5'UTRs")
save(sim_f_utr_sites, file="output/sim_f_utr_sites_raw.RData") 