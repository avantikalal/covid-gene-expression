
# Import requirements
library(data.table)
library(seqinr)
library(Biostrings)
source("covid-gene-expression/project/human-RBP-analysis/rbp_functions.R")

##############################################

# Load binding sites from real sequences
print("Loading binding sites on the reference genome")
load("output/sites.RData")
load("output/annotated_sites.RData")

# Load binding sites from simulated sequences
print("Loading binding sites on simulated sequences")
load("output/sim_genome_sites_raw.RData")
load("output/sim_f_utr_sites_raw.RData")
load("output/sim_t_utr_sites_raw.RData")

##############################################

# RBPs
print("Loading RBPs")
load("output/filtered_rbp.RData")

# Select unique RBP-PWM mappings
print("Mapping RBPs to PWMs")
rbp_pwm = unique(rbp[, .(Gene_name, Matrix_id)])

##############################################

# Match binding sites with protein name
print("Matching binding sites on simulated genomes to RBPs")
sim_genome_sites = merge(sim_genome_sites, rbp_pwm, by="Matrix_id", allow.cartesian=TRUE)

# Find duplicate binding sites for the same protein and choose the one with the highest score
print("Removing duplicate sites")
initial_n = nrow(sim_genome_sites)
sim_genome_sites = sim_genome_sites[sim_genome_sites[, .I[score == max(score)], by=.(seqname, start, end, strand, Gene_name)]$V1]
sim_genome_sites = sim_genome_sites[!duplicated(sim_genome_sites[, .(seqname, start, end, strand, Gene_name)]),]
print(paste0("Reduced number of sites from ", initial_n, " to ", nrow(sim_genome_sites)))

# Enrichment
print("calculating enrichment in the whole genome")
site_count = enrich_rbps_real(sites, sim_genome_sites)

# Save
print("Saving results")
save(site_count, file="output/site_count.RData")

############################################

# 3'UTR

# Match binding sites with protein name
print("Matching binding sites on simulated 3'UTRs to RBPs")
sim_t_utr_sites = merge(sim_t_utr_sites, rbp_pwm, by="Matrix_id")

# Find duplicate binding sites for the same protein and choose the one with the highest score
print("Removing duplicate sites")
initial_n = nrow(sim_t_utr_sites)
sim_t_utr_sites = sim_t_utr_sites[sim_t_utr_sites[, .I[score == max(score)], by=.(seqname, start, end, strand, Gene_name)]$V1]
sim_t_utr_sites = sim_t_utr_sites[!duplicated(sim_t_utr_sites[, .(seqname, start, end, strand, Gene_name)]),]
print(paste0("Reduced number of sites from ", initial_n, " to ", nrow(sim_t_utr_sites)))

# Get predicted binding sites in the real 3' UTR
print("Selecting binding sites on the real 3'UTR")
t_utr_sites = as.data.table(annotated_sites[annotated_sites$type == "three_prime_UTR",])

# Enrichment
print("calculating enrichment in 3'UTR")
t_utr_site_count = enrich_rbps_real(t_utr_sites, sim_t_utr_sites)

# Save
print("Saving results")
save(t_utr_site_count, file="output/t_utr_site_count.RData")

################################################

# Match binding sites with protein name
print("Matching binding sites on simulated 5'UTRs to RBPs")
sim_f_utr_sites = merge(sim_f_utr_sites, rbp_pwm, by="Matrix_id")

# Find duplicate binding sites for the same protein and choose the one with the highest score
print("Removing duplicate sites")
initial_n = nrow(sim_f_utr_sites)
sim_f_utr_sites = sim_f_utr_sites[sim_f_utr_sites[, .I[score == max(score)], by=.(seqname, start, end, strand, Gene_name)]$V1]
sim_f_utr_sites = sim_f_utr_sites[!duplicated(sim_f_utr_sites[, .(seqname, start, end, strand, Gene_name)]),]
print(paste0("Reduced number of sites from ", initial_n, " to ", nrow(sim_f_utr_sites)))

# Get predicted binding sites in the real 5' UTR
print("Selecting binding sites on the real 5'UTR")
f_utr_sites = as.data.table(annotated_sites[annotated_sites$type == "five_prime_UTR",])

# Enrichment
print("calculating enrichment in 5'UTR")
f_utr_site_count = enrich_rbps_real(f_utr_sites, sim_f_utr_sites)

# Save
print("Saving results")
save(f_utr_site_count, file="output/f_utr_site_count.RData")
