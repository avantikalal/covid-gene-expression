### Modify this section for your local environment!!### 

# Set working directory
setwd("/covid-omics")

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

# RBPs
print("Loading RBPs")
rbp = fread("ATtRACT/ATtRACT_db.txt")

# Load PWMs
print("Loading PWMs")
pwm = readPWMsFromFasta("ATtRACT/pwm.txt")
##############################################

# Select only human RBPs
print("Filtering human RBPs")
initial_n = nrow(rbp)
rbp = rbp[Organism == "Homo_sapiens",]
print(paste0("Reduced number of RBPs from ", initial_n, " to ", nrow(rbp)))

# Select unique RBP-PWM mappings
print("Mapping RBPs to PWMs")
rbp_pwm = unique(rbp[, .(Gene_name, Matrix_id)])

##############################################

# Scan simulated genomes

# List candidate pwms - PWMs whose binding sites were found
print("Filtering PWMs that were found in the genome")
cand_pwm = pwm[names(pwm) %in% sites[, Matrix_id]]
print(paste0("Reduced number of PWMs from ", length(pwm), " to ", length(cand_pwm)))

# Scan the simulated genomes with these PWMs
print("Scanning simulated genomes with candidate PWMs")
sim_genome_sites = ScanSeqWithPWMs(sim_genomes, cand_pwm)
print(paste0("obtained ", nrow(sim_genome_sites), " binding sites from 1000 simulated genomes"))

# Merge sim_sites with protein name
print("Matching binding sites on simulated genomes to RBPs")
sim_genome_sites = merge(sim_genome_sites, rbp_pwm, by="Matrix_id")

# Find duplicate binding sites for the same protein and choose the one with the highest score
print("Removing duplicate sites")
initial_n = nrow(sim_genome_sites)
sim_genome_sites = sim_genome_sites[sim_genome_sites[, .I[score == max(score)], by=.(start, end, strand, Gene_name)]$V1]
sim_genome_sites = sim_genome_sites[!duplicated(sim_genome_sites[, .(start, end, strand, Gene_name)]),]
print(paste0("Reduced number of sites from ", initial_n, " to ", nrow(sim_genome_sites)))

# sort by position
print("sorting sites")
sim_genome_sites = sim_genome_sites[order(start),]

# Save sites
print("saving binding sites on simulated genomes")
save(sim_genome_sites, file="output/sim_genome_sites.RData") 

##############################################

# Enrichment calculation for simulated genomes

# Count sites per protein, by strand
print("Counting binding sites per protein per strand, on the real genome")
site_count = sites[, .N, by=.(seqname, Gene_name, strand)]

print("Counting binding sites per protein per strand, on the simulated genomes")
sim_genome_site_count = sim_genome_sites[, .N, by=.(seqname, Gene_name, strand)]
sim_genome_site_count = sim_genome_site_count[, .(mean_count=mean(N), sd_count=sd(N)), by=.(Gene_name, strand)]

print("Comparing binding site counts")
site_count = merge(site_count, sim_genome_site_count, by = c("Gene_name", "strand"))

# Calculate z-score
print("Calculating z-scores")
site_count[, z:=(N-mean_count)/sd_count]
site_count[, pval:=pnorm(-z)]

# Multiple hypothesis correction
print("FDR correction")
site_count[, qval:=p.adjust(pval, "fdr"), by=strand]

# Save
print("Saving results")
save(site_count, file="output/site_count.RData")

# Select enriched
print(site_count[qval<0.01, ][order(qval),])

################################################

# Scan simulated 3' UTRs
print("Extracting binding sites in the 3'UTR")
t_utr_sites = annotated_sites[type == "three_prime_UTR"]
print(paste0("Found ", nrow(t_utr_sites), " sites on the 3'UTR."))

# List candidate pwms - PWMs whose binding sites were found
print("Filtering PWMs that were found in the 3'UTR")
cand_pwm = pwm[names(pwm) %in% t_utr_sites[, Matrix_id]]
print(paste0("Reduced number of PWMs from ", length(pwm), " to ", length(cand_pwm)))

# Scan the simulated genomes with these PWMs
print("Scanning simulated 3'UTRs with candidate PWMs")
sim_t_utr_sites = ScanSeqWithPWMs(sim_t_utrs, cand_pwm)
print(paste0("obtained ", nrow(sim_t_utr_sites), " binding sites from 1000 simulated 3'UTRs"))

# Merge sim_sites with protein name
print("Matching binding sites on simulated 3'UTRs to RBPs")
sim_t_utr_sites = merge(sim_t_utr_sites, rbp_pwm, by="Matrix_id")

# Find duplicate binding sites for the same protein and choose the one with the highest score
print("Removing duplicate sites")
initial_n = nrow(sim_t_utr_sites)
sim_t_utr_sites = sim_t_utr_sites[sim_t_utr_sites[, .I[score == max(score)], by=.(start, end, strand, Gene_name)]$V1]
sim_t_utr_sites = sim_t_utr_sites[!duplicated(sim_t_utr_sites[, .(start, end, strand, Gene_name)]),]
print(paste0("Reduced number of sites from ", initial_n, " to ", nrow(sim_t_utr_sites)))

# sort by position
print("sorting sites")
sim_t_utr_sites = sim_t_utr_sites[order(start),]

# Save sites
print("saving binding sites on simulated 3'UTRs")
save(sim_t_utr_sites, file="output/sim_t_utr_sites.RData") 

################################################

# Enrichment calculation separately for 3'UTR

# Count sites per protein, by strand
print("Counting binding sites per protein per strand, on the real 3'UTR")
t_utr_site_count = t_utr_sites[, .N, by=.(seqname, Gene_name, strand)]

print("Counting binding sites per protein per strand, on the simulated 3'UTRs")
sim_t_utr_site_count = sim_t_utr_sites[, .N, by=.(seqname, Gene_name, strand)]
sim_t_utr_site_count = sim_t_utr_site_count[, .(mean_count=mean(N), sd_count=sd(N)), by=.(Gene_name, strand)]

print("Comparing binding site counts")
t_utr_site_count = merge(site_count, sim_t_utr_site_count, by = c("Gene_name", "strand"))

# Calculate z-score
print("Calculating z-scores")
t_utr_site_count[, z:=(N-mean_count)/sd_count]
t_utr_site_count[, pval:=pnorm(-z)]

# Multiple hypothesis correction
print("FDR correction")
t_utr_site_count[, qval:=p.adjust(pval, "fdr"), by=strand]

# Save
print("Saving results")
save(t_utr_site_count, file="output/t_utr_site_count.RData")

# Select enriched
print(t_utr_site_count[qval<0.01, ][order(qval),])
################################################

# Scan simulated 5' UTRs
print("Extracting binding sites in the 5'UTR")
f_utr_sites = annotated_sites[type == "five_prime_UTR"]

# List candidate pwms - PWMs whose binding sites were found
print("Filtering PWMs that were found in the 5'UTR")
cand_pwm = pwm[names(pwm) %in% f_utr_sites[, Matrix_id]]
print(paste0("Reduced number of PWMs from ", length(pwm), " to ", length(cand_pwm)))

# Scan the simulated genomes with these PWMs
print("Scanning simulated 5'UTRs with candidate PWMs")
sim_f_utr_sites = ScanSeqWithPWMs(sim_f_utrs, cand_pwm)
print(paste0("obtained ", nrow(sim_f_utr_sites), " binding sites from 1000 simulated 5'UTRs"))

# Merge sim_sites with protein name
print("Matching binding sites on simulated 5'UTRs to RBPs")
sim_f_utr_sites = merge(sim_f_utr_sites, rbp_pwm, by="Matrix_id")

# Find duplicate binding sites for the same protein and choose the one with the highest score
print("Removing duplicate sites")
initial_n = nrow(sim_f_utr_sites)
sim_f_utr_sites = sim_f_utr_sites[sim_f_utr_sites[, .I[score == max(score)], by=.(start, end, strand, Gene_name)]$V1]
sim_f_utr_sites = sim_f_utr_sites[!duplicated(sim_f_utr_sites[, .(start, end, strand, Gene_name)]),]
print(paste0("Reduced number of sites from ", initial_n, " to ", nrow(sim_f_utr_sites)))

# sort by position
print("sorting sites")
sim_f_utr_sites = sim_f_utr_sites[order(start),]

# Save sites
print("Saving binding sites on simulated 5'UTRs")
save(sim_f_utr_sites, file="output/sim_f_utr_sites.RData") 

################################################

# Enrichment calculation separately for 5'UTR
# Count sites per protein, by strand
print("Counting binding sites per protein per strand, on the real 5'UTR")
f_utr_site_count = f_utr_sites[, .N, by=.(seqname, Gene_name, strand)]

print("Counting binding sites per protein per strand, on the simulated 5'UTRs")
sim_f_utr_site_count = sim_f_utr_sites[, .N, by=.(seqname, Gene_name, strand)]
sim_f_utr_site_count = sim_f_utr_site_count[, .(mean_count=mean(N), sd_count=sd(N)), by=.(Gene_name, strand)]

print("Comparing binding site counts")
f_utr_site_count = merge(site_count, sim_f_utr_site_count, by = c("Gene_name", "strand"))

# Calculate z-score
print("Calculating z-scores")
f_utr_site_count[, z:=(N-mean_count)/sd_count]
f_utr_site_count[, pval:=pnorm(-z)]

# Multiple hypothesis correction
print("FDR correction")
f_utr_site_count[, qval:=p.adjust(pval, "fdr"), by=strand]

# Save
print("Saving results")
save(f_utr_site_count, file="output/f_utr_site_count.RData")

# Select enriched
print(f_utr_site_count[qval<0.01, ][order(qval),])