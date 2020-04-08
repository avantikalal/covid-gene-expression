### Modify this section for your local environment!!### 

# Set working directory
setwd("proj_2_rbp_data")

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

load("output/sim_genomes.RData")
load("output/sim_t_utrs.RData")
load("output/sim_f_utrs.RData")

# Load binding sites
load("output/sites.RData")
load("output/annotated_sites.RData")

# RBPs
rbp = fread("ATtRACT/ATtRACT_db.txt")

# Load PWMs
pwm = readPWMsFromFasta("ATtRACT/pwm.txt")
##############################################

# Select only human RBPs
rbp = rbp[Organism == "Homo_sapiens",]

# Select unique RBP-PWM mappings
rbp_pwm = unique(rbp[, .(Gene_name, Matrix_id)])

##############################################

# Scan simulated genomes

# Start with 100 genomes
#sim_genomes = sim_genomes[1:100]

# List candidate pwms - PWMs whose binding sites were found
cand_pwm = pwm[names(pwm) %in% sites[, Matrix_id]]

# Scan the simulated genomes with these PWMs
sim_genome_sites = ScanSeqWithPWMs(sim_genomes, cand_pwm)

# Merge sim_sites with protein name
sim_genome_sites = merge(sim_genome_sites, rbp_pwm, by="Matrix_id")

# Find duplicate binding sites for the same protein and choose the one with the highest score
sim_genome_sites = sim_genome_sites[sim_genome_sites[, .I[score == max(score)], by=.(start, end, strand, Gene_name)]$V1]
sim_genome_sites = sim_genome_sites[!duplicated(sim_genome_sites[, .(start, end, strand, Gene_name)]),]

# sort by position
sim_genome_sites = sim_genome_sites[order(start),]

# Save sites
save(sim_genome_sites, file="output/sim_genome_sites.RData") 

##############################################

# Enrichment calculation for simulated genomes

# Count sites per protein, by strand
site_count = sites[, .N, by=.(seqname, Gene_name, strand)]
sim_genome_site_count = sim_genome_sites[, .N, by=.(seqname, Gene_name, strand)]
sim_genome_site_count = sim_genome_site_count[, .(mean_count=mean(N), sd_count=sd(N)), by=.(Gene_name, strand)]
site_count = merge(site_count, sim_genome_site_count, by = c("Gene_name", "strand"))

# Calculate z-score
site_count[, z:=(N-mean_count)/sd_count]
site_count[, pval:=pnorm(-z)]

# Multiple hypothesis correction
site_count[, qval:=p.adjust(pval, "fdr"), by=strand]

# Save
save(site_count, file="output/site_count.RData")

# Select enriched
site_count[qval<0.01, ][order(qval),]

################################################

# Scan simulated 3' UTRs

t_utr_sites = annotated_sites[type == "three_prime_UTR"]

# List candidate pwms - PWMs whose binding sites were found
cand_pwm = pwm[names(pwm) %in% t_utr_sites[, Matrix_id]]

# Scan the simulated genomes with these PWMs
sim_t_utr_sites = ScanSeqWithPWMs(sim_t_utrs, cand_pwm)

# Merge sim_sites with protein name
sim_t_utr_sites = merge(sim_t_utr_sites, rbp_pwm, by="Matrix_id")

# Find duplicate binding sites for the same protein and choose the one with the highest score
sim_t_utr_sites = sim_t_utr_sites[sim_t_utr_sites[, .I[score == max(score)], by=.(start, end, strand, Gene_name)]$V1]
sim_t_utr_sites = sim_t_utr_sites[!duplicated(sim_t_utr_sites[, .(start, end, strand, Gene_name)]),]

# sort by position
sim_t_utr_sites = sim_t_utr_sites[order(start),]

# Save sites
save(sim_t_utr_sites, file="output/sim_t_utr_sites.RData") 

################################################

# Enrichment calculation separately for 3'UTR

# Count sites per protein, by strand
t_utr_site_count = t_utr_sites[, .N, by=.(seqname, Gene_name, strand)]
sim_t_utr_site_count = sim_t_utr_sites[, .N, by=.(seqname, Gene_name, strand)]
sim_t_utr_site_count = sim_t_utr_site_count[, .(mean_count=mean(N), sd_count=sd(N)), by=.(Gene_name, strand)]
t_utr_site_count = merge(site_count, sim_t_utr_site_count, by = c("Gene_name", "strand"))

# Calculate z-score
t_utr_site_count[, z:=(N-mean_count)/sd_count]
t_utr_site_count[, pval:=pnorm(-z)]

# Multiple hypothesis correction
t_utr_site_count[, qval:=p.adjust(pval, "fdr"), by=strand]

# Save
save(t_utr_site_count, file="output/t_utr_site_count.RData")

# Select enriched
t_utr_site_count[qval<0.01, ][order(qval),]
################################################

# Scan simulated 5' UTRs

f_utr_sites = annotated_sites[type == "five_prime_UTR"]

# List candidate pwms - PWMs whose binding sites were found
cand_pwm = pwm[names(pwm) %in% f_utr_sites[, Matrix_id]]

# Scan the simulated genomes with these PWMs
sim_f_utr_sites = ScanSeqWithPWMs(sim_f_utrs, cand_pwm)

# Merge sim_sites with protein name
sim_f_utr_sites = merge(sim_f_utr_sites, rbp_pwm, by="Matrix_id")

# Find duplicate binding sites for the same protein and choose the one with the highest score
sim_f_utr_sites = sim_f_utr_sites[sim_f_utr_sites[, .I[score == max(score)], by=.(start, end, strand, Gene_name)]$V1]
sim_f_utr_sites = sim_f_utr_sites[!duplicated(sim_f_utr_sites[, .(start, end, strand, Gene_name)]),]

# sort by position
sim_f_utr_sites = sim_f_utr_sites[order(start),]

# Save sites
save(sim_f_utr_sites, file="output/sim_f_utr_sites.RData") 

################################################

# Enrichment calculation separately for 5'UTR
