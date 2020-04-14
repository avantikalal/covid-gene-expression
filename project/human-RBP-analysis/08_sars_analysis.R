##############################################

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

# List input data files

## NCBI genome for bat coronavirus
ref_file="reference/SARS1/sequence.fasta.txt"

#GFF file for SARS
gff_file="reference/SARS1/sequence.gff3"

n_sims = 5000 # Number of simulated genomes

################################################

# Read input data

# Reference genome
print("Reading reference genome")
ref = readFasta(ref_file)

#GFF
print("Reading GFF")
gff = readGFF(gff_file)

# RBPs
print("Reading filtered RBPs")
load("output/filtered_rbp.RData")

# PWMs
print("Reading filtered PWMs")
load("output/filtered_pwm.RData")

##############################################

# Scan the genome with these PWMs
print(paste0("Scanning the reference genome ", ref_file, " for PWMs"))
sites = ScanSeqWithPWMs(ref[[1]], pwm, names(ref))
print(paste0("Found ", nrow(sites), " sites."))

# Match binding sites with protein name
print("Matching sites to RBPs")
sites = merge(sites, unique(rbp[, .(Gene_name, Matrix_id)]), by="Matrix_id")

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
save(sites, file=paste0("output/", names(ref), "/sites.RData"))

##############################################################

# Annotate binding sites
print("overlapping sites with genomic features")
sites_ranges = GRanges(sites)
site_overlaps = findOverlaps(sites_ranges, gff, ignore.strand=F, type="any")

print("Annotating sites")
annotated_sites = sites_ranges[queryHits(site_overlaps)]
annotated_sites$type = gff$type[subjectHits(site_overlaps)]
annotated_sites$gene = gff$gene[subjectHits(site_overlaps)]

# Save sites
print("Saving annotated sites")
save(annotated_sites, file=paste0("output/", names(ref), "/annotated_sites.RData"))

##############################################################

# Simulate whole genomes

# Create simulated genomes
print(paste0("Simulating ", n_sims, " genomes"))
freqs = alphabetFrequency(ref[[1]])[1:4]
sim_genomes = simulateSeq(freqs, n_sims)
print(sim_genomes)

#Save
print("Saving simulated genomes")
save(sim_genomes, file=paste0("output/", names(ref), "/sim_genomes.RData"))

# Simulate 3' UTRs

# Extract 3' UTR sequence
print("Extracting 3'UTR")
t_utr_range = ranges(gff[gff$gbkey=="3'UTR"])
t_utr_seq = ref[[1]][start(t_utr_range):end(t_utr_range)]

# Create simulated 3'UTRs
print(paste0("Simulating ", n_sims, " 3'UTRs"))
freqs = alphabetFrequency(t_utr_seq)[1:4]
sim_t_utrs = simulateSeq(freqs, n_sims)
print(sim_t_utrs)

print("Saving simulated 3'UTRs")
save(sim_t_utrs, file=paste0("output/", names(ref), "/sim_t_utrs.RData"))

###############################################

# Simulate 5' UTRs

# Extract 5' UTR sequence
print("Extracting 5'UTR")
f_utr_range = ranges(gff[gff$gbkey=="5'UTR"])
f_utr_seq = ref[[1]][start(f_utr_range):end(f_utr_range)]

# Create simulated 3'UTRs
print(paste0("Simulating ", n_sims, " 5'UTRs"))
freqs = alphabetFrequency(f_utr_seq)[1:4]
sim_f_utrs = simulateSeq(freqs, n_sims)
print(sim_f_utrs)

print("Saving simulated 5'UTRs")
save(sim_f_utrs, file=paste0("output/", names(ref), "/sim_f_utrs.RData"))

####################################################################

# Scan simulated genomes

# List all PWMs for RBPs that were found in the genome
print("Filtering PWMs for RBPs that were found in the genome")
cand_pwm_ids = rbp[Gene_name %in% sites[, Gene_name], unique(Matrix_id)]
cand_pwm = pwm[names(pwm) %in% cand_pwm_ids]
print(paste0("Reduced number of PWMs from ", length(pwm), " to ", length(cand_pwm)))

# Scan the simulated genomes with these PWMs
print("Scanning simulated genomes with candidate PWMs")
sim_genome_sites = ScanSeqWithPWMs(sim_genomes, cand_pwm)
print(paste0("obtained ", nrow(sim_genome_sites), " binding sites from ", n_sims, "simulated genomes"))

# Save sites
print("saving binding sites on simulated genomes")
save(sim_genome_sites, file=paste0("output/", names(ref), "/sim_genome_sites_raw.RData"))

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
print(paste0("obtained ", nrow(sim_t_utr_sites), " binding sites from ", n_sims, "simulated 3'UTRs"))

# Save sites
print("saving binding sites on simulated 3'UTRs")
save(sim_t_utr_sites, file=paste0("output/", names(ref), "/sim_t_utr_sites_raw.RData"))

################################################

# Scan simulated 5' UTRs
print("Extracting binding sites in the 5'UTR")
f_utr_sites = annotated_sites[annotated_sites$type == "five_prime_UTR",]
print(paste0("Found ", nrow(f_utr_sites), " sites on the 5'UTR."))

# List candidate pwms - PWMs for RBPs whose binding sites were found
print("Filtering PWMs for RBPs that were found in the 5'UTR")
cand_pwm_ids = rbp[Gene_name %in% f_utr_sites$Gene_name, unique(Matrix_id)]
cand_pwm = pwm[names(pwm) %in% cand_pwm_ids]
print(paste0("Reduced number of PWMs from ", length(pwm), " to ", length(cand_pwm)))

# Scan the simulated genomes with these PWMs
print("Scanning simulated 5'UTRs with candidate PWMs")
sim_f_utr_sites = ScanSeqWithPWMs(sim_f_utrs, cand_pwm)
print(paste0("obtained ", nrow(sim_f_utr_sites), " binding sites from ", n_sims, " simulated 5'UTRs"))

# Save sites
print("Saving binding sites on simulated 5'UTRs")
save(sim_f_utr_sites, file=paste0("output/", names(ref), "/sim_f_utr_sites_raw.RData"))

###############################################

# Match binding sites with protein name
print("Matching binding sites on simulated genomes to RBPs")
sim_genome_sites = merge(sim_genome_sites, unique(rbp[, .(Gene_name, Matrix_id)]), by="Matrix_id", allow.cartesian=TRUE)

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
save(site_count, file=paste0("output/", names(ref), "/site_count.RData"))

############################################

# 3'UTR

# Match binding sites with protein name
print("Matching binding sites on simulated 3'UTRs to RBPs")
sim_t_utr_sites = merge(sim_t_utr_sites, unique(rbp[, .(Gene_name, Matrix_id)]), by="Matrix_id")

# Find duplicate binding sites for the same protein and choose the one with the highest score
print("Removing duplicate sites")
initial_n = nrow(sim_t_utr_sites)
sim_t_utr_sites = sim_t_utr_sites[sim_t_utr_sites[, .I[score == max(score)], by=.(seqname, start, end, strand, Gene_name)]$V1]
sim_t_utr_sites = sim_t_utr_sites[!duplicated(sim_t_utr_sites[, .(seqname, start, end, strand, Gene_name)]),]
print(paste0("Reduced number of sites from ", initial_n, " to ", nrow(sim_t_utr_sites)))

# Enrichment
print("calculating enrichment in 3'UTR")
t_utr_site_count = enrich_rbps_real(as.data.table(t_utr_sites), sim_t_utr_sites)

# Save
print("Saving results")
save(t_utr_site_count, file=paste0("output/", names(ref), "/t_utr_site_count.RData"))

################################################

# Match binding sites with protein name
print("Matching binding sites on simulated 5'UTRs to RBPs")
sim_f_utr_sites = merge(sim_f_utr_sites, unique(rbp[, .(Gene_name, Matrix_id)]), by="Matrix_id")

# Find duplicate binding sites for the same protein and choose the one with the highest score
print("Removing duplicate sites")
initial_n = nrow(sim_f_utr_sites)
sim_f_utr_sites = sim_f_utr_sites[sim_f_utr_sites[, .I[score == max(score)], by=.(seqname, start, end, strand, Gene_name)]$V1]
sim_f_utr_sites = sim_f_utr_sites[!duplicated(sim_f_utr_sites[, .(seqname, start, end, strand, Gene_name)]),]
print(paste0("Reduced number of sites from ", initial_n, " to ", nrow(sim_f_utr_sites)))

# Enrichment
print("calculating enrichment in 5'UTR")
f_utr_site_count = enrich_rbps_real(as.data.table(f_utr_sites), sim_f_utr_sites)

# Save
print("Saving results")
save(f_utr_site_count, file=paste0("output/", names(ref), "/f_utr_site_count.RData"))

#########################################################################

genome_hits = site_count[qval<0.01][strand=="+"][order(qval)][N>2]
t_utr_hits = t_utr_site_count[qval<0.01][strand=="+"][order(qval)][N>2]
f_utr_hits = f_utr_site_count[qval<0.01][strand=="+"][order(qval)][N>2]
genome_neg_hits = site_count[qval<0.01][strand=="-"][order(qval)][N>2]

save(genome_hits, file=paste0("output/", names(ref), "/genome_hits.RData"))
save(t_utr_hits, file=paste0("output/", names(ref), "/t_utr_hits.RData"))
save(f_utr_hits, file=paste0("output/", names(ref), "/f_utr_hits.RData"))
save(genome_neg_hits, file=paste0("output/", names(ref), "/genome_neg_hits.RData"))

