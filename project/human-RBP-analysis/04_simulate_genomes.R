##############################################

n_sims = 5000

# Import requirements
library(Biostrings)
library(rtracklayer)
source("covid-gene-expression/project/human-RBP-analysis/rbp_functions.R")

#######################

# load reference genome
print("Reading reference genome")
ref = readFasta("reference/sequence.fasta.txt")

# load gff
print("Reading GFF")
gff = readGFF("reference/GCF_009858895.2_ASM985889v3_genomic.gff")

######################

# Simulate whole genomes

# Create simulated genomes
print(paste0("Simulating ", n_sims, " genomes"))
freqs = alphabetFrequency(ref[[1]])[1:4]
sim_genomes = simulateSeq(freqs, n_sims)
print(sim_genomes)

#Save
print("Saving simulated genomes")
save(sim_genomes, file="output/sim_genomes.RData")

######################

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
save(sim_t_utrs, file="output/sim_t_utrs.RData")

######################

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
save(sim_f_utrs, file="output/sim_f_utrs.RData")
