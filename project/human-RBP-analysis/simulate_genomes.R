##############################################

### Modify this section for your local environment!!### 

# Set working directory
setwd("proj_2_rbp_data")

# Import requirements
library(Biostrings)
library(rtracklayer)
source("covid-gene-expression/project/human-RBP-analysis/rbp_functions.R")

#######################

# load reference genome
ref = readFasta("reference/sequence.fasta.txt")

# load gff
gff = readGFF("reference/GCF_009858895.2_ASM985889v3_genomic.gff")

######################

# Simulate 1000 whole genomes

# Create simulated genomes
freqs = alphabetFrequency(ref)[1:4]
sim_genomes = simulateSeq(freqs, 1000)

#Save
save(sim_genomes, file="sim_genomes.RData")

######################

# Simulate 3' UTRs

# Extract 3' UTR sequence
t_utr_range = ranges(gff[gff$gbkey=="3'UTR"])
t_utr_seq = ref[[1]][start(t_utr_range):end(t_utr_range)]

# Create simulated 3'UTRs
freqs = alphabetFrequency(t_utr_seq)[1:4]
sim_t_utrs = simulateSeq(freqs, 1000)

save(sim_t_utrs, file="sim_t_utrs.RData")

######################

# Simulate 5' UTRs

# Extract 5' UTR sequence
f_utr_range = ranges(gff[gff$gbkey=="5'UTR"])
f_utr_seq = ref[[1]][start(f_utr_range):end(f_utr_range)]

# Create simulated 3'UTRs
freqs = alphabetFrequency(f_utr_seq)[1:4]
sim_f_utrs = simulateSeq(freqs, 1000)

save(sim_f_utrs, file="sim_f_utrs.RData")
