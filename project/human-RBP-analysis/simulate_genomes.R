##############################################

### Modify this section for your local environment!!### 

# Set working directory
setwd("/covid-omics")

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

# Simulate 1000 whole genomes

# Create simulated genomes
print("Simulating 1000 genomes")
freqs = alphabetFrequency(ref)[1:4]
sim_genomes = simulateSeq(freqs, 1000)

#Save
print("Saving simulated genomes")
save(sim_genomes, file="sim_genomes.RData")

######################

# Simulate 3' UTRs

# Extract 3' UTR sequence
print("Extracting 3'UTR")
t_utr_range = ranges(gff[gff$gbkey=="3'UTR"])
t_utr_seq = ref[[1]][start(t_utr_range):end(t_utr_range)]

# Create simulated 3'UTRs
print("Simulating 1000 3'UTRs")
freqs = alphabetFrequency(t_utr_seq)[1:4]
sim_t_utrs = simulateSeq(freqs, 1000)

print("Saving simulated 3'UTRs")
save(sim_t_utrs, file="sim_t_utrs.RData")

######################

# Simulate 5' UTRs

# Extract 5' UTR sequence
print("Extracting 5'UTR")
f_utr_range = ranges(gff[gff$gbkey=="5'UTR"])
f_utr_seq = ref[[1]][start(f_utr_range):end(f_utr_range)]

# Create simulated 3'UTRs
print("Simulating 1000 5'UTRs")
freqs = alphabetFrequency(f_utr_seq)[1:4]
sim_f_utrs = simulateSeq(freqs, 1000)

print("Saving simulated 5'UTRs")
save(sim_f_utrs, file="sim_f_utrs.RData")
