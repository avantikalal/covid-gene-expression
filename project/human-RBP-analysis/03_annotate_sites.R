
# List input data files

## NCBI # Downloaded from https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.
gff_file = "reference/GCF_009858895.2_ASM985889v3_genomic.gff"

################################################

# Import requirements
library(data.table)
library(seqinr)
library(Biostrings)
library(TFBSTools)
library(foreach)
library(doParallel)
source("covid-gene-expression/project/human-RBP-analysis/rbp_functions.R")

################################################

# Read GFF
print("Reading GFF")
gff = readGFF(gff_file)

# Load sites
load("output/sites.RData")

##############################################

# Annotate sites

print("overlapping sites with genomic features")
sites_ranges = GRanges(sites)
site_overlaps = findOverlaps(sites_ranges, gff, ignore.strand=F, type="any")

print("Annotating sites")
annotated_sites = sites_ranges[queryHits(site_overlaps)]
annotated_sites$type = gff$type[subjectHits(site_overlaps)]
annotated_sites$gene = gff$gene[subjectHits(site_overlaps)]

# Save sites
print("Saving annotated sites")
save(annotated_sites, file="output/annotated_sites.RData")
