# Import requirements

library(data.table)
source("covid-gene-expression/project/human-RBP-analysis/rbp_functions.R")

################################################

# Compare hits in whole genome

load("output/genome_hits.RData")
sars2 = genome_hits

load("output/MN996532.1/genome_hits.RData")
bat = genome_hits

load("output/NC_004718.3/genome_hits.RData")
sars = genome_hits

setdiff(sars2[, Gene_name], bat[, Gene_name]) #character(0)
setdiff(sars2[, Gene_name], sars[, Gene_name]) #"PABPC1" "ZNF638" "PABPC5" "CELF2"  "PABPC3" "YBX1" 
################################################

# Compare hits in 3'UTR
load("output/t_utr_hits.RData")
sars2 = genome_hits

load("output/MN996532.1/t_utr_hits.RData")
bat = genome_hits

load("output/NC_004718.3/t_utr_hits.RData")
sars = genome_hits

setdiff(sars2[, Gene_name], bat[, Gene_name]) #character(0)
setdiff(sars2[, Gene_name], sars[, Gene_name]) #character(0)

###############################################
# Compare hits in 5'UTR

load("output/f_utr_hits.RData")
sars2 = genome_hits

load("output/MN996532.1/f_utr_hits.RData")
bat = genome_hits

load("output/NC_004718.3/f_utr_hits.RData")
sars = genome_hits

setdiff(sars2[, Gene_name], bat[, Gene_name]) #character(0)
setdiff(sars2[, Gene_name], sars[, Gene_name]) #character(0)

###############################################

# Compare hits in - strand
load("output/genome_neg_hits.RData")
sars2 = genome_neg_hits

load("output/MN996532.1/genome_neg_hits.RData")
bat = genome_neg_hits

load("output/NC_004718.3/genome_neg_hits.RData")
sars = genome_neg_hits

setdiff(sars2[, Gene_name], bat[, Gene_name]) #"ELAVL1"  "HNRNPDL" "RBFOX1" 
setdiff(sars2[, Gene_name], sars[, Gene_name]) #"SRSF3"  "TIA1"   "ELAVL1" "FUS"    

