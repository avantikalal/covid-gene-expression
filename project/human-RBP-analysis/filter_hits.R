# Set working directory
setwd("/covid-omics")

library(data.table)

load("output/site_count.RData")
load("output/t_utr_site_count.RData")
load("output/f_utr_site_count.RData")

genome_hits = site_count[qval<0.01][strand=="+"][order(qval)][N>2]
t_utr_hits = t_utr_site_count[qval<0.01][strand=="+"][order(qval)][N>2]
f_utr_hits = f_utr_site_count[qval<0.01][strand=="+"][order(qval)][N>2]

save(genome_hits, file="output/genome_hits.RData")
save(t_utr_hits, file="output/t_utr_hits.RData")
save(f_utr_hits, file="output/f_utr_hits.RData")
