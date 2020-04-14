# Import requirements

library(data.table)

# Load site counts
load("output/site_count.RData")
load("output/t_utr_site_count.RData")
load("output/f_utr_site_count.RData")

# Filter significant hits
genome_hits = site_count[qval<0.01][strand=="+"][order(qval)][N>2]
t_utr_hits = t_utr_site_count[qval<0.01][strand=="+"][order(qval)][N>2]
f_utr_hits = f_utr_site_count[qval<0.01][strand=="+"][order(qval)][N>2]
genome_neg_hits = site_count[qval<0.01][strand=="-"][order(qval)][N>2]

# Save RData
save(genome_hits, file="output/genome_hits.RData")
save(genome_neg_hits, file="output/genome_neg_hits.RData")
save(t_utr_hits, file="output/t_utr_hits.RData")
save(f_utr_hits, file="output/f_utr_hits.RData")

# Save tsv
write.table(genome_hits, file="output/genome_hits.tsv", sep="\t", quote=F)
write.table(genome_neg_hits, file="output/genome_neg_hits.tsv", sep="\t", quote=F)
write.table(t_utr_hits, file="output/t_utr_hits.tsv", sep="\t", quote=F)
write.table(f_utr_hits, file="output/f_utr_hits.tsv", sep="\t", quote=F)