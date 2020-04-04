# covid-gene-expression
Hackathon team: gene expression analysis for Covid-19

[Tasks](https://github.com/avantikalal/covid-gene-expression/projects/1)

# workflow
0. Check literature to select interesting genes to study
1. Gathering datasets and downloading RNAseq from SRA/GEO 
2. Pipeline to clean reads and map against human genome
3. SNP calling: [kissplice]http://kissplice.prabi.fr/)
4. Splicing events
5. Differential gene expression analysis (for all genes and SNPs)
6. Pathway analysis

# software-tools

source("http://bioconductor.org/biocLite.R")  
biocLite("DESeq2", dep=T)

# data-sources
- SARS-MERS: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56192  
- COVID: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507  
