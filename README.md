# covid-gene-expression
Hackathon team: gene expression analysis for Covid-19

# Main Objective
We want to perform RNAseq-based analyses on published datasets in order to better understand the interaction between human host and virus.

## Deliverables
_Biological:_ Perform a global RNA-Seq analysis with SARS-CoV-2 infected datasets to search for new candidate genes for testing experimentally
_Methodological:_ Create a packaged reproducible pipeline in Docker to help scientists to easily treat their RNA-Seq data and for us if any new dataset comes out

[Tasks](https://github.com/avantikalal/covid-gene-expression/projects/1)

# workflow
0. Check literature to select interesting genes/datasets to study
1. Downloading RNAseqs from SRA/GEO 
2. Pipeline to clean reads
3. Map against viral genome (+ viral DBs)
4. Map against human genome (genes and isoforms pipelines)
5. Check shared reads between both genomes
6. Map reads against transposable elements
7. Check the existence of chimeric reads
8. Perform differential expression analysis on genes, transcripts and TEs
9. Functional enrichment analysis
10. SNP/Splicing on risk factor datasets for selected genes

# software-tools

source("http://bioconductor.org/biocLite.R")  
biocLite("DESeq2", dep=T)

## code repositories

# data-sources
- SARS-MERS: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56192  
- COVID: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507 
- Murine coronavirus (M-CoV): https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-4111/

# references
- https://bioconductor.org/packages/release/bioc/html/DESeq2.html
