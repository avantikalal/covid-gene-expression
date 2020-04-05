# Covid-19 Gene Expression Work Group (Private project)
Hackathon team: gene expression analysis for Covid-19 Virtual Biohackathon (vBH)

### Link to vBH Github
https://github.com/virtual-biohackathons/covid-19-bh20
### Official Gene Expression Work Group Page (Public page)
https://github.com/virtual-biohackathons/covid-19-bh20/wiki/GeneExpression

## Main Objective
We want to perform RNAseq-based analyses on published datasets in order to better understand the interaction between human host and virus.

## Deliverables
_Biological:_ Perform a global RNA-Seq analysis with SARS-CoV-2 infected datasets to search for new candidate genes for testing experimentally

_Methodological:_ Create a packaged reproducible pipeline in Docker to help scientists to easily treat their RNA-Seq data and for us if any new dataset comes out

[Tasks](https://github.com/avantikalal/covid-gene-expression/projects/1)

## Workflow
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

# Software Used during this project

## Command line tools
- FastQC (https://github.com/s-andrews/FastQC)
- Fastq-screen (https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
- trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic)
- STAR (https://github.com/alexdobin/STAR)
- Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- TEtools (https://github.com/l-modolo/TEtools)
- TEtranscripts (http://hammelllab.labsites.cshl.edu/software)
- LIONS (www.github.com/ababaian/LIONS)
- samtools / picard (http://samtools.sourceforge.net/; https://broadinstitute.github.io/picard/)
- featureCounts (http://subread.sourceforge.net)
- MultiQC (https://github.com/ewels/MultiQC)
- KissSplice (http://kissplice.prabi.fr/)

## R packages

### Direct install
   _ggplot2_
```
install.packages("ggplot2")
```

### Biocmanager install
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```
_DESEq2_
```
BiocManager::install("DESeq2")
```
_EdgeR_
```
BiocManager::install("EdgeR")
```
_Limma-Voom_
```
BiocManager::install("limma")
```

### Devtools install
install.packages("devtools")
_SARTools_
install_github("PF2-pasteur-fr/SARTools", build_opts="--no-resave-data")

## Code repositories

## Data-sources

### Virus infection studies
- SARS-MERS: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56192  
- COVID: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507 
- Murine coronavirus (M-CoV): https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-4111/
### Increased risk factors studies


## References

