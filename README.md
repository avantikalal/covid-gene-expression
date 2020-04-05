# Covid-19 Gene Expression Work Group (Private project)
Hackathon team: gene expression analysis for Covid-19 Virtual Biohackathon (vBH)

### Link to vBH Github
https://github.com/virtual-biohackathons/covid-19-bh20
### Official Gene Expression Work Group Page (Public page)
https://github.com/virtual-biohackathons/covid-19-bh20/wiki/GeneExpression

## Main Objective
We want to perform RNAseq-based analyses on published datasets in order to better understand the interaction between human host and virus.

<figure>
  <img src="https://github.com/avantikalal/covid-gene-expression/blob/master/Diagram.png" width="200">
  <figcaption>Fig.1 - We want to focus on already known genes (such as ACE2 and TMPRSS2) but also on new candidate ones that may play a punctual or a general role in the interaction between host and virus. To this end, we will perform extensive RNAseq analyses as described in the workflow section below</figcaption>
</figure>

## Deliverables
_Biological:_ Perform a global RNA-Seq analysis with SARS-CoV-2 infected datasets to search for new candidate genes for testing experimentally

_Methodological:_ Create a packaged reproducible pipeline in Docker to help scientists to easily treat their RNA-Seq data and for us if any new dataset comes out

We have four main projects going on. To summarize them :

1. SARS-CoV-2 infection global analyses: will include global gene expression, functional and regulation of gene expression analyses.
2. Human-virus interaction analyses : will try to search RNA-binding proteins that might be key in the interaction between human and SARS-CoV-2.
3. Increased risk factors analyses: will look for gene expression data in other datasets to compare especially with selected genes from previous analyses.
4. Subtyping of expression response to drugs after COVID infection: will focus on trying to search for potential drugs that could impact in important genes for the interaction of human and virus.
5. Reporting findings to electronic medical records : will try to make these findings arrive at the hands of clinicians.

## Progress tracking
We are tracking progress on project-specific boards here: https://github.com/avantikalal/covid-gene-expression/projects

# Software Used during this project

## Command line tools
- FastQC (https://github.com/s-andrews/FastQC)
- Fastq-screen (https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
- trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic)
- STAR (https://github.com/alexdobin/STAR)
- Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- HISAT2 (https://daehwankimlab.github.io/hisat2/)
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

### BiocManager install
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
```
install.packages("devtools")
```
_SARTools_
```
library(devtools)
install_github("PF2-pasteur-fr/SARTools", build_opts="--no-resave-data")
```

## Code repositories


## Data-sources


### Virus infection studies
- SARS-MERS: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56192  
- COVID: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507 
- Murine coronavirus (M-CoV): https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-4111/
### Increased risk factors studies


## References

