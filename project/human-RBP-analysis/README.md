# Project 2: Interaction of the SARS-CoV-2 genome with human RNA-binding proteins

## Workflow:
1. Download a database of human RNA-binding proteins and their PWMs
2. Download the viral reference genome
3. Scan the viral genome for these PWMs.
4. Calculate enrichment of PWMs in the viral genome.
5. Annotate potential RBP binding sites on the viral genome
6. Check for expression/activity of these RBPs in relevant cell types
7. Check whether any of the RBP binding sites are modified in different isolates of SARS-CoV-2, or in closely related viruses, e.g. the bat coronavirus sequences similar to SARS-CoV-2.

## Task tracker: 
https://github.com/avantikalal/covid-gene-expression/projects/2

## Data sources:
1. ATtRACT, a database of RNA-binding proteins and PWMs: https://attract.cnic.es
2. SARS-CoV-2 reference genome: https://www.ncbi.nlm.nih.gov/nuccore/NC_045512

## Data:
See the wiki for link to the latest data. This contains the reference genome and data on RNA-binding proteins and PWMs from ATtRACT.

## Identifying binding sites
See `find_binding_sites.R` for a basic script to find RBP-binding sites on the SARS-CoV-2 genome. 

## Requirements
The script is in R and uses the following libraries:
1. data.table
2. TFBSTools
3. seqinr
4. Biostrings

### Dockerfile
A docker file with all of the R dependencies is available at [hpobiolab/rbp-pwm-r](https://hub.docker.com/orgs/hpobiolab/repositories)

## Output
The data folder linked in the wiki page also contains a file `sites.RData`. This is the table of potential binding sites produced by `find_binding_sites.R`. To open this, type in an R terminal:
```
load("sites.RData")
```
