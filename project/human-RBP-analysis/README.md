# Project 2: Interaction of the SARS-CoV-2 genome with human RNA-binding proteins

## Workflow:
1. Download a database of human RNA-binding proteins and their PWMs
2. Download the viral reference genome
3. Scan the viral genome for these PWMs.
4. Annotate potential RBP binding sites on the viral genome
5. Check for expression/activity of these RBPs in relevant cell types

## Task tracker: 
https://github.com/avantikalal/covid-gene-expression/projects/2

## Data sources:
1. RBPDB, a database of RNA-binding proteins and PWMs: http://rbpdb.ccbr.utoronto.ca/index.php
2. SARS-CoV-2 reference genome: https://www.ncbi.nlm.nih.gov/nuccore/NC_045512

## Data download:
1. Download the reference genome for SARS-CoV-2 from NCBI: https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.
2. Download data on human RNA-binding proteins from RBPDB: http://rbpdb.ccbr.utoronto.ca/downloads/RBPDB_v1.3.1_human_2012-11-21_TDT.zip. Unzip the downloaded directory to get the data tables.
3. Download PWMs (Position Weight Matrices) for these human RNA-binding proteins from RBPDB: http://rbpdb.ccbr.utoronto.ca/downloads/matrices_human.zip. This download contains both PWMs and PFMs. PWM files are in the `PWMDir` folder.

## Identifying binding sites
See `find_binding_sites.R` for a basic script to find RBP-binding sites on the SARS-CoV-2 genome. To run this script you need to modify the first section with paths to your downloaded data files.

## Requirements
The script is in R and uses the following libraries:
1. data.table
2. TFBSTools
3. seqinr
4. Biostrings

## Output
`sitesets.RData` is the table of potential binding sites produced by `find_binding_sites.R`. To open this, type in an R terminal:
```
load("sitesets.RData")
```
