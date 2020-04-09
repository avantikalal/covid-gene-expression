# Project 2: Interaction of the SARS-CoV-2 genome with human RNA-binding proteins

## Workflow:
1. Download PWMs for human RNA-binding proteins
2. Download the SARS-CoV-2 reference genome
3. Scan the viral genome for these PWMs.
4. Calculate enrichment of PWMs in the viral genome and its 3' and 5' UTRs.
5. Annotate potential RBP binding sites on the viral genome
6. Check for expression/activity of these RBPs in relevant cell types
7. Check whether the RBP binding sites are modified in different isolates of SARS-CoV-2, or in related viruses.

## Task tracker: 
https://github.com/avantikalal/covid-gene-expression/projects/2

## Data sources:
1. ATtRACT, a database of RNA-binding proteins and PWMs: https://attract.cnic.es
2. SARS-CoV-2 reference genome: https://www.ncbi.nlm.nih.gov/nuccore/NC_045512

## Requirements
The R scripts for this folder use the following libraries:
1. data.table
2. TFBSTools
3. seqinr
4. Biostrings
5. foreach
6. doParallel
7. rtracklayer

## Scripts (to be run in the order below):
1. `find_binding_sites.R`: Identifies RBP-binding sites on the SARS-CoV-2 genome. 
2. `annotate_sites.R`: Annotates the RBP-binding sites identified in step 1.
3. `simulate_genomes.R`: Generates scrambled viral genome and UTR sequences for a background model.
4. `find_sim_binding_sites.R`: Identifies RBP-binding sites on the scrambled sequences.
5. `enrichment_analysis.R`: Calculates enrichment of binding sites for each RBP in the viral genome and UTRs.

## Output files
Running all the scripts in order produces the following output files:

### Predicted binding sites
1. `sites.RData`: Predicted binding sites on the SARS-CoV-2 genome
2. `annotated_sites.RData`: Predicted binding sites on the SARS-CoV-2 genome, annotated with genetic element (UTR, gene)

### Enrichment files
1. `site_count.RData`: counts and p-values for RBP binding to the SARS-CoV-2 genome
2. `f_utr_site_count.RData`: counts and p-values for RBP binding to the SARS-CoV-2 genome 5'UTR
3. `t_utr_site_count.RData`: counts and p-values for RBP binding to the SARS-CoV-2 genome 3'UTR

### Dockerfile
A docker file with all of the R dependencies is available at [hpobiolab/rbp-pwm-r](https://hub.docker.com/orgs/hpobiolab/repositories)

## Data and results
Input data (ATtRACT db and reference genome) can be found at `/home/shared/proj_2_rbp_data`.
Current output files (the 5 files listed above) can be found at `/home/shared/proj_2_rbp_output`
