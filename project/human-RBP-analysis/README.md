# Project 2: Interaction of the SARS-CoV-2 genome with human RNA-binding proteins

## Workflow:
1. Download PWMs for human RNA-binding proteins
2. Download the SARS-CoV-2 reference genome
3. Scan the viral genome for these PWMs.
4. Annotate potential RBP binding sites on the viral genome
5. Calculate enrichment of PWMs in the viral genome and its 3' and 5' UTRs.
6. Filter significantly enriched RBPs.
7. Repeat steps 2-6 for related coronaviruses.
8. Compare significant hits across species.
9. Check for expression/activity of these RBPs in relevant cell types
10. Check whether the RBP binding sites are modified in different isolates of SARS-CoV-2, or in related viruses.

## Task tracker: 
https://github.com/avantikalal/covid-gene-expression/projects/2

## Data sources:
1. ATtRACT, a database of RNA-binding proteins and PWMs: https://attract.cnic.es
2. SARS-CoV-2 reference genome: https://www.ncbi.nlm.nih.gov/nuccore/NC_045512
3. SARS-CoV-1 reference genome: https://www.ncbi.nlm.nih.gov/nuccore/NC_004718.3
4. RaTG13 reference genome: https://www.ncbi.nlm.nih.gov/nuccore/MN996532

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
1. `01_filter_ATtRACT.R`: Filters the entries in the ATtRACT database to high-quality RBPs and PWMs
2. `02_find_binding_sites.R`: Identifies RBP-binding sites on the SARS-CoV-2 genome. 
3. `03_annotate_sites.R`: Annotates the RBP-binding sites identified in step 1.
4. `04_simulate_genomes.R`: Generates scrambled viral genome and UTR sequences for a background model.
5. `05_find_sim_binding_sites.R`: Identifies RBP-binding sites on the scrambled sequences.
6. `06_enrichment_analysis.R`: Calculates enrichment of binding sites for each RBP in the viral genome and UTRs.
7. `07_filter_hits.R`: Filters significantly enriched RBPs.
8. `08_sars_analysis.R`: Analysis of SARS-CoV-1 genome
9. `09_bat_coronavirus_analysis.R`: Analysis of RaTG13 genome
10. `10_compare_hits.R`: Cross-species comparison

## Output files
Running all the scripts in order produces the following output files in folder `output`:

### Predicted binding sites
1. `sites.RData`: Predicted binding sites on the SARS-CoV-2 genome
2. `annotated_sites.RData`: Predicted binding sites on the SARS-CoV-2 genome, annotated with genetic element (UTR, gene)

### Enrichment files
1. `site_count.RData`: counts and p-values for RBP binding to the SARS-CoV-2 genome
2. `f_utr_site_count.RData`: counts and p-values for RBP binding to the SARS-CoV-2 genome 5'UTR
3. `t_utr_site_count.RData`: counts and p-values for RBP binding to the SARS-CoV-2 genome 3'UTR

### Filtered significant hits
1. `genome_hits.tsv`: RBPs whose binding sites are significantly enriched in the SARS-CoV-2 genome
2. `f_utr_hits.tsv`: RBPs whose binding sites are significantly enriched in the 5'UTR
3. `t_utr_hits.tsv`: RBPs whose binding sites are significantly enriched in the 3'UTR
4. `genome_neg_hits.tsv`: RBPs whose binding sites are significantly enriched in the negative strand of the SARS-CoV-2 genome

Similarly named results files for SARS-CoV-1 and RaTG13 are produced in folders `output/NC_004718.3` and `output/MN996532`.

### Dockerfile
A docker file with all of the R dependencies is available at [hpobiolab/rbp-pwm-r](https://hub.docker.com/orgs/hpobiolab/repositories)

## Data and results
Input data are at the following path:
```
/home/shared/proj_2_rbp_data
```
Current output files are at the following path:
```
/home/shared/proj_2_rbp_output_0410
```
