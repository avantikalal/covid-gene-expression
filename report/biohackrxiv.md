---
title: 'BioHackrXiv contribution template'
tags:
  - BioHackathon
  - Preprint
  - Gene expression
  - RNA-seq
  - COVID-19
authors:
  - name:Mariana G. Ferrarini
    orcid: https://orcid.org/0000-0002-9574-9991
    affiliation: 1, *
  - name: Vanessa Aguiar-Pulido
    affiliation: 2, #
  - name: Eric T. Dawson
    affiliation: 3, #
  - name: Andrea Guarracino
    orcid: https://orcid.org/0000-0001-9744-131X
    affiliation: 4, #
  - name: Andreas Gruber
    orcid: https://orcid.org/0000-0001-7664-4257
    affiliation: 5, #
  - name: Lukas Heumos
    orcid: https://orcid.org/0000-0002-8937-3457
    affiliation: 6, #
  - name: Alexander Kanitz
    orcid: https://orcid.org/0000-0002-3468-0652
    affiliation: 7, #
  - name: Avantika Lal
    orcid: https://orcid.org/0000-0002-5827-0826
    affiliation: 8, #
  - name: Brett E. Pickett
    orcid: https://orcid.org/0000-0001-7930-8160
    affiliation: 9, #
  - name: Rita Rebollo
    orcid: https://orcid.org/0000-0002-8138-5082
    affiliation: 1, #
  - name: Carlos Ruiz-Arenas
    affiliation: 10, #
  - name: Olaitan Awe
    affiliation: 11
  - name: Suhana Bedi
    affiliation: 12
  - name: Ben Busby
    affiliation: 13
  - name: Marielena Georgaki
    affiliation: 14
  - name: Chela James
    affiliation: 15
  - name: Itziar Martinez Gonzalez
    affiliation: 16
  - name: Birgit Meldal
    affiliation: 17
  - name: Scheila G. Mucha
    affiliation: 18
  - name: Jakke Neiro
    affiliation: 19
  - name: Nuria Queralt Rosinach
    affiliation: 14
  - name: Philippe Rocca-Serra
    affiliation: 20
  - name: Daniel Siqueira de Oliveira
    affiliation: 21
  - name: Maria Tsagiopoulou
    affiliation: 22
affiliations:
 - name: Group coordinator
   index: *
 - name: Project leader
   index: #
 - name: University of Lyon, INSA-Lyon, INRA, BF2I, Villeurbanne, France
   index: 1
 - name: Center for Neurogenetics, Weill Cornell Medicine, Cornell University, New York, NY, USA
   index: 2
 - name: Division of Cancer Epidemiology and Genetics, National Cancer Institute, Rockville Maryland, USA
   index: 3
 - name: Centre for Molecular Bioinformatics, Department of Biology, University Of Rome Tor Vergata, Rome, Italy
   index: 4
 - name: Oxford Big Data Institute, Nuffield Department of Medicine, University of Oxford, Oxford, UK
   index: 5
 - name: Quantitative Biology Center (QBiC), University of Tübingen, Germany
   index: 6
 - name: Swiss Institute of Bioinformatics, ELIXIR Switzerland
   index: 7
 - name: NVIDIA Corporation, Santa Clara, CA, USA
   index: 8
 - name: Brigham Young University, Provo, UT, USA
   index: 9
 - name: Centro de Investigacion Biomedica en Red de Enfermedades Raras, Universitat Pompeu Fabra, Barcelona, Spain
   index: 10
 - name: University of Ibadan, Ibadan, Nigeria
   index: 11
 - name: University of Texas at Dallas, TX, USA
   index: 12
 - name: DNAnexus, Mountain View, CA, USA
   index: 13
 - name: Virtual BioHackathon participant
   index: 14
 - name: Institute of Cancer Research, London, UK
   index: 15
 - name: Amsterdam UMC, Amsterdam, The Netherlands
   index: 16
 - name: European Bioinformatics Institute (EMBL-EBI), European Molecular Biology Laboratory, Wellcome Genome Campus, Hinxton, UK
   index: 17
 - name: INRIA Grenoble Rhone-Alpes, Montbonnot-Saint-Martin, France
   index: 18
 - name: Department of Zoology, University of Oxford, Oxford, UK
   index: 19
 - name: Oxford e-Research Centre, Department of Engineering Science, University of Oxford, Oxford, UK
   index: 20
 - name: Universite Claude Bernard Lyon 1, Villeurbanne, France
   index: 21
 - name: Institute of Applied Biosciences, Centre for Research and Technology Hellas, Thessaloniki, Greece
   index: 22
date: 12 April 2020
bibliography: paper.bib
---

# Introduction

As part of the virtual BioHackathon 2020, we formed a working group that focused on the analysis of gene expression in the context of COVID-19. More specifically, we performed transcriptome analyses on published datasets in order to better understand the interaction between the human host and the SARS-CoV-2 virus.

The ideas proposed during this hackathon were divided into five projects (Fig. 1):

1. SARS-CoV-2 infection global analyses: Understanding how global gene expression in human cells responds to infection by the SARS-CoV-2 virus, including changes in gene regulatory networks.
2. Human-virus interaction analyses: Identification of human RNA-binding proteins that might be key in the interaction between human cells and the RNA genome of SARS-CoV-2.
3. Increased risk factors analyses: Investigating gene expression in other datasets with the goal of identifying commonalities and differences with the two previous analyses, focusing on specific genes.
4. Identification of potential pharmacological treatments: Searching for potential drugs that could impact the expression of human genes that are important for the interaction of human and virus.
5. Workflows for reproducibility of analysis: Packaging the workflows devised within the Gene Expression group to enable seamless integration and approach reproducibility. 

Projects 1 and 2 aim to identify human genes that are important in the process of viral infection of human cells. Projects 3 and 4 aim to take the candidate genes identified in projects 1 and 2, as well as by independent studies, and relate them to clinical information and to possible therapeutic interventions. All data analyzed during this study are fully available and meet the FAIR principles of Findability, Accessibility, Interoperability, and Reusability. Finally, Project 5 aims to package and containerize software and workflows used and generated here in a reusable manner, ultimately providing scalable and reproducible workflows.

<figure>
  <img src="https://github.com/avantikalal/covid-gene-expression/blob/mariferrarini-patch-11/Diagram_projects.png" width="800" align="middle">
  <figcaption>Fig. 1. Project structure and interaction. Project 1 and 2 along with literature research will provide a list of candidate genes for Project 3 and 4 that will take into account external factors (comorbidities, and potential drug treatments). All data analyzed during this project are fully available to the medical community and meet the FAIR principles. Finally, Project 5 allows the efforts of all of the previous projects to be clearly detailed into workflows for increased reproducibility.</figcaption>
</figure>


# Background


# Results


# Discussion


# Conclusions

This working group has focused on establishing methodologies that could be followed by other researchers to perform similar analyses. The aim of this work was to lay the foundation for further research in COVID-19 taking advantage of existing RNA-Seq datasets. Standard protocols have been developed and packaged in containers, allowing scientific reproducibility and scalability. Everything has been made publicly available and linked through the group’s main GitHub, as listed in the pertinent section.


# Future work

Although a lot has been accomplished during the virtual BioHackathon, many avenues can be pursued in the future and several of the proposed projects could be taken to completion. Below is listed some potential future work based on the methodologies and workflows described here.
* Projects 1, 2 and 3 would benefit from the inclusion of additional datasets as these become available. For this purpose, a thorough research of the literature and public databases should be performed on a regular basis. More specifically, for projects 1 and 2, the list of potentially interesting genes obtained could be further refined by including these, increasing the robustness of the final results. 
* A deeper analysis to interpret all the results generated here is required, especially with the help of experts in the fields of virology, immunology, molecular and cellular biology, and clinicians.
* Finally, packaging and developing additional reproducible workflows using containers such as those proposed here would be extremely beneficial for the research community, as many of these can be used to analyze datasets in other domains.


# Methods


## Datasets

Raw sequencing data were downloaded from the Gene Expression Omnibus (GEO) database (see Table XX). These studies compare the transcriptional response of several cell types to different viruses, with different time points and concentrations.

Table **XX**. Summary of the datasets used in this work

|  Study  |  Virus  |  Strain  |  Cell line  |  MOI  |  Time (hpi)  |
|  -----  |  -------  |  -----  |  -----  |  -----  |  -----    |
|  [GSE147507](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507)  |  SARS-CoV-2  |  USA-WA1/2020  |  NHBE  |  2  |  24 h  |
|  [GSE147507](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507)  |  SARS-CoV-2  |  USA-WA1/2020  |  A549  |  0.02  |  24 h  |
|  [GSE147507](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507)  |  RSV  |  A2  |  A549  |  15  |  24 h  |
|  [GSE147507](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507)  |  IAV  |  Puerto Rico/8/1934 (H1N1)  |  A549  |  5  |  9 h  |
|  [GSE122876](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122876)  |  MERS-CoV  |  HCoV-EMC/2012  |  Calu-3  |  2  |  24 h  |
|  [GSE56192](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56192)  |  MERS-CoV  |  HCoV-EMC/2012  |  MRC5  |  0.1  |  24 h  |
|  [GSE56192](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56192)  |  MERS-CoV  |  HCoV-EMC/2012  |  MRC5  |  3  |  24 h  |
|  [GSE56192](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56192)  |  MERS-CoV  |  HCoV-EMC/2012  |  MRC5  |  0.1  |  48 h  |
|  [GSE56192](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56192)  |  MERS-CoV  |  HCoV-EMC/2012  |  MRC5  |  3  |  48 h  |
|  [GSE56192](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56192)  |  SARS  |  Urbani strain  |  MRC5  |  0.1  |  24 h  |
|  [GSE56192](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56192)  |  SARS  |  Urbani strain  |  MRC5  |  3  |  24 h  |
|  [GSE56192](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56192)  |  SARS  |  Urbani strain  |  MRC5  |  0.1  |  48 h  |
|  [GSE56192](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56192)  |  SARS  |  Urbani strain  |  MRC5  |  3  |  48 h |
| [GSE57148](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57148) | N/A | N/A | Primary Lung | N/A | N/A |


# Data, GitHub repositories and reproducible workflows

More detailed information, code, software and containers created as part of this effort is available at: https://github.com/avantikalal/covid-gene-expression. Data generated within this effort is available at: zenodo url. The software developed is registered under the MIT license and the data generated under a CC0 license.


# Acknowledgement

We would like to thank the clinical personnel who are currently making an extraordinary effort to take care of patients in such a challenging time, everyone who has contributed to research in COVID-19 and those who participated in the BioHackathon, especially within the Gene Expression supergroup. We would also like to acknowledge the computational resources shared with us during the hackathon: a CentOS node with shared storage and a virtual Slurm cluster was provided by the IT Center for Science (CSC), ELIXIR Finland, and an AWS S3 bucket was provided by Artem Babaian, University of British Columbia.


# References

(.bib)
