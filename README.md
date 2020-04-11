# Covid-19 Gene Expression Working Group

[![code_license](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![data_license](https://img.shields.io/badge/License-CC0%201.0-lightgrey.svg)](http://creativecommons.org/publicdomain/zero/1.0/)

Hackathon team: gene expression analysis for Covid-19 Virtual Biohackathon (vBH)

### Link to vBH Github
https://github.com/virtual-biohackathons/covid-19-bh20
### Official Gene Expression Work Group Page
https://github.com/virtual-biohackathons/covid-19-bh20/wiki/GeneExpression

## Main Objective
As part of the virtual BioHackathon 2020, we formed a working group that focused on the analysis of gene expression in the context of COVID-19. More specifically, we performed transcriptome analyses on published datasets in order to better understand the interaction between the human host and the SARS-CoV-2 virus.

<figure>
  <img src="https://github.com/avantikalal/covid-gene-expression/blob/mariferrarini-patch-11/Diagram_projects.png" width="800">
  <figcaption>Fig. 1. Project structure and interaction. Projects 1 and 2 along with literature research will provide a list of candidate genes for Project 3 and 4 that will take into account external factors (comorbidities, and potential drug treatments). Project 5 will implement the findings into electronic medical records. Finally, Project 6 allows the efforts of all of the previous projects to be clearly detailed into workflows for increased reproducibility.</figcaption>
</figure>

## Deliverables
_Biological:_ Perform a global RNA-Seq analysis with SARS-CoV-2 infected datasets to search for new candidate genes for testing experimentally

_Methodological:_ Create a packaged reproducible pipeline in Docker to help scientists to easily treat their RNA-Seq data and for us if any new dataset comes out

## Projects
The ideas proposed during this hackathon were divided into six projects (Fig. 1):

1. SARS-CoV-2 infection global analyses: Understanding how global gene expression in human cells responds to infection by the SARS-CoV-2 virus, including changes in gene regulatory networks.

2. Human-virus interaction analyses: Identification of human RNA-binding proteins that might be key in the interaction between human cells and the RNA genome of SARS-CoV-2.

3. Increased risk factors analyses: Investigating gene expression in other datasets with the goal of identifying commonalities and differences with the two previous analyses, focusing on specific genes.

4. Subtyping of expression response to drugs after COVID infection: Searching for potential drugs that could impact the expression of human genes that are important for the interaction of human and virus.

5. Reporting findings to electronic medical records: Determining how to bring the results obtained from the previous analyses to clinical practice.

6. Workflows for reproducibility of analysis: Packaging the workflows devised within the Gene Expression group to enable seamless integration and approach reproducibility.

Please see the [project](project) folder for details on each individual project.

## Contributors
See the [contributors table](contributors.md) for a full list of the amazing people who have 
contributed to the project.

## Progress tracking
We are tracking progress on project-specific boards here: https://github.com/avantikalal/covid-gene-expression/projects

## Licensing
This working group is dedicated to Open Science. All code in this repository is licensed under the
[MIT license](LICENSE.md). Data generated during the course of this project is licensed under the
[CC0 license](DATA_LICENSE.md).
