# Single-plant-omics reveals the cascade of transcriptional changes during the vegetative-to-reproductive transition - R code

This repository contains the code to run the analysis and create figures for this project. 

For additional details see: `data/data_summary.md` for the sources and explanation of the data; and `R/script_summary.md` for an explanation of the R scripts.

## Installation

This project uses [`renv`](https://rstudio.github.io/renv/articles/renv.html). `renv` should automatically set up and download required packages when you first load this project.

Follow [this guide](https://rstudio.github.io/renv/articles/collaborating.html) for more instructions.

## Folder structure

-   `R` contains the code to perform analysis
-   `data` contains the input data required for the analysis (including quantifications from Salmon)
-   `plots` contains all the plots produced in the analysis
-   `outputs` contains all outputs from the code that aren't plots
-   `renv` contains the `activate.R` script required for install
-   `external_scripts` contains published code that cannot be installed using renv
