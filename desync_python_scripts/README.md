# Python scripts for *Single-plant-omics reveals the cascade of transcriptional changes during the vegetative-to-reproductive transition*

This repository contains the code to run the analysis and create figures for the elastic net section of the project.

## Installation

This folder uses `conda` to ensure reproducibility with Python packages. The `conda` environment can be reproduced from `flowering_variability_env.yml`.

## Folder structure

The top folder contains 3 files:

-   `flowering_variability_env.yml`, a YAML file which can be used to [recreate the Conda environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
-   `train_elastic_nets.py`, a Python script. This produces the outputs in `elastic_net_results` from the `data` folder. Note that this needs to be run on a high-performance cluster due to the tiered cross validation.
-   `analysis_elastic_nets.ipynb`, a Jupyter notebook. This takes results in `elastic_net_results` and data in `data` as inputs.

There are also the following subfolders:

-   `data` contains the input data required for the analysis (including outputs from R code)
-   `plots` contains all the plots produced in the analysis, with separate subfolders for categories of plots
-   `elastic_net_results` contains the outputs from `train_elastic_nets.py`, which was run on a high performance compute cluster (Viking, University of York)
-   `gene_regulator_lists` contains the coefficients from the trained elastic net models
-   `hyperparam_lists` contains the hyperparameters chosen cross-validation, per sample
-   `external_scripts` contains published code that cannot be installed using renv
-   `images` contains images used in the Jupyter notebooks