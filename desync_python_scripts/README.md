# Single-plant-omics reveals the cascade of transcriptional changes during the vegetative-to-reproductive transition - Python scripts

This repository contains the code used to train elastic net models on the gene expression, variant, and physiology data sets.

## Folder structure

The top folder includes this README and a Conda environment file `flowering_variability_env.yml`, which can be used to reproduce the conda environment used for this project.

- `pheno_prediction` contains: the two Python files used to train the models on gene expression data and variant data; the input data files; and outputs from the training procedure. (Note that some of the outputs contain the suffix `_smaller_save`, since they have been reduced in size by the analysis jupyter notebook).
- `analysis` contains: a Jupyter notebook to create figures; and, all of the figures and supplemental data produced by the analysis.

## Reproducing the Conda environment

Please see [this guide](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#activating-an-environment) for how to use the `.yml` file.

## Training the models 

Data and scripts are included to train the models. They can be found in the `pheno_prediction` folder. Make sure your working directory is set to the `pheno_prediction` folder before running. 

To train models on the gene expression data, *with varying values for `l1_ratio`*, you should run the command `python train_elastic_nets.py`. To train models with *fixed* `l1_ratio` values, you should run `python train_elastic_nets.py -f`. Note, the fixed values for l1_ratio have been hard coded, so they will need to be edited in the code if you want to change them. The training parallelizes the different LOOCV folds.

To train models on variants, run `python train_elastic_nets_variants.py`.

For the paper, these commands were run on the [Viking computing cluster](https://vikingdocs.york.ac.uk/). For the most time consuming stage (running models on gene expression with varying `l1_ratio`), we used 16 cores and approximately 350GB of RAM, for a toal run time of about 3 hours.

## Analysing the outputs

A Jupyter notebook has been included to create the plots for the paper. After installing the conda environment, you should be able to run the command `jupyter notebook` and open this folder in your web browser. Make sure the environment has access to both the `pheno_prediction` and `analysis` folders. 

*Note: the logistic regression outputs have been reduced in size by running `remove_coefs_from_dict` in `analysis_elastic_nets.ipynb`, which removes the large `coefs_paths_` attribute.*
