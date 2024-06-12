# Summary of scripts

Scripts in this folder are numbered by steps in the analysis.

### List of scripts scripts

-   `0_load_packages_and_functions.R` loads the packages and custom functions used throughout
-   `1_basic_data_processing_and_filtering.R` loads in the RNA-seq and phenotypic data, produces plots explaining the filtering of genes / plants, procceses data from `biomaRt`, saves neccesary outputs of the above
-   `2_DGE.R` performs differential gene expression, GO term over-representation and produces heatmaps of relevant results
-   `3_physiology_analysis.R` produces summary plots of the physiology data, and compares DGE results to Elastic Net models
-   `4_linear_pseudotime.R` performs the bootstrapping to calculate consensus pseudotime
-   `5_smoothing_pseudotime.R` fits B-spline curves to genes over pseudotime and filters for appropriately smooth genes
-   `6_pseudotime_analysis.R` calculates AUC values and related results
-   `7_network_analysis.R` calculates GRNs and appropriate benchmarks
-   `8_variant_calling_analysis.R` processes the variants and relates them to pseudotime and physiology
-   `9_ape_phylogeny.R` uses the `ape` package to produce a bootstrapped neighbour-joining tree based on variants
-   `10_dtwclustering.R` uses the `DTWClust` package to perform clustering on the gene expression data over pseudotime
-   `11_SNPrelate_phylogeny.R` uses `SNPRelate` to test for variant subgroups using a different method
