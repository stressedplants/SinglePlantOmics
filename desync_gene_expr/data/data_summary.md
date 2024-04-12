# Summary of data

## Subfolders

-   `pheno` contains the physiological data (i.e. biomass and leaf size).
-   `tair10_salmon_quants` contains the quantifications from Salmon
    against the TAIR10 transcriptome.
-   `dap_seq` contains the processed DAP-seq data from (O'Malley, Huang,
    et. al., Cell, 2016). Specifically, the putative regulatory targets
    (with FRiP \> 5%) were downloaded for all sequenced TF binding
    sites. This was retrieved from
    <http://neomorph.salk.edu/dap_data_v4/fullset/dap_download_may2016_genes.zip>
    on 01/06/2023, via
    <http://neomorph.salk.edu/dap_web/pages/browse_table_aj.php>.
    and contains information on the TF families.
-   `biomaRt_data` contains data about GO term annotations and gene names,
    downloaded using the R package `biomaRt`. Note: these files were downloaded
    in February 2024.
-   `gProfiler_data` contains data about GO term over-representation, downloaded
    using the R package `gprofiler2` in February 2024. This called the g:Profiler 
    server (version e111_eg58_p18_30541362), which utilises the g:SCS multiple 
    testing correction method, and we then applied a significance threshold of 
    1e-4 (Kolberg et al., 2023).

## Individual files

-   `Ath_TF_list.txt` was downloaded in 2023 from 
    <http://planttfdb.gao-lab.org/download/TF_list/Ath_TF_list.txt.gz> and contains
    information on the TF families.
-   `corrected_gene_names.csv` includes manually changed gene names for some
    regulators in the gene regulatory network.
-   `go_size_thresholds.csv` contains the 'medium' and 'large' thresholds for sizes
    of GO terms, after pseudotime filtering.
-   `\PP2015-01929R1_Supplemental_Table_S3.xls` was reproduced from 
    <https://doi.org/10.1104/pp.15.01929>. It is used to compare the pseudotime 
    analysis with the mature-to-senscent phase of Woo et al.'s leaf development
    RNA-seq time series.

