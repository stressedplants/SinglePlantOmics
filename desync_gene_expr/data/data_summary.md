---
editor_options: 
  markdown: 
    wrap: 72
---

# Summary of data files

-   `pheno` contains the physiological data (i.e. biomass and leaf size)
    for wild type and mutants.
-   `tair10_salmon_quants` contains the quantifications from Salmon
    against the TAIR10 transcriptome.
-   `dap_seq` contains the processed DAP-seq data from (O'Malley, Huang,
    et. al., Cell, 2016). Specifically, the putative regulatory targets
    (with FRiP \> 5%) were downloaded for all sequenced TF binding
    sites. This was retrieved from
    <http://neomorph.salk.edu/dap_data_v4/fullset/dap_download_may2016_genes.zip>
    on 01/06/2023, via
    <http://neomorph.salk.edu/dap_web/pages/browse_table_aj.php>.
-   `Ath_TF_list.txt` was downloaded in 2023 from
    <http://planttfdb.gao-lab.org/download/TF_list/Ath_TF_list.txt.gz>
    and contains information on the TF families.
-   `biomaRt_data` contains data about GO term annotations and gene names,
    downloaded using the R package `biomaRt`. Note: the regulatory genes were
    downloaded in February 2023 and the other files were downloaded in June 2023.
-   `gProfiler_data` contains data about GO term over-representation, downloaded
    using the R package `gprofiler2` in August 2023.
