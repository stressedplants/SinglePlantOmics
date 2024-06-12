# Load physiological data =================

print('Loading phsyiological data')

wt_data <- read_excel("./data/pheno/wt_2021_cleaned.xlsx")

colnames(wt_data) <- c("id", "leaf_num", "bolting", "biomass", "leaf_1", 
                       "leaf_2", "notes", "rna", "date_of_extraction", 
                       "conc", "260_280", "260_230", "tray_num")


# Wild type cleaning ===============

# Fix data types - wild type data
# NAs will be produced, but that is expected
wt_data$leaf_1 <- as.numeric(wt_data$leaf_1)
wt_data$leaf_2 <- as.numeric(wt_data$leaf_2)
wt_data$tray_num <- as.factor(wt_data$tray_num)

# Remove samples from wt_data that weren't sequenced
wt_data_seq <- subset(wt_data, id != "P083" & id != "P218" & id != "P115" & 
                      id != "P190" & id != "P124")

rm(wt_data) # only work with sequenced data so this is no problem

wt_data_seq <- wt_data_seq %>% rowwise() %>%
  mutate(combined_pheno = paste(leaf_num, bolting, sep = ""), .after = bolting)
wt_data_seq$combined_pheno <- as.factor(wt_data_seq$combined_pheno)

wt_data_seq <- wt_data_seq %>% rowwise() %>%
  mutate(leaf_avg = mean(c(leaf_1, leaf_2)), .after = leaf_2)

# Order by id
wt_data_seq <- wt_data_seq[order(wt_data_seq$id),]

# Get conditions ONLY df for use below
conds <- as.data.frame(wt_data_seq[,c("id", "leaf_num", "bolting")])
rownames(conds) <- conds$id
conds$leaf_num <- as.factor(conds$leaf_num)
conds$bolting <- as.factor(conds$bolting)

# Input RNA-seq data =============

print('Inputting RNA-seq data')

folder_list <- list.dirs(path = "./data/tair10_salmon_quants",
                         recursive = FALSE)

# Note: currently keep and P128 and P158 initial load, but these are removed
# in the initial processing step

# initalise data table with P002 - just use this for row names
first_fold <- "./data/tair10_salmon_quants/P002_quant/"

# get the first count matrix -> in order to access transcript names
count_raw_init <- read.table(paste(first_fold, "/quant.sf", sep = ""),
                             header = TRUE)

# bit of a hack - treats .d as file extension
gene_only <- tools::file_path_sans_ext(count_raw_init$Name) 

# use this data frame for input to tximport, with tx2gene
trans_to_gene <- data.frame(count_raw_init$Name, gene_only)
colnames(trans_to_gene) <- c("TXNAME", "GENEID")

# create list of files for tximport
# ALSO have P158_rerun
file_list <- sapply(folder_list, function(i){paste0(i, "/quant.sf")})

# Set names so that tximport works correctly 
# Chop off  "./data/tair10_salmon_quants/" from beginning and "_quant" from end

temp_name_list <- str_remove(names(file_list),
                             fixed("./data/tair10_salmon_quants/"))
temp_name_list <- str_remove(temp_name_list, fixed("_quant"))
names(file_list) <- temp_name_list

txi_data_genes <- tximport(file_list, type = "salmon", tx2gene = trans_to_gene)
# txi_data_tx <- tximport(file_list, type = "salmon", txOut = T)

# these were only used in this section
rm(file_list, first_fold, folder_list, count_raw_init, trans_to_gene,
   temp_name_list, gene_only)

# Removal of P128 ========

print('Making plots about P128')

TPM_pre_removal <- as.data.frame(txi_data_genes$abundance)

# Testing correlation of each col - and making a heatmap
correlation_matrix_with_P128 <- apply(log2(TPM_pre_removal + 1), 2, function(i){
  apply(log2(TPM_pre_removal + 1), 2, function(j){
    cor(i,j)
  })
})

# This justifies P128 being removed - since correlation with every other plant
# is very low. 

# Note that P158 and P158_rerun are within the same cluster. Justifies keeping 
# P158_rerun (which has higher read depth so better to use that than P158)
svg("./plots/initial_processing/plant_removal_heatmap.svg",
    width = 16,
    height = 16)
heatmap_plants_with_P128 <- pheatmap(
  correlation_matrix_with_P128,
  color = inferno(10),
  # show_rownames = FALSE, 
  # show_colnames = FALSE,
  cellwidth = 13, cellheight = 13,
  # treeheight_row = 30, treeheight_col = 30,
  fontsize = 15,
  # filename = "./plots/initial_processing/plant_removal_heatmap.png"
)
draw(heatmap_plants_with_P128)
dev.off()

TPM_master <- TPM_pre_removal[,!(names(TPM_pre_removal) == "P128")]

rm(correlation_matrix_with_P128, heatmap_plants_with_P128)

# Thresholds for minimum gene expression =========

# First filtering = remove all genes where <= 10 samples had a TPM > 0.
# Justified since all true 'step' genes will have at least this many non-zero values.
number_of_samples <- ncol(TPM_master)
low_count_boolean_matrix <- TPM_master == 0
low_count_samples <- rowCounts(low_count_boolean_matrix)

# get a boolean vector of whether to REMOVE each row (ie gene)
low_count_positions <- low_count_samples >= number_of_samples - 10

first_low_filtering_genes <- row.names(TPM_master[low_count_positions, ])

# Can do additional filtering by an absolute threshold 
# this plot justifies the threshold
if (!file.exists("./plots/initial_processing/low_plot.svg")) {
  print('Plotting low threshold limits')
  
  low_threshold_lims <- seq(from = 0, to = 5, by = 0.1)
  num_removed_low <- sapply(low_threshold_lims, function(i){
    second_low_filtering_genes <- row.names(TPM_master %>%
                                              filter_all(all_vars(. <= i)))
    return(length(union(second_low_filtering_genes, first_low_filtering_genes)))
  })
  low_df <- data.frame(limits = low_threshold_lims,
                       nums = num_removed_low)
  
  low_plot <- ggplot(low_df, aes(x = limits, y = nums)) +
    geom_point() +
    geom_vline(xintercept = 0.5, color = 'red') +
    theme_minimal(base_size = 18) +
    theme(axis.line = element_line(linewidth = 1, colour = "black")) +
    scale_y_continuous(limits = c(1e4, 2e4)) +
    xlab("Lower threshold limit") + 
    ylab("Number of genes removed from further analysis")
  
  ggsave("./plots/initial_processing/low_plot.svg", low_plot,
         width = 10, height = 7)
}

final_low_removed_genes <- row.names(TPM_master %>%
                                     filter_all(all_vars(. <= 0.5)))
final_low_removed_genes <- union(final_low_removed_genes,
                                 first_low_filtering_genes)

print(paste0("Total number of genes removed by the low expression criteria: ",
             length(final_low_removed_genes)))

# Save lists of the genes which have been removed...
TPM_low_removed <- TPM_master[final_low_removed_genes, ]

# Basic filtering based on these threshold values =====

print('Producing final filtered expression values')

genes_to_keep <- setdiff(row.names(TPM_master), row.names(TPM_low_removed))

TPM_filtered <- TPM_master[genes_to_keep, ]

print(paste0("Number of remaining genes after filtering: ",
             nrow(TPM_filtered)))

# # Hierarchical clustering plot =========
# 
# plant_hclust <- hclust(dist(t(TPM_filtered)))
# 
# svg(filename = './plots/initial_processing/hclust_all.svg', 
#     width = 12, height = 7)
# plot(plant_hclust, xlab = 'Plant ID', sub = "")
# dev.off()
# 
# rm(plant_hclust)

# Remove variables related to low / high gene expression ======

rm(number_of_samples,
   low_count_boolean_matrix,
   low_count_samples,
   low_count_positions,
   first_low_filtering_genes)

# PCA plot - justification to remove P196 and P158 =========

pca_pre_P196 <- prcomp(t(log2(TPM_filtered + 1)))
plant_pca_summary <- summary(pca_pre_P196)$importance

first_pc <- 1
second_pc <- 2

pca_plant_for_plot <- data.frame(pca_pre_P196$x[, c(first_pc, second_pc)])

# modify rownames because of P158_rerun
modified_rownames <- row.names(pca_plant_for_plot)
modified_rownames[modified_rownames == "P158_rerun"] <- "P158"

pca_plant_for_plot$bolting <- conds[modified_rownames, "bolting"]

pca_filtering <- ggplot(pca_plant_for_plot, aes(x = PC1, y = PC2,
                                                color = bolting)) + 
  geom_point() +
  scale_color_manual(values = c("N" = no_col,
                                "Y" = yes_col)) +
  xlab(paste0("PC", first_pc, " (", 
              sprintf(plant_pca_summary[2, first_pc]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC", second_pc," (",
              sprintf(plant_pca_summary[2, second_pc]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  geom_text_repel(label = rownames(pca_plant_for_plot),
                  max.overlaps = 20,
                  seed = 1234) +
  theme_classic(base_size = 12)

ggsave("./plots/initial_processing/pca_pre_196.svg", pca_filtering,
       width = 10, height = 7)

# Remove P196 and (original) P158 ====

TPM_filtered <- TPM_filtered[
  , !(colnames(TPM_filtered) %in% c("P158", "P196"))
]

colnames(TPM_filtered)[which(colnames(TPM_filtered) == "P158_rerun")] <- "P158"

# Ensure genes meet the filtering criteria after removal of samples ====

# First filtering = remove all genes where <= 10 samples had a TPM > 0.
# Justified since all true 'step' genes will have at least this many non-zero values.
number_of_samples <- ncol(TPM_filtered)
low_count_boolean_matrix <- TPM_filtered == 0
low_count_samples <- rowCounts(low_count_boolean_matrix)

# get a boolean vector of whether to REMOVE each row (ie gene)
low_count_positions <- low_count_samples >= number_of_samples - 10
first_low_filtering_genes <- row.names(TPM_master[low_count_positions, ])

final_low_removed_genes <- row.names(TPM_filtered %>%
                                       filter_all(all_vars(. <= 0.5)))
final_low_removed_genes <- union(final_low_removed_genes,
                                 first_low_filtering_genes)
TPM_filtered <- TPM_filtered[
  !(row.names(TPM_filtered) %in% final_low_removed_genes), 
]

print(paste0("Final number of genes after removing P128 and P196: ",
             nrow(TPM_filtered)))

# Repeat PCA after this filtering ====

pca_post_P196 <- prcomp(t(log2(TPM_filtered + 1)))
plant_pca_summary <- summary(pca_post_P196)$importance

first_pc <- 1
second_pc <- 2

pca_plant_for_plot <- data.frame(pca_post_P196$x[, c(first_pc, second_pc)])
pca_plant_for_plot$bolting <- conds[row.names(pca_plant_for_plot), "bolting"]

pca_filtering <- ggplot(pca_plant_for_plot, aes(x = PC1, y = PC2,
                                                color = bolting)) + 
  geom_point() +
  scale_color_manual(values = c("N" = no_col,
                                "Y" = yes_col)) +
  xlab(paste0("PC", first_pc, " (", 
              sprintf(plant_pca_summary[2, first_pc]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC", second_pc," (",
              sprintf(plant_pca_summary[2, second_pc]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  geom_text_repel(label = rownames(pca_plant_for_plot),
                  max.overlaps = 20,
                  seed = 1234) +
  theme_classic(base_size = 12)

ggsave("./plots/initial_processing/pca_post_196.svg", pca_filtering,
       width = 10, height = 7)

# Save the values for prediction =============

print("Outputting data related to phenotype prediction")

phenos_to_predict <- as.data.frame(
  wt_data_seq[wt_data_seq$id %in% colnames(TPM_filtered),
              c("id", "leaf_num", "bolting",
                "leaf_avg", "biomass")])
rownames(phenos_to_predict) <- phenos_to_predict$id
phenos_to_predict <- within(phenos_to_predict, rm(id))

write.csv(phenos_to_predict, file = "outputs/tpm_tables/phenos_to_predict.csv",
          quote = FALSE)

TPMs_filtered_z_scored <- get_z_score_matrix(TPM_filtered)

write.csv(TPMs_filtered_z_scored, file = "outputs/tpm_tables/TPM_z_scored.csv",
          quote = FALSE)

# Write other values for supplementary figures =========

write.csv(TPM_master,
          "outputs/tpm_tables/TPM_master.csv")

write.csv(TPM_filtered,
          "outputs/tpm_tables/TPM_filtered.csv")

# Plant similarity heatmap =======

print('Plotting similarity of plants AFTER all filtering')

conds_filtered <- data.frame(conds[(conds$id %in% colnames(TPM_filtered)),])
id_temp <- conds_filtered$id

anno_col <- data.frame(conds_filtered[,"bolting"])
colnames(anno_col) <- c("Bolting status")
rownames(anno_col) <- id_temp

anno_colours <- list("Bolting status" = c(N = no_col, Y = yes_col))

correlation_matrix <- apply(log2(TPM_filtered + 1), 2, function(i){
  apply(log2(TPM_filtered + 1), 2, function(j){
    cor(i,j)
  })
})

svg("./plots/initial_processing/plant_heatmap.svg",
    width = 16,
    height = 16)
heatmap_plants_cor <- pheatmap(
  correlation_matrix,
  color = inferno(10),
  annotation_col = anno_col,
  annotation_colors = anno_colours,
  annotation_row = anno_col,
  # annotation_legend = FALSE,
  # show_rownames = FALSE,
  cellwidth = 13, cellheight = 13,
  # treeheight_row = 30, treeheight_col = 30,
  fontsize = 15,
  # filename = "./plots/initial_processing/plant_heatmap.png"
)
draw(heatmap_plants_cor)
dev.off()

rm(id_temp, correlation_matrix, heatmap_plants_cor)

# Download the gene names of all genes included after filtering =========

# Load the file if it already exists - to be consistent across runs
biomart_data_file <- "./data/biomaRt_data/biomart_data.Rdata"

if (!file.exists(biomart_data_file)) {
  print("Downloading gene names and GO terms from TAIR")
  
  # BiomaRt set up - only need to run once and then save
  ensembl_arabidopsis <- useEnsemblGenomes(biomart = "plants_mart",
                                           dataset = "athaliana_eg_gene")
  
  
  ## Use these functions to find what attributes and filters you can find
  # test_atr <- listAttributes(ensembl_arabidopsis)
  # test_filt <- listFilters(ensembl_arabidopsis)
  attributes_to_retrieve <- c('external_gene_name', 'external_gene_source',
                              'external_transcript_name', 'external_transcript_source_name',
                              'external_synonym',
                              'tair_locus', 
                              'go_id', 'name_1006', 'definition_1006',
                              'go_linkage_type', 'namespace_1003')
  
  genes_after_filtering <- row.names(TPM_filtered)
  
  BM_after_filtering <- getBM(attributes = attributes_to_retrieve,
                              filters = 'tair_locus',
                              values = genes_after_filtering,
                              mart = ensembl_arabidopsis)
  
  save(BM_after_filtering,
       file = biomart_data_file)
} else {
  print("Loading biomart data from file")
  load(biomart_data_file)
}

# Produce a named vector of go term names and namespaces, labelled by ID ====

# some GO terms do not have names attached - just remove these for now
# but could manually add these later
GO_names_df <- unique(
  BM_after_filtering[, c("go_id", "name_1006", "namespace_1003")]
) %>% filter(go_id != "") %>% filter(name_1006 != "")
GO_names_df$namespace_1003 <- as.factor(GO_names_df$namespace_1003)

row.names(GO_names_df) <- GO_names_df$go_id

# Produce a data frame linking TAIR ID to gene names ========

# Load if already exists

id_name_file <- "data/biomaRt_data/id_to_name_df.csv"

if (!file.exists(id_name_file)) {
  
  print("Creating gene ID to name data frame")
  id_to_name_df <- data.frame(matrix(nrow = length(genes_after_filtering),
                                     ncol = 2))
  row.names(id_to_name_df) <- genes_after_filtering
  names(id_to_name_df) <- c('gene_name', 'combined_name_id')
  
  for (gene_id in genes_after_filtering) {
    BM_only_one_gene <- BM_after_filtering[
      BM_after_filtering$tair_locus == gene_id,
    ]
    
    if (dim(BM_only_one_gene)[1] == 0) {
      # then we can only keep the gene id, since gene is not in BM table
      id_to_name_df[gene_id ,1] <- gene_id
      id_to_name_df[gene_id ,2] <- gene_id
    } else {
      # if there is a non-empty gene name, use this
      poss_names <- unique(BM_only_one_gene$external_gene_name)
      if ((length(poss_names) == 1) && (poss_names == "")) {
        id_to_name_df[gene_id ,1] <- gene_id
        id_to_name_df[gene_id ,2] <- gene_id
      } else {
        # cycle through poss_names and choose first non empty one
        temp_names <- poss_names[poss_names != ""]
        if (length(temp_names) > 0 ) {
          id_to_name_df[gene_id ,1] <- temp_names[1]
          id_to_name_df[gene_id ,2] <- paste0(temp_names[1], ' (', gene_id, ')')
        } else {
          id_to_name_df[gene_id ,1] <- gene_id
          id_to_name_df[gene_id ,2] <- gene_id
        }
      }
    }
  }
  
  # save the id to name data frame
  write.csv(id_to_name_df,
            file = id_name_file,
            quote = FALSE)
} else {
  print("Loading gene ID to name data frame")
  id_to_name_df <- read.csv(
    file = id_name_file,
    row.names = 1
  )
}

# Filter for regulators =========

# Load regulatory genes from a file, to ensure consistency with the Feb 23 database

regulator_genes_file <- "data/biomaRt_data/biomaRt_regulatory_genes.txt"

if (!file.exists(regulator_genes_file)) {
  list_of_regulatory_GO_terms <- c("GO:0003700", "GO:0004871", "GO:0006355")
  
  BM_for_regulators <- BM_after_filtering[
    BM_after_filtering$go_id %in% list_of_regulatory_GO_terms, 
  ]
  
  regulator_genes <- unique(BM_for_regulators$tair_locus)
  regulator_genes <- regulator_genes[regulator_genes %in% row.names(TPM_filtered)]
  
  cat(sapply(regulator_genes, toString), 
      file = regulator_genes_file, sep = "\n")
} else {
  regulator_genes <- readLines(regulator_genes_file)
}

print(paste0("Number of regulatory genes included after filtering: ",
             length(regulator_genes)))

TPMs_z_scored_only_regulators <- TPMs_filtered_z_scored[
  row.names(TPMs_filtered_z_scored) %in% regulator_genes,
]

write.csv(TPMs_z_scored_only_regulators,
          file = "outputs/tpm_tables/TPM_z_scored_only_regs.csv",
          quote = FALSE)

# Load the transcription factor family data frame ==============

TF_family_preproc_df <- read.delim("data/Ath_TF_list.txt")
TF_family_preproc_df <- TF_family_preproc_df[, c("Gene_ID", "Family")]

TF_family_df <- unique(TF_family_preproc_df)

# for now remove duplicates (ie one gene associated to two families)
# but we may need to add them back in if needed
dup_TFs <- TF_family_df$Gene_ID[duplicated(TF_family_df$Gene_ID)]
TF_family_df <- TF_family_df[!(TF_family_df$Gene_ID %in% dup_TFs), ]
rownames(TF_family_df) <- TF_family_df$Gene_ID

rm(TF_family_preproc_df)

# Remove any leftover variables that will not be used again ============

rm(TPM_low_removed,
   conds,
   genes_to_keep,
   final_low_removed_genes)

gc()
