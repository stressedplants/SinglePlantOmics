# Data prep ===========

TPM_clipped_and_01 <- as.data.frame(t(apply(TPM_clipped, 1, function(i){
  temp <- as.numeric(unlist(i))
  return((temp - min(temp)) / max((temp - min(temp))))
})))
colnames(TPM_clipped_and_01) <- colnames(TPM_clipped)

UP_list <- intersect(DE_UP, row.names(TPM_clipped_and_01))
DOWN_list <- intersect(DE_DOWN, row.names(TPM_clipped_and_01))

TPM_DE_filtered <- TPM_clipped_and_01[UP_list, ]

# if genes are DOWN regulated, then reverse the min and max
TPM_DE_filtered <- rbind(TPM_DE_filtered, 1 - TPM_clipped_and_01[DOWN_list, ])

# Order plants with resampling ===========

return_order_from_genes <- function(gene_ids,
                                    TPM_table = TPM_DE_filtered) {
  TPM_temp <- TPM_table[intersect(gene_ids, row.names(TPM_table)), ]
  
  ordered_samples <- names(sort(colSums(TPM_temp)))
  return(ordered_samples)
}

# first predict lots of orderings based on a partition of the genes

partition_crossvalidation <- function(num_bins,
                                TPM_table = TPM_DE_filtered,
                                seed = 123) {
  
  set.seed(seed)  # for reproducibility
  
  num_samples <- ncol(TPM_table)
  
  chunk <- function(x,n) split(x, cut(sample(seq_along(x)), n, labels = FALSE))
  gene_lists <- chunk(row.names(TPM_table), num_bins)
  
  ptimes <- lapply(gene_lists, function(i){
    return_order_from_genes(i, TPM_table = TPM_table)
  })
  
  ptimes_df <- data.frame(matrix(NA, 
                                 nrow = num_bins, 
                                 ncol = num_samples))
  colnames(ptimes_df) <- colnames(TPM_table)
  
  for (vec_i in 1:length(ptimes)) {
    temp_vec <- 1:num_samples
    names(temp_vec) <- ptimes[[vec_i]]
    ptimes_df[vec_i, names(temp_vec)] <- temp_vec
  }
  
  set.seed(Sys.time())  # unset
  
  return(ptimes_df)
}

cv_output <- partition_crossvalidation(100)

# order by the mode than add a heatmap
mode_val <- function(x) {
  which.max(tabulate(x))
}

set.seed(123)  # ensure that decisions about pseudotime as reproducible

mode_values <- apply(cv_output, 2, mode_val)
cv_output <- cv_output[, order(mode_values)]
rownames(cv_output) <- as.character(rownames(cv_output))

pseudotime_output <- 1:ncol(cv_output)
names(pseudotime_output) <- colnames(cv_output)

reordered_bolting <- anno_col[colnames(cv_output), ]
reordered_anno_col <- data.frame("bolting" = reordered_bolting)
rownames(reordered_anno_col) <- colnames(cv_output)

renamed_anno_colours <- anno_colours
names(renamed_anno_colours) <- c("bolting")

svg(filename = "plots/pseudotime/cv_heatmap.svg",
    width = 10, height = 11.5, pointsize = 7.5)
pheatmap(cv_output,
         color = inferno(10),
         cluster_cols = FALSE,
         show_rownames = TRUE,
         annotation_col = reordered_anno_col,
         annotation_colors = renamed_anno_colours)
dev.off()


# Make a metastable state heatmap =========

# metastable states are clusters of 'uncertainty' along the diagonal

metastable_mat <- data.frame(matrix(
  NA, nrow = ncol(cv_output), ncol = ncol(cv_output)
))

rownames(metastable_mat) <- colnames(cv_output)
colnames(metastable_mat) <- colnames(cv_output)

for (sample_i in colnames(cv_output)) {
  for (sample_j in colnames(cv_output)) {
    simple_df <- cv_output[, c(sample_i, sample_j)]
    
    # sample probability of sample i < sample j in pseudotime
    sample_proportion <- (sum(simple_df[,1] < simple_df[,2]) / nrow(simple_df))
    metastable_mat[sample_i, sample_j] <- sample_proportion
  }
}

svg(filename = "plots/pseudotime/metastable_heatmap.svg",
    width = 10, height = 10, pointsize = 7.5)
pheatmap(metastable_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE, 
         color = plasma(10))
dev.off()

# Create final pseudotime TPMs ======

TPM_pseudotime <- TPM_clipped_and_01[, colnames(cv_output)]

write.csv(TPM_pseudotime, "outputs/tpm_tables/TPM_pseudotime.csv")

# Plot pseudotime along PCA of plants =====

pca_on_plants <- prcomp(t(TPM_pseudotime))
plant_pca_summary <- summary(pca_on_plants)$importance

first_pc <- 1
second_pc <- 2

# correct order for conditions
conds_pseudotime <- conds_filtered[colnames(cv_output), ]  

pca_plant_for_plot <- data.frame(pca_on_plants$x[, c(first_pc, second_pc)])
pca_plant_for_plot$bolting <- conds_pseudotime$bolting
pca_plant_for_plot$pseudotime <- pseudotime_output

pca_bolting <- ggplot(pca_plant_for_plot, aes(x = PC1, y = PC2, 
                                              color = bolting)) +
  geom_point() +
  # geom_text() +
  theme_classic(base_size = 12) +
  xlab(paste0("PC", first_pc, " (", 
              sprintf(plant_pca_summary[2, first_pc]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC", second_pc," (",
              sprintf(plant_pca_summary[2, second_pc]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  scale_color_manual(values = c("N" = no_col,
                                "Y" = yes_col))

pca_bolting <- pca_bolting + labs(color = "Bolting") + theme(
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12),
)

ggsave("plots/pseudotime/pca_bolting_scatter.svg",
       plot = pca_bolting,
       width = 7, height = 5)

# Compare pseudotime to dimensionality reduction ======

pca_pseudotime <- ggplot(pca_plant_for_plot, aes(x = PC1, y = PC2, 
                                                 color = pseudotime)) +
  geom_point() +
  scale_color_gradient(low = "#FBB03B", high = "#D4145A") +
  # geom_text() +
  theme_classic(base_size = 12) +
  xlab(paste0("PC", first_pc, " (", 
              sprintf(plant_pca_summary[2, first_pc]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC", second_pc," (",
              sprintf(plant_pca_summary[2, second_pc]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
  )

ggsave("plots/pseudotime/pca_pseudotime_scatter.svg",
       plot = pca_pseudotime,
       width = 7, height = 5)

# Replot physiology data but with pseudotime ========

phenos_with_pseudotime <- phenos_to_predict
phenos_with_pseudotime$pseudotime <- pseudotime_output[
  row.names(phenos_with_pseudotime)
]

pheno_pseudotime <- ggplot(phenos_with_pseudotime, aes(x = leaf_avg, 
                                                       y = biomass, 
                                                       color = pseudotime)) +
  geom_point() +
  scale_color_gradient(low = "#FBB03B", high = "#D4145A") +
  # geom_text() +
  geom_text_repel(label = rownames(phenos_with_pseudotime),
                  max.overlaps = 20,
                  max.time = 5,
                  max.iter = 1e7,
                  seed = 1234) +
  theme_classic(base_size = 12) +
  xlab("Average leaf size (mm^2)") +
  ylab("Biomass (mg)") +
  labs(color = "Pseudotime") + 
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
  )

ggsave("plots/physiology_analysis/bolting_leaf_pseudotime_scatter.svg",
       plot = pheno_pseudotime,
       width = 10, height = 7)

# Visualise whether key genes look sensible over ptime ====

# a random sample of 400 genes
# 
# pheatmap(TPM_pseudotime[sample(union(DE_DOWN, DE_UP), size = 400), ],
#          cluster_cols = FALSE,
#          annotation_col = reordered_anno_col,
#          annotation_colors = renamed_anno_colours)

GO_heatmap_over_pseudotime <- function(
    GO_term_list,
    filter_list = union(DE_UP, DE_DOWN),
    TPM_table = TPM_pseudotime,
    new_group_labels = NULL,  # use this to provide custom gene group names
    rotate_labels = 0,
    row_dend_thickness = 10, # use to change the thickness of row dendogram in mm
    use_raster = FALSE, # set to true to get the heatmaps themselves to be images - useful for large svgs and removing grey borders
    raster_quality = 5
) {
  
  gene_lists <- get_genes_associated_GO_terms(GO_term_list, filter_list)
  specific_TPM <- TPM_table[unique(unlist(gene_lists)), ]
  
  # create heatmaps
  heatmap_list = vector("list", length(GO_term_list))
  
  shared_col_fun <- colorRamp2(seq(from = 0, to = 1, length.out = 10),
                               inferno(10))
  
  shared_lgd_list <- list(
    col_fun = shared_col_fun,
    title = "Scaled TPM",
    at = c(0, 0.5, 1),
    labels = c("0", "0.5", "1")
  )
  
  # this annotation will only be used for the top heatmap
  bolting_status <- HeatmapAnnotation(
    bolting = conds_pseudotime$bolting,
    annotation_label = "Bolted",
    annotation_name_side = "left",
    col = list(bolting = c("Y" = yes_col, "N" = no_col))
  )
  
  for (term_i in 1:length(GO_term_list)) {
    
    row_label <- GO_term_list[term_i]
    if (!is.null(new_group_labels)) {
      row_label <- new_group_labels[term_i]
    }
    
    if (term_i == 1) {
      single_heatmap <- Heatmap(
        TPM_table[gene_lists[[term_i]], ],
        name = GO_term_list[term_i],
        rect_gp = gpar(col = 0),
        col = shared_col_fun,
        row_dend_gp = gpar(lwd = row_dend_thickness),
        row_title = row_label,
        row_title_rot = rotate_labels,
        cluster_columns = FALSE,  # don't cluster the columns because pseudotime
        top_annotation = bolting_status,
        heatmap_legend_param = shared_lgd_list,
        show_row_names = FALSE,
        use_raster = use_raster,
        raster_quality = raster_quality
      )
    } else {
      
      single_heatmap <- Heatmap(
        TPM_table[gene_lists[[term_i]], ],
        name = GO_term_list[term_i],
        rect_gp = gpar(col = 0),
        col = shared_col_fun,
        row_title = row_label,
        row_title_rot = rotate_labels,
        cluster_columns = FALSE,
        row_dend_gp = gpar(lwd = row_dend_thickness),
        show_heatmap_legend = FALSE,
        show_row_names = FALSE,
        use_raster = use_raster,
        raster_quality = raster_quality
      )
    }
    
    heatmap_list[[term_i]] <- single_heatmap
  }
  
  return(heatmap_list)
}

GO_of_interest <- c("GO:0012501", "GO:0051726", "GO:0015979", "GO:0009753")
short_names <- c("A", "B", "C", "D")

hlist_ptime <- GO_heatmap_over_pseudotime(
  GO_of_interest,
  new_group_labels = short_names,
  rotate_labels = 0,
  row_dend_thickness = 0.4,
  use_raster = TRUE,
  raster_quality = 10
)

total_ptime <- hlist_ptime[[1]] %v% hlist_ptime[[2]] %v% hlist_ptime[[3]] %v% hlist_ptime[[4]]

svg(filename = "plots/pseudotime/pseudotime_heatmap.svg",
    width = 10.5, height = 7, pointsize = 7.5)
draw(total_ptime)
dev.off()

tiff(filename = "plots/pseudotime/pseudotime_heatmap.tiff",
     width = 10.5, height = 7, units = "in", res = 600)
draw(total_ptime)
dev.off()
