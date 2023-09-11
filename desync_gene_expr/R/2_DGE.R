# Input data ready for DE =====

# Note should remove outliers from previous section and replace P158 by P158_rerun 

txi_for_DE <- txi_data_genes

coerce_samples <- function(matrix) {
  # used specifically to delete P128 and singletons and replace P158 by P158_rerun
  keep_cols = which(!(colnames(matrix) %in% c("P128", "P169", "P196", "P098", "P079", "P158")))
  out_matrix <- matrix[, keep_cols]
  P158_rerun_loc <- which(colnames(out_matrix) == "P158_rerun")
  colnames(out_matrix)[P158_rerun_loc] <- "P158"
  return(out_matrix)
}

txi_for_DE$abundance <- coerce_samples(txi_for_DE$abundance)
txi_for_DE$counts <- coerce_samples(txi_for_DE$counts)
txi_for_DE$length <- coerce_samples(txi_for_DE$length)

# Perform Wilcoxon test on CPMs ===========

# separate out the bolting and non-bolting TPMs in order to run test

bolting_samples <- conds_filtered$id[
  which(conds_filtered$bolting == "Y")
]
non_bolting_samples <- conds_filtered$id[
  which(conds_filtered$bolting == "N")
]

# Replace TPM with CPM - as suggested in ref. (Li et al. Genome Biol. 2022)
CPM_DE <- txi_for_DE$counts
keep_genes_cpm <- which(
  row.names(CPM_DE) %in% row.names(TPM_filtered)
)
CPM_DE <- CPM_DE[keep_genes_cpm, ]
sample_sums <- colSums(CPM_DE)

for (col_i in 1:ncol(CPM_DE)) {
  CPM_DE[, col_i] <- (CPM_DE[, col_i] / sample_sums[col_i]) * 1e6
}

CPMs_bolting <- CPM_DE[, bolting_samples]
CPMs_non_bolting <- CPM_DE[, non_bolting_samples]

wilcox_output <- sapply(1:nrow(CPMs_bolting), function(gene_num) {
  wilcox.test(as.numeric(CPMs_bolting[gene_num, ]), 
              as.numeric(CPMs_non_bolting[gene_num, ]))
})

# apply correction
adjusted_p_vals <- p.adjust(wilcox_output["p.value", ], 
                            method = "BH")
# plot(wilcox_output["p.value", ], adjusted_p_vals)

# find the log fold changes
log_fold_changes <- log2(rowMeans(as.matrix(CPMs_bolting)) /
                           rowMeans(as.matrix(CPMs_non_bolting)))

# Make a custom volcano plot
volcano_data <- data.frame(id = row.names(TPM_filtered),
                           log2FoldChange = log_fold_changes,
                           p_adj = adjusted_p_vals)

volcano_data$diffexpressed <- "NO"
volcano_data$diffexpressed[volcano_data$log2FoldChange > 0.1 & 
                             volcano_data$p_adj < 0.05] <- "UP"
volcano_data$diffexpressed[volcano_data$log2FoldChange < -0.1 & 
                             volcano_data$p_adj < 0.05] <- "DOWN"

print(paste0("Number of up-regulated genes ",
             "(log2FoldChange > 0.1 and p_adj < 0.05): ",
             sum(volcano_data$diffexpressed == "UP")))

print(paste0("Number of down-regulated genes ",
             "(log2FoldChange < -0.1 and p_adj < 0.05): ",
             sum(volcano_data$diffexpressed == "DOWN")))

list_total_DEs <- volcano_data$id[
  volcano_data$diffexpressed == "UP" | volcano_data$diffexpressed == "DOWN"
]

# plot adding up all layers we have seen so far
volcano_plot <- ggplot(
    data = volcano_data, 
    aes(x = log2FoldChange,
        y = -log10(p_adj), 
        col = diffexpressed)
  ) +
  geom_point() + 
  theme_minimal(base_size = 14) +
  # geom_text_repel() +
  scale_color_manual(values = c("red", "black", "blue")) +
  geom_vline(xintercept = c(-0.1, 0.1), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  xlim(c(-3, 3))

ggsave("plots/clustering_and_DE/overall_volcano_plot.svg", volcano_plot,
       width = 10, height = 7)

volcano_data$gene_name <- id_to_name_df[volcano_data$id, "gene_name"]

write.xlsx(
  volcano_data[volcano_data$diffexpressed != "NO", ],
  "outputs/clustering_and_DE/DEGs.xlsx",
  row.names = FALSE,
  sheetName = "DEGs"
)

# Save results to upload to gProfiler for GO terms =========

# this is still very useful for uploading to website - a better interface
save_lists_gprofiler <- function(list_of_gene_lists, file_name) {

  list_to_write <- c()
  for (lvl in 1:length(list_of_gene_lists)) {
    list_to_write <- append(list_to_write,
                            paste0(">Cluster ", lvl))

    # select all genes in that cluster
    list_to_write <- append(list_to_write,
                            list_of_gene_lists[[lvl]])
  }
  fileConn <- file(file_name)
  writeLines(list_to_write, file_name)
  close(fileConn)
}

save_gene_lists_GOrilla <- function(
      list_of_vecs,
      folder_name  # should end in /
  ) {

    for (lvl in 1:length(list_of_vecs)) {
      # output to a separate file for each cluster

      file_out <- paste0(folder_name, "GOrilla_cluster_", lvl, ".txt")
      # select all genes in that cluster
      genes_in_lvl <- list_of_vecs[[lvl]]
      
      fileConn <- file(file_out)
      writeLines(genes_in_lvl, fileConn)
      close(fileConn)
    }

}

DE_UP <- as.vector(volcano_data[volcano_data$diffexpressed == "UP", ]$id)
DE_DOWN <- as.vector(volcano_data[volcano_data$diffexpressed == "DOWN", ]$id)

# save all DEs
save_lists_gprofiler(
  list(DE_UP, DE_DOWN),
  "outputs/clustering_and_DE/differentially_expressed.txt"
)

save_gene_lists_GOrilla(
  list(DE_UP, DE_DOWN),
  "outputs/clustering_and_DE/DE_for_GOrilla/"
)

# also save a background set of all genes in TAIR10 from TPM_master
fileConn <- file("outputs/clustering_and_DE/DE_for_GOrilla/background.txt")
writeLines(row.names(TPM_master), fileConn)
close(fileConn)

# Visualisation of TPMs with specific GO / KEGG terms ===========

# create TPMs where all values are clipped if z-score is > 3 or < -3

TPM_clipped <- apply(TPMs_filtered_z_scored, c(1,2), function(val) {
  if (val > 3) {
    val <- 3
  }
  if (val < -3) {
    val <- -3
  }
  return(val)
})

# Load gprofiler query if possible 

gprofiler_query_file <- "data/gProfiler_data/gprofiler_query.Rdata"

if (!file.exists(gprofiler_query_file)) {
  gprofiler_query_DEGs <- gost(
    query = list(DE_UP, DE_DOWN), 
    organism = "athaliana",
    multi_query = TRUE,
    user_threshold = 1e-4
  )
  
  save(gprofiler_query_DEGs, file = gprofiler_query_file)
} else {
  load(gprofiler_query_file)
}

# Save a version of the results for use in a supplementary fig.
# Note that many columns will have to be separated out into DOWN and UP

# This function splits a data frame with *vector values* into separate columns.
separate_columns <- function(
    df_combined,
    cols_to_split,  # give the NAMES of the columns to be split
    col_suffixes  # give the desired suffixes for each new column
){
  
  df_output <- df_combined[, !(colnames(df_combined) %in% cols_to_split)]
  
  for (col_name in cols_to_split) {
    for (col_suf_i in 1:length(col_suffixes)) {
      # get a vector of all 'col_suf_i'th elements from 'col_name' column
      new_col <- unlist(lapply(df_combined[[col_name]], function(vec_in){
        return(vec_in[col_suf_i])
      }))
      
      new_name <- paste0(col_name, col_suffixes[col_suf_i])
      # add this to the output data frame
      df_output[[new_name]] <- new_col
    }
  }
  
  return(df_output)
}

gprofiler_for_excel <- separate_columns(
  gprofiler_query_DEGs$result,
  c("p_values", "significant", "query_sizes", "intersection_sizes"),
  c("_up", "_down")
)

write.xlsx2(gprofiler_for_excel,
            "outputs/clustering_and_DE/gprofiler_DE_res.xlsx",
            row.names = FALSE)

# Make a heatmap based on interesting terms =======

get_genes_single_GO <- function(
  GO_term    
) {
  BM_single_GO <- BM_after_filtering[
    BM_after_filtering$go_id == GO_term,
  ]
  
  return(unique(BM_single_GO$tair_locus))
}

# retrieves a list of gene vectors, based on the GO terms provided.
# also filters these lists based on filter_list
get_genes_associated_GO_terms <- function(
  GO_term_list,
  filter_list
) {

  gene_lists_out <- vector("list", length(GO_term_list))

  for (term_i in 1:length(GO_term_list)) {
  
    unfiltered_remaining <- get_genes_single_GO(GO_term_list[term_i])
    
    filtered_remaining <- unfiltered_remaining[
      unfiltered_remaining %in% filter_list
    ]

    gene_lists_out[[term_i]] <- filtered_remaining

  }

  return(gene_lists_out)

}

make_heatmap_based_on_GO <- function(
  GO_term_list,
  TPM_table = TPM_clipped,
  new_group_labels = NULL,  # use this to provide custom gene group names
  rotate_labels = 0,
  row_dend_thickness = 10, # use to change the thickness of row dendogram in mm
  use_dendsort = FALSE,
  rotate_top_branches_point = NULL,
  reverse_dend = FALSE,
  use_single_clustering = FALSE,
  use_raster = FALSE, # set to true to get the heatmaps themselves to be images - useful for large svgs and removing grey borders
  raster_quality = 5
) {
  
  # get associated genes to terms
  possible_genes <- c(DE_DOWN, DE_UP)
  gene_lists <- get_genes_associated_GO_terms(GO_term_list, possible_genes)
  
  # find the dendogram overall for all genes
  specific_TPM <- TPM_table[unique(unlist(gene_lists)), ]
  
  if (use_single_clustering) {
    dend_out <- as.dendrogram(hclust(dist(t(specific_TPM)), method = "single"))
  } else {
    dend_out <- as.dendrogram(hclust(dist(t(specific_TPM))))
  }
  
  if (use_dendsort) {
    dend_out <- dendsort(dend_out)
  }
  
  if (!(is.null(rotate_top_branches_point))) {
    dend_out <- dendextend::rotate(
      dend_out, 
      c((rotate_top_branches_point + 1):ncol(specific_TPM),
        1:rotate_top_branches_point)
    )
  }
  
  if (reverse_dend) {
    dend_out <- rev(dend_out)
  }
  
  # create heatmaps
  heatmap_list = vector("list", length(GO_term_list))
  
  shared_col_fun <- colorRamp2(seq(from = -3, to = 3, length.out = 10),
                               inferno(10))
  
  shared_lgd_list <- list(
    col_fun = shared_col_fun,
    title = "Scaled TPM",
    at = c(-3, 0, 3),
    labels = c("-3", "0", "3")
  )
  
  # this annotation will only be used for the top heatmap
  bolting_status <- HeatmapAnnotation(
    bolting = conds_filtered$bolting,
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
        TPM_table[gene_lists[[term_i]], conds_filtered$id],
        name = GO_term_list[term_i],
        rect_gp = gpar(col = 0),
        col = shared_col_fun,
        row_dend_gp = gpar(lwd = row_dend_thickness),
        # rect_gp = gpar(col = "black", lwd = 0),
        row_title = row_label,
        row_title_rot = rotate_labels,
        cluster_columns = dend_out,
        top_annotation = bolting_status,
        heatmap_legend_param = shared_lgd_list,
        show_row_names = FALSE,
        use_raster = use_raster,
        raster_quality = raster_quality
      )
    } else {
      
      single_heatmap <- Heatmap(
        TPM_table[gene_lists[[term_i]], conds_filtered$id],
        name = GO_term_list[term_i],
        rect_gp = gpar(col = 0),
        col = shared_col_fun,
        # rect_gp = gpar(col = "black", lwd = 0),
        row_title = row_label,
        row_title_rot = rotate_labels,
        cluster_columns = dend_out,
        row_dend_gp = gpar(lwd = row_dend_thickness),
        # top_annotation = bolting_status,
        # heatmap_legend_param = shared_lgd_list,
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

# Plot the heatmap =============

GO_of_interest <- c("GO:0012501", "GO:0051726", "GO:0015979", "GO:0009753")

# spit out p-vals of these ...
sapply(GO_of_interest, function(go_term){
  select_row <- which(gprofiler_query_DEGs$result$term_id == go_term)
  p_vals <- unlist(gprofiler_query_DEGs$result[select_row, 2])
  cat(paste0(
    "GO term: ", go_term,
    "\nP-val of up-regulation: ", p_vals[1],
    "\nP-val of down-regulation: ", p_vals[2], "\n"
  ))
  return(p_vals)
})

short_names <- c("A", "B", "C", "D")

hlist <- make_heatmap_based_on_GO(
  GO_of_interest,
  new_group_labels = short_names,
  rotate_labels = 0,
  row_dend_thickness = 0.4,
  use_dendsort = TRUE,
  use_raster = TRUE,
  raster_quality = 10
)

total <- hlist[[1]] %v% hlist[[2]] %v% hlist[[3]] %v% hlist[[4]]

pdf(file = "plots/clustering_and_DE/overall_heatmap.pdf",
    width = 10, height = 7, pointsize = 7.5)
draw(total)
dev.off()

svg(filename = "plots/clustering_and_DE/overall_heatmap.svg",
    width = 10, height = 7, pointsize = 7.5)
draw(total)
dev.off()

tiff(filename = "plots/clustering_and_DE/overall_heatmap.tiff",
     width = 10, height = 7, units = "in", res = 600)
draw(total)
dev.off()

# Explore using box plots to see how genes depend on bolting =====

single_gene_bolting <- function(gene_id,
                                new_title = NULL,
                                p_val_to_display = NULL) {

  reduced_df <- data.frame(
    c(t(TPM_filtered[gene_id, ])),
    c(conds_filtered$bolting)
  )

  colnames(reduced_df) <- c("expr", "bolting")
  
  # print(test_output)
  
  if (is.null(new_title)) {
    new_title <- gene_id
  }
  
  if (is.null(p_val_to_display)) {
    test_output <- wilcox.test(expr ~ bolting, data = reduced_df)
    p_val_to_display <- test_output[["p.value"]]
  }
  
  output <- ggplot(reduced_df, aes(x = bolting, y = expr)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.1) +
    labs(title = new_title) +
    annotation_compass(label = paste0("Adjusted p value: ", 
                          formatC(signif(p_val_to_display, digits = 3), 
                                  digits = 3,
                                  format = "fg", 
                                  flag = "#")),
                       position = "SE")
  
  return(output)
}
