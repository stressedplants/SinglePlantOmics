# Load data with the gene ID info included ========

gene_id_vcf <- read.delim("data/variant_calling/vcf_with_gene_names.txt", header = FALSE)

# add the column names from the vcf only file
temp_vcf_file <- read.delim("data/variant_calling/final_filtered.recode.vcf", skip = 8)
colnames(gene_id_vcf)[1:ncol(temp_vcf_file)] <- colnames(temp_vcf_file)

# add the standard gtf column names
colnames(gene_id_vcf)[75:83] <- c("seqname", "source", "feature", "start", 
                                  "end", "score", "strand", "frame",
                                  "attribute")

# filter for only protein coding and non-duplicated values
gene_id_vcf_protein_only <- gene_id_vcf[
  grep("protein_coding", gene_id_vcf[,"attribute"]), 
]

combined_chrom_position_protein <- sapply(1:nrow(gene_id_vcf_protein_only), function(i){
  paste0(gene_id_vcf[i,1], "_", gene_id_vcf[i,2])
})

gene_id_vcf_protein_only <- gene_id_vcf_protein_only[
  !duplicated(combined_chrom_position_protein), 
]

# Basic PCA of variants to see whether bolting can be determined from variants =========

only_geno <- gene_id_vcf_protein_only[, 10:74]

only_geno_numeric <- apply(only_geno, c(1,2), function(geno_in){
  if (geno_in == "0|0") return(0)
  if (geno_in == "0|1" | geno_in == "1|0") return(1)
  if (geno_in == "1|1") return(2)
})

# remove any where all rows are 1 - i.e. all samples heterogeneous
all_1_rows <- apply(only_geno_numeric, 1, function(row_in) {
  all(row_in == 1)
})
only_geno_numeric <- only_geno_numeric[!all_1_rows, ]

# are there any rows with a 0, 1, and 2 in?
any_012_rows <- apply(only_geno_numeric, 1, function(row_in){
  (0 %in% row_in) & (1 %in% row_in) & (2 %in% row_in)
})
print(paste0(sum(any_012_rows), " variants include all possible genotypes in population"))

pca_out <- prcomp(t(only_geno_numeric))
pca_summary <- summary(pca_out)
sum_of_importances <- sum(pca_summary$importance[2, ])

# plot with pseudotime
pca_snps_for_plot <- data.frame(pca_out$x[, c(1, 2)])
pca_snps_for_plot$bolting <- conds_filtered$bolting
pca_snps_for_plot$pseudotime <- pseudotime_output[row.names(pca_snps_for_plot)]

pca_snps_pseudotime <- ggplot(pca_snps_for_plot, aes(x = PC1, y = PC2, 
                                              color = pseudotime)) +
  geom_point() +
  # geom_text() +
  theme_classic(base_size = 12) +
  ggtitle("Pseudotime compared to PCA of variants") +
  xlab(paste0("PC", first_pc, " (", 
              sprintf(pca_summary$importance[2, 1]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC", second_pc," (",
              sprintf(pca_summary$importance[2, 2]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  scale_color_gradient(low = "#FBB03B", high = "#D4145A")

ggsave("plots/variant_calling/pca_snps_pseudotime.svg", pca_snps_pseudotime)

# bolting ggplot

pca_snps_bolting <- ggplot(pca_snps_for_plot, aes(x = PC1, y = PC2, 
                                                  color = bolting)) +
  geom_point() +
  theme_classic(base_size = 12) +
  ggtitle("Bolting status compared to PCA of variants") +
  xlab(paste0("PC", first_pc, " (", 
              sprintf(pca_summary$importance[2, 1]*100, fmt = "%#.1f"),
              "% of variance)")) +
  ylab(paste0("PC", second_pc," (",
              sprintf(pca_summary$importance[2, 2]*100, fmt = "%#.1f"),
              "% of variance)")) + 
  scale_color_manual(values = c("N" = no_col,
                                "Y" = yes_col))

ggsave("plots/variant_calling/pca_snps_bolting.svg", pca_snps_bolting)

# Make a simple pseudotime plot visualising the variants =========

# first need to filter for the most important variants according to PC1 (quick and easy)
pc1_comps <- abs(pca_out$rotation[,1])
which_snps <- order(pc1_comps, decreasing = TRUE)
bolting_important_variants <- pc1_comps[pc1_comps > 0.06]
variant_ready_for_plot <- only_geno_numeric[names(bolting_important_variants), ]

# get gene names too
bolting_variant_genes <- gene_id_vcf_protein_only[names(bolting_important_variants), ]

gene_names_only <- unname(sapply(bolting_variant_genes[, "attribute"], 
  function(str_i) {
    first_match <- str_match(str_i, "ID=gene:[^;]*;")
    second_match <- str_match(first_match, "AT[0-9]G[0-9]*")
    return(second_match)
  }
))

bolting_variant_genes$gene <- gene_names_only

# # reorder by pseudotime
# variant_ready_for_plot <- variant_ready_for_plot[, names(pseudotime_output)]

# # TODO: replace with ComplexHeatmap and add annotations using rowAnnotation and anno_text
# pheatmap(variant_ready_for_plot, cluster_cols = FALSE, 
#          annotation_row = data.frame(gene = gene_names_only),
#          annotation_col = conds_filtered)

# Find correlations between pseudotime / physiology and variants =========

ptime_variant_corrs <- apply(only_geno_numeric, 1, function(row_in) {
  # ensure that pseudotime is reordered to consider the correlation
  cor(row_in, pseudotime_output[names(row_in)])
})

# repeat for biomass
biomass_variant_corrs <- apply(only_geno_numeric, 1, function(row_in) {
  # ensure that pseudotime is reordered to consider the correlation
  cor(row_in, phenos_to_predict[names(row_in), "biomass"])
})


# repeat for leaf size
leaf_variant_corrs <- apply(only_geno_numeric, 1, function(row_in) {
  # ensure that pseudotime is reordered to consider the correlation
  cor(row_in, phenos_to_predict[names(row_in), "leaf_avg"])
})

# create pseudotime SNPs heatmap ==========

anno_col_ptime_snps <- data.frame(bolting = conds_pseudotime$bolting)
row.names(anno_col_ptime_snps) <- row.names(conds_pseudotime)

pseudotime_snps <- pheatmap(
  only_geno_numeric[which(abs(ptime_variant_corrs) >= 0.5),
                    names(pseudotime_output)],
  color = c("#45BBB2", "#F0F0F0", "#EC8254"),
  cluster_cols = FALSE,
  cutree_rows = 6,
  show_rownames = FALSE,
  annotation_legend = FALSE,
  legend_breaks = c(0, 1, 2), 
  legend_labels = c("0|0", "0|1", "1|1"),
  annotation_col = anno_col_ptime_snps,
  annotation_colors = list(bolting = c("N" = no_col, "Y" = yes_col))
)

svg(filename = "plots/variant_calling/pseudotime_snps_heatmap.svg",
    width = 10, height = 7, pointsize = 7.5)
draw(pseudotime_snps)
dev.off()

# bolting_status_anno <- HeatmapAnnotation(
#   bolting = conds_pseudotime$bolting,
#   annotation_label = "Bolted",
#   annotation_name_side = "left",
#   col = list(bolting = c("Y" = yes_col, "N" = no_col))
# )

# pseudotime_snps <- Heatmap(
#   only_geno_numeric[which(abs(ptime_variant_corrs) >= 0.5),
#                     names(pseudotime_output)],
#   # cluster_cols = FALSE,
#   # show_rownames = FALSE,
#   top_annotation = bolting_status_anno, 
#   # annotation_legend = FALSE,
#   heatmap_legend_param = list(
#     title = "variant",
#     at = c(0,1,2),
#     labels = c("0|0", "0|1", "1|1")
#   )
# )


# Make a plot of correlations by chromosome and position ======

library(qqman)

# # this is an example data frame - shows the required format
# View(gwasResults)
# manhattan(gwasResults, logp = FALSE)  # use logp = FALSE for other scores

chromosome_info <- gene_id_vcf_protein_only[!all_1_rows, "X.CHROM"]
position_info <- gene_id_vcf_protein_only[!all_1_rows, "POS"]

plot_manhattan_style_correlation <- function(correlation_vec,
                                             suggestiveline = 0.5,
                                             ylab = "Correlation",
                                             ...) {
  manhattan_plot_df <- data.frame(SNP = rep("SNP", length(chromosome_info)),
                                   CHR = chromosome_info,
                                   BP = position_info,
                                   P = correlation_vec)
  
  # remove rows which are not on chromosomes 1 through 5
  manhattan_plot_df <- manhattan_plot_df[
    manhattan_plot_df$CHR %in% c("1", "2", "3", "4", "5"),
  ]
  manhattan_plot_df$CHR <- as.numeric(manhattan_plot_df$CHR)
  
  manhattan(manhattan_plot_df, 
            logp = FALSE, 
            suggestiveline = suggestiveline, 
            genomewideline = FALSE,
            ylab = ylab,
            ...)
}

svg(filename = "plots/variant_calling/stacked_manhattan_plot.svg", width = 5,
    height = 7)
old_pars <- par(mfrow = c(3,1), mar = c(4, 4, 1.5, 1.5))
plot_manhattan_style_correlation(ptime_variant_corrs, main = "Pseudotime", 
                                 xlab = "",
                                 cex.lab = 1.5,
                                 cex.axis = 1.5,
                                 cex.main = 2,
                                 cex.sub = 1.5)
plot_manhattan_style_correlation(biomass_variant_corrs, main = "Biomass", 
                                 xlab = "",
                                 cex.lab = 1.5,
                                 cex.axis = 1.5,
                                 cex.main = 2,
                                 cex.sub = 1.5)
plot_manhattan_style_correlation(leaf_variant_corrs, main = "Leaf size",
                                 cex.lab = 1.5,
                                 cex.axis = 1.5,
                                 cex.main = 2,
                                 cex.sub = 1.5)
par(old_pars)
dev.off()

# Supplemental: correlation of variants to one another ==========

variant_corr_mat <- cor1(t(
  only_geno_numeric[which(abs(ptime_variant_corrs) >= 0.5),]
))
pheatmap(abs(variant_corr_mat))
