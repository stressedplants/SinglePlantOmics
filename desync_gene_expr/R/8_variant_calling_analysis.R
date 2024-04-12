library(vcfR)

# Load data with the gene ID info included ========

# add the column names from the vcf only file
temp_vcf_file <- read.vcfR("data/variant_calling/final_filtered.recode.vcf")

# Basic PCA of variants to see whether bolting can be determined from variants =========

only_geno <- temp_vcf_file@gt[, 2:ncol(temp_vcf_file@gt)]

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

# save the variants as numbers
write.csv(only_geno_numeric, 
          file = "outputs/variant_calling/removed_all_heterozygous.csv")

# save as a vcf too - useful for further analysis
vcf_filtered <- temp_vcf_file
vcf_filtered@fix <- vcf_filtered@fix[!all_1_rows, ]
vcf_filtered@gt <- vcf_filtered@gt[!all_1_rows, ]
write.vcf(vcf_filtered, 
          file = "outputs/variant_calling/removed_all_heterozygous.vcf")

# are there any rows with a 0, 1, and 2 in?
any_012_rows <- apply(only_geno_numeric, 1, function(row_in){
  (0 %in% row_in) & (1 %in% row_in) & (2 %in% row_in)
})
print(paste0(sum(any_012_rows), " variants include all possible genotypes in population"))

pca_out <- prcomp(t(only_geno_numeric))
pca_summary <- summary(pca_out)
# sum_of_importances <- sum(pca_summary$importance[2, ])

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

# Make a plot of correlations by chromosome and position ======

library(qqman)

chromosome_info <- vcf_filtered@fix[, "CHROM"]
position_info <- as.numeric(vcf_filtered@fix[, "POS"])

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
