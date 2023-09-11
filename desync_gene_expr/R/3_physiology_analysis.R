# Scatter plot of leaf size vs bolting =======

pheno_scatter <- ggplot(phenos_to_predict, aes(x = leaf_avg, y = biomass, 
                                               color = bolting)) +
                    geom_point() +
                    scale_color_manual(values = c("N" = no_col,
                                                  "Y" = yes_col)) +
                    # geom_text() +
                    geom_text_repel(label = rownames(phenos_to_predict),
                                    max.overlaps = 20,
                                    max.time = 5,
                                    max.iter = 1e7,
                                    seed = 1234) +
                    theme_classic(base_size = 12) +
                    xlab("Average leaf size (mm^2)") +
                    ylab("Biomass (mg)")
                    
pheno_scatter <- pheno_scatter + labs(color = "Bolting") + theme(
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12),
)

ggsave("plots/physiology_analysis/bolting_leaf_scatter.svg",
       pheno_scatter,
       width = 10, height = 7)

# Make a copy of the phenotype scatter plot with bigger text / symbols ======

pheno_scatter_big <- ggplot(phenos_to_predict, aes(x = leaf_avg, y = biomass)) +
  geom_point(aes(color = bolting), size = 3.5) +
  scale_color_manual(values = c("N" = no_col,
                                "Y" = yes_col)) +
  theme_classic(base_size = 18) +
  xlab("Average leaf size (mm^2)") +
  ylab("Biomass (mg)") + 
  labs(color = "Bolting status") + 
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.position = c(0.9, 0.2)
  ) +
  # can add a smooth line if this helps the explanation of points which are 
  # unexpectedly bigger biomass
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2, col = "black")

phys_model <- lm(biomass ~ leaf_avg, data = phenos_to_predict)
summary(phys_model)

ggsave("plots/physiology_analysis/bolting_leaf_scatter_big.svg",
       pheno_scatter_big,
       width = 7*1.2, height = 5*1.3)

# Comparison of DGEs with physiology prediction ==========

biomass_regs_enet_params <- as.data.frame(read_excel(
  path = "data/physiology_pred/biomass_regs_coefs.xlsx"
))
row.names(biomass_regs_enet_params) <- biomass_regs_enet_params[, 1]
biomass_regs_enet_params <- biomass_regs_enet_params[, -1]

non_zero_biomass_genes <- colnames(biomass_regs_enet_params)[
  which(abs(biomass_regs_enet_params["median", ]) > 0)
]

leaf_regs_enet_params <- as.data.frame(read_excel(
  path = "data/physiology_pred/leaf_regs_coefs.xlsx"
))
row.names(leaf_regs_enet_params) <- leaf_regs_enet_params[, 1]
leaf_regs_enet_params <- leaf_regs_enet_params[, -1]

non_zero_leaf_genes <- colnames(leaf_regs_enet_params)[
  which(abs(leaf_regs_enet_params["median", ]) > 0)
]

# How many of these are DEGs?

all_DGEs <- row.names(volcano_data)[volcano_data$diffexpressed != "NO"]

print(paste0(sum(non_zero_biomass_genes %in% all_DGEs), " out of ",
             length(non_zero_biomass_genes), " possible non-zero median biomass ",
             "predictors are differentially expressed"))

print(paste0(sum(non_zero_leaf_genes %in% all_DGEs), " out of ",
             length(non_zero_leaf_genes), " possible non-zero median leaf ",
             "predictors are differentially expressed"))
