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
                                    max.iter = 1e3,
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

# Test to compare leaf size and biomass based on bolting status ====

wilcox.test(biomass ~ bolting, data = phenos_to_predict)
wilcox.test(leaf_avg ~ bolting, data = phenos_to_predict)
