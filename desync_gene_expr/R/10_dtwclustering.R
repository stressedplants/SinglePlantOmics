library(dtwclust)
library(cluster)
library(factoextra)

# providing the values of k simulataneously allows tsclust to precompute 
# the distance matrix. In testing, this is most time consuming part.
multiple_k <- 2:50
seed_id_vec <- 1:5

evaluated_matrix_for_clustering <- t(sapply(monotonic_genes, function(i){
  eval.fd(seq(from = range_of_timeseries[1], 
              to = range_of_timeseries[2],
              by = 1), rescaled_fd[i])
}))

clustering_output_file <- "outputs/pseudotime_analysis/multiple_dtw_clustering.RData"

# avoid running if the file already exists
# save the distance matrix and clustering outputs for each run - NOT everything!
if (!file.exists(clustering_output_file)) {
  # run this for 5 different seeds. this will have to recompute the dist matrix
  # but the argument nrep didn't seem to work with multiple k values.
  
  cluster_df <- data.frame(matrix(nrow = length(multiple_k),
                                  ncol = length(seed_id_vec)))
  colnames(cluster_df) <- seed_id_vec
  rownames(cluster_df) <- multiple_k
  
  clusinfo_df <- data.frame(matrix(nrow = length(multiple_k),
                                   ncol = length(seed_id_vec)))
  colnames(clusinfo_df) <- seed_id_vec
  rownames(clusinfo_df) <- multiple_k
  
  for (seed_id in seed_id_vec) {
    
    # only required to precompute the distance matrix once
    if (seed_id == 1) {
      single_output <- tsclust(
        evaluated_matrix_for_clustering,
        distance = "sbd",
        k = multiple_k,
        seed = seed_id,
        preproc = NULL,
        trace = TRUE
      )
      
      sbd_dist <- single_output[[1]]@distmat
    } else {
      predefined_distmat <- partitional_control(
        distmat = sbd_dist
      )
      
      single_output <- tsclust(
        evaluated_matrix_for_clustering,
        distance = "sbd",
        k = multiple_k,
        seed = seed_id,
        preproc = NULL,
        trace = TRUE,
        control = predefined_distmat
      )
    }
    
    # set all the clusters for one seed
    cluster_temp_list <- vector("list", length(multiple_k))
    for (i in 1:length(multiple_k)) {
      cluster_temp_list[[i]] <- single_output[[i]]@cluster
    }
    # USE CHARS TO INDEX THESE DATA FRAMES
    cluster_df[[as.character(seed_id)]] <- cluster_temp_list
    
    clusinfo_temp_list <- vector("list", length(multiple_k))
    for (i in 1:length(multiple_k)) {
      clusinfo_temp_list[[i]] <- single_output[[i]]@cldist
    }
    clusinfo_df[[as.character(seed_id)]] <- clusinfo_temp_list
    
    rm(single_output, cluster_temp_list, clusinfo_temp_list)
    
  }
  save(sbd_dist, cluster_df, clusinfo_df, file = clustering_output_file)
} else {
  cat("Loading clustering outputs\n")
  load(clustering_output_file)
}

# PLAN
# 1) choose an optimal k via the elbow method
# 2) make a plot of what the clusters look like

# 1) choosing k ===========

calculate_WCSS <- function(cluster_membership_vec,
                           cluster_distance) {
  clusters <- unique(cluster_membership_vec)
  total_sum <- 0
  for (id in clusters) {
    # select scores only for this cluster
    cluster_distance_specific <- cluster_distance[
      which(cluster_membership_vec == id)
    ]
    total_sum <- total_sum + sum(cluster_distance_specific**2)
  }
  return(total_sum)
}

# make a summary data frame
WCSS_df <- data.frame("k" = multiple_k,
                      "mean" = multiple_k,
                      "min" = multiple_k,
                      "max" = multiple_k)
row.names(WCSS_df) <- as.character(multiple_k)

for (k_temp in as.character(multiple_k)) {
  temp_outputs_per_k <- rep(0, length(seed_id_vec))
  names(temp_outputs_per_k) <- as.character(seed_id_vec)
  for (seed_id_temp in as.character(seed_id_vec)) {
    WCSS_single <- calculate_WCSS(
      cluster_df[k_temp, seed_id_temp][[1]], 
      clusinfo_df[k_temp, seed_id_temp][[1]]
    )
    temp_outputs_per_k[seed_id_temp] <- WCSS_single
  }
  WCSS_df[k_temp, "mean"] <- mean(temp_outputs_per_k)
  WCSS_df[k_temp, "min"] <- min(temp_outputs_per_k)
  WCSS_df[k_temp, "max"] <- max(temp_outputs_per_k)
}

opt_num_clusters <- 25

# Overall plot of all clusters
p_wcss <- ggplot(WCSS_df, aes(x = k, y = mean)) + 
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = min, ymax = max), width = .2) +
  scale_y_log10() +
  geom_vline(xintercept = opt_num_clusters, linetype = "dotted", 
             color = "blue", linewidth = 1.5)

ggsave("plots/pseudotime_analysis/wcss_on_all.pdf", plot = p_wcss, device = "pdf")

# 2) visualise the optimal number of clusters ======

predefined_distmat <- partitional_control(
  distmat = sbd_dist
)

opt_clustering <- tsclust(
  evaluated_matrix_for_clustering,
  distance = "sbd",
  k = opt_num_clusters,
  seed = 1,
  preproc = NULL,
  trace = TRUE,
  control = predefined_distmat
)@cluster

# prepare a list of the genes in each cluster
genes_per_cluster <- vector(mode = 'list', length = opt_num_clusters)

for (clust_i in 1:opt_num_clusters) {
  genes_per_cluster[[clust_i]] <- row.names(evaluated_matrix_for_clustering)[
    opt_clustering == clust_i
  ]
}

# prepare a list of mean AUC values per cluster
AUC_per_cluster <- sapply(genes_per_cluster, function(gene_vec) {
  mean(AUC_total_df[gene_vec, "AUC_value"])
})

# manually reordered to put decreasing clusters before increasing clusters
decr_clusters <- c(
  1, 2, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 21, 23, 24, 25
)
incr_clusters <- c(
  3, 6, 12, 18, 19, 20, 22
)
order_for_opt_clusters <- c(
  decr_clusters[order(AUC_per_cluster[decr_clusters])],
  incr_clusters[order(AUC_per_cluster[incr_clusters])]
)

plot_time_series_clusters <- function(
    time_series_df,
    cluster_assignments,
    AUC_values,
    order_of_plots,  # should be a reordering of 1:length(unique(cluster_assignments))
    total_ncol = 3,
    alpha_level = 0.2
) {
  
  # Number of clusters
  num_clusters <- length(unique(cluster_assignments))
  
  # Create a list to store plots
  plots <- list()
  
  # Iterate over each cluster
  for (clust_selector in 1:num_clusters) {
    
    clust_i <- order_of_plots[clust_selector]
    
    cluster_data <- as.data.frame(
      time_series_df[, cluster_assignments == clust_i]
    )
    
    plotting_df <- data.frame(
      matrix(NA, nrow = (nrow(cluster_data) * ncol(cluster_data)), ncol = 3)
    )
    names(plotting_df) <- c("time", "value", "id")
    
    plotting_df$time <- rep(1:nrow(cluster_data), times = ncol(cluster_data))
    plotting_df$id <- rep(1:ncol(cluster_data), each = nrow(cluster_data))
    
    for (sample_i in 1:ncol(cluster_data)) {
      plotting_df$value[plotting_df$id == sample_i] <- cluster_data[, sample_i]
    }
    
    # make title including AUC value
    formatted_AUC_value <- formatC(AUC_values[clust_i], digits = 3)
    
    plot_title <- paste0("cluster ", clust_i, 
                         " | mean AUC : ", formatted_AUC_value)
    
    # plot time series for the current cluster
    plot <- ggplot(plotting_df, aes(x = time, y = value, group = id)) +
      geom_line(alpha = alpha_level) +
      labs(title = plot_title) +
      xlab("Pseudotime") +
      ylab("Smoothed gene expression (scaled)") +
      theme_minimal()
    
    # Add plot to the list
    plots[[clust_selector]] <- plot  
    
  }
  
  # Combine plots into a grid
  return(grid.arrange(grobs = plots, ncol = total_ncol))
}

cluster_viz <- plot_time_series_clusters(
  t(evaluated_matrix_for_clustering),
  opt_clustering,
  AUC_per_cluster,
  order_for_opt_clusters,
  total_ncol = 5
)

ggsave(
  "plots/pseudotime_analysis/clustering_output.pdf",
  cluster_viz,
  width = 18.75,
  height = 15,
  device = 'pdf'
)

set.seed(123) 

# for some extra context, plot a random gene from each of these clusters
pdf("plots/pseudotime_analysis/cluster_examples.pdf",
    width = 18.75,
    height = 18.75)
old_pars <- par(mfrow = c(5,5))
for (clust_i in order_for_opt_clusters) {
  sample_gene <- sample(genes_per_cluster[[clust_i]], 1)
  plot_title <- paste0(
    "Cluster ", clust_i, 
    " | ", sample_gene, " | ",
    "AUC: ", 
    formatC(AUC_total_df[sample_gene, "AUC_value"], digits = 3)
  )
  
  plot_fitted_values(
    sample_gene,
    plot_title =  plot_title
  )
}
dev.off()
par(old_pars)

set.seed(Sys.time())