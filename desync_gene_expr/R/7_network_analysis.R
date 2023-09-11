# DAP-seq data loading ==============

# First, need to produce this network by looping through all folders and reading
# in the separate link lists. Combine these lists.

DAP_seq_network_file <- "outputs/network_analysis/DAP_seq_link_list.csv"

if (!file.exists(DAP_seq_network_file)) {
  produce_DAP_seq_df <- function(verbose = FALSE){
    
    DAP_seq_files <- list.files('data/dap_seq/',recursive = TRUE)
    DAP_seq_files <- sapply(DAP_seq_files, function(in_file){
      paste0('data/dap_seq/', in_file)
    })
    
    out_df <- data.frame(regulator = 'E.g.',
                         target = 'E.g.')
    # could search to find duplicates / use a different structure
    # if this is inefficient
    
    for (file_name in DAP_seq_files) {
      
      if (verbose) print(paste0('Loading ', file_name))
      
      temp_tbl <- read.delim(file_name)
      names(temp_tbl) <- c('regulator', 'target')
      
      out_df <- rbind(out_df, temp_tbl)
    }
    
    if (verbose) print('Removing duplicates')
    
    before_len <- dim(out_df)[1]
    out_df <- unique(out_df)
    after_len <- dim(out_df)[1]
    
    if (verbose) print(paste0('Removed ', (before_len - after_len), ' rows'))
    
    # also remove e.g.
    out_df <- out_df[-1,]
    
    return(out_df)
    
  }
  
  DAP_ll <- produce_DAP_seq_df(verbose = TRUE)
  write.csv(DAP_ll, DAP_seq_network_file)
  
  gc()
} else {
  DAP_ll <- read.csv(DAP_seq_network_file, row.names = 1)
}

TFs_in_DAP <- unique(DAP_ll$regulator)
g_DAP_seq <- graph.data.frame(DAP_ll)
rm(DAP_ll)  # use the iGraph object instead

# Prepare for DynGENIE3 ======

time_vec_dynGENIE3 <- 1:number_of_timepoints
genes_in_dynGENIE3 <- rescaled_fd$fdnames$reps
regs_list <- TFs_in_DAP[TFs_in_DAP %in% genes_in_dynGENIE3]

source("external_scripts/dynGENIE3_compiled/dynGENIE3.R")

link_list_filepath_ns <- 'outputs/network_analysis/link_list_all_dynGENIE3_nonsmoothed.csv'

TPM_dynGENIE3_ns <- TPM_for_smoothing[pseudotime_output, genes_in_dynGENIE3]

# these variables need to be contained in lists for DynGENIE3
time_points_dynGENIE3_ns <- list(time_vec_dynGENIE3)
expr_data_dynGENIE3_ns <- list(t(as.matrix(TPM_dynGENIE3_ns)))

if (!file.exists(link_list_filepath_ns)) {
  print("Running DynGENIE3 on unsmoothed data")
  # this is intensive (take ~20 mins on my laptop)
  dynGENIE3_output_ns <- dynGENIE3(expr_data_dynGENIE3_ns, time_points_dynGENIE3_ns,
                                   regulators = regs_list,
                                   verbose = TRUE, seed = 123)
  
  # this is not the final threshold! just a low threshold in order to not save a 
  # massive data frame
  ll_dynGENIE3_ns <- get.link.list(dynGENIE3_output_ns$weight.matrix, 
                                   threshold = 1e-10)
  
  write.csv(ll_dynGENIE3_ns, link_list_filepath_ns, row.names = F)
  
  gc()
} else {
  ll_dynGENIE3_ns <- read.csv(link_list_filepath_ns)
}

# set the number of edges threshold for other methods - to compare
n_edges_ns <- length(ll_dynGENIE3_ns$weight)
n_edges_threshold_ns <- ceiling(5/100 * n_edges_ns)
weight_threshold_ns <- ll_dynGENIE3_ns$weight[n_edges_threshold_ns]

# Plot a sample of the link list
sample_indicies <- seq(from = 1, to = n_edges_ns, by = 200)

svg("plots/network/DynGENIE3_context.svg",
    width = 7, height = 5)
plot(sample_indicies, ll_dynGENIE3_ns$weight[sample_indicies],
     xlab = 'Index of link list', ylab = 'Weight', type = "l")
abline(h = weight_threshold_ns, col = "red")
title("Cut off for edge weight of DynGENIE3 in context")
dev.off()

# Filter edges by threshold and then create an iGraph object for further use
top_links_dynGENIE3_ns <- ll_dynGENIE3_ns[
  (ll_dynGENIE3_ns$weight >= weight_threshold_ns),
]
g_dynGENIE3_ns <- graph.data.frame(top_links_dynGENIE3_ns)

print(paste0(
  'Keep ', 
  dim(top_links_dynGENIE3_ns)[1], 
  ' edges in the network (no smoothing) before DAP-seq verification'
))

print(paste0(
  'Keep ', 
  length(unique(top_links_dynGENIE3_ns[, 1])),
  ' regulators in the network (no smoothing) before DAP-seq verification'
))

print(paste0(
  'Keep ',
  length(unique(top_links_dynGENIE3_ns[, 2])),
  ' targets in the network (no smoothing) before DAP-seq verification'
))

# Clean up after loading data =====

gc()

# Output a graphml file for plotting of VERIFIED GRN ======

# intersects graph_1 with graph_2
# returns the largest component of graph_1 that remains
intersect_and_reduce_graph <- function(graph_1, graph_2) {
  graph_out <- intersection(graph_1, graph_2,
                            keep.all.vertices = FALSE)
  components <- clusters(graph_out, mode = "weak")
  biggest_cluster_id <- which.max(components$csize)
  vert_ids <- V(graph_out)[components$membership == biggest_cluster_id]
  graph_out <- induced_subgraph(graph_out, vert_ids)
  return(graph_out)
}

# Create a subsetted version with intersection of DAP-seq and dynGENIE3
g_dynGENIE3_verified <- intersect_and_reduce_graph(g_dynGENIE3_ns, g_DAP_seq)

# save a copy of this as a csv
save_dynGENIE3_df <- as_long_data_frame(g_dynGENIE3_verified)[, 3:5]
colnames(save_dynGENIE3_df) <- c("weight", "regulatory gene", "target gene")
write.csv(save_dynGENIE3_df, 
          file = "outputs/network_analysis/g_dynGENIE3_verified.csv",
          row.names = FALSE, quote = FALSE)

# Manually correct unlabelled regulatory genes - load from file
# (where biomaRt didn't provide a short gene name)

relabelled_genes <- read.csv("data/corrected_gene_names.csv")
colnames(relabelled_genes) <- c("id", "name")

# include a combined name / id column
combined_unlabelled <- rep("", nrow(relabelled_genes))
for (i in 1:nrow(relabelled_genes)) {
  combined_unlabelled[i] <- paste0(relabelled_genes[i, "name"], " (",
                                relabelled_genes[i, "id"], ")")
}
relabelled_genes$combined <- combined_unlabelled

# use this to edit the id_to_name_df for further use
id_to_name_df_network <- id_to_name_df
id_to_name_df_network[relabelled_genes$id, 1] <- relabelled_genes$name
id_to_name_df_network[relabelled_genes$id, 2] <- relabelled_genes$combined

# Add additional info to iGraph =======

vertex_list <- as.vector(V(g_dynGENIE3_verified)$name)
auc_exact_dynGENIE3 <- AUC_total_df[vertex_list, "AUC_value"]

# first add the AUC values - with scaling within up or down

# TRUE if beginning value is larger than end value, FALSE otherwise
classify_down <- sapply(vertex_list, function(gname) {
  print(gname)
  temp_evals <- eval.fd(range_of_timeseries, rescaled_fd[gname])
  return(temp_evals[1] > temp_evals[2])
})

down_simple <- names(classify_down[classify_down])
up_simple <- names(classify_down[!classify_down])

# create the z-scored AUC value for up and down regs
down_z_scores_df <- AUC_total_df[AUC_total_df$gene_id %in% down_simple, ]
down_z_scores <- as.vector(scale(down_z_scores_df$AUC_value))
names(down_z_scores) <- down_z_scores_df$gene_id

up_z_scores_df <- AUC_total_df[AUC_total_df$gene_id %in% up_simple, ]
up_z_scores <- as.vector(scale(up_z_scores_df$AUC_value))
names(up_z_scores) <- up_z_scores_df$gene_id

scaled_auc_score <- sapply(vertex_list, function(vname){
  if (vname %in% down_simple) {
    return(down_z_scores[vname])
  } else {
    return(up_z_scores[vname])
  }
})

down_rank <- rank(down_z_scores_df$AUC_value)
names(down_rank) <- down_z_scores_df$gene_id

up_rank <- rank(up_z_scores_df$AUC_value)
names(up_rank) <- up_z_scores_df$gene_id

rank_auc_score <- sapply(vertex_list, function(vname){
  if (vname %in% down_simple) {
    return(down_rank[vname])
  } else {
    return(up_rank[vname])
  }
})

# add an extra rank variable that only for decreasing or increasing TFs
# also scale so that down and up are comparable (between 0 and 1)

down_TFs <- intersect(regs_list, down_simple)
up_TFs <- intersect(regs_list, up_simple)

rank_only_down_TFs <- AUC_total_df[down_TFs, "AUC_value"]
rank_only_down_TFs <- rank(rank_only_down_TFs)
rank_only_down_TFs <- rank_only_down_TFs / length(rank_only_down_TFs)
names(rank_only_down_TFs) <- down_TFs

rank_only_up_TFs <- AUC_total_df[up_TFs, "AUC_value"]
rank_only_up_TFs <- rank(rank_only_up_TFs)
rank_only_up_TFs <- rank_only_up_TFs / length(rank_only_up_TFs)
names(rank_only_up_TFs) <- up_TFs

rank_only_TFs_AUC_score <- sapply(vertex_list, function(vname){
  if (vname %in% down_TFs) {
    return(rank_only_down_TFs[vname])
  } else {
    if (vname %in% up_TFs) {
      return(rank_only_up_TFs[vname])
    } else {
      return(NA)
    }
  }
})

# dummy variables for color, position, and size of vertices

vertices_dummy_col <- sapply(vertex_list, function(vname){
  if (vname %in% regs_list) {  # regulators should be in -1 or 1
    if (vname %in% down_simple) {
      return(-1)
    } else {
      return(1)
    }
  } else {
    return(0)
  }
})

set.seed(123)
jitter_amt <- 0.2
vertices_dummy_position <- sapply(vertex_list, function(vname){
  if (vname %in% regs_list) {  # regulators should be in -1 or 1
    if (vname %in% down_simple) {
      return(-1)
    } else {
      return(1)
    }
  } else {
    return(runif(1, min = -jitter_amt, max = jitter_amt))
  }
})
set.seed(Sys.time())

vertices_dummy_size <- sapply(vertex_list, function(vname){
  if (vname %in% regs_list) {  # regulators should be in -1 or 1
    return(1)
  } else {
    return(0)
  }
})

#  change the label so that only TFs have a label - this will be useful in gephi
vertex_labels <- sapply(vertex_list, function(vname){
  if (vname %in% regs_list) {
    return(id_to_name_df_network[vname, "gene_name"])
  } else {
    return(" ")
  }
})

# get a table of the different TF families
regs_TF_families <- TF_family_df[regs_list, ]
rownames(regs_TF_families) <- regs_list

# check all are included in the df
no_TF_regs <- regs_list[!(regs_list %in% unique(TF_family_df$Gene_ID))]

for (un_reg in no_TF_regs) {
  regs_TF_families[un_reg, "Gene_ID"] <- un_reg
  regs_TF_families[un_reg, "Family"] <- "Unknown"
}

regs_TF_families$Family <- as.factor(regs_TF_families$Family)

# add these to the edges of the graph, by first assigning families
TF_family_by_vertex <- sapply(
  as.vector(names(V(g_dynGENIE3_verified))),
  function(gene_id) {
    if (gene_id %in% rownames(regs_TF_families)) {
      return(as.character(regs_TF_families[gene_id, "Family"]))
    } else {
      return("NA")
    }
  }
)

g_dynGENIE3_labelled <- g_dynGENIE3_verified %>%
  set_vertex_attr("gene_name", value = id_to_name_df_network[vertex_list, 1]) %>%
  set_vertex_attr("id", value = vertex_list) %>%
  set_vertex_attr("label", value = vertex_labels) %>%
  set_vertex_attr("scaled_auc", value = scaled_auc_score) %>%
  set_vertex_attr("rank_auc", value = rank_auc_score) %>%
  set_vertex_attr("rank_only_TFs_auc", value = rank_only_TFs_AUC_score) %>%
  set_vertex_attr("dummy_pos", value = vertices_dummy_position) %>%
  set_vertex_attr("dummy_col", value = vertices_dummy_col) %>%
  set_vertex_attr("dummy_size", value = vertices_dummy_size) %>%
  set_edge_attr("weight", value = rep(1, length(E(g_dynGENIE3_verified)))) %>%
  set_vertex_attr("TF_family", value = TF_family_by_vertex)

write_graph(g_dynGENIE3_labelled,
            'outputs/network_analysis/g_dynGENIE3_verified.graphml',
            format = 'graphml')

# Summary stats of the verified GRN ==========

# how many genes regulated by only decreasing TFs / increasing / both ?

down_targs <- adjacent_vertices(g_dynGENIE3_verified, 
                                intersect(down_simple, regs_list),
                                mode = "out")
down_targs <- unname(unlist(down_targs))

up_targs <- adjacent_vertices(g_dynGENIE3_verified, 
                              intersect(up_simple, regs_list),
                              mode = "out")
up_targs <- unname(unlist(up_targs))


print(paste0("Only decreasing TF targets: ",
             length(unique(setdiff(down_targs, up_targs)))))

print(paste0("Only increasing TF targets: ",
             length(unique(setdiff(up_targs, down_targs)))))

print(paste0("Both increasing and decreasing TF targets: ",
             length(unique(intersect(up_targs, down_targs)))))

# Plots of each of the up- and down- TFs over pseudotime ==========

cex_scale = 2

svg(filename = "plots/network/up_TFs_pseudotime.svg",
    width = 11, height = 7)
plot(rescaled_fd[up_TFs], lty = 1, col = "#5733FF", lwd = 2, 
     cex.lab = cex_scale, cex.axis = cex_scale, cex.main = cex_scale,
     family = "Arial")
dev.off()

svg(filename = "plots/network/down_TFs_pseudotime.svg",
    width = 11, height = 7)
plot(rescaled_fd[down_TFs], lty = 1, col = "#FF5733", lwd = 2, 
     cex.lab = cex_scale, cex.axis = cex_scale, cex.main = cex_scale,
     family = "Arial")
dev.off()
