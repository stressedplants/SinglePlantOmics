
# Find optimal fits AFTER IQR FILTERING =============

range_of_timeseries <- c(1, number_of_timepoints)
penalty_basis <- create.bspline.basis(rangeval = range_of_timeseries,
                                      nbasis = optimal_nbasis)

# Define a linear differential operator

second_deriv_operator = int2Lfd(2)

fdParobj = fdPar(penalty_basis, second_deriv_operator, optimal_lambda)

filtered_fd <- smooth.basis(1:number_of_timepoints,
                            TPM_for_smoothing[, keep_after_IQR],
                            fdParobj)$fd

# Functions for scaling FDA objects =======

# This function takes an fd_obj and rescales each of the fits, so that the 
# smooth functions have a min of 0 and max of 1.
rescale_fd_to_01 <- function(fd_obj,
                             original_data,
                             fdParobj_in,
                             timepoints = timepoint_vec,
                             make_comparison_plot = FALSE) {
  
  curve_evals <- eval.fd(seq(from = timepoints[1],
                             to = tail(timepoints, 1),
                             by = 0.05),
                         fd_obj)
  
  min_max_fda <- apply(curve_evals, 2, function(in_col) {
    return(c(min(in_col), max(in_col) - min(in_col)))
  })
  
  # rescale data - then make a new FDA fit
  # this was easier than scaling the functions directly in FDA package
  scaled_data <- original_data
  
  for (col_num in 1:ncol(scaled_data)) {
    scaled_data[, col_num] <- (original_data[, col_num] - 
                                 min_max_fda[1, col_num]) / 
      (min_max_fda[2, col_num])
  }
  
  scaledfd <- smooth.basis(timepoints,  # time values
                           scaled_data,
                           fdParobj_in)$fd
  
  if (make_comparison_plot) {
    par(mfrow = c(1,2))
    indices_to_plot = sample(1:ncol(scaled_data), 10)
    plot(fd_obj[indices_to_plot], main = "original")
    plot(scaledfd[indices_to_plot], main = "scaled to 0-1")
    par(mfrow = c(1,1))
  }
  
  return(scaledfd)
}

rescaled_fd <- rescale_fd_to_01(
  filtered_fd,
  TPM_for_smoothing[, keep_after_IQR],
  fdParobj,
  1:number_of_timepoints
)

# How many genes are non-monotonic (have a turning point)? =========

find_if_monotonic <- function(fd_obj,
                              soft_threshold = 1e-3,
                              rng = range_of_timeseries,  # might want to trim
                              tp_interval = 0.5) {
  # find the x values where the y value is 1/2
  fine_timepoints <- seq(from = rng[1], to = rng[2], by = tp_interval)
  d_eval <- eval.fd(fine_timepoints, deriv.fd(fd_obj))
  
  d_eval[(d_eval < soft_threshold) & (d_eval > -soft_threshold)] <- 0
  d_eval <- d_eval[!(d_eval == 0)]
  
  # find the crossing points by finding differences in *sign*
  idx <- which(abs(diff(sign(d_eval))) == 2)
  monotonic <- (length(idx) == 0)
  
  return(monotonic)
}

test_all_monotonic <- sapply(1:length(rescaled_fd$fdnames$reps), function(i){
  find_if_monotonic(rescaled_fd[i, ], 
                    rng = c(5,60))
})

print(paste0(sum(test_all_monotonic), " out of ", num_after_IQR,
             " genes are monotonic over pseudotime"))

monotonic_genes <- rescaled_fd$fdnames$reps[test_all_monotonic]

# Find AUC statistics for monotonic genes - early / late changing genes ======

integrate_fd_flip <- function(fd_obj,
                              rng = range_of_timeseries,  # c(1, 65)
                              make_plot = FALSE) {
  
  # can't pass fd_obj directly to integrate, need another fn
  eval_fn <- function(eval_pts) {
    return(eval.fd(eval_pts, fd_obj))
  }
  
  integral_value <- integrate(eval_fn, 
                              lower = rng[1],
                              upper = rng[2])$value
  
  end_vals <- eval.fd(rng, fd_obj)
  end_minus_beginning <- tail(end_vals, 1) - head(end_vals, 1)
  
  # if increasing, return area OVER the curve!
  if (end_minus_beginning > 0) {
     integral_value <- (rng[2] - rng[1]) - integral_value
  } 
  
  if (make_plot) {
    plot(fd_obj, sub = paste0("Integral: ", sprintf("%.2f", integral_value)))
  }
  
  return(integral_value)
}


# Pre-processing for AUC values and monotonicity =========

# create increasing / decreasing + monotonic lists -
# based on whether the initial value is less than the end value

monotonic_dir_classifier <- sapply(monotonic_genes, function(gname) {
  temp_evals <- eval.fd(range_of_timeseries, rescaled_fd[gname])
  return(temp_evals[1] > temp_evals[2])
})

decreasing_monotonic <- names(monotonic_dir_classifier[monotonic_dir_classifier])
increasing_monotonic <- names(monotonic_dir_classifier[!monotonic_dir_classifier])

# E.g. to check that the AUC calculation is working and accounting for 
# increasing / decreasing genes correctly ...
# integrate_fd_flip(rescaled_fd[sample(decreasing_monotonic, 1), ],
#                   make_plot = TRUE)

all_increasing_AUC <- sapply(increasing_monotonic, function(gene_id){
  integrate_fd_flip(rescaled_fd[gene_id])
})
all_increasing_AUC_df <- data.frame(gene_id = increasing_monotonic,
                                    AUC_value = all_increasing_AUC)


all_decreasing_AUC <- sapply(decreasing_monotonic, function(gene_id){
  integrate_fd_flip(rescaled_fd[gene_id])
})
all_decreasing_AUC_df <- data.frame(gene_id = decreasing_monotonic,
                                    AUC_value = all_decreasing_AUC)


nonmonotonic_genes <- setdiff(keep_after_IQR, monotonic_genes)
all_nonmonotonic <- sapply(nonmonotonic_genes, 
  function(gene_id){
    integrate_fd_flip(rescaled_fd[gene_id])
  }
)
all_nonmonotonic_AUC_df <- data.frame(gene_id = nonmonotonic_genes,
                                      AUC_value = all_nonmonotonic)

AUC_total_df <- rbind(all_decreasing_AUC_df, all_increasing_AUC_df, all_nonmonotonic_AUC_df)
row.names(AUC_total_df) <- AUC_total_df$gene_id

# Save output of AUC for supplemental data =====

AUC_output_df <- AUC_total_df
AUC_output_df <- AUC_output_df[order(AUC_total_df$gene_id), ]

# add a column to say either "increasing", "decreasing" or "nonmonotonic"
AUC_output_df$monotonicity <- rep("", nrow(AUC_output_df)) 
AUC_output_df[increasing_monotonic, "monotonicity"] <- "increasing"
AUC_output_df[decreasing_monotonic, "monotonicity"] <- "decreasing"
AUC_output_df[nonmonotonic_genes, "monotonicity"] <- "nonmonotonic"

write.xlsx(AUC_output_df, file = "outputs/pseudotime_analysis/AUC_values.xlsx",
           row.names = FALSE)
# Plots of AUC by GO term ========

produce_AUC_df_GO <- function(
    go_ids,
    filter_list = monotonic_genes,  # by default, filter for monotonic genes
    verbose = FALSE
) {
  
  num_go <- length(go_ids)
  
  stat_df <- data.frame(matrix(ncol = 3, nrow = 0))
  stat_colnames <- c("GO", "gene", "statistic")
  colnames(stat_df) <- stat_colnames
  
  genes_associated <- get_genes_associated_GO_terms(
    go_ids,
    filter_list
  )
  
  # go through each GO term, and add the statistic related to each gene
  
  for (i_go in 1:num_go) {
    
    gene_list <- genes_associated[[i_go]]
    
    # don't update the large dataframe if no genes are appropriate
    if (!is_empty(gene_list)) {
      stat_values <- sapply(gene_list, function(gname){
        return(AUC_total_df[gname, "AUC_value"])
      })
      
      for (i_gene in 1:length(gene_list)) {
        stat_df[nrow(stat_df) + 1, ] <- c(
          go_ids[i_go],
          gene_list[i_gene],
          stat_values[i_gene]
        )  
      }
    }
    
    # verbose if needed
    if (verbose & ((i_go %% 10 == 0) || (i_go == num_go))) {
      print(paste0("Completed ", i_go, " out of ", num_go, " GO terms"))
    }
    
  }
  
  stat_df$GO <- as.factor(stat_df$GO)
  stat_df$statistic <- as.numeric(stat_df$statistic)
  
  return(stat_df)
}

# this does not include nonmonotonic genes!
produce_grouped_AUC_df <- function(
    go_ids,
    do_coord_flip = TRUE,
    verbose = FALSE
) {

  inc_df <- produce_AUC_df_GO(
    go_ids,
    filter_list = increasing_monotonic,
    verbose = verbose
  )
  
  inc_df$dir <- rep("increasing", nrow(inc_df))
  
  dec_df <- produce_AUC_df_GO(
    go_ids,
    filter_list = decreasing_monotonic,
    verbose = verbose
  )
  
  dec_df$dir <- rep("decreasing", nrow(dec_df))
  
  auc_w_direction <- rbind(inc_df, dec_df)
  
  return(auc_w_direction)
  
}

# SNIPPET TO FILTER GO BY NUMBER OF GENES
#  %>% group_by(GO) %>% filter(n_distinct(gene) >= minimum_number_genes)

# Characterise AUC values for all GO highlighted in DGE
possible_GO_terms <- unique(gprofiler_query_DEGs$result$term_id)

# Filter these by GO:MF, GO:CC, and GO:BP
mol_func_terms <- possible_GO_terms[
  possible_GO_terms %in% GO_names_df$go_id[
    GO_names_df$namespace_1003 == "molecular_function"
  ]
]

mol_func_df <- produce_grouped_AUC_df(
  mol_func_terms,
  verbose = TRUE
)

cell_comp_terms <- possible_GO_terms[
  possible_GO_terms %in% GO_names_df$go_id[
    GO_names_df$namespace_1003 == "cellular_component"
  ]
]

cell_comp_df <- produce_grouped_AUC_df(
  cell_comp_terms,
  verbose = TRUE
)

bio_proc_terms <- possible_GO_terms[
  possible_GO_terms %in% GO_names_df$go_id[
    GO_names_df$namespace_1003 == "biological_process"
  ]
]

bio_proc_df <- produce_grouped_AUC_df(
  bio_proc_terms,
  verbose = TRUE
)

filter_GO_df <- function(
    input_GO_df,
    min_genes = 100,
    max_genes = 4000,
    dir_filter = "decreasing",
    replace_GO_by_name = FALSE    
) {
  
 # filter by decreasing and minimum number of distinct genes
  filtered_GO_data <- input_GO_df %>% 
    filter(dir == dir_filter) %>% 
    group_by(GO) %>% 
    filter(n_distinct(gene) >=  min_genes) %>% 
    filter(n_distinct(gene) <  max_genes) %>% 
    ungroup()
  
  filtered_GO_data$GO <- as.character(filtered_GO_data$GO)
  
  # this relies on GO_named linking the GO term IDs and names
  if (replace_GO_by_name) {
    new_names <- sapply(filtered_GO_data$GO, function(go_id){
      return(GO_names_df[go_id, "name_1006"])
    })
    filtered_GO_data$GO <- new_names
  }
  
  return(filtered_GO_data)
}
  
plot_GO_violin_filtering <- function(
    input_GO_df,
    min_genes = 100,
    max_genes = 4000,
    dir_filter = "decreasing",
    replace_GO_by_name = TRUE,
    do_kruskal_wallis = TRUE,
    kw_signif_value = 0.01,
    do_post_hoc = TRUE,
    use_chisquare = TRUE
){

  filtered_GO_data <- filter_GO_df(
    input_GO_df,
    min_genes = min_genes,
    max_genes = max_genes,
    dir_filter = dir_filter,
    replace_GO_by_name = replace_GO_by_name
  )
  
  # now order so that median value is increasing
  combined_GO_AUC_df <- filtered_GO_data %>%
    group_by(GO) %>%
    mutate(median_statistic = median(statistic)) %>%
    ungroup() %>%
    arrange(factor(GO, levels = unique(GO[order(median_statistic)])))
  
  kw_output <- NA
  ph_output <- NA
  # if appropriate, set these values to something real
  if (do_kruskal_wallis) {
    kw_output <- kruskal.test(filtered_GO_data$statistic,
                              as.factor(filtered_GO_data$GO))
    if (do_post_hoc & (kw_output$p.value < kw_signif_value)) {
      ph_output <- kwAllPairsNemenyiTest(
        filtered_GO_data$statistic,
        as.factor(filtered_GO_data$GO),
        dist = ifelse(use_chisquare, "Chisquare", "Tukey")
      )
    }
  }

  combined_GO_AUC_df$GO <- factor(
    combined_GO_AUC_df$GO,
    levels = unique(combined_GO_AUC_df$GO[order(combined_GO_AUC_df$median_statistic)])
  )

  GO_violin <- ggplot(data = combined_GO_AUC_df,
                      mapping = aes(GO, statistic)) +
    geom_violin() +
    geom_jitter(position = position_dodge(.9)) +
    coord_flip() +
    theme_minimal() +
    ylab("Integral") +
    xlab("GO term name") +
    theme(plot.title = element_text(hjust = 0.5))
  
  output <- list("plot" = GO_violin, 
                 "kw" = kw_output,
                 "ph" = ph_output)
  
  return(output)
}

# Load the thresholds of GO terms from file, then produce plots ========

make_GO_size_table <- function(go_auc_df,
                               dir_filter = "decreasing",
                               keep_zero = FALSE) {
  temp_df <- go_auc_df %>% filter(dir == dir_filter)
  out_tbl <- table(temp_df$GO)
  
  if (!keep_zero) {
    out_tbl <- out_tbl[out_tbl >= 1]  
  }
  
  return(out_tbl)
}

go_size_thresholds_df <- read.csv("data/go_size_thresholds.csv")
go_category_to_df <- list(
  "molecular_function" = mol_func_df,
  "biological_process" = bio_proc_df,
  "cellular_component" = cell_comp_df
)

# cycle through possible directions and GO categories and make appropriate plots
for (direction in unique(go_size_thresholds_df$direction)) {
  
  svg(filename = paste0("plots/pseudotime_analysis/",
                        direction,
                        "_GO_thresholds.svg"),
      width = 10, height = 7, pointsize = 18)
  
  old_plotting_par <- par(mfrow = c(1, 3), xaxt = 'n', 
                          cex.lab = 1.5, cex.main = 1.5)
  
  first_plot <- TRUE
  
  for (go_cat in unique(go_size_thresholds_df$go_type)) {

    go_cat_tbl <- make_GO_size_table(go_category_to_df[[go_cat]], 
                                     dir_filter = direction)
    plot(rep(10, length(go_cat_tbl)), as.numeric(go_cat_tbl), 
         log = 'y', 
         xlab = "",
         ylab = ifelse(first_plot, "Size of GO term", ""),
         main = go_cat)
    
    high_thresh <- go_size_thresholds_df[
      (go_size_thresholds_df$go_type == go_cat) & 
      (go_size_thresholds_df$direction == direction),
      "high"
    ]
    med_thresh <- go_size_thresholds_df[
      (go_size_thresholds_df$go_type == go_cat) & 
        (go_size_thresholds_df$direction == direction),
      "medium"
    ]
    
    abline(h = high_thresh, col = "red")
    abline(h = med_thresh, col = "purple")
    
    first_plot <- FALSE
  }
  
  par(old_plotting_par)
  dev.off()
}

# Produce plots and xlsx outputs from GO term AUC violins ==========

# cycle through possible directions and GO categories and make appropriate plots
for (direction in unique(go_size_thresholds_df$direction)) {
  
  stats_wb = createWorkbook()

  for (go_cat in unique(go_size_thresholds_df$go_type)) {
    short_name = paste0(go_cat, "-", direction)
    high_thresh <- go_size_thresholds_df[
      (go_size_thresholds_df$go_type == go_cat) & 
        (go_size_thresholds_df$direction == direction),
      "high"
    ]
    med_thresh <- go_size_thresholds_df[
      (go_size_thresholds_df$go_type == go_cat) & 
        (go_size_thresholds_df$direction == direction),
      "medium"
    ]
    
    # make the high threshold plot
    GO_high_plot_etc <- plot_GO_violin_filtering(
      go_category_to_df[[go_cat]],
      min_genes = high_thresh,
      dir_filter = direction
    )
    ggsave(
      paste0("plots/pseudotime_analysis/", short_name,"-high.svg"),
      plot = GO_high_plot_etc$plot + ggtitle(paste0(go_cat, " - ", direction, " - high")),
      width = 7, height = 7
    )
    
    
    
    # add the high threshold p values
    stats_sheet = createSheet(stats_wb, sheetName = paste0(short_name, "-high"))
    addDataFrame(data.frame(kw_p_value = GO_high_plot_etc[["kw"]]$p.value),
                 sheet = stats_sheet,
                 startColumn = 1,
                 row.names = FALSE)
    if (class(GO_high_plot_etc[["ph"]]) == "PMCMR") {
      addDataFrame(GO_high_plot_etc[["ph"]]$p.value,
                   sheet = stats_sheet,
                   startColumn = 3)
    }
    
    # make the medium threshold plot
    GO_med_plot_etc <- plot_GO_violin_filtering(
      go_category_to_df[[go_cat]],
      max_genes = high_thresh,
      min_genes = med_thresh,
      dir_filter = direction
    )
    ggsave(
      paste0("plots/pseudotime_analysis/", go_cat, "-", direction, "-medium.svg"),
      plot = GO_med_plot_etc$plot + ggtitle(paste0(go_cat, " - ", direction, " - medium")),
      width = 7, height = 7
    )
    
    # add the medium threshold p values
    stats_sheet = createSheet(stats_wb, sheetName = paste0(short_name, "-medium"))
    addDataFrame(data.frame(kw_p_value = GO_med_plot_etc[["kw"]]$p.value),
                 sheet = stats_sheet,
                 startColumn = 1,
                 row.names = FALSE)
    if (class(GO_med_plot_etc[["ph"]]) == "PMCMR") {
      addDataFrame(GO_med_plot_etc[["ph"]]$p.value,
                   sheet = stats_sheet,
                   startColumn = 3)
    }
  }
  saveWorkbook(stats_wb, paste0("outputs/pseudotime_analysis/", direction,".xlsx"))
}

# Make some finalised plots for the 'interesting' GO terms =====

mf_of_interest <- c("GO:0005515", "GO:0003735")
cc_of_interest <- c("GO:0005730", "GO:0005737", "GO:0009507", "GO:0009506")
bp_of_interest <- c("GO:0042254", "GO:0006397", "GO:0009793", "GO:0016310")

mf_temp_df <- filter_GO_df(mol_func_df[mol_func_df$GO %in% mf_of_interest, ],
                           min_genes = 1, replace_GO_by_name = TRUE)
mf_temp_df$category <- "molecular function"
bp_temp_df <- filter_GO_df(bio_proc_df[bio_proc_df$GO %in% bp_of_interest, ],
                           min_genes = 1, replace_GO_by_name = TRUE)
bp_temp_df$category <- "biological process"
cc_temp_df <- filter_GO_df(cell_comp_df[cell_comp_df$GO %in% cc_of_interest, ],
                           min_genes = 1, replace_GO_by_name = TRUE)
cc_temp_df$category <- "cellular component"

combined_plot_df <- rbind(mf_temp_df,
                          bp_temp_df,
                          cc_temp_df)
combined_plot_df <- combined_plot_df[, c("GO", "statistic", "category")]

combined_plot_df <- combined_plot_df %>%
  group_by(GO) %>%
  mutate(mean_statistic = mean(statistic)) %>%
  ungroup() %>%
  arrange(factor(GO, levels = unique(GO[order(mean_statistic)])))

combined_plot_df$GO <- factor(
  combined_plot_df$GO,
  levels = unique(combined_plot_df$GO[order(combined_plot_df$mean_statistic)])
)

combined_plot <- ggplot(data = combined_plot_df,
       mapping = aes(GO, statistic, color = category)) +
  facet_grid(rows = vars(category), scales = 'free') +
  geom_violin() +
  geom_beeswarm(size = 0.5, cex = 0.8) +
  coord_flip() +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none") +
  ylab("Area under the curve") +
  xlab("GO term name")
  
ggsave("plots/pseudotime_analysis/combinded_violin.svg",
       plot = combined_plot,
       device = "svg",
       width = 10,
       height = 7,
       pointsize = 18)

