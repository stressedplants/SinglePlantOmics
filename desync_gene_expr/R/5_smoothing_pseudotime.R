######################################
# Find optimal values for number of basis functions and smoothing parameters
######################################

# Functions used to produce a smoothed function ===============

set.seed(123) # for reproducibility

# - Separate timeseries should be stored in columns - used by FDA this way
#   sample with very high / low expression in one gene)
# - Use 0 to 1 normalisation in order for GCV to be consistent across samples
TPM_for_smoothing <- t(TPM_pseudotime)
number_of_timepoints <- nrow(TPM_for_smoothing)

#' Fit a spline to all genes, ordered over pseudotime
#'
#' @param nbasis_in The number of basis functions to use
#' @param lambda_in The value of the smoothing parameter, lambda
#'
#' @return n fdSmooth object, with fits to every gene and summary parameters.
#'  Usually, this output would be immediately used to find sum(output$gcv).
spline_from_params <- function(nbasis_in, lambda_in) {
  range_of_timeseries <- c(1, number_of_timepoints)
  penalty_basis <- create.bspline.basis(rangeval = range_of_timeseries,
                                        nbasis = nbasis_in)
  
  # Define a linear differential operator
  second_deriv_operator = int2Lfd(2)
  
  fdParobj = fdPar(penalty_basis, second_deriv_operator, lambda_in)
  optimal_outputs = smooth.basis(1:number_of_timepoints,
                             TPM_for_smoothing, fdParobj)
  return(optimal_outputs)

}

# Find optimal number of basis functions and smoothing parameters =======

test_nbasis_lambda <- function(nbasis_vector, lambda_vector) {
  cat("Running B-spline smoothing with different seeds\n")
  gcv_for_params <- expand.grid(nbasis_vector, lambda_vector)
  names(gcv_for_params) <- c("nbasis", "lambda")
  gcv_for_params$gcv_sum <- 0
  
  for (param_num in 1:(nrow(gcv_for_params))) {
    nbasis_param <- gcv_for_params[param_num, "nbasis"]
    lambda_param <- gcv_for_params[param_num, "lambda"]
    
    gcv_for_params[param_num, "gcv_sum"] <- sum(
      spline_from_params(nbasis_param, lambda_param)$gcv
    )
    
    # progress tracker
    if (param_num %% 30 == 0) cat(paste0("Completed: ", param_num, "\n"))
  }
  
  return(gcv_for_params)
}


different_bases_file <- "outputs/pseudotime/GCV_over_different_bases.RData"

nbasis_vec <- seq(from = 5, to = 10, by = 1)
lambda_vec <- 10**seq(from = 0, to = 3.5, length.out = 50)

# avoid running if the file already exists
if (!file.exists(different_bases_file)) {
  cat("Testing different smoothing hyperparameters\n")
  gcv_for_params = test_nbasis_lambda(nbasis_vec, lambda_vec)
  save(gcv_for_params, file = different_bases_file)
} else {
  cat("Loading B-spline gcv outputs\n")
  load(different_bases_file)
}

# Plots related to different GCV values ========

smoothing_folder <- "plots/pseudotime/"

create_gcv_matrix <- function(gcv_param_df) {
  gcv_outputs_wide <- gcv_param_df %>% spread("lambda", "gcv_sum")
  row.names(gcv_outputs_wide) <- gcv_outputs_wide$nbasis
  gcv_outputs_wide <- within(gcv_outputs_wide, rm("nbasis"))

  return(gcv_outputs_wide)
}

gcv_first_matrix <- create_gcv_matrix(gcv_for_params)

write.csv(gcv_first_matrix,
          file = "outputs/pseudotime/gcv_sums_different_params.csv")

gcv_lambda_plot <- function(gcv_results_mat,
                            lambda_vec_in) {
  gcv_for_plot <- data.frame(t(gcv_results_mat))
  gcv_for_plot$lambda <- lambda_vec_in
  
  gcv_for_plot <- melt(gcv_for_plot, id.vars = "lambda")
  
  gcv_for_plot$variable <- as.factor(sapply(gcv_for_plot$variable, 
                                             function(ch_in) {
                                               as.numeric(gsub("X", "", ch_in))
                                             }))
  colnames(gcv_for_plot)[2] <- "n_basis"
  out_plot <- ggplot(gcv_for_plot,
         aes(x = lambda, y = value, colour = n_basis)) + 
    scale_x_log10() +
    theme_minimal(base_size = 18) +
    geom_line() +
    xlab("Smoothing penalty (lambda)") +
    ylab("Sum of GCV")
  
  return(out_plot)
}

gcv_first_lambdas <- gcv_lambda_plot(gcv_first_matrix,
                                     lambda_vec)

ggsave("plots/pseudotime/lambda_search.svg", gcv_first_lambdas,
       width = 9, height = 5, device = "svg")


# Repeat the above with only optimal hyperparameters ======

optimal_params <- gcv_for_params[
  which.min(gcv_for_params[, "gcv_sum"]), 
]

optimal_nbasis <- as.numeric(optimal_params["nbasis"])
optimal_lambda <- as.numeric(optimal_params["lambda"])

optimal_outputs <- spline_from_params(optimal_nbasis, optimal_lambda)

# Filter genes by GCV values ==========

# INCREASE GCV THRESHOLD FOR NOW
threshold_to_keep <- 4000

svg("plots/pseudotime/GCV_parameter_order.svg",
    width = 9, height = 5.5, pointsize = 18)
plot(optimal_outputs$gcv[order(optimal_outputs$gcv)],
     xlab = "Position in order of GCV values",
     ylab = "GCV value",
     type = "l",
     lwd = 2)
abline(v = threshold_to_keep)
dev.off()

# keep all genes below GCV threshold
ordered_GCV_genes <- names(optimal_outputs$gcv)[order(optimal_outputs$gcv)]

keep_after_GCV <- ordered_GCV_genes[1:threshold_to_keep]
discard_after_GCV <- ordered_GCV_genes[
  (threshold_to_keep + 1):length(ordered_GCV_genes)
]

gcv_threshold <- optimal_outputs$gcv[tail(keep_after_GCV, 1)]
print(paste0("Keeping all genes with a GCV no more than ", 
             sprintf(gcv_threshold, fmt = "%#.4f")))

# Comparison of this to differentially expressed genes ======

pseudotime_DE_filtered <- volcano_data[keep_after_GCV, ]

# useful to see that there is a good balance between down- and up-regulated genes
pseudotime_DE_counts <- table(pseudotime_DE_filtered$diffexpressed)

cat("Counts for differentially expressed genes which passed GCV filtering\n")
print(pseudotime_DE_counts)

# Visualise pseudotime fits compared to real gene expression ============

plot_fitted_values <- function(sample_id,
                               plot_title = NA,
                               verbose = FALSE,
                               scale_to_01 = FALSE) {
  
  if (verbose) {  # useful if selecting random genes
    print(sample_id)
  }
  
  if (is.na(plot_title)) {
    plot_title <- paste0(sample_id, " over pseudotime")
  }
  
  if (scale_to_01) {
    expression_values <- TPM_for_smoothing[, sample_id]
    
    plot(1:nrow(TPM_for_smoothing), expression_values,
         xlab = "Order of sample over pseudotime", ylab = "Normalised expression")
    lines(optimal_outputs$fd[sample_id], col = "purple", lwd = 2)
    title(main = plot_title,
          sub = paste0("GCV = ",
                       sprintf(optimal_outputs$gcv[sample_id], fmt = "%#.5f")))
  } else {
    expression_values <- as.numeric(TPM_filtered[sample_id, rownames(TPM_for_smoothing)])
    
    # calculate CV^2 too ...
    cv2 <- var(expression_values)/(mean(expression_values)**2)
    if (verbose) {
      print(cv2)
    }
    
    # rescale by multiplying by range, then adding the minimum
    rng = max(expression_values) - min(expression_values)
    
    rescaled_single <- (optimal_outputs$fd[sample_id] * rng) + min(expression_values)
    
    plot(1:ncol(TPM_filtered), expression_values,
         xlab = "Order of sample over pseudotime", ylab = "Expression (TPM)",
         ylim = c(0, max(expression_values)))
    lines(rescaled_single, col = "purple", lwd = 2)
    title(main = plot_title,
          sub = paste0("GCV = ",
                       sprintf(optimal_outputs$gcv[sample_id], fmt = "%#.5f"),
                       "  ||  CV^2 = ",
                       sprintf(cv2, fmt = "%.3f")))
  }

}

# Filter out genes which are mostly 0 but with a few jumps =======

# What happens if I filter genes out based on ratios between the total range 
# and the inner quartile?

ratios_vec <- sapply(keep_after_GCV, function(i){
  vec_in <- TPM_for_smoothing[, i]
  return((max(vec_in) - min(vec_in))/stats::IQR(vec_in))
})

# set infinite values to slightly higher than non-Infinite max
ratios_vec[is.infinite(ratios_vec)] <- max(
  ratios_vec[!is.infinite(ratios_vec)]
) + 10

# plot the quantiles to find an appropriate cut off point 
evaluated_probs <- seq(0, 1, by = 0.0001)
quantiles_ratios <- quantile(ratios_vec, probs = evaluated_probs)

IQR_threshold <- 6.75

svg(filename = "plots/pseudotime/IQR_filtering.svg",
    width = 9, height = 5.5, pointsize = 18)
plot(quantiles_ratios, 1 - evaluated_probs, log = "x",
     xlab = "Cutoff point", ylab = "Proportion of samples removed",
     type = "l")
abline(v = IQR_threshold, col = "red", lty = "dashed", lwd = 2)
abline(v = ratios_vec["AT5G04380"], col = "red")
dev.off()

# sample through the bad ones if interested
svg(filename = "plots/pseudotime/IQR_filtering_example.svg",
    width = 9, height = 5.5, pointsize = 18)
plot_fitted_values("AT5G04380")
dev.off()

# produce multiple examples of genes with bad GCV scores ========

indexes_of_bad_genes <- c(10000, 15000, 19250)

svg("plots/pseudotime/GCV_bad_egs.svg",
    width = 9, height = 5.5, pointsize = 10)
old_pars <- par(mfrow = c(2,2))

plot(optimal_outputs$gcv[order(optimal_outputs$gcv)],
     xlab = "Position in order of GCV values",
     ylab = "GCV value",
     type = "l",
     lwd = 2)
abline(v = threshold_to_keep)
abline(v = indexes_of_bad_genes[1], col = "red")
abline(v = indexes_of_bad_genes[2], col = "red")
abline(v = indexes_of_bad_genes[3], col = "red")

for (bad_index in indexes_of_bad_genes) {
  bad_gene <- ordered_GCV_genes[bad_index]
  
  plot_fitted_values(
    bad_gene, 
    plot_title = bad_gene,
    scale_to_01 = TRUE
  )
}
dev.off()

par(old_pars)

# find genes within these ranges and remove them =====

remove_after_IQR <- names(which(ratios_vec > IQR_threshold))
print(paste0("Removing ", length(remove_after_IQR), " genes after IQR filtering"))

keep_after_IQR <- setdiff(keep_after_GCV, remove_after_IQR)
num_after_IQR <- length(keep_after_IQR)
print(paste0("Leaves ", num_after_IQR, " genes after GCV and IQR filtering"))
