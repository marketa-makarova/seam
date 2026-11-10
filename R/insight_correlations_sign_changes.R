#' Visualize variable-trait correlations and negative SPCA weights
#'
#' Selects module-trait pairs exceeding the absolute SPCA correlation
#' threshold in at least one group. For each pair, plots the density of
#' unweighted variable-trait correlations among variables active in SPCA.
#' The title reports the percentage of active variables with negative
#' SPCA weights in each group.
#'
#' @param module_trait_correlations Result of compute_module_trait_correlations()
#'   with mode = "spca" or "both".
#' @param spca_data Result of spca_all_modules() or spca_all_modules_grid().
#' @param labels Group labels, if present, in the same order as the expression data.
#' @param col Optional colors for the groups.
#' @param interest_threshold Miminal absolute module-trait correlation threshold.
#' @param cex Base text size.
#' @param xlab,ylab Axis labels.
#' 
#' @return A named list of ggplot objects is returned invisibly. Plots are printed separately.
#' 
#' @examples
#' \donttest{
#' data("seam_demo", package = "seam")
#' fit = spca_all_modules(
#'  seam_demo$expression_data,
#'  seam_demo$module_assignment,
#'  lambda = 0.05,
#'  bootstrap = FALSE
#' )
#' mtc = compute_module_trait_correlations(
#'  seam_demo$trait_data, 
#'  fit, 
#'  seam_demo$module_assignment,
#'  mode = "both")
#' insight_correlations_sign_changes(
#'  mtc, 
#'  fit,
#'  labels = seam_demo$labels,
#'  interest_threshold = 0.2
#' )
#' }
#' 
#' @export
insight_correlations_sign_changes = function(
    module_trait_correlations,
    spca_data,
    labels,
    col = NULL,
    interest_threshold = 0.6,
    cex = 20,
    xlab = "Correlations",
    ylab = ""){

  if(!is.list(spca_data) || is.null(module_trait_correlations$mtc_spca)){
    stop("spca_data and SPCA module-trait correlations are required.")
  }
  if(!is.numeric(interest_threshold) || length(interest_threshold) != 1 ||
     !is.finite(interest_threshold) || interest_threshold <= 0 || interest_threshold >= 1){
    stop("interest_threshold must be between 0 and 1.")
  }
  nsets = length(spca_data$vs)
  if(length(labels) != nsets || length(module_trait_correlations$mtc_spca) != nsets ||
     length(module_trait_correlations$trait_data) != nsets){
    stop("Group labels and correlation results must match spca_data.")
  }

  # prepare necessary variables
  reference = module_trait_correlations$mtc_spca[[1]]
  module_colors = spca_data$module_colors
  if(is.null(colnames(reference))) stop("Trait columns must have names.")
  interesting = matrix(FALSE, nrow(reference), ncol(reference))
  for(set in seq_len(nsets)){
    current = module_trait_correlations$mtc_spca[[set]]
    interesting = interesting | (is.finite(current) & abs(current) > interest_threshold)
  }
  positions = which(interesting, arr.ind = TRUE)
  if(nrow(positions) == 0){
    message("No interesting correlations; lower interest_threshold.")
    return(invisible(list()))
  }
  colors = if(is.null(col)) grDevices::hcl.colors(nsets, "Dark 3") else rep(col, length.out = nsets)
  names(colors) = labels
  plots = list()

  # prepare plot data and plot
  for(i in seq_len(nrow(positions))){
    module_index = positions[i, 1]
    trait_index = positions[i, 2]
    module_name = rownames(reference)[module_index]
    trait_name = colnames(reference)[trait_index]
    if(is.null(module_name)) module_name = spca_data$module_colors[module_index]
    if(is.null(trait_name)) trait_name = as.character(trait_index)

    densities = list()
    percentages = rep(NA_real_, nsets)

    for(set in seq_len(nsets)){
      x = spca_data$module_matrices[[set]][[module_index]]
      weights = spca_data$vs[[set]][[module_index]]
      trait = module_trait_correlations$trait_data[[set]][, trait_index]
      if(length(weights) != ncol(x) || length(trait) != nrow(x)) stop("Weights, expression matrix, and trait data dimensions do not match.")

      # filter out dead variables
      active = is.finite(weights) & weights != 0
      if(any(active)){
        percentages[set] = mean(weights[active] < 0) * 100
        cors = suppressWarnings(stats::cor(x[, active, drop = FALSE], trait, use = "pairwise.complete.obs"))
        values = as.vector(cors)
        values = values[is.finite(values)]
        if(length(values)){
          densities[[length(densities) + 1]] = data.frame(cor = values, label = labels[set], stringsAsFactors = FALSE)
        }
      }
    }

    # prepare heading
    if(!length(densities)) next
    plot_data = do.call(rbind, densities)
    plot_data$label = factor(plot_data$label, levels = labels)
    percentage_text = paste0(
      labels, ": ",
      ifelse(is.na(percentages), "NA", sprintf("%.0f", percentages)),
      "% negative weights"
    )
    heading = if(length(spca_data$module_colors) == 1 && spca_data$module_colors[1] == "all"){
      paste0("Correlation with ", trait_name)
    } else {
      paste0(module_name, " module correlated with ", trait_name)
    }
    heading = paste0(heading, "\n", paste(percentage_text, collapse = "\n"))

    # plot
    p = ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$cor, fill = .data$label)) +
      ggplot2::geom_density(alpha = 0.4, na.rm = TRUE) +
      ggplot2::scale_fill_manual(values = colors, name = "Group", drop = FALSE) +
      ggplot2::theme_classic(base_size = cex) +
      ggplot2::labs(x = xlab, y = ylab, title = heading)

    key = paste(module_name, trait_name, sep = "__")
    plots[[key]] = p
    print(p)
  }

  invisible(plots)
}
