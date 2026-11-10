#' Reconstruct module-trait correlations from top contributing variables
#'
#' Ranks variables within each group by the absolute product of their
#' SPCA weight and variable-trait correlation. Reconstructs the
#' SPCA score from the top variables and compares its trait correlation
#' with the original full-module correlation. Optional bootstrap fits do
#' not change the reconstruction; only original estimates are used.
#'
#' @param module_trait_correlations Result of compute_module_trait_correlations()
#'   calculated with mode = "spca" or "both".
#' @param spca_data Result of spca_all_modules() or spca_all_modules_grid().
#' @param labels Group labels, if present, in the same order as the expression data.
#' @param col Optional colors for the groups.
#' @param cex Base text size.
#' @param top_n_upper Maximum number of selected variables.
#' @param top_n_lower Minimum number of selected variables.
#' @param interest_threshold Miminal absolute module-trait correlation threshold.
#' @param return_information `TRUE` for refined undead classifications,
#'   `FALSE` for raw selection, or `"both"` to print both plots.
#' @param threshold Absolute within-group correlation threshold for reclassification to undead status.
#' @param ylab,xlab Axis labels.
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
#'  mode = "both"
#' )
#' insight_correlations_reconstruct(
#'  mtc, 
#'  fit,
#'  labels = seam_demo$labels,
#'  top_n_lower = 1, 
#'  top_n_upper = 8,
#'  interest_threshold = 0.2
#' )
#' }
#' 
#' @export
insight_correlations_reconstruct = function(
    module_trait_correlations,
    spca_data,
    labels,
    col = NULL,
    cex = 20,
    top_n_upper = 50,
    top_n_lower = 10,
    interest_threshold = 0.6,
    threshold = 0.8,
    return_information = TRUE,
    ylab = "Correlation",
    xlab = "Top n variables"){

  if(!is.list(spca_data) || is.null(module_trait_correlations$mtc_spca)){
    stop("spca_data and SPCA module-trait correlations are required.")
  }
  if(!is.numeric(interest_threshold) || length(interest_threshold) != 1 ||
     !is.finite(interest_threshold) || interest_threshold <= 0 || interest_threshold >= 1){
    stop("interest_threshold must be between 0 and 1.")
  }
  if(threshold <= 0 || threshold >= 1 || !is.numeric(threshold) || length(threshold) != 1 || is.na(threshold)){
    stop("threshold argument value must be a number between 0 and 1.")
  }
  if(length(top_n_lower) != 1 || length(top_n_upper) != 1 ||
     !is.finite(top_n_lower) || !is.finite(top_n_upper) ||
     top_n_lower < 1 || top_n_upper < top_n_lower ||
     top_n_lower != as.integer(top_n_lower) || top_n_upper != as.integer(top_n_upper)){
    stop("top_n_lower and top_n_upper must be positive integers in increasing order.")
  }
  nsets = length(spca_data$vs)
  if(length(labels) != nsets || length(module_trait_correlations$mtc_spca) != nsets){
    stop("Group labels and correlation results must match spca_data.")
  }
  if(anyDuplicated(labels)) stop("Group labels must be unique.")

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
  activity = compute_activity(
    spca_data = spca_data,
    labels = labels,
    return_information = return_information,
    threshold = threshold
  )
  colors = if(is.null(col)){
    stats::setNames(grDevices::hcl.colors(nsets, "Dark 3"), labels)
  } else {
    stats::setNames(rep(col, length.out = nsets), labels)
  }

  # prepare plot data
  plots = list()
  for(i in seq_len(nrow(positions))){
    module_index = positions[i, 1]
    trait_index = positions[i, 2]
    module_name = module_colors[module_index]
    trait_name = colnames(reference)[trait_index]
    top_values = seq.int(top_n_lower, top_n_upper)
    reconstruction = matrix(NA_real_, length(top_values), nsets, dimnames = list(NULL, labels))

    for(set in seq_len(nsets)){
      x = spca_data$module_matrices[[set]][[module_index]]
      variables = colnames(x)
      if(is.null(variables)) stop("Module expression matrices need variable names.")
      trait = module_trait_correlations$trait_data[[set]][, trait_index]
      weights = spca_data$vs[[set]][[module_index]]
      if(!is.null(names(weights))) weights = weights[variables]
      if(length(weights) != ncol(x)) stop("Weights do not match module variables.")

      # filter out dead variables
      active = activity$module == module_name & activity$variable_id %in% variables & activity$activity != "dead"
      eligible = activity$variable_id[active]
      variable_cor = as.vector(stats::cor(x, trait, use = "pairwise.complete.obs"))
      contribution = abs(weights * variable_cor)
      contribution[!is.finite(contribution) | !(variables %in% eligible)] = -Inf
      ranking = order(contribution, decreasing = TRUE)
      ranking = ranking[is.finite(contribution[ranking])]

      # reconstruct
      for(k in seq_along(top_values)){
        j = top_values[k]
        selected = head(ranking, j)
        partial_weights = numeric(length(weights))
        partial_weights[selected] = weights[selected]
        score = as.vector(x %*% partial_weights)
        if(stats::sd(score) == 0) next
        # weights have already been oriented in spca_all_modules() - do not flip again
        reconstruction[k, set] = suppressWarnings(stats::cor(score, trait, use = "pairwise.complete.obs"))
      }
    }

    plot_data = data.frame(
      top_n = rep(top_values, times = nsets),
      group = factor(rep(labels, each = length(top_values)), levels = labels),
      correlation = as.vector(reconstruction)
    )

    hline_data = data.frame(
      group = factor(labels, levels = labels),
      correlation = sapply(seq_len(nsets), function(set){ # pull original correlations
        module_trait_correlations$mtc_spca[[set]][module_index, trait_index]
      })
    )

    title = if(module_name == "all"){
      paste0("Correlation with ", trait_name)
    } else {
      paste0(module_name, " module correlated with ", trait_name)
    }
    title = paste0(title, " using ", top_n_lower, " to ", top_n_upper, " variables")

    p = ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$top_n, y = .data$correlation, color = .data$group)) +
      ggplot2::geom_line(linewidth = 1, alpha = 0.65, na.rm = TRUE) +
      ggplot2::geom_hline(data = hline_data, ggplot2::aes(yintercept = .data$correlation, color = .data$group),
        linewidth = 1.5, show.legend = FALSE
      ) +
      ggplot2::geom_hline(yintercept = 0, color = 8) +
      ggplot2::scale_color_manual(values = colors, name = "Group") +
      ggplot2::theme_classic(base_size = cex) +
      ggplot2::labs(title = title, x = xlab, y = ylab)

    key = paste(module_name, trait_name, sep = "__")
    plots[[key]] = p
    print(p)
  }

  invisible(plots)
}
