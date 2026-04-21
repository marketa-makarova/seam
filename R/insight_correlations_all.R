#' Combined visualization of SPCA–trait correlation insights
#'
#' This wrapper function sequentially prints all diagnostic plots for SPCA-based module–trait correlations:
#' lollipop plots, reconstruction plots, and sign-change density plots. It applies a unified interest_threshold and aesthetic 
#' settings across all visualizations, allowing comprehensive interpretation of SPCA module behavior across conditions.
#'
#' @param module_trait_correlations A list of module trait correlations
#' computed using the \code{compute_module_trait_correlations} function.
#' @param spca_data A list containing SPCA decomposition results from the 
#' function \code{spca_all_modules} or \code{spca_all_modules_grid}.
#' @param labels Character vector. Provides group labels. 
#' These labels are used for activity assignment.
#' @param col Color or vector of colors if groups are present. 
#' Used for plotting points.
#' @param cex Numeric. Font size for plot elements. 
#' Passed to \code{theme_classic()} argument. 
#' Defaults to 20.
#' @param top_n_upper Integer. 
#' Maximum number of variables to include in the reconstruction. 
#' Defaults to 50.
#' @param top_n_lower Integer. 
#' Minimum number ofvariables to include. 
#' Defaults to 10.
#' @param interest_threshold Numeric. Only module–trait correlations greater
#' than this value (in absolute terms) are considered interesting and plotted. 
#' Default is 0.6.
#' @param return_information Logical. If \code{FALSE}, returns only raw activity 
#' assignments. If \code{TRUE} performs additional correlation-based checks to 
#' detect **undead** variables and updates the activity status accordingly. 
#' Defaults to \code{TRUE}.
#' @param threshold Numeric. The correlation threshold used to determine 
#' whether a **dead** variable highly correlated with an active variable, 
#' suggesting it should be assigned the status **undead**. 
#' Default value is 0.8.
#' @param col_connective Color of a horizontal line at 0. 
#' Defaults to 8 (grey).
#' @param col_active Color for variables that are active in both groups. 
#' Default is 3 (green).
#' @param y_dodge Integer. Used to spread out y-axis labels to avoid label overlap. 
#' Defaults to 2.
#' @param n Integer. The number of top variables 
#' (by absolute weighted correlation) to label on the y-axis. 
#' Defaults to 10.
#'
#' @details
#' This function is a convenience wrapper that:
#' \itemize{
#'   \item Calls \code{\link{insight_correlations_lollipop}} to show per-variable trait correlation contributions.
#'   \item Calls \code{\link{insight_correlations_reconstruct}} to visualize how well trait correlations can be reconstructed using top variables.
#'   \item Calls \code{\link{insight_correlations_sign_changes}} to show the distribution and polarity of correlation signs pre and post weighing.
#' }
#'
#' All three visualizations are printed sequentially using the same \code{interest_threshold}, \code{labels}, and coloring scheme, ensuring coherence across analyses.
#'
#' @return One or more ggplot2 plots are printed to the active graphics device.
#' If no correlations exceed the interest_threshold, a message is returned instead.
#'
#' @seealso \code{\link{insight_correlations_lollipop}}, \code{\link{insight_correlations_reconstruct}}, \code{\link{insight_correlations_sign_changes}}
#' 
#' @export
insight_correlations_all = function(
    module_trait_correlations, 
    spca_data, 
    labels,
    col,
    cex = 20, 
    top_n_upper = 50, 
    top_n_lower = 10, 
    interest_threshold = 0.6,
    return_information = TRUE,
    threshold = 0.8,
    col_connective = 8,
    col_active = 3,
    y_dodge = 1,
    n = 10){
  print(
    insight_correlations_lollipop(
      module_trait_correlations = module_trait_correlations, 
      spca_data = spca_data, 
      cex = cex, 
      labels = labels, 
      col = col,
      n = n, 
      interest_threshold = interest_threshold, 
      xlab = "",
      ylab = "Weighted correlations",
      col_connective = col_connective,
      col_active = col_active,
      y_dodge = y_dodge
    )
  )
  print(
    insight_correlations_reconstruct(
      module_trait_correlations = module_trait_correlations, 
      spca_data = spca_data, 
      cex = cex, 
      labels = labels, 
      col = col,
      top_n_upper = top_n_upper, 
      top_n_lower = top_n_lower, 
      interest_threshold = interest_threshold,
      threshold = threshold,
      return_information = return_information, 
      ylab = "Correlation",
      xlab = "Top n variables",
      col_connective = col_connective
    )
  )
  print(
    insight_correlations_sign_changes(
      module_trait_correlations = module_trait_correlations, 
      spca_data = spca_data, 
      col = col,
      interest_threshold = interest_threshold, 
      cex = cex, 
      labels = labels,
      xlab = "Correlations",
      ylab = ""
    )
  )
}