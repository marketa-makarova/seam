#' Plot percentage of inactive variables (sparsity) across modules
#'
#' Visualizes the proportion of zero weights ("dead variables") 
#' in each module for groups, if present. 
#' The resulting heatmap helps assess sparsity induced by SPCA.
#'
#' @param spca_data A list containing SPCA decomposition results from the 
#' function \code{spca_all_modules} or \code{spca_all_modules_grid}.
#' @param labels Character vector. Provides group labels. 
#' These labels are used for activity assignment.
#' @param modules Character or character vector. Allows for module selection. Defaults to \code{"all"}.
#' @param cex Numeric. Scaling factor for the heatmap text labels. Defaults to 1.2.
#' @param main Character. Title of the heatmap.
#' @param color_palette Color palette for the heatmap. Default is \code{blueWhiteRed(50)}.
#' @param show_legend Logical. Defaults to TRUE.
#' @param show_size Logical. Defaults to TRUE.
#' @param y_axis_labels A character vector for optional y-axis labels. Defaults to NULL.
#'
#' @details
#' For each module and each group, 
#' the percentage of zero entries in the module’s SPCA loading
#' vector is calculated.
#' These percentages indicate how sparse the component is. 
#' This is useful for understanding how many variables
#' are “inactive” in each module, per group.
#'
#' The function requires the existence of a helper function \code{zeroPercent()}
#' that computes the proportion of zero entries in a vector or matrix.
#'
#' @return A heatmap produced by \code{labeledHeatmap()} showing the percentage of zero weights per module and group.
#' @importFrom WGCNA labeledHeatmap
#' @export

plot_sparsity = function(
    spca_data,
    modules = "all",
    labels = NA,
    y_axis_labels = NA,
    cex = 1.2,
    main = "Dead variables [in %]",
    color_palette = blueWhiteRed(50),
    show_legend = TRUE,
    show_size = TRUE){
  
  if(!is.list(spca_data)){
    stop("Invalid spca_data argument, must be a list.")
  }
  
  nsets = length(spca_data$PCs)
  
  if(length(labels) != nsets){
    stop("Length of labels does not coincide with number of groups in spca_data.")
  }
  
  if(sum(modules == "all") == 1){
    modules = 1:length(spca_data$module_colors)
    module_colors = spca_data$module_colors
  } else {
    module_colors = modules
    modules = which(spca_data$module_colors %in% modules) 
  }
  
  percentage_of_zeroes_matrix = matrix(NA, ncol = nsets, nrow = length(modules))
  module_sizes = spca_data$module_sizes
  
  for(set in 1:nsets){
    for(i in modules){
      percentage_of_zeroes_matrix[which(modules == i), set] = spca_data$vs[[set]][[i]] |> zeroPercent()
    }
  }

  
  y_axis_labels = paste0(
    ifelse(is.na(module_colors), yes = "", no = module_colors),
    if(show_size) paste0(" (", module_sizes, ")") else ""
  )[modules]
  
  labeledHeatmap(Matrix = percentage_of_zeroes_matrix,
                 xLabels = labels,
                 yLabels = y_axis_labels,
                 ySymbols = module_colors,
                 colors = color_palette,
                 plotLegend = show_legend,
                 textMatrix = percentage_of_zeroes_matrix,
                 cex.text = cex,
                 zlim = c(0, 100),
                 main = main)
}