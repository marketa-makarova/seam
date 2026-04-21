#' Plot explained variance in modules as a heatmap
#'
#' This function visualizes the proportion of variance explained 
#' by the sparse eigengene for each module in both groups, using a heatmap.
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
#' This function builds a heatmap of explained variance 
#' (from 0 to 100, converted to percentage display) 
#' for each module and group, if present. 
#' The heatmap uses the \code{labeledHeatmap()} function from WGCNA.
#'
#' @return A heatmap showing module-wise explained variance per groups.
#'
#' @export
plot_explained_variance = function(
    spca_data, 
    labels, 
    modules = "all", 
    y_axis_labels = NA,
    cex = 1.2,
    main = "Approximate explained variance [in %]",
    color_palette = blueWhiteRed(50),
    show_legend = TRUE,
    show_size = TRUE){
  
  if(!is.list(spca_data)){
    stop("Invalid spca_data argument, must be a list.")
  }
  
  if(length(labels) != length(spca_data$PCs)){
    stop("Length of labels does not coincide with number of groups in spca_data.")
  }
  
  nsets = length(spca_data$PCs)
  module_sizes = lapply(spca_data$vs[[1]], length)
  if(sum(modules == "all") == 1){
    modules = 1:length(module_sizes)
  } else {
    modules = which(spca_data$module_colors %in% modules) 
  }
  
  module_colors = spca_data$module_colors[modules]
  explained_variance_matrix = matrix(NA, ncol = nsets, nrow = length(module_colors))
  
  for(set in 1:nsets){
    explained_variance_matrix[, set] = spca_data$var_exp[[set]][modules] * 100
  }

  y_axis_labels = paste0(
    ifelse(is.na(module_colors), yes = "", no = module_colors),
    if(show_size) paste0(" (", module_sizes[modules], ")") else ""
  )
  
  labeledHeatmap(Matrix = explained_variance_matrix,
                 xLabels = labels,
                 yLabels = y_axis_labels,
                 ySymbols = module_colors,
                 colors = color_palette,
                 plotLegend = show_legend,
                 textMatrix = round(explained_variance_matrix, 2),
                 cex.text = cex,
                 zlim = c(0, 100),
                 main = main)
}