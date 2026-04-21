#' Plot module–trait correlations using sparse and classic eigengene
#'
#' This function visualizes Pearson's module–trait correlation matrices 
#' using both sparse and classic eigengenes
#' across experimental groups. 
#' It also optionally plots root-mean-square errors (RMSE) 
#' between methods and group correlations.
#'
#' @param module_trait_correlations A list of module trait correlations
#' computed using the \code{compute_module_trait_correlations} function.
#' @param labels Character vector. Provides group labels. 
#' @param modules Character or character vector. Allows for module selection. Defaults to \code{"all"}.
#' These labels are used for activity assignment.
#' @param mfrow Layout of plots, default is \code{c(1, 2)}.
#' @param main Optional title for plots. Default is \code{NA}.
#' @param plot_avg_differences Logical. Indicates whether to plot RMSE between sparse eigengene-based and classic module–trait correlations. Defaults to \code{TRUE}.
#' @param mode Character string, either \code{"spca"} or \code{"wgcna"}, indicating which method to use for eigengene computation. Default is \code{"spca"}.
#' @param color_palette Color palette to use for heatmaps. Default is \code{blueWhiteRed(50)}.
#' @param cex Numeric. Scaling factor for the heatmap text labels. Defaults to 1.2.
#' @param show_legend Logical. Default is TRUE.
#'
#' @details
#' This function produces several heatmaps:
#' \itemize{
#'   \item Heatmaps for each method (sparse and classic) and each group.
#'   \item RMSE heatmap between sparse eigengene-based and classic correlations (optional).
#'   \item Consensus correlation heatmap for each module across groups (optional).
#' }
#'
#' Consensus correlations reflect agreement in sign and strength between groups:
#' the weaker (more conservative) correlation is retained if both signs agree.
#'
#' @return Produces multiple heatmaps using \code{labeledHeatmap()}.
#' @export

plot_module_trait_correlations = function(
    module_trait_correlations, 
    labels,
    modules = "all", 
    mfrow = c(1, 2),
    main = NA, 
    plot_avg_differences = FALSE, 
    mode = "spca",
    color_palette = blueWhiteRed(50),
    cex = 1.2,
    show_legend = TRUE){
  
  if(!(mode %in% c("classic", "spca", "both"))){
    stop("Type can be either classic, spca, or both.")
  }
  
  if(plot_avg_differences == TRUE & mode != "both"){
    stop("Cannot plot average differences unless type is set to both.")
  }
  
  if(length(mfrow) != 2){
    stop("Argument mfrow should be a vector of length 2.")
  }
  
  nsets = length(module_trait_correlations[[1]])
  
  if(length(labels) != nsets){
    stop("Length of labels does not coincide with number of groups in module_trait_correlations.")
  }
  
  if(((length(modules) == 1) & (modules[1] != "all")) | ncol(module_trait_correlations$mtc_spca[[1]]) < 2){
    stop("Select at least two modules and traits.")
  }
  
  module_colors = rownames(module_trait_correlations$mtc_spca[[1]])
  
  # error, must have at least 2 colors
  
  if(sum(modules == "all") == 1){
    modules = 1:length(module_colors)
  } else {
    modules = which(module_colors %in% modules) 
  }
  
  module_sizes = module_trait_correlations$module_sizes[modules]
  reduced_module_sizes = module_trait_correlations$reduced_module_sizes
  sizes = list(list(module_sizes, module_sizes), list(reduced_module_sizes, reduced_module_sizes))
  
  if(length(module_trait_correlations) == 5){ # type == both
    types = c("WGCNA", "SPCA")
    for(type in 1:2){
      for(set in 1:nsets){
        labeledHeatmap(Matrix = module_trait_correlations[[type]][[set]][modules, ],
                       xLabels = colnames(module_trait_correlations[[type]][[set]][modules, ]),
                       yLabels = paste(
                         module_colors[modules],
                         " (",
                         if(type == 1){module_sizes} else {reduced_module_sizes[[set]][modules]},
                         ")",
                         sep = ""),
                       colors = color_palette,
                       textMatrix = round(module_trait_correlations[[type]][[set]][modules, ], 2),
                       setStdMargins = FALSE,
                       plotLegend = show_legend,
                       cex.text = cex,
                       main = paste("Module trait correlations by", types[type], "for", labels[set]))
      }
    }
    
    if(plot_avg_differences == TRUE){ # change to cycles, prepare matrix of rmses to store 
      diff_rmse = matrix(NA, nrow = length(modules), ncol = nsets)
      for(set in 1:nsets){
        diff_rmse[, set] = apply(module_trait_correlations[[1]][[set]][modules, ] - module_trait_correlations[[2]][[set]][modules, ], 1, function(x){ sqrt(mean(x^2)) })
      }
      
      labeledHeatmap(Matrix = diff_rmse,
                     xLabels = labels,
                     yLabels = paste(module_colors[modules], " (", module_sizes, ")", sep = ""),
                     ySymbols = module_colors[modules],
                     colors = color_palette,
                     textMatrix = round(diff_rmse, 2),
                     cex.text = cex,
                     plotLegend = show_legend,
                     zlim = c(-1, 1),
                     main = "RMSE of correlation in modules between COR by WGCNA and SPCA") 
    }
  }
}