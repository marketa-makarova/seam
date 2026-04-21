#' Plot full module heatmaps, highlights dimension reduction and plots SPCA reduced module heatmaps
#'
#' This function creates correlation-based heatmaps of expression values.
#' If modules and groups are present, this is done for each module and group 
#' separately.
#' For each instnace, three heatmaps are shown:
#' full, sparse (showing which variables were assigned zero weights by SPCA),¨
#' and reduced (with the aforementioned variables removed).
#'
#' @param spca_data A list containing SPCA decomposition results from the 
#' function \code{spca_all_modules} or \code{spca_all_modules_grid}.
#' @param labels Character vector. Provides group labels. 
#' computed using the \code{compute_module_trait_correlations} function.
#' @param modules Character or character vector. Allows for module selection. Defaults to \code{"all"}.
#' @param cex Numeric. Font size for plot elements. 
#' Passed to \code{theme_classic()} argument. 
#' Defaults to 20.
#' @param return_information Logical. If \code{FALSE}, returns only raw activity 
#' assignments. If \code{TRUE} performs additional correlation-based checks to 
#' detect **undead** variables and updates the activity status accordingly. 
#' Defaults to \code{TRUE}.
#' @param threshold Numeric. The correlation threshold used to determine 
#' whether a **dead** variable highly correlated with an active variable, 
#' suggesting it should be assigned the status **undead**. 
#' Default value is 0.8.
#' @param midpoint Numeric. Midpoint for the heatmap color gradient. Default is 0.5.
#' @param col_low Color for low values in the heatmap. Default is \code{"darkred"}.
#' @param col_mid Color for mid values (at midpoint) in the heatmap. Default is \code{"orange"}.
#' @param col_high Color for high values in the heatmap. Default is \code{"white"}.
#'
#' @details
#' For each module and each group, 
#' the function computes a distance matrix (1 - |correlation|) and displays three versions:
#' \itemize{
#'   \item \strong{Full:} all variable correlations sorted using a dendrogram.
#'   \item \strong{Sparse:} with variable assigned zero weighs entries shown as grey.
#'   \item \strong{Reduced:} subset of correlations assigned non-zero SPCA weights variables.
#' }
#'
#' Clustering is based on average linkage hierarchical clustering using the full matrix. 
#' 
#' A procedure refining variable selection may be implemented using 
#' the \code{return_information} parameter. For more information, refer to 
#' \code{\link{compute_activity}} function documentation.
#'
#' @return Plots heatmaps (one for each module) using \code{ggplot2}.
#'
#' @export

plot_module_heatmaps = function(
    spca_data,
    labels,
    modules = "all",
    return_information = TRUE,
    threshold = 0.8,
    cex = 20,
    midpoint = 0.5, 
    col_low = "darkred", 
    col_mid = "orange", 
    col_high = "white"){
  
  if(!is.list(spca_data)){
    stop("Invalid spca_data argument, must be a list.")
  }
  
  if(length(labels) != length(spca_data$PCs)){
    stop("Length of labels does not coincide with number of groups in spca_data.")
  }
  
  if(!is.numeric(threshold)){
    stop("Invalid threshold argument value.")
  }
  
  if(threshold <= 0 | threshold >= 1){
    warning("Threshold argument value should be between 0 and 1.")
  }
  
  if(midpoint <= -1 | midpoint >= 1){
    warning("Midpoint argument value should be between -1 and 1.")
  }
  
  if(sum(modules == "all") == 1){
    modules = 1:length(spca_data$module_colors)
  } else {
    modules = which(spca_data$module_colors %in% modules) 
  }
  
  n_sets = length(spca_data$module_matrices)
  module_sizes = spca_data$module_sizes
  module_colors = spca_data$module_colors
  module_activity_table = compute_activity(spca_data = spca_data,
                                           labels = labels,
                                           return_information = return_information,
                                           threshold = threshold)
  
  for(set in 1:n_sets){
    # load data
    module_matrices = spca_data$module_matrices[[set]]
    n_modules = length(module_matrices)
    
    # iterate through modules
    for(module in modules){
      
      # create heatmap data
      diss = 1 - abs(cor(as.matrix(module_matrices[[module]]))) 
      correlation_data = data.frame(
        corr = as.vector(diss),
        x = colnames(diss), 
        y = rep(colnames(diss), each = ncol(diss))) 
      
      hc = hclust(as.dist(diss), method = "average") 
      
      ordered_labels = colnames(diss)[hc$order]
      
      correlation_data$x = factor(correlation_data$x, levels = ordered_labels)
      correlation_data$y = factor(correlation_data$y, levels = ordered_labels)
      

      sparse_correlation_data = correlation_data |> 
        left_join(module_activity_table, by = join_by("x" == "variable_id")) |>
        left_join(module_activity_table, by = join_by("y" == "variable_id")) |>
        mutate(
          corr = if_else(.data[[paste0("w_", labels[set], ".x")]] == 0, NA, .data$corr)
        ) |>
        mutate(
          corr = if_else(.data[[paste0("w_", labels[set], ".y")]] == 0, NA, .data$corr)
        )
      
        
      n_returned = sum(grepl("undead", sparse_correlation_data$activity.x)) / spca_data$module_sizes[module]
      sparse_correlation_data = sparse_correlation_data |> select(c("corr", "x", "y"))

      reduced_correlation_data = na.omit(sparse_correlation_data)
      
      if(return_information == FALSE){
        heatmap_data = rbind(correlation_data, sparse_correlation_data, reduced_correlation_data) |> 
          mutate(
            type = factor(
              c(rep("full", nrow(correlation_data)), 
                rep("sparse", nrow(sparse_correlation_data)), 
                rep("reduced", nrow(reduced_correlation_data))), 
              levels = c("full", "sparse", "reduced")
            )
          )
      } 
      
      returned_correlation_data = correlation_data |> 
        left_join(module_activity_table, by = join_by("x" == "variable_id")) |>
        left_join(module_activity_table, by = join_by("y" == "variable_id")) |>
        mutate(corr = ifelse(
          test = .data$activity.x %in% c("undead", paste0("undead_", labels[set]), paste0("active_", labels[set]), "active") &
            .data$activity.y %in% c("undead", paste0("undead_", labels[set]), paste0("active_", labels[set]), "active"),
          yes = .data$corr,
          no = NA)
        ) |>
        select(c("corr", "x", "y"))
      
      returned_reduced_correlation_data = na.omit(returned_correlation_data)
      
      if(return_information == TRUE){
        heatmap_data = rbind(correlation_data, sparse_correlation_data, returned_correlation_data, returned_reduced_correlation_data) |> 
          mutate(
            type = factor(
              c(rep("full", nrow(correlation_data)), 
                rep("sparse", nrow(sparse_correlation_data)), 
                rep("returned", nrow(returned_correlation_data)),
                rep("reduced", nrow(returned_reduced_correlation_data))), 
              levels = c("full", "sparse", "returned", "reduced")
            )
          )
      }
      
      if(sum(is.na(spca_data$module_colors)) == 0){
        title = paste(
          labels[set],
          " ",
          module_colors[module],
          " module (", 
          module_sizes[module],
          ") removing ",
          sum(spca_data$vs[[set]][[module]] == 0),
          " variables,",
          ifelse(
            return_information == TRUE, 
            yes = paste(
              " returning ",
              n_returned,
              " variables, ",
              sep = ""
            ), 
            no = " "),
          "final size ",
          ifelse(
            return_information == TRUE,
            yes = module_sizes[module] - sum(spca_data$vs[[set]][[module]] == 0) + n_returned,
            no = module_sizes[module] - sum(spca_data$vs[[set]][[module]] == 0)),
          sep = ""
        )
      } else {
        title = paste(
          labels[set],
          " (", 
          module_sizes[module],
          ") removing ",
          sum(spca_data$vs[[set]][[module]] == 0),
          " variables, final size ",
          module_sizes[module] - sum(spca_data$vs[[set]][[module]] == 0),
          sep = ""
        )
      }
      
      # plot heatmap
      print(
        ggplot(
          heatmap_data, aes(x = .data$x, y = .data$y, fill = .data$corr)) + 
          geom_tile() + 
          facet_wrap(~ type, scales = "free", ncol = ifelse(return_information == TRUE, yes = 4, no = 3)) +
          theme_classic(base_size = cex) +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank()
          ) +
          scale_fill_gradient2(midpoint = 0.5, low = col_low, mid = col_mid, high = col_high) +
          labs(fill = "Distance") + 
          ggtitle(title)
      )
    } 
  }
}