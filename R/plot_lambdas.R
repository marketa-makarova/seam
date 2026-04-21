#' Plot SPCA sparsity and explained variance for different values of lambda
#'
#' Visualizes the effect of different \eqn{\lambda} values (sparsity penalty) 
#' on sparsity across modules and groups, if present.
#'
#' @param expression_data A list containing expression matrices for each group.
#' @param module_assignment Character vector. Assigns
#'  module/cluster to each variable by color. Cannot
#'  be set to \code{NA}.
#' @param modules Character or character vector. Allows for module selection. Defaults to \code{"all"}.
#' @param lambda_grid A numeric vector of lambda values to iterate over.
#' @param cex Numeric. Font size for plot elements. 
#' Passed to \code{theme_classic()} argument. 
#' Defaults to 20.
#' @param labels Character vector. Provides group labels. 
#' @param legend_title Character. Optional title for a legend. 
#' @param legend_position Character. Position of a legend to be passed to \code{ggplot} code. Default is \code{"right"} for a case with modules and \code{"none"} for a case with no modules.
#' @param linewidth Numeric. Width of plotted lines. Default is 2.
#' @param xlab Character. X-axis label.
#' @param ylab Character. Y-axis label.
#'
#' @details
#' This function builds a heatmap of the percentages of variables
#' which were assigned zero weights by SPCA 
#' for each module and group, if present. 
#' The heatmap uses the \code{labeledHeatmap()} function from WGCNA.
#'
#' @return Two printed ggplot objects. No values are returned.
#' @export

plot_lambdas = function(
    expression_data, 
    module_assignment, 
    modules = "all", 
    lambda_grid, 
    cex = 17,
    labels,
    legend_title = "Module",
    legend_position = "right",
    linewidth = 2,
    xlab = "Lambda", 
    ylab = "Percentage of dead variables"){
  
  if(!is.list(expression_data)){
    stop("Invalid expression_data argument, must be a list.")
  }
  
  if(any(lambda_grid < 0)) { stop("Select non-negative lambdas.") }
  
  nsets = length(expression_data)

  for(i in 1:nsets){
    if(length(module_assignment) != ncol(expression_data[[i]])){
      stop("Number of variables in expression_data does not coincide with length of module_assignment.")
    }
  }
  
  if(length(labels) != length(expression_data)){
    stop("Length of labels does not coincide with number of groups in expression_data.")
  }

  
  exp_var = numeric(0) # percentage of explained variability
  percentage_null = numeric(0) # percentage of variables with weight pushed to zero

  
  if(sum(modules == "all") == 1){
    module_colors = unique(module_assignment)
  } else {
    module_colors = modules
  }
  
  for(i in 1:length(lambda_grid)){
    module_spca = spca_all_modules(
      expression_data = expression_data, 
      module_assignment = module_assignment, 
      modules = module_colors, 
      lambda = lambda_grid[i]
    )
    for(set in 1:nsets){
      percentage_null = c(percentage_null, lapply(module_spca$vs[[set]], zeroPercent) |> as.numeric())
      exp_var = c(exp_var, module_spca$var_exp[[set]][, 1])
    }
  }
  
  module_sizes = module_spca$module_sizes
  
  spca_data = data.frame(lambda = rep(rep(1:length(lambda_grid), each = nsets), each = length(module_colors)),
                         percentage_null = percentage_null * 100,
                         exp_var = exp_var * 100,
                         module = rep(rep(paste(module_colors, ' (', module_sizes, ')', sep = ""), nsets), length(lambda_grid)),
                         group = rep(rep(labels, each = length(module_colors)), length(lambda_grid))) # fix labels
  
  color_values = if(any(is.na(module_colors))){
    rep("black", length(unique(spca_data$module)))  
  } else {
    sort(module_colors)
  }
  
  a = ggplot(spca_data, aes(x = .data$lambda, y = .data$exp_var, group = .data$module, color = .data$module)) +
    geom_line(linewidth = linewidth) +
    scale_color_manual(values = color_values) +
    scale_x_continuous(breaks = 1:length(lambda_grid), labels = lambda_grid) + 
    facet_grid(~ .data$group) + 
    theme_minimal(base_size = cex) + 
    labs(x = "Lambda", y = "Explained variance [in %]", color = legend_title) +
    theme(legend.position = ifelse(sum(is.na(module_colors)) != 0, yes = "none", no = legend_position))
  
  print(a)
  
  b = ggplot(spca_data, aes(x = .data$lambda, y = .data$percentage_null, group = .data$module, color = .data$module)) +
    geom_line(linewidth = linewidth) +
    scale_color_manual(values = color_values) +
    scale_x_continuous(breaks = 1:length(lambda_grid), labels = lambda_grid) +
    facet_grid(~ .data$group) +
    theme_minimal(base_size = cex) +
    labs(x = xlab, y = ylab, color = legend_title) 
  
  print(b)
  
}