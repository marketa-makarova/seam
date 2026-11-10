#' Plot SPCA sparsity and explained variance across different lambda values
#'
#' A connected scatter plot show the percentage of explained variance and 
#' variables assigned zero weights. 
#'
#' @param expression_data List of numeric expression matrices (samples in rows).
#' @param module_assignment Vector of modules assigned to variables. One module per variable.
#' @param modules Module labels to plot, or `"all"` (default).
#' @param lambda_grid Numeric vector of non-negative sparsity penalties.
#' @param cex Base text size.
#' @param labels Group labels, if present, in the same order as the expression data.
#' @param legend_title Legend title.
#' @param legend_position ggplot legend position.
#' @param linewidth Line width.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param separate_modules Logical. Whether to facet by module as well as group.
#' 
#' @return A named list of ggplot objects is returned invisibly. Plots are printed separately.
#' 
#' @examples
#' \donttest{
#' data("seam_demo", package = "seam")
#' plot_lambdas(
#'  seam_demo$expression_data,
#'  seam_demo$module_assignment,
#'  lambda_grid = c(0.02, 0.05, 0.1),
#'  labels = seam_demo$labels
#' )
#' }
#' 
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
    ylab = "Percentage of dead variables",
    separate_modules = FALSE){

  if(!is.list(expression_data) || length(expression_data) == 0){
    stop("expression_data must be a non-empty list.")
  }
  if(length(labels) != length(expression_data) || anyDuplicated(labels) || any(grepl("-", labels))){
    stop("Provide unique labels, one per group, without '-'.")
  }
  if(!is.numeric(lambda_grid) || length(lambda_grid) == 0 || anyNA(lambda_grid) || any(!is.finite(lambda_grid)) || any(lambda_grid < 0)){
    stop("lambda_grid must contain finite non-negative numbers.")
  }
  no_modules = length(module_assignment) == 1 && is.na(module_assignment)
  module_colors = if(no_modules){
    "all"
  } else if(length(modules) == 1 && modules == "all"){
    unique(as.character(module_assignment))
  } else {
    as.character(modules)
  }

  # prepare plot data
  results = vector("list", length(lambda_grid) * length(expression_data) * length(module_colors))
  k = 0

  for(j in seq_along(lambda_grid)){
    fit = spca_all_modules(
      expression_data = expression_data,
      module_assignment = module_assignment,
      modules = module_colors,
      lambda = lambda_grid[j],
      bootstrap = FALSE
    )

    for(set in seq_along(expression_data)){
      for(i in seq_along(fit$module_colors)){
        k = k + 1
        weights = fit$vs[[set]][[i]]
        sparsity = mean(weights == 0) * 100
        explained = fit$var_exp[[set]][i, 1] * 100

        row = data.frame(
          lambda = lambda_grid[j],
          group = labels[set],
          module = fit$module_colors[i],
          module_size = fit$module_sizes[i],
          percentage_null = sparsity,
          exp_var = explained
        )

        results[[k]] = row
      }
    }
  }

  spca_data = do.call(rbind, results)
  module_labels = paste0(module_colors, " (", fit$module_sizes, ")")
  spca_data$module = factor(paste0(spca_data$module, " (", spca_data$module_size, ")"), levels = module_labels)
  spca_data$group = factor(spca_data$group, levels = labels)

  if(no_modules){
    color_values = setNames("black", module_labels)
  } else {
    # use module names as colors when they are valid R colors, otherwise assign a palette
    valid_colors = vapply(module_colors, function(z){
      tryCatch({ grDevices::col2rgb(z); TRUE }, error = function(e) FALSE)
    }, logical(1))
    color_values = if(all(valid_colors)) module_colors else
      grDevices::hcl.colors(length(module_colors), "Dark 3")
    names(color_values) = module_labels
  }

  base_plot = function(y, lower, upper, y_label){
    p = ggplot2::ggplot(spca_data, ggplot2::aes(x = .data$lambda, y = .data[[y]], group = .data$module, color = .data$module)) +
      ggplot2::geom_line(linewidth = linewidth) +
      ggplot2::scale_color_manual(values = color_values) +
      ggplot2::scale_x_continuous(breaks = unique(lambda_grid)) +
      ggplot2::theme_minimal(base_size = cex) +
      ggplot2::labs(x = xlab, y = y_label, color = legend_title) +
      ggplot2::theme(legend.position = if(no_modules) "none" else legend_position)

    if(separate_modules){
      p = p + ggplot2::facet_grid(rows = ggplot2::vars(.data$group), cols = ggplot2::vars(.data$module)
      )
    } else {
      p = p + ggplot2::facet_grid(rows = ggplot2::vars(.data$group)
      )
    }
    p
  }

  a = base_plot("exp_var", "exp_var_lower", "exp_var_upper", "Explained variance [in %]")
  b = base_plot("percentage_null", "percentage_null_lower", "percentage_null_upper", ylab)

  print(a)
  print(b)
  invisible(list(explained_variance = a, sparsity = b, data = spca_data))
}
