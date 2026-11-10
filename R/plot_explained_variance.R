#' Plot explained variance in modules as a heatmap
#'
#' Cell colors show percentages of variance explained by SPCA.
#'
#' @param spca_data Result of spca_all_modules() or spca_all_modules_grid().
#' @param labels Group labels, if present, in the same order as the expression data.
#' @param modules Module labels to plot, or `"all"` (default).
#' @param y_axis_labels Optional custom module labels.
#' @param cex Base text size.
#' @param main Plot title.
#' @param color_palette Colors for the continuous heatmap scale.
#' @param show_legend Logical. Whether to display the color legend.
#' @param show_size Logical. Whether to include module sizes in row labels.
#' 
#' @return A ggplot object, printed and returned invisibly.
#' 
#' @examples
#' data("seam_demo", package = "seam")
#' fit = spca_all_modules(
#'  seam_demo$expression_data,
#'  seam_demo$module_assignment,
#'  lambda = 0.05
#' )
#' plot_explained_variance(
#'  fit, 
#'  labels = seam_demo$labels
#' )
#' 
#' @export
plot_explained_variance = function(
  spca_data,
  labels,
  modules = "all",
  y_axis_labels = NA,
  cex = 1.2,
  main = "Approximate explained variance [in %]",
  color_palette = grDevices::colorRampPalette(c("blue", "white", "red"))(50),
  show_legend = TRUE,
  show_size = TRUE){

  if(!is.list(spca_data)) stop("Invalid spca_data.")
  nsets = length(spca_data$PCs)
  if(length(labels) != nsets) stop("Length of labels does not match the number of groups.")
  all_modules = as.character(spca_data$module_colors)
  if(!length(modules) || anyNA(modules) || anyDuplicated(modules)){
    stop("Unknown module in modules.")
  }
  if(length(modules) == 1 && identical(modules, "all")){
    modules = seq_along(all_modules)
  } else {
    if(!all(modules %in% all_modules)) stop("Unknown module in modules.")
    modules = match(modules, all_modules)
  }

  # prepare plot data
  module_colors = spca_data$module_colors[modules]
  module_sizes = spca_data$module_sizes[modules]
  explained_variance_matrix = matrix(NA_real_, nrow = length(modules), ncol = nsets)
  text_matrix = matrix(NA_character_, nrow = length(modules), ncol = nsets)

  for(set in seq_len(nsets)){
    for(i in seq_along(modules)){
      explained_variance = spca_data$var_exp[[set]][modules[i], 1] * 100
      explained_variance_matrix[i, set] = explained_variance
      text_matrix[i, set] = as.character(round(explained_variance, 2))
    }
  }

  if(length(y_axis_labels) == 1 && is.na(y_axis_labels)){
    y_axis_labels = paste0(
      ifelse(is.na(module_colors), "", module_colors),
      if(show_size) paste0(" (", module_sizes, ")") else ""
    )
  }

  if(length(y_axis_labels) != length(modules)){
    stop("Length of y_axis_labels must match the number of selected modules.")
  }

  plot_data = data.frame(
    group = rep(seq_len(nsets), each = length(modules)),
    module = rep(seq_along(modules), times = nsets),
    percentage = as.vector(explained_variance_matrix),
    text = as.vector(text_matrix)
  )
  plot_data$group = factor(plot_data$group, levels = seq_len(nsets), labels = labels)
  plot_data$module = factor(plot_data$module, levels = rev(seq_along(modules)), labels = rev(y_axis_labels))

  p = ggplot(plot_data, aes(x = .data$group, y = .data$module, fill = .data$percentage)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = .data$text), size = cex * 3.5) +
    scale_fill_gradientn(colors = color_palette, limits = c(0, 100), name = "Explained variance [%]") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL, title = main) +
    theme_minimal(base_size = 14) +
    theme(
      aspect.ratio = 0.8,
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = if(show_legend) "right" else "none"
    )

  print(p)
  invisible(p)
}
