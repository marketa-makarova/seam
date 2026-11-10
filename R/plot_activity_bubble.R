#' Plot variable activity statuses as a bubble chart
#'
#' Bubble size describes counts of variables assigned a given status.
#' When bootstrap results are available, a second number optionally displays
#' an average status count across bootstrap replicates.
#'
#' @param spca_data Result of spca_all_modules() or spca_all_modules_grid().
#' @param labels Group labels, if present, in the same order as the expression data.
#' @param return_information `TRUE` for refined undead classifications,
#'   `FALSE` for raw selection, or `"both"` to print both plots.
#' @param modules Module labels to plot, or `"all"` (default).
#' @param threshold Absolute within-group correlation threshold for reclassification to undead status.
#' @param cex Base text size.
#' @param scale_range Minimum and maximum bubble sizes.
#' @param number_color Color of numbers inside bubbles.
#' @param bubble_colors Optional colors for modules; NA uses module colors when they are valid.
#' @param show_bootstrap Show average bootstrap counts below original counts.
#' 
#' @return A ggplot object, printed and returned invisibly.
#' 
#' @examples
#' data("seam_demo", package = "seam")
#' fit = spca_all_modules(
#'  seam_demo$expression_data,
#'  seam_demo$module_assignment,
#'  lambda = 0.05, 
#'  bootstrap = FALSE
#' )
#' plot_activity_bubble(
#'  fit, 
#'  labels = seam_demo$labels
#' )
#' 
#' @export
plot_activity_bubble = function(
    spca_data,
    labels,
    modules = "all",
    threshold = 0.8,
    return_information = TRUE,
    cex = 20,
    scale_range = c(5, 30),
    number_color = "white",
    bubble_colors = NA,
    show_bootstrap = TRUE){

  if(!is.list(spca_data) || is.null(spca_data$module_colors) || is.null(spca_data$PCs)) stop("Invalid spca_data.")
  if(length(labels) != length(spca_data$PCs)){
    stop("Length of labels does not match the number of groups.")
  }
  if(length(scale_range) != 2 || any(!is.finite(scale_range))){
    stop("scale_range must contain two finite numbers.")
  }
  if(!is.numeric(threshold) || length(threshold) != 1 || is.na(threshold) || threshold <= 0 || threshold >= 1){
    stop("threshold argument value must be a number between 0 and 1.")
  }
  module_colors = as.character(spca_data$module_colors)
  all_modules = as.character(spca_data$module_colors)
  selected_modules = if(length(modules) == 1 && identical(modules, "all")){
    module_colors
  } else {
    as.character(modules)
  }
  if(!length(selected_modules) || anyNA(selected_modules) || anyDuplicated(selected_modules) || !all(selected_modules %in% all_modules)){
    stop("Unknown module in modules.")
  } 

  # prepare plot data
  activity_table = compute_activity(
    spca_data = spca_data,
    labels = labels,
    threshold = threshold,
    return_information = return_information
  )
  activity_table$module = as.character(activity_table$module)
  activity_table = activity_table[activity_table$module %in% selected_modules, , drop = FALSE]
  pct_cols = grep("^pct_", names(activity_table), value = TRUE)
  statuses = unique(c(as.character(activity_table$activity), sub("^pct_", "", pct_cols)))
  # keep the same activity order regardless of which categories are observed
  canonical = c("dead", "active", paste0("active_", labels), "undead", paste0("undead_", labels))
  statuses = c(intersect(canonical, statuses), setdiff(statuses, canonical))

  data = expand.grid(
    module = selected_modules,
    activity = statuses,
    stringsAsFactors = FALSE
  )
  data$counts = integer(nrow(data))
  data$bootstrap_mean = rep(NA_real_, nrow(data))

  for(i in seq_len(nrow(data))){
    rows = activity_table$module == data$module[i]
    data$counts[i] = sum(rows & activity_table$activity == data$activity[i])
    col = paste0("pct_", data$activity[i])
    if(show_bootstrap && col %in% pct_cols){
      data$bootstrap_mean[i] = sum(activity_table[[col]][rows]) / 100
    }
  }

  has_bootstrap = show_bootstrap && any(is.finite(data$bootstrap_mean))
  data$text = if(has_bootstrap){
    paste0(data$counts, "\n(", format(round(data$bootstrap_mean, 1), nsmall = 1), ")")
  } else {
    as.character(data$counts)
  }
  data$module = factor(data$module, levels = selected_modules)
  data$activity = factor(data$activity, levels = rev(statuses))

  # selecting colours when none are supplied
  if(length(bubble_colors) == 1 && is.na(bubble_colors)){
    valid_colors = vapply(selected_modules, function(x){
      !inherits(try(grDevices::col2rgb(x), silent = TRUE), "try-error")
    }, logical(1))
    colors = if(all(valid_colors)) selected_modules else
      grDevices::hcl.colors(length(selected_modules), "Dark 3")
  } else {
    if(length(bubble_colors) != length(selected_modules)){
      stop("bubble_colors must have one color per selected module.")
    }
    colors = bubble_colors
  }
  names(colors) = selected_modules

  # plot
  p = ggplot2::ggplot(data, ggplot2::aes(x = .data$module, y = .data$activity, color = .data$module, size = .data$counts)) +
    ggplot2::geom_point() +
    ggplot2::geom_text(ggplot2::aes(label = .data$text), color = number_color, size = 4, lineheight = 0.9) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_size(range = scale_range) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_classic(base_size = cex) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  if(all(is.na(spca_data$module_colors))){
    p = p + ggplot2::theme(axis.text.x.top = ggplot2::element_blank())
  }

  print(p)
  invisible(p)
}
