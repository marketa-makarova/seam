#' Plot variable activity statuses as a stacked barplot
#'
#' Displays counts of original SPCA activity classifications across modules
#' as a stacked barplot. Statuses are differentiated by colour.
#'
#' @param spca_data Result of spca_all_modules() or spca_all_modules_grid().
#' @param return_information `TRUE` for refined undead classifications,
#'   `FALSE` for raw selection, or `"both"` to print both plots.
#' @param modules Module labels to plot, or `"all"` (default).
#' @param threshold Absolute within-group correlation threshold for reclassification to undead status.
#' @param labels Group labels, if present, in the same order as the expression data.
#' @param cex Base text size.
#' @param col_both_active Color for variables selected in all groups.
#' @param col_active Colors for variables selected in exactly one group;
#'   one color per group, in `labels` order.
#' @param col_dead Color for variables excluded in every group.
#' @param col_both_undead Color for variables undead in all groups.
#' @param col_undead Colors for variables undead in exactly one group;
#'   one color per group, in `labels` order.
#' @param xlab,ylab Axis labels.
#' 
#' @return A named list of ggplot objects is returned invisibly. Multiple plots are printed separately.
#' 
#' @examples
#' data("seam_demo", package = "seam")
#' fit = spca_all_modules(
#'  seam_demo$expression_data,
#'  seam_demo$module_assignment,
#'  lambda = 0.05, 
#'  bootstrap = FALSE
#' )
#' plot_activity_bar(
#'  fit, 
#'  labels = seam_demo$labels
#' )
#' 
#' @export
plot_activity_bar = function(
    spca_data,
    return_information = TRUE,
    modules = "all",
    threshold = 0.6,
    labels,
    cex = 20,
    col_both_active = "forestgreen",
    col_active = NULL,
    col_dead = "grey55",
    col_both_undead = "lightgrey",
    col_undead = NULL,
    xlab = "Module",
    ylab = "Count"){

  if(!is.list(spca_data) || is.null(spca_data$module_colors) || is.null(spca_data$PCs)) stop("Invalid spca_data.")
  nsets = length(spca_data$PCs)
  if(length(labels) != nsets || anyNA(labels) || anyDuplicated(labels) || any(!nzchar(labels))){
    stop("Provide unique labels, one per group, without '-'.")
  }
  if(!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold) || threshold <= 0 || threshold >= 1){
    stop("threshold argument value must be a number between 0 and 1.")
  }
  if(!(return_information == TRUE || return_information == FALSE || return_information == "both")){
    stop('return_information must be TRUE, FALSE, or "both".')
  }
  all_modules = as.character(spca_data$module_colors)
  selected = if(length(modules) == 1 && identical(modules, "all"))
    all_modules else as.character(modules)
  if(!length(selected) || anyNA(selected) || anyDuplicated(selected) || !all(selected %in% all_modules)){
    stop("Unknown module in modules.")
  } 
  if(is.null(col_active)) col_active = grDevices::hcl.colors(nsets, "Dark 3")
  if(is.null(col_undead)) col_undead = grDevices::adjustcolor(col_active, alpha.f = 0.5)
  if(length(col_active) != nsets || length(col_undead) != nsets){
    stop("col_active and col_undead must contain one color per group.")
  }

  # prepare colours
  palette = c(
    dead = col_dead,
    active = col_both_active,
    stats::setNames(col_active, paste0("active_", labels)),
    undead = col_both_undead,
    stats::setNames(col_undead, paste0("undead_", labels))
  )

  # prepare plot data
  sizes = spca_data$module_sizes[match(selected, all_modules)]
  tick_labels = stats::setNames(paste0(selected, " (", sizes, ")"), selected)

  # plot
  make_plot = function(refined){ # as a function to avoid duplicating code

    activity = compute_activity(
      spca_data = spca_data,
      return_information = refined,
      threshold = threshold,
      labels = labels
    )

    activity = activity[as.character(activity$module) %in% selected, , drop = FALSE]
    activity$module = factor(as.character(activity$module), levels = selected)
    observed = unique(as.character(activity$activity))
    status_order = c("active", paste0("active_", labels), "dead", "undead", paste0("undead_", labels))
    activity$activity = factor(as.character(activity$activity), levels = status_order)

    p = ggplot2::ggplot(activity, ggplot2::aes(x = .data$module, fill = .data$activity)) +
      ggplot2::geom_bar() +
      ggplot2::scale_fill_manual(values = palette, drop = TRUE) +
      ggplot2::scale_x_discrete(labels = tick_labels, drop = FALSE) +
      ggplot2::theme_classic(base_size = cex) +
      ggplot2::labs(x = xlab, y = ylab, fill = "Activity")
    if(length(all_modules) == 1 && identical(all_modules, "all")){
      p = p + ggplot2::theme(axis.ticks.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank()) +
        ggplot2::labs(x = NULL)
    }

    print(p)
    p
  }

  if(identical(return_information, "both")){
    plots = list(raw = make_plot(FALSE), refined = make_plot(TRUE))
    return(invisible(plots))
  }

  invisible(make_plot(return_information))
}
