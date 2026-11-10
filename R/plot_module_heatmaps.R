#' Plot module correlation heatmaps and sparsified heatmaps
#'
#' The SPCA fit determines which variables are excluded. When bootstrap
#' results are present, masked cells are shaded according to how often either
#' variable received a zero weight across bootstrap replicates.
#'
#' @param spca_data Result of spca_all_modules() or spca_all_modules_grid().
#' @param labels Group labels, if present, in the same order as the expression data.
#' @param modules Module labels to plot, or `"all"` (default).
#' @param return_information `TRUE` for refined undead classifications,
#'   `FALSE` for raw selection, or `"both"` to print both plots.
#' @param threshold Absolute within-group correlation threshold for reclassification to undead status.
#' @param cex Base text size.
#' @param midpoint Midpoint of the distance color scale.
#' @param col_low,col_mid,col_high Colors of the distance scale.
#' 
#' @return A named list of ggplot objects is returned invisibly. Plots are printed separately.
#' 
#' @examples
#' \donttest{
#' data("seam_demo", package = "seam")
#' fit = spca_all_modules(
#'  seam_demo$expression_data,
#'  seam_demo$module_assignment,
#'  lambda = 0.05, 
#'  bootstrap = FALSE
#' )
#' plot_module_heatmaps(
#'  fit, 
#'  labels = seam_demo$labels,
#'  modules = "blue",
#'  return_information = FALSE
#' )
#' }
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

  if(!is.list(spca_data)) stop("Invalid spca_data.")
  if(anyDuplicated(labels) || any(grepl("-", labels))){
    stop("Provide unique labels, one per group, without '-'.")
  }
  if(length(modules) == 1 && modules == "all"){
    selected = seq_along(spca_data$module_colors)
  } else {
    selected = match(modules, spca_data$module_colors)
    if(anyNA(selected)) stop("Unknown module in modules.")
  }

  # prepare plot data
  activity = NULL
  if(return_information){
    activity = compute_activity(spca_data, labels = labels, return_information = TRUE, threshold = threshold)
  }

  boot = spca_data$bootstrap
  plots = list()

  for(set in seq_along(labels)){
    for(module in selected){
      x = as.matrix(spca_data$module_matrices[[set]][[module]])
      variables = colnames(x)
      if(is.null(variables)) variables = as.character(seq_len(ncol(x)))
      weights = spca_data$vs[[set]][[module]]
      original_dead = weights == 0
      names(original_dead) = variables

      diss = 1 - abs(stats::cor(x))
      diss[!is.finite(diss)] = NA_real_
      if(ncol(x) > 1 && !anyNA(diss)){
        order = stats::hclust(stats::as.dist(diss), method = "average")$order
      } else {
        order = seq_along(variables)
      }
      ordered = variables[order]

      dead_frequency = as.numeric(original_dead)
      names(dead_frequency) = variables
      if(!is.null(boot) && !is.null(boot$vs)){
        bootstrap_weights = boot$vs[[set]][[module]]
        dead_frequency = colMeans(bootstrap_weights == 0, na.rm = TRUE)
        names(dead_frequency) = colnames(bootstrap_weights)
        dead_frequency = dead_frequency[variables]
      }

      returned = rep(FALSE, length(variables))
      if(return_information){
        current = activity[activity$module == spca_data$module_colors[module], , drop = FALSE]
        current = current[match(variables, current$variable_id), , drop = FALSE]
        returned = original_dead & (current$activity == "undead" | current$activity == paste0("undead_", labels[set]))
        returned[is.na(returned)] = FALSE
      }

      # distance scale is always based on the full expression matrix
      # no correlations are recomputed after SPCA
      make_panel = function(panel, keep){
        vars = ordered[keep[match(ordered, variables)]]
        if(!length(vars)) return(NULL)
        grid = expand.grid(x = vars, y = vars, stringsAsFactors = FALSE)
        grid$distance = diss[cbind(match(grid$x, variables), match(grid$y, variables))]
        grid$panel = panel
        grid$masked = FALSE
        grid$frequency = NA_real_
        if(panel %in% c("sparse", "returned")){
          excluded = if(panel == "sparse") original_dead else original_dead & !returned
          grid$masked = excluded[grid$x] | excluded[grid$y]
          grid$frequency = pmax(dead_frequency[grid$x], dead_frequency[grid$y])
        }
        grid
      }

      panels = list(
        make_panel("full", rep(TRUE, length(variables))),
        make_panel("sparse", rep(TRUE, length(variables)))
      )
      if(return_information){
        panels[[length(panels) + 1]] = make_panel("returned", rep(TRUE, length(variables)))
      }
      keep = !original_dead | returned
      panels[[length(panels) + 1]] = make_panel("reduced", keep)
      plot_data = do.call(rbind, Filter(Negate(is.null), panels))
      plot_data$panel = factor(plot_data$panel, levels = c("full", "sparse", "returned", "reduced"))
      plot_data$x = factor(plot_data$x, levels = ordered)
      plot_data$y = factor(plot_data$y, levels = rev(ordered))

      masked_data = plot_data[plot_data$masked, , drop = FALSE]
      p = ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$x, y = .data$y)) +
        ggplot2::geom_tile(ggplot2::aes(fill = .data$distance)) +
        ggplot2::scale_fill_gradient2(
          low = col_low, mid = col_mid, high = col_high,
          midpoint = midpoint, limits = c(0, 1),
          name = "Distance")

      if(nrow(masked_data)){
        if(!is.null(boot) && !is.null(boot$vs)){
          p = p + ggplot2::geom_tile(data = masked_data, ggplot2::aes(alpha = .data$frequency), fill = "grey15", inherit.aes = TRUE) +
            ggplot2::scale_alpha_continuous(
              range = c(0.15, 0.85), limits = c(0, 1),
              name = "Dead frequency",
              labels = function(x) paste0(round(100 * x), "%")
            )
        } else {
          p = p + ggplot2::geom_tile(data = masked_data, fill = "grey50", alpha = 0.65)
        }
      }

      n_dead = sum(original_dead)
      n_returned = sum(returned)
      p = p +
        ggplot2::facet_wrap(~ panel, scales = "free", nrow = 1) +
        ggplot2::scale_x_discrete(drop = TRUE, expand = ggplot2::expansion(add = 0)) +
        ggplot2::scale_y_discrete(drop = TRUE, expand = ggplot2::expansion(add = 0)) +
        ggplot2::theme_classic(base_size = cex) +
        ggplot2::theme(axis.text = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.title = ggplot2::element_blank(),
                       panel.spacing = grid::unit(0.5, "lines")) +
        ggplot2::labs(title = paste0(labels[set], " ", spca_data$module_colors[module],
                                      " module (", ncol(x), "): ", n_dead,
                                      " excluded", if(return_information)
                                      paste0(", ", n_returned, " restored") else "",
                                      ", final size ", sum(keep)))

      key = paste(labels[set], spca_data$module_colors[module], sep = "_")
      plots[[key]] = p
      print(p)
    }
  }
  invisible(plots)
}
