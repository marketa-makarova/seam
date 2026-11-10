#' Compare SPCA variable weights between two groups
#'
#' Draws a lollipop plot of the original SPCA weight for each
#' variable in two groups. Variables are ranked by between-group weight
#' difference, weight magnitude, or bootstrap selection instability. When
#' bootstrap results are available, optional 95% percentile intervals and
#' selection frequencies describe the variability of each loading.
#'
#' @param spca_data Result of spca_all_modules() or spca_all_modules_grid().
#' @param labels Group labels, if present, in the same order as the expression data.
#' @param modules Module labels to plot, or `"all"` (default).
#' @param n Maximum number of variables shown per module. Default is 20.
#' @param bootstrap Logical. Whether to use available bootstrap results.
#'   Default is `TRUE`. If results are absent, only original weights are shown.
#' @param show_intervals Logical. Whether to draw 95% percentile intervals for bootstrap
#'   weights. Default is `TRUE`.
#' @param show_stability Logical. Whether to encode bootstrap selection frequency as
#'   interval opacity (darker means more frequently selected). Default is `TRUE`.
#' @param rank_by One of `"difference"` (absolute difference in aligned
#'   original weights), `"magnitude"` (largest absolute weight in either
#'   group), or `"instability"` (lowest mean bootstrap selection frequency).
#'   `"instability"` requires bootstrap results.
#' @param align_signs Logical. Whether to align the second group's original and bootstrap
#'   loading vectors to the first group by the sign of their dot product.
#'   Default is `TRUE`.
#' @param cex Base text size.
#' @param col Optional colors for the groups.
#' @param xlab Label for the weight axis.
#'
#' @details
#' The same variables must be present in both groups for each module.
#' Group B is sign-aligned to group A *for visualization only*; the input
#' object is not modified. Bootstrap loading vectors are assumed to have
#' already been aligned to their respective original solutions by the SPCA
#' fitting function. The sign of an SPCA component is arbitrary, and the
#' dot-product alignment can be ambiguous for nearly orthogonal solutions.
#'
#' Selection frequency is the percentage of bootstrap fits with a nonzero
#' loading. Intervals are the 2.5th and 97.5th percentiles of the full
#' bootstrap loading distribution, including zero values. These percentile
#' ranges are descriptive and may not have nominal confidence coverage after
#' sparse variable selection.
#'
#' @return A named list of ggplot objects is returned invisibly. Plots are printed separately.
#' 
#' @examples
#' data("seam_demo", package = "seam")
#' fit = spca_all_modules(
#'  seam_demo$expression_data,
#'  seam_demo$module_assignment,
#'  lambda = 0.05, 
#'  bootstrap = FALSE
#' )
#' insight_weights_lollipop(
#'  fit, labels = seam_demo$labels,
#'  modules = "blue", 
#'  n = 8,
#'  bootstrap = FALSE
#' )
#' 
#' @export
insight_weights_lollipop = function(
    spca_data,
    labels,
    modules = "all",
    n = 20,
    bootstrap = TRUE,
    show_intervals = TRUE,
    show_stability = TRUE,
    rank_by = c("difference", "magnitude", "instability"),
    align_signs = TRUE,
    cex = 14,
    col = c("#377EB8", "#E69F00"),
    xlab = "SPCA weight"){

  rank_by = match.arg(rank_by)
  if(length(spca_data$vs) != 2 || length(labels) != 2 || anyDuplicated(labels)){
    stop("This plot requires exactly two groups with distinct labels.")
  }
  if(length(col) != 2 || n < 1 || length(n) != 1){
    stop("Provide two colors and a positive value of n.")
  }
  available = as.character(spca_data$module_colors)
  selected = if(length(modules) == 1 && modules == "all") available else as.character(modules)
  if(!all(selected %in% available)){
    stop("Unknown module in modules.")
  }
  has_boot = isTRUE(bootstrap) && !is.null(spca_data$bootstrap) &&
    !is.null(spca_data$bootstrap$vs)
  if(rank_by == "instability" && !has_boot){
    stop("rank_by = 'instability' requires bootstrap weights.")
  }

  # draw both group colors in every key of the stability legend
  stability_key = function(data, params, size){
    opacity = if(is.null(data$alpha) || is.na(data$alpha[1])) 1 else data$alpha[1]
    grid::grobTree(
      grid::segmentsGrob(
        x0 = grid::unit(0.05, "npc"), x1 = grid::unit(0.95, "npc"),
        y0 = grid::unit(0.72, "npc"), y1 = grid::unit(0.72, "npc"),
        gp = grid::gpar(col = grDevices::adjustcolor(col[1], alpha.f = opacity), lwd = 2.5)
      ),
      grid::segmentsGrob(
        x0 = grid::unit(0.05, "npc"), x1 = grid::unit(0.95, "npc"),
        y0 = grid::unit(0.28, "npc"), y1 = grid::unit(0.28, "npc"),
        gp = grid::gpar(col = grDevices::adjustcolor(col[2], alpha.f = opacity), lwd = 2.5)
      )
    )
  }

  # prepare plot data
  plots = list()
  for(module in selected){
    w_a = spca_data$vs[[1]][[module]]
    w_b = spca_data$vs[[2]][[module]]
    ids_a = names(w_a)
    ids_b = names(w_b)
    if(is.null(ids_a)) ids_a = colnames(spca_data$module_matrices[[1]][[module]])
    if(is.null(ids_b)) ids_b = colnames(spca_data$module_matrices[[2]][[module]])
    if(is.null(ids_a) || is.null(ids_b) || anyDuplicated(ids_a) || anyDuplicated(ids_b) || !setequal(ids_a, ids_b)){
      stop(paste("Variable identifiers do not match between groups for module", module))
    }
    w_a = stats::setNames(as.numeric(w_a), ids_a)
    w_b = stats::setNames(as.numeric(w_b), ids_b)[ids_a] # make sure variables match

    # align
    dot = sum(w_a * w_b)
    direction = if(isTRUE(align_signs) && is.finite(dot) && dot < 0) -1 else 1
    w_b = direction * w_b

    # compute intervals
    wide = data.frame(variable_id = ids_a, A = unname(w_a), B = unname(w_b))
    if(has_boot){
      boot_a = spca_data$bootstrap$vs[[1]][[module]]
      boot_b = spca_data$bootstrap$vs[[2]][[module]]
      if(is.null(colnames(boot_a))) colnames(boot_a) = ids_a
      if(is.null(colnames(boot_b))) colnames(boot_b) = ids_b
      if(!setequal(colnames(boot_a), ids_a) || !setequal(colnames(boot_b), ids_a)){
        stop(paste("Bootstrap variable identifiers do not match for module", module))
      }
      boot_a = boot_a[, ids_a, drop = FALSE]
      boot_b = direction * boot_b[, ids_a, drop = FALSE]
      wide$stability_A = colMeans(boot_a != 0, na.rm = TRUE) * 100
      wide$stability_B = colMeans(boot_b != 0, na.rm = TRUE) * 100
      wide$lower_A = apply(boot_a, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
      wide$upper_A = apply(boot_a, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
      wide$lower_B = apply(boot_b, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
      wide$upper_B = apply(boot_b, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
    }

    # rank
    wide$ranking = switch(rank_by,
      difference = abs(wide$A - wide$B),
      magnitude = pmax(abs(wide$A), abs(wide$B)),
      instability = 100 - (wide$stability_A + wide$stability_B) / 2
    )
    wide = wide[order(-wide$ranking, wide$variable_id), , drop = FALSE]
    wide = utils::head(wide, n)
    wide$variable_id = factor(wide$variable_id, levels = rev(wide$variable_id))

    long = data.frame(
      variable_id = rep(wide$variable_id, 2),
      group = factor(rep(labels, each = nrow(wide)), levels = labels),
      weight = c(wide$A, wide$B)
    )
    if(has_boot){
      long$stability = c(wide$stability_A, wide$stability_B)
      long$lower = c(wide$lower_A, wide$lower_B)
      long$upper = c(wide$upper_A, wide$upper_B)
    }

    # plot
    p = ggplot2::ggplot() +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60")
    if(has_boot && isTRUE(show_intervals)){
      if(isTRUE(show_stability)){
        p = p + ggplot2::geom_errorbar(
          data = long, ggplot2::aes(
            y = .data$variable_id, 
            xmin = .data$lower,
            xmax = .data$upper, 
            color = .data$group,
            alpha = .data$stability),
          orientation = "y", width = 0.22, linewidth = 0.9,
          position = ggplot2::position_dodge(width = 0.35),
          key_glyph = stability_key,
          show.legend = c(colour = FALSE, alpha = TRUE)
        ) + ggplot2::scale_alpha_continuous(
          name = "Selected in bootstrap (%)", limits = c(0, 100),
          range = c(0.2, 1),
          breaks = c(0, 25, 50, 75, 100)
        )
      } else {
        p = p + ggplot2::geom_errorbar(data = long, ggplot2::aes(
            y = .data$variable_id, 
            xmin = .data$lower,
            xmax = .data$upper, 
            color = .data$group),
          orientation = "y", width = 0.22, linewidth = 0.9,
          position = ggplot2::position_dodge(width = 0.35)
        )
      }
    }

    p = p + ggplot2::geom_point(data = long, ggplot2::aes(x = .data$weight, y = .data$variable_id, color = .data$group),
      size = 3, position = ggplot2::position_dodge(width = 0.35)
    )

    p = p +
      ggplot2::scale_color_manual(values = stats::setNames(col, labels), name = "Group") +
      ggplot2::labs(
        title = paste0(module, " module: comparative SPCA weights"),
        subtitle = if(isTRUE(align_signs)) "Group 2 aligned to group 1 by loading-vector dot product" else NULL, x = xlab, y = NULL) +
      ggplot2::theme_classic(base_size = cex) +
      ggplot2::theme(legend.position = "right")

    print(p)
    plots[[module]] = p
  }
  invisible(plots)
}
