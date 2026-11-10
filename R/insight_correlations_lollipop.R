#' Lollipop plot of variable-trait correlations and SPCA-weighted correlations
#'
#' Selects module-trait pairs whose SPCA eigengene correlation
#' exceeds interest_threshold in at least one group. For each pair, plots
#' variable-trait correlations and/or correlations multiplied by SPCA weights.
#' Original estimates are used whether or not spca_data contains bootstrap fits.
#'
#' @param module_trait_correlations Result of compute_module_trait_correlations()
#'   calculated with mode = "spca" or "both".
#' @param spca_data Result of spca_all_modules() or spca_all_modules_grid().
#' @param type One of "classic" (unweighted), "weighted", or "both".
#' @param interest_threshold Miminal absolute module-trait correlation threshold.
#' @param cex Base text size.
#' @param cex_point Point size.
#' @param labels Group labels, if present, in the same order as the expression data.
#' @param n Number of top variable labels to display per plot.
#' @param xlab,ylab Axis labels.
#' @param col Optional colors for the groups.
#' @param col_active Color of segments for variables active in multiple groups.
#' @param y_dodge Number of rows used to dodge variable labels.
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
#' mtc = compute_module_trait_correlations(
#'  seam_demo$trait_data, 
#'  fit,
#'  seam_demo$module_assignment,
#'  mode = "both")
#' insight_correlations_lollipop(
#'  mtc, 
#'  fit, 
#'  labels = seam_demo$labels,
#'  interest_threshold = 0.2, 
#'  n = 5
#' )
#' }
#' 
#' @export

insight_correlations_lollipop <- function(
    module_trait_correlations,
    spca_data,
    type = "both",
    interest_threshold = 0.6,
    cex = 20,
    cex_point = 5,
    labels,
    n = 10,
    xlab = "",
    ylab = "Weighted correlations",
    col = NULL,
    col_active = 3,
    y_dodge = 2){

  if(!type %in% c("classic", "weighted", "both")){
    stop("Type must be 'classic', 'weighted', or 'both'.")
  }
  if(!is.list(spca_data) || is.null(module_trait_correlations$mtc_spca)){
    stop("spca_data and SPCA module-trait correlations are required.")
  }
  if(!is.numeric(interest_threshold) || length(interest_threshold) != 1 ||
     !is.finite(interest_threshold) || interest_threshold <= 0 || interest_threshold >= 1){
    stop("interest_threshold must be between 0 and 1.")
  }
  nsets = length(spca_data$vs)
  if(length(labels) != nsets || length(module_trait_correlations$mtc_spca) != nsets){
    stop("Group labels and correlation results must match spca_data.")
  }
  module_colors = spca_data$module_colors
  reference = module_trait_correlations$mtc_spca[[1]]
  if(is.null(colnames(reference))){
    stop("Module-trait correlation matrices need trait column names.")
  }

  # select module-trait combinations above interest_threshold
  interesting = matrix(FALSE, nrow(reference), ncol(reference), dimnames = dimnames(reference))
  for(set in seq_len(nsets)){
    current = module_trait_correlations$mtc_spca[[set]]
    interesting = interesting | (is.finite(current) & abs(current) > interest_threshold)
  }
  positions = which(interesting, arr.ind = TRUE)
  if(nrow(positions) == 0){
    message("No interesting correlations; lower interest_threshold.")
    return(invisible(list()))
  }

  # prepare necessary variables
  activity_table = compute_activity(spca_data = spca_data,
                                    labels = labels,
                                    return_information = FALSE)
  plots = list()
  group_colors = if(is.null(col)){
    stats::setNames(grDevices::hcl.colors(nsets, "Dark 3"), labels)
  } else {
    stats::setNames(rep(col, length.out = nsets), labels)
  }

  # prepare data for plotting
  for(k in seq_len(nrow(positions))){
    module_index = positions[k, 1]
    trait_index = positions[k, 2]
    module_name = module_colors[module_index]
    trait_name = colnames(reference)[trait_index]
    variables = colnames(spca_data$module_matrices[[1]][[module_index]])
    if(is.null(variables)){ stop("Expression matrices must have variable column names.") }
    cors = data.frame(variable_id = variables, module = module_name,
                      stringsAsFactors = FALSE)

    for(set in seq_len(nsets)){
      x = spca_data$module_matrices[[set]][[module_index]]
      trait = module_trait_correlations$trait_data[[set]][, trait_index]
      if(!identical(colnames(x), variables)){ stop("Variable order must be identical across groups within each module.") }
      cors[[paste0("cor_", labels[set])]] = as.vector(
        stats::cor(x, trait, use = "pairwise.complete.obs"))
    }

    activity = activity_table[
      activity_table$module == module_name & activity_table$variable_id %in% variables,
      c("variable_id", "module", "activity", paste0("w_", labels)),
      drop = FALSE
    ]
    cors = dplyr::left_join(cors, activity, by = c("variable_id", "module"))
    cors = cors[cors$activity != "dead" & !is.na(cors$activity), , drop = FALSE]
    if(nrow(cors) == 0) next

    for(set in seq_len(nsets)){
      cors[[paste0("w_cor_", labels[set])]] = cors[[paste0("w_", labels[set])]] * cors[[paste0("cor_", labels[set])]]
    }

    # rank variable labels by their largest absolute weighted correlation
    scores = apply(abs(as.matrix(cors[, paste0("w_cor_", labels), drop = FALSE])), 1, max, na.rm = TRUE)
    scores[!is.finite(scores)] = -Inf
    top = head(order(scores, decreasing = TRUE), n)
    top_ids = cors$variable_id[top]
    cors = cors[order(cors$activity, -scores), , drop = FALSE]
    variable_order = cors$variable_id
    label_map = stats::setNames(ifelse(variable_order %in% top_ids, variable_order, ""), variable_order)

    plot_data = list()
    for(set in seq_len(nsets)){
      for(method in c("unweighted", "weighted")){
        value_col = if(method == "weighted") paste0("w_cor_", labels[set]) else paste0("cor_", labels[set])
        plot_data[[length(plot_data) + 1]] = data.frame(
          variable_id = cors$variable_id,
          activity = cors$activity,
          group = labels[set],
          method = method,
          value = cors[[value_col]],
          stringsAsFactors = FALSE
        )
      }
    }
    plot_data = do.call(rbind, plot_data)
    plot_data$variable_id = factor(plot_data$variable_id, levels = rev(variable_order))
    plot_data$group = factor(plot_data$group, levels = labels)

    if(type == "classic") plot_data = plot_data[plot_data$method == "unweighted", ]
    if(type == "weighted") plot_data = plot_data[plot_data$method == "weighted", ]

    correlations = vapply(seq_len(nsets), function(set){
      sprintf("%s = %.2f", labels[set], module_trait_correlations$mtc_spca[[set]][module_index, trait_index])
    }, character(1))
    title = if(length(spca_data$module_colors) == 1 && spca_data$module_colors[1] == "all"){
      paste0("Correlation with ", trait_name)
    } else {
      paste0(module_name, " module correlated with ", trait_name)
    }
    title = paste0(title, "\nSPCA correlations: ", paste(correlations, collapse = ", "))

    p = ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$variable_id, y = .data$value, color = .data$group)) +
      ggplot2::geom_segment(ggplot2::aes(xend = .data$variable_id, yend = 0), na.rm = TRUE) +
      ggplot2::geom_point(size = cex_point, na.rm = TRUE) +
      ggplot2::coord_flip() +
      ggplot2::geom_hline(yintercept = 0, color = 8, linewidth = 1) +
      ggplot2::scale_color_manual(values = group_colors, name = "Group") +
      ggplot2::scale_x_discrete(labels = label_map, guide = ggplot2::guide_axis(n.dodge = y_dodge)) +
      ggplot2::theme_classic(base_size = cex) +
      ggplot2::labs(x = xlab, y = if(type == "classic") "Correlations" else ylab, title = title)

    if(type == "both"){
      p = p + ggplot2::facet_wrap(~ method, scales = "free")
    }
    plots[[paste(module_name, trait_name, sep = "__")]] = p
    print(p)
  }

  invisible(plots)
}
