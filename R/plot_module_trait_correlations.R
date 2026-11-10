#' Plot module-trait correlations for sparse and classic eigengenes
#'
#' Cell colors represent strength of correlations. 
#' When bootstrap results exist, 95% percentile intervals are also included.
#'
#' @param module_trait_correlations Result of compute_module_trait_correlations().
#' @param labels Group labels, if present, in the same order as the expression data.
#' @param modules Module labels to plot, or `"all"` (default).
#' @param mfrow Number of rows and columns per page (grid layout).
#' @param main Optional title. If NA, titles identify method and group.
#' @param plot_avg_differences Show per-module RMSE across traits between
#'   classic and SPCA correlations; requires mode = "both".
#' @param mode "spca", "classic", or "both".
#' @param color_palette Colors for the diverging correlation scale.
#' @param cex Base text size.
#' @param show_legend Logical. Whether to display color legends.
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
#' mtc = compute_module_trait_correlations(
#'  seam_demo$trait_data, 
#'  fit, 
#'  seam_demo$module_assignment,
#'  mode = "both"
#' )
#' plot_module_trait_correlations(
#'  mtc, 
#'  labels = seam_demo$labels,
#'  mode = "spca"
#' )
#' 
#' @export
plot_module_trait_correlations = function(
    module_trait_correlations,
    labels,
    modules = "all",
    mfrow = c(1, 2),
    main = NA,
    plot_avg_differences = FALSE,
    mode = "spca",
    color_palette = grDevices::colorRampPalette(c("blue", "white", "red"))(50),
    cex = 1.2,
    show_legend = TRUE){

  mode = tolower(mode)
  if(!mode %in% c("spca", "classic", "both")) stop("Mode must be 'spca', 'classic', or 'both'.")
  if(plot_avg_differences && mode != "both") stop("plot_avg_differences requires mode = 'both'.")
  if(length(mfrow) != 2 || any(mfrow < 1) || any(mfrow != as.integer(mfrow))) stop("mfrow must contain two positive integers.")
  methods = if(mode == "both") c("classic", "spca") else mode
  for(method in methods){
    if(is.null(module_trait_correlations[[paste0("mtc_", method)]])){
      stop(paste("Missing correlations for", method, "mode."))
    }
  }
  first = module_trait_correlations[[paste0("mtc_", methods[1])]]
  n_sets = length(first)
  if(length(labels) != n_sets || anyDuplicated(labels)) stop("Provide unique labels, one per group, without '-'.")
  all_modules = rownames(first[[1]])
  if(is.null(all_modules)) all_modules = as.character(seq_len(nrow(first[[1]])))
  selected = if(length(modules) == 1 && modules == "all"){
    seq_along(all_modules)
  } else {
    match(modules, all_modules)
  }
  if(anyNA(selected) || !length(selected)) stop("Unknown or empty module selection.")

  # prepare plot data
  module_names = all_modules[selected]
  boot = module_trait_correlations$bootstrap
  sizes = unlist(module_trait_correlations$module_sizes)
  reduced_sizes = module_trait_correlations$reduced_module_sizes
  plots = list()

  make_plot = function(values, annotations, title, legend_title, limits){
    nr = nrow(values)
    nc = ncol(values)
    d = data.frame(
      module = factor(rep(rownames(values), times = nc), levels = rev(rownames(values))),
      trait = factor(rep(colnames(values), each = nr), levels = colnames(values)),
      value = as.vector(values),
      text = as.vector(annotations)
    )
    ggplot2::ggplot(d, ggplot2::aes(x = .data$trait, y = .data$module, fill = .data$value)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::geom_text(ggplot2::aes(label = .data$text), size = cex * 3.5) +
      ggplot2::scale_fill_gradientn(colors = color_palette, limits = limits, name = legend_title, na.value = "grey90") +
      ggplot2::labs(x = NULL, y = NULL, title = title) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::scale_x_discrete(position = "top", expand = expansion(add = 0)) +
      ggplot2::scale_y_discrete(expand = expansion(add = 0)) +
      ggplot2::theme(
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x.top = element_text(
          angle = 45,
          hjust = 0,
          vjust = 0
        ),
        aspect.ratio = NULL
      )
  }

  for(method in methods){
    original = module_trait_correlations[[paste0("mtc_", method)]]
    boot_method = if(!is.null(boot)) boot[[paste0("mtc_", method)]] else NULL

    for(set in seq_len(n_sets)){
      mat = original[[set]][selected, , drop = FALSE]
      if(method == "spca"){
        counts = unlist(reduced_sizes[[set]])[selected]
      } else {
        counts = sizes[selected]
      }
      rownames(mat) = paste0(module_names, " (", counts, ")")
      annotations = matrix(sprintf("%.2f", mat), nrow = nrow(mat), dimnames = dimnames(mat))
      annotations[is.na(mat)] = "NA"

      if(!is.null(boot_method)){
        arr = boot_method[[set]][selected, , , drop = FALSE]
        for(i in seq_len(nrow(mat))){
          for(j in seq_len(ncol(mat))){
            v = as.numeric(arr[i, j, ])
            v = v[is.finite(v)]
            if(length(v)){
              q = stats::quantile(v, c(0.025, 0.975), names = FALSE)
              annotations[i, j] = paste0(annotations[i, j], "\n[", sprintf("%.2f", q[1]), ", ", sprintf("%.2f", q[2]), "]")
            }
          }
        }
      }

      title = if(length(main) == 1 && is.na(main)){
        paste("Module-trait correlations by", toupper(method), "for", labels[set])
      } else {
        paste(main, "-", toupper(method), "-", labels[set])
      }
      key = paste(method, labels[set], sep = "_")
      plots[[key]] = make_plot(mat, annotations, title, "Correlation", c(-1, 1))
    }
  }

  if(plot_avg_differences){
    rmse = matrix(NA_real_, length(selected), n_sets,
                  dimnames = list(module_names, labels))
    rmse_text = matrix(NA_character_, length(selected), n_sets,
                       dimnames = dimnames(rmse))
    for(set in seq_len(n_sets)){
      classic = module_trait_correlations$mtc_classic[[set]][selected, , drop = FALSE]
      sparse = module_trait_correlations$mtc_spca[[set]][selected, , drop = FALSE]
      for(i in seq_along(selected)){
        differences = classic[i, ] - sparse[i, ]
        rmse[i, set] = if(all(is.na(differences))) NA_real_ else
          sqrt(mean(differences^2, na.rm = TRUE))
        rmse_text[i, set] = sprintf("%.2f", rmse[i, set])
      }

      if(!is.null(boot) && !is.null(boot$mtc_classic) && !is.null(boot$mtc_spca)){
        a = boot$mtc_classic[[set]][selected, , , drop = FALSE]
        b = boot$mtc_spca[[set]][selected, , , drop = FALSE]
        for(i in seq_along(selected)){
          values = vapply(seq_len(dim(a)[3]), function(k){
            diff = a[i, , k] - b[i, , k]
            if(all(is.na(diff))) return(NA_real_)
            sqrt(mean(diff^2, na.rm = TRUE))
          }, numeric(1))
          values = values[is.finite(values)]
          if(length(values)){
            q = stats::quantile(values, c(0.025, 0.975), names = FALSE)
            rmse_text[i, set] = paste0(rmse_text[i, set], "\n[", sprintf("%.2f", q[1]), ", ", sprintf("%.2f", q[2]), "]")
          }
        }
      }
    }
    rownames(rmse) = paste0(module_names, " (", sizes[selected], ")")
    rownames(rmse_text) = rownames(rmse)
    plots$rmse = make_plot(rmse, rmse_text, "RMSE between classic and SPCA correlations", "RMSE", c(0, 2))
  }

  for(p in plots){
    print(p)
  }

  invisible(plots)
}
