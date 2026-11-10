#' Compute variable–trait correlations
#'
#' Computes Pearson correlations between expression variables and clinical traits,
#' alogside SPCA-weighted correlations (correlation multiplied by corresponding SPCA weight).
#' When bootstrap fits are available, adds bootstrap confidence intervals using
#' combined resampling of traits and samples.
#'
#' @param trait_data List of numeric trait data frames, one per group, with rows
#'   in the same order as the corresponding expression matrix.
#' @param spca_data Result of spca_all_modules() or spca_all_modules_grid().
#' @param return_information `TRUE` for refined undead classifications,
#'   `FALSE` for raw selection, or `"both"` to print both plots.
#' @param threshold Absolute within-group correlation threshold for reclassification to undead status.
#' @param labels Group labels, if present, in the same order as the expression data.
#'
#' @return List of data frames, one per group. Each contains variable_id,
#'   module, original SPCA weights, activity, and `<trait>_cor` and `<trait>_w_cor` columns. 
#' 
#' @examples
#' data("seam_demo", package = "seam")
#' fit = spca_all_modules(
#'  seam_demo$expression_data,
#'  seam_demo$module_assignment,
#'  lambda = 0.05, 
#'  bootstrap = FALSE
#' )
#' vtc = compute_variable_trait_correlations(
#'  seam_demo$trait_data, 
#'  fit, 
#'  labels = seam_demo$labels)
#' head(vtc[[1]])
#' 
#' @export
compute_variable_trait_correlations = function(
    trait_data,
    spca_data,
    return_information = TRUE,
    threshold = 0.8,
    labels){

  if(!is.list(spca_data) || !is.list(trait_data)){
    stop("spca_data and trait_data must be lists.")
  }
  n_sets = length(spca_data$module_matrices)
  if(length(trait_data) != n_sets || length(labels) != n_sets){
    stop("Provide one trait dataset and one label per group.")
  }

  # prepare necessary variables
  activity = compute_activity(
    spca_data = spca_data,
    return_information = return_information,
    threshold = threshold,
    labels = labels
  )
  boot = spca_data$bootstrap
  results = vector("list", n_sets)
  names(results) = labels

  # compute correlations
  for(set in seq_len(n_sets)){
    modules = spca_data$module_matrices[[set]]
    x = do.call(cbind, modules)
    traits = as.data.frame(trait_data[[set]])

    # additional checks
    if(nrow(traits) != nrow(x)){
      stop("Trait and expression data must have the same number of rows for each group.")
    }
    if(!is.null(rownames(traits)) && !is.null(rownames(x)) && !identical(rownames(traits), rownames(x))){
      stop("Trait and expression sample row names must be in the same order.")
    }
    if(anyDuplicated(colnames(traits)) || is.null(colnames(traits))){
      stop("Trait columns must have unique names.")
    }
    if(!all(vapply(traits, is.numeric, logical(1)))){
      stop("Trait columns must be numeric.")
    }

    # align expression columns to the activity table by module and variable ID
    # supports variable names that repeat in different modules
    keys = paste(activity$module, activity$variable_id, sep = "\r")
    x_keys = unlist(lapply(seq_along(modules), function(m){
      paste(spca_data$module_colors[m], colnames(modules[[m]]), sep = "\r")
    }), use.names = FALSE)
    order_index = match(keys, x_keys)
    if(anyNA(order_index) || anyDuplicated(x_keys)){
      stop("Cannot match expression columns to activity rows by module and variable_id.")
    }
    x = as.matrix(x[, order_index, drop = FALSE])
    out = activity
    original_weights = out[[paste0("w_", labels[set])]]

    # compute correlations
    for(t in seq_along(traits)){
      trait_name = names(traits)[t]
      y = traits[[t]]
      raw_cor = as.vector(suppressWarnings(stats::cor(x, y, use = "pairwise.complete.obs")))
      out[[paste0(trait_name, "_cor")]] = raw_cor
      out[[paste0(trait_name, "_w_cor")]] = raw_cor * original_weights

      # bootstrap correlations, if bootstrap is present
      if(!is.null(boot) && !is.null(boot$vs) && !is.null(boot$samples)){
        bootstrap_n = nrow(boot$samples[[set]])
        raw_boot = matrix(NA_real_, bootstrap_n, ncol(x))
        weighted_boot = matrix(NA_real_, bootstrap_n, ncol(x))
        for(b in seq_len(bootstrap_n)){
          sampled_rows = boot$samples[[set]][b, ]
          r = as.vector(suppressWarnings(stats::cor(x[sampled_rows, , drop = FALSE], y[sampled_rows], use = "pairwise.complete.obs")))
          weights_b = numeric(ncol(x))
          for(m in seq_along(modules)){
            ix = which(out$module == spca_data$module_colors[m])
            w = boot$vs[[set]][[m]][b, ]
            if(!is.null(colnames(boot$vs[[set]][[m]]))){
              w = w[match(out$variable_id[ix], colnames(boot$vs[[set]][[m]]))]
            }
            weights_b[ix] = w
          }
          raw_boot[b, ] = r
          weighted_boot[b, ] = r * weights_b
        }
        
        # make sure quantiles do not break with missing data
        quantile2 = function(z, prob){
          apply(z, 2, function(v){
            if(all(is.na(v))) return(NA_real_)
            as.numeric(stats::quantile(v, prob, na.rm = TRUE))
          })
        }

        out[[paste0(trait_name, "_cor_q025")]] = quantile2(raw_boot, 0.025)
        out[[paste0(trait_name, "_cor_q975")]] = quantile2(raw_boot, 0.975)
        out[[paste0(trait_name, "_w_cor_q025")]] = quantile2(weighted_boot, 0.025)
        out[[paste0(trait_name, "_w_cor_q975")]] = quantile2(weighted_boot, 0.975)
      }
    }
    results[[set]] = out
  }
  return(results)
}
