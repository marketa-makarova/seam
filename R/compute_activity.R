#' Compute activity statuses of variables using SPCA
#' 
#' The activity column describes the original SPCA fit. When bootstrap fits are
#' available, pct_* columns give the percentage of bootstrap replicates assigned
#' to each status. Percentages sum to 100 per variable (up to rounding).
#'
#' @param spca_data Result of spca_all_modules() or spca_all_modules_grid().
#' @param return_information `TRUE` for refined undead classifications,
#'   `FALSE` for raw selection, or `"both"` to print both plots.
#' @param threshold Absolute within-group correlation threshold for reclassification to undead status.
#' @param labels Group labels, if present, in the same order as the expression data.
#' 
#' @return Data frame with variable_id, module, `w_<label>` columns, activity,
#'   and optional `pct_<status>` bootstrap percentages.
#'
#' @details
#' A variable is active in a group when its SPCA loading is nonzero. A variable
#' excluded in all groups is classified as undead in a given group only when it
#' has absolute correlation above `threshold` with at least one variable
#' selected by SPCA in that same group. Variables selected in other groups
#' do not contribute to this classification. 
#' 
#' Bootstrap resampling is used to assess the stability of variable selection
#' under repeated sampling. Selection frequencies describe the reproducibility 
#' of nonzero SPCA loadings and should not be interpreted as 
#' statistical significance or evidence of differential population support.
#'
#' With the presence of groups, activity labels are `dead`, `active`, `active_<label>`,
#' `undead`, and `undead_<label>`. The labels `active` and `undead` indicate
#' all groups. Bootstrap percentages use the selected variables from each
#' corresponding group and bootstrap replicate, with within-group correlations
#' calculated from the original expression matrices.
#' 
#' @examples
#' data("seam_demo", package = "seam")
#' fit = spca_all_modules(
#'  seam_demo$expression_data,
#'  seam_demo$module_assignment,
#'  lambda = 0.05, 
#'  bootstrap = FALSE
#' )
#' activity = compute_activity(
#'  fit, 
#'  labels = seam_demo$labels
#' )
#' table(
#'  activity$module, 
#'  activity$activity
#' )
#'
#' @export

compute_activity = function(spca_data, return_information = TRUE, threshold = 0.8, labels){

  if(!is.list(spca_data)){
    stop("spca_data must be a list.")
  } 
  if(threshold <= 0 || threshold >= 1 || !is.numeric(threshold) || length(threshold) != 1 || is.na(threshold)){
    stop("threshold argument value must be a number between 0 and 1.")
  }
  nsets = length(spca_data$vs)
  if(length(labels) != nsets || anyDuplicated(labels) || any(grepl("-", labels))){
    stop("Provide unique labels, one per group, without '-'.")
  }

  # prepare necessary variables
  modules = spca_data$module_colors
  boot = spca_data$bootstrap
  results = vector("list", length(modules))
  statuses = if(nsets == 1){
    c("dead", "active", if(isTRUE(return_information) || return_information == "both") "undead")
  } else {
    c("dead", "active", paste0("active_", labels),
      if(isTRUE(return_information) || return_information == "both") c("undead", paste0("undead_", labels)))
  }

  # prepare classification function
  classify = function(selected, adjacency){
    n_active = rowSums(selected)
    activity = rep("dead", nrow(selected))
    activity[n_active > 1] = "active"
    if(nsets == 1){
      activity[n_active == 1] = "active"
    } else {
      single = which(n_active == 1)
      if(length(single)){
        activity[single] = paste0("active_", labels[max.col(selected[single, , drop = FALSE])])
      }
    }

    if((isTRUE(return_information) || return_information == "both") && any(n_active == 0)){
      dead = which(n_active == 0)
      undead = matrix(FALSE, nrow(selected), nsets)
      for(set in seq_len(nsets)){
        # only variables selected in this same group can restore information
        undead[, set] = as.vector(adjacency[[set]] %*% as.integer(selected[, set])) > 0
      }
      undead[!seq_len(nrow(selected)) %in% dead, ] = FALSE
      n_undead = rowSums(undead)
      if(nsets == 1){
        activity[dead[n_undead[dead] > 0]] = "undead"
      } else {
        activity[dead[n_undead[dead] > 1]] = "undead"
        single_undead = dead[n_undead[dead] == 1]
        if(length(single_undead)){
          activity[single_undead] = paste0(
            "undead_", labels[max.col(undead[single_undead, , drop = FALSE])]
          )
        }
      }
    } 
    return(activity)
  }

  # classify activity
  for(module in seq_along(modules)){
    variable_ids = colnames(spca_data$module_matrices[[1]][[module]])
    if(is.null(variable_ids)){
      variable_ids = as.character(seq_len(ncol(spca_data$module_matrices[[1]][[module]])))
    }
    p = length(variable_ids)
    weights = matrix(NA_real_, p, nsets, dimnames = list(variable_ids, labels))
    adjacency = vector("list", nsets)

    for(set in seq_len(nsets)){
      x = spca_data$module_matrices[[set]][[module]]
      w = spca_data$vs[[set]][[module]]
      if(!is.null(names(w))){
        w = w[variable_ids]
      } else if(!identical(colnames(x), variable_ids)){
        stop("Module variable order differs across groups.")
      }
      weights[, set] = w
      if(isTRUE(return_information) || return_information == "both"){
        x = as.matrix(x[, variable_ids, drop = FALSE])
        r = suppressWarnings(abs(stats::cor(x)))
        adjacency[[set]] = r > threshold & is.finite(r)
        diag(adjacency[[set]]) = FALSE
      }
    }

    out = data.frame(variable_id = variable_ids,
                     module = rep(modules[module], p),
                     stringsAsFactors = FALSE)
    for(set in seq_len(nsets)) out[[paste0("w_", labels[set])]] = weights[, set]
    out$activity = classify(weights != 0, adjacency)

    if(!is.null(boot) && !is.null(boot$vs)){
      bootstrap_n = nrow(boot$vs[[1]][[module]])
      counts = matrix(0L, p, length(statuses), dimnames = list(NULL, statuses))
      for(b in seq_len(bootstrap_n)){
        selected = matrix(FALSE, p, nsets)
        for(set in seq_len(nsets)){
          boot_weights = boot$vs[[set]][[module]]
          if(!is.null(colnames(boot_weights))){
            selected[, set] = boot_weights[b, variable_ids] != 0
          } else {
            selected[, set] = boot_weights[b, ] != 0
          }
        }
        current = classify(selected, adjacency)
        for(k in seq_along(statuses)){
          counts[, k] = counts[, k] + as.integer(current == statuses[k])
        }
      }
      for(k in seq_along(statuses)){
        out[[paste0("pct_", statuses[k])]] = 100 * counts[, k] / bootstrap_n
      }
    }
    results[[module]] = out
  }

  activity_table = do.call(rbind, results)
  rownames(activity_table) = NULL
  return(activity_table)

}
