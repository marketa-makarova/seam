#' Reconstruct correlations with clinical traits using top contributing variables
#'
#' This function visualizes how well module–trait correlations can be
#' reconstructed using only the top \code{n} variables 
#' with the highest absolute values of SPCA weights. It generates line plots showing
#' how correlation strength evolves as more top variables 
#' are included in the reconstruction. If modules/clusters are present, 
#' each module will have its own plot created.
#'
#' @param module_trait_correlations A list of module trait correlations
#' computed using the \code{compute_module_trait_correlations} function.
#' @param spca_data A list containing SPCA decomposition results from the 
#' function \code{spca_all_modules} or \code{spca_all_modules_grid}.
#' @param labels Character vector. Provides group labels. 
#' These labels are used for activity assignment.
#' @param col Color or vector of colors if groups are present. 
#' Used for plotting points.
#' @param cex Numeric. Font size for plot elements. 
#' Passed to \code{theme_classic()} argument. 
#' Defaults to 20.
#' @param top_n_upper Integer. 
#' Maximum number of variables to include in the reconstruction. 
#' Defaults to 50.
#' @param top_n_lower Integer. 
#' Minimum number of variables to include. 
#' Defaults to 10.
#' @param interest_threshold Numeric. Only module–trait correlations greater
#' than this value (in absolute terms) are considered interesting and plotted. 
#' Default is 0.6.
#' @param return_information Logical. If \code{FALSE}, returns only raw activity 
#' assignments. If \code{TRUE} performs additional correlation-based checks to 
#' detect **undead** variables and updates the activity status accordingly. 
#' Defaults to \code{TRUE}.
#' @param threshold Numeric. The correlation threshold used to determine 
#' whether a **dead** variable highly correlated with an active variable, 
#' suggesting it should be assigned the status **undead**. 
#' Default value is 0.8.
#' @param xlab Character. X-axis label.
#' @param ylab Character. Y-axis label.
#' @param col_connective Color of a horizontal line at 0. 
#' Defaults to 8 (grey).
#'
#' @details
#' For each correlation above the threshold:
#' \itemize{
#'   \item The top contributing variables (based on absolute weighted SPCA correlation) are identified.
#'   \item The correlation is re-calculated for using a subset of these \code{n}
#'   top contributing variables. It ranges from \code{top_n_lower} to \code{top_n_upper}.
#'   \item A line plot shows the reconstruction performance across \code{n}, 
#'   compared to the original module–trait correlation.
#' }
#'
#' This helps assess how few key variables are needed to approximate the full 
#' module-trait correlation.
#'
#' @return One or more ggplot2 plots are printed to the active graphics device.
#' If no correlations exceed the interest_threshold, a message is returned instead.
#'
#' @import ggplot2
#' @export

insight_correlations_reconstruct = function(
    module_trait_correlations,
    spca_data, 
    labels,
    col,
    cex = 20, 
    top_n_upper = 50, 
    top_n_lower = 10, 
    interest_threshold = 0.6,
    threshold = 0.8,
    return_information = TRUE, 
    ylab = "Correlation",
    xlab = "Top n variables",
    col_connective = 8){
  
  if(!is.numeric(interest_threshold)){
    stop("Invalid interest_threshold argument value.")
  }
  
  if(!is.list(spca_data)){
    stop("Invalid spca_data argument, must be a list.")
  }
  
  if(!is.numeric(threshold)){
    stop("Invalid threshold argument value.")
  }
  
  if(interest_threshold <= 0 | interest_threshold >= 1){
    warning("Threshold argument value should be between 0 and 1.")
  }
  
  if(threshold <= 0 | threshold >= 1){
    warning("Threshold argument value should be between 0 and 1.")
  }
  
  if(length(labels) != length(spca_data$PCs)){
    stop("Length of labels does not coincide with number of groups in spca_data.")
  }
  
  # top n uppper & lower oesnt excee moule sizes !!
  
  nsets = length(spca_data$module_matrices)
  
  # select positions of interesting modules & traits
  positions = matrix(NA, ncol = 2, nrow = 0)
  for(set in 1:nsets){
    positions = rbind(positions, which(abs(module_trait_correlations[[1]][[set]]) > interest_threshold, arr.ind = TRUE) |> as.data.frame() |> mutate(set = set))
  }
  positions$trait = colnames(module_trait_correlations[[1]][[1]])[positions$col]
  positions$module = rownames(module_trait_correlations[[1]][[1]])[positions$row]
  
  if(nrow(positions) == 0){return(paste("No interesting correlations, lower interest_threshold."))}
  
  for(i in 1:nrow(positions)){
    
    sets = positions$set[i]
    module_variable_names = colnames(spca_data$module_matrices[[1]][[positions$row[i]]])
    cors = numeric(0)
    
    # compute protein-trait unweighted correlations for all groups
    for(set in 1:nsets){
      module_variables = spca_data$module_matrices[[set]][[positions$row[i]]]
      given_trait = module_trait_correlations$trait_data[[set]][positions[i, 2]]
      cors = c(cors, cor(module_variables, given_trait, use = 'pairwise.complete.obs') |> as.vector())
    }
    
    cors = data.frame(
      variable_id = rep(module_variable_names, nsets),
      num_id = rep(1:length(module_variable_names), nsets),
      cors = cors
    )
    
    # compute activity and bind it to unweighted correlations 
    module_activity_table = compute_activity(
      spca_data = spca_data, 
      return_information = return_information, 
      threshold = threshold,
      labels = labels
    ) |>
      tidyr::pivot_longer(
        cols = starts_with("w_"),   # all columns that begin with "cor_"
        names_to = "label",           # new column for names
        values_to = "weights"             # new column for values
      )
    
    # compute weighted correlations
    cors = cors |>
      left_join(module_activity_table, by = join_by("variable_id"), relationship = "many-to-many") |>
      filter(.data$activity != "dead") |> # ignore dead
      arrange(.data$activity) |> #                           # sort by activity status
      mutate(variable_id = factor(.data$variable_id),
             w_cor = .data$weights * .data$cors) |>       # change variable id to factor
      mutate(variable_id = reorder(.data$variable_id, .data$num_id, mean)) 
    
    row_index = 1
    reconstruct_results = matrix(NA, nrow = length(top_n_lower:top_n_upper), ncol = nsets)
    
    for(top_n in top_n_lower:top_n_upper){
      
      # identify top n variables with highest effect on correlation, save their ids and group label
      top_n_variable_id = cors |>
        filter(.data$module == positions$module[i]) |>
        group_by(.data$label) |>
        arrange(desc(abs(.data$w_cor))) |>
        slice_head(n = top_n) |>
        ungroup() |>
        select(.data$label, .data$variable_id)
      
      # reconstructed correlation based on top n variables
      for(set in 1:nsets){
        # select variables and trait
        module_variables = spca_data$module_matrices[[set]][[positions$row[i]]]
        given_trait = module_trait_correlations$trait_data[[set]][positions[i, 2]]
        module_variable_names = colnames(module_variables)
        top_variable_ids = top_n_variable_id |> filter(grepl(labels[set], .data$label)) |> pull(.data$variable_id)
        
        # reconstruct eigengene and correlation, fix its sign
        reconstruct_SPCA = module_variables %*% (spca_data$vs[[set]][positions$row[i]][[1]] * (module_variable_names %in% top_variable_ids))
        reconstruct_cor = suppressWarnings(cor(reconstruct_SPCA, given_trait, use = "pairwise.complete.obs"))
        sign_change = spca_data$module_signs[[set]][positions$row[i]]
        reconstruct_results[row_index, set] = reconstruct_cor * sign_change
      }
      row_index = row_index + 1
    }
    
    reconstruct_results = reconstruct_results |> 
      data.frame() |> 
      mutate(top_n = top_n_lower:top_n_upper) 
    names(reconstruct_results)[1:nsets] = labels
    reconstruct_results = reconstruct_results |>
      tidyr::pivot_longer(
        cols = 1:nsets,
        names_to = "label",
        values_to = "cor"
      )
    
    
    # title
    if(sum(is.na(spca_data$module_colors)) == 0){
      title = paste(
        positions$module[i],
        'module correlated with',
        positions$trait[i],
        'using',
        top_n_lower,
        'to',
        top_n_upper,
        'variables'
      )
    } else {
      title = paste(
        'correlation with',
        positions$trait[i],
        'using',
        top_n_lower,
        'to',
        top_n_upper,
        'variables'
      )
    }
    
    # plot top_n results
    hline_data = data.frame(
      label = labels,
      yintercept = vapply(
        seq_len(nsets),
        \(x) module_trait_correlations[[1]][[x]][positions$row[i], positions$col[i]],
        numeric(1)
      )
    )
    
    print(
      ggplot(data = reconstruct_results,
             aes(x = .data$top_n, y = .data$cor, group = .data$label, color = .data$label)) +
        geom_line(alpha = 0.4) +
        geom_hline(
          data = hline_data,
          aes(yintercept = .data$yintercept, color = .data$label),
          linewidth = 2
        ) +
        geom_hline(yintercept = 0, color = col_connective) +
        scale_color_manual(labels = labels, values = col, name = "Group") +
        theme_classic(base_size = cex) +
        ggtitle(title) +
        labs(y = ylab, x = xlab)
    )
  }
}