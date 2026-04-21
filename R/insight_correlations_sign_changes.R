#' Visualize sign changes in variable-trait correlations due to SPCA weights
#'
#' This function explores how weights between variables differ in signs between 
#' classic eigengenes and SPCA-derived eigengenes. 
#' It visualizes the distribution 
#' of raw and weighted correlations and 
#' reports the proportion of variables whose weighted correlation signs 
#' have changed due to sign change in weights. If modules/clusters are present, 
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
#' @param interest_threshold Numeric. Only module–trait correlations greater
#' than this value (in absolute terms) are considered interesting and plotted. 
#' Default is 0.6.
#' @param cex Numeric. Font size for plot elements. 
#' Passed to \code{theme_classic()} argument. 
#' Defaults to 20.
#' @param xlab Character. X-axis label.
#' @param ylab Character. Y-axis label.
#'
#' @return One or more ggplot2 plots are printed to the active graphics device.
#' If no correlations exceed the interest_threshold, a message is returned instead.
#'
#' @export

insight_correlations_sign_changes = function(
    module_trait_correlations, 
    spca_data, 
    labels,
    col,
    interest_threshold = 0.6, 
    cex = 20, 
    xlab = "Correlations",
    ylab = ""){
  
  if(!is.list(spca_data)){
    stop("Invalid spca_data argument, must be a list.")
  }
  
  if(!is.numeric(interest_threshold)){
    stop("Invalid interest_threshold argument value.")
  }
  
  if(interest_threshold <= 0 | interest_threshold >= 1){
    warning("Threshold argument value should be between 0 and 1.")
  }
  
  if(length(labels) != length(spca_data$PCs)){
    stop("Length of labels does not coincide with number of groups in spca_data.")
  }
  
  nsets = length(spca_data$vs)
  
  # select positions of interesting modules & traits
  positions = matrix(NA, ncol = 2, nrow = 0)
  for(set in 1:nsets){
    positions = rbind(positions, which(abs(module_trait_correlations[[1]][[set]]) > interest_threshold, arr.ind = TRUE) |> as.data.frame() |> mutate(set = set))
  }
  positions$trait = colnames(module_trait_correlations[[1]][[1]])[positions$col]
  positions$module = rownames(module_trait_correlations[[1]][[1]])[positions$row]
  
  # error
  if(nrow(positions) == 0){return(paste("No interesting correlations, lower interest_threshold."))}
  
  # compute correlations 
  for(i in 1:nrow(positions)){
    set = positions$set[i]
    # create correlations
    cors = data.frame(variable_id = colnames(spca_data$module_matrices[[set]][[positions[i, 1]]]))
    for(set in 1:nsets){
      new_column = paste("cor", labels[set], sep = "_")
      cors = cors |>
        mutate(
          !!new_column := cor(
            spca_data$module_matrices[[set]][[positions[i, 1]]], 
            module_trait_correlations$trait_data[[set]][positions[i, 2]],
            use = "pairwise.complete.obs"
          ) |> as.vector()
        )
    }
    
    # create table of correlations, weights, correlations * weights
    activity_table = compute_activity(spca_data = spca_data, labels = labels, return_information = FALSE)
    cors =  cors |>
      left_join(activity_table, by = join_by("variable_id")) |>
      arrange(.data$activity) |>
      mutate(
        num = 1:nrow(cors),                      
        variable_id = factor(.data$variable_id),
        variable_id = reorder(.data$variable_id, .data$num, mean)) |>
      filter(!.data$activity %in% c("dead", "undead"))
    
    new_cols = setNames(
      lapply(labels, function(lbl) {
        cors[[paste0("w_", lbl)]] * cors[[paste0("cor_", lbl)]]
      }),
      paste0("w_cor_", labels)
    )
    
    cors = bind_cols(cors, as.data.frame(new_cols))
    
    # % of variables switching signs of correlation
    sign_changes = cors |>
      summarize(
        across(
          starts_with("w_"),
          ~ mean(sign(.x) == -1) * 100,
          .names = "percent_neg_{.col}"
        )
      ) |> select(1:length(labels)) # select only active, how many % are -1
    
    # for plotting
    cors_density = cors |> 
      tidyr::pivot_longer(
        cols = matches("^(cor_|w_cor_)"),
        names_to = c(".value", "label"),
        names_pattern = "(w?_?cor)_(.*)"
      )    
    
    # title
    if (sum(is.na(spca_data$module_colors)) == 0) {
      title <- paste(
        positions$module[i],
        " module correlated with ",
        positions$trait[i],
        "\n",
        paste(
          paste0(labels, ": ", round(sign_changes), " % changed signs"),
          collapse = "\n"
        ),
        sep = ""
      )
    } else {
      title <- paste(
        "correlation with ",
        positions$trait[i],
        "\n",
        paste(
          paste0(labels, ": ", round(sign_changes), " % changed signs"),
          collapse = "\n"
        ),
        sep = ""
      )
    }
    
    
    # density plot
    print(ggplot(data = cors_density, aes(x = cor, group = .data$label, fill = .data$label)) +
            scale_fill_manual(values = col, name = "Group") +
            geom_density(alpha = 0.4) +
            theme_classic(base_size = cex) + 
            labs(x = xlab, y = ylab) + 
            ggtitle(title))
  }
}