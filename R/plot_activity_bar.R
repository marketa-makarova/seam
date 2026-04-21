#' Plot variable activity statuses as a stacked bar chart
#'
#' This function visualizes the distribution of variable activity statuses 
#' using stacked bar charts. 
#' Activity statuses include three main states **active**, **dead**, and **undead**.
#' Function is able to incorporate both groups and modules.
#' 
#' @param spca_data A list containing SPCA decomposition results from the 
#' function \code{spca_all_modules} or \code{spca_all_modules_grid}.
#' @param labels Character vector. Provides group labels. 
#' These labels are used for activity assignment.
#' @param modules Character or character vector. Allows for module selection. Defaults to \code{"all"}.
#' @param return_information Logical. If \code{FALSE}, returns only raw activity 
#' assignments. If \code{TRUE} performs additional correlation-based checks to 
#' detect **undead** variables and updates the activity status accordingly. 
#' Defaults to \code{TRUE}.
#' @param threshold Numeric. The correlation threshold used to determine 
#' whether a **dead** variable highly correlated with an active variable, 
#' suggesting it should be assigned the status **undead**. 
#' Default value is 0.8.
#' @param cex Numeric. Font size for plot elements. 
#' Passed to \code{theme_classic()} argument. 
#' Defaults to 20.
#' @param col_both_active Color for variables that are active in both groups. 
#' Default is 3 (green).
#' @param col_active Color or vector of colors if groups are present. 
#' Used for plotting points of **active** variables. 
#' @param col_dead Color for **dead** variables. 
#' Defaults to 8 (grey).
#' @param col_both_undead Color for **undead** variables. 
#' Default is \code{"lightgrey"}.
#' @param col_undead Color or vector of colors if groups are present. 
#' Used for plotting points of **undead** variables.
#' @param xlab Character. X-axis label.
#' @param ylab Character. Y-axis label.
#'
#' @details
#' The function determines activity states based on whether the SPCA weights
#' for a variables are zero or non-zero.
#' Variables labeled **both_dead** but showing strong intra-module correlations 
#' with active variables are reclassified as **undead** 
#' for interpretability if option \code{return_information = TRUE}.
#' 
#' For details about activity statuses refer to \code{\link{compute_activity}}.
#'
#' @return A (\code{ggplot}) showing the activity distribution per module. 
#'
#' @export

plot_activity_bar = function(
    spca_data, 
    return_information = TRUE, 
    modules = "all",
    threshold = 0.6, 
    labels,
    cex = 20, 
    col_both_active = 3,
    col_active,
    col_dead = 8,
    col_both_undead = "lightgrey",
    col_undead,
    xlab = "Module", 
    ylab = "Count"){
  
  if(!is.list(spca_data)){
    stop("Invalid spca_data argument, must be a list.")
  }
  
  if(!is.numeric(threshold)){
    stop("Invalid threshold argument value.")
  }
  
  if(threshold <= 0 | threshold >= 1){
    warning("Threshold argument value should be between 0 and 1.")
  }
  
  if(length(labels) != length(spca_data$PCs)){
    stop("Length of labels does not coincide with number of groups in spca_data.")
  }
  
  if(any(is.na(spca_data$module_colors))) {
    ticks_theme = element_blank()
    axis_labels = element_blank()
    xlab = element_blank()
  } else {
    ticks_theme = element_line()
    axis_labels = element_text(size = cex)
    xlab = xlab
  }
  
  current_module_colors = if(modules[1] == "all"){spca_data$module_colors}else{modules} 
  current_module_sizes = spca_data$module_sizes[match(current_module_colors, spca_data$module_colors)]
  
  if(return_information == FALSE | return_information == "both"){
    
    activity_table = compute_activity(
      spca_data = spca_data, 
      return_information = FALSE,
      threshold = threshold,
      labels = labels
    ) |>
      filter(.data$module %in% current_module_colors) |> 
      mutate(module = factor(.data$module, levels = current_module_colors)) # force correct order in barplot
    
    p = ggplot(activity_table, aes(fill = .data$activity, x = .data$module)) + 
      geom_bar() +
      scale_fill_manual(values = c(col_both_active, col_active, col_dead, "black")) + 
      scale_x_discrete(labels = paste(current_module_colors, ' (', current_module_sizes, ')', sep = '')) + 
      theme_classic(base_size = cex) + 
      labs(x = xlab, y = ylab) + 
      theme(axis.ticks = ticks_theme,
            axis.text.x = axis_labels)
    
    print(p)
  }
  
  if(return_information == TRUE | return_information == "both"){
    
    activity_table = compute_activity(
      spca_data = spca_data, 
      return_information = TRUE,
      threshold = threshold,
      labels = labels
    ) |>
      filter(.data$module %in% current_module_colors) |> 
      mutate(module = factor(.data$module, levels = current_module_colors)) # force correct order in barplot
    
    p = ggplot(activity_table, aes(fill = .data$activity, x = .data$module)) + 
      geom_bar() +
      scale_fill_manual(values = c(col_both_active, col_active, col_dead, col_both_undead, col_undead)) +
      scale_x_discrete(labels = paste(current_module_colors, ' (', current_module_sizes, ')', sep = '')) + 
      theme_classic(base_size = cex) + 
      labs(x = xlab, y = ylab) + 
      theme(axis.ticks = ticks_theme,
            axis.text.x = axis_labels)
    
    print(p)
  }
}