#' Plot SPCA-based activity across modules as a bubble chart
#'
#' This function visualizes the distribution of variable activity statuses 
#' using a bubble chart. 
#' Each bubble represents the count of variable in a given activity 
#' category. The bubble size reflects 
#' the count, and colors correspond to module assignments, if present.
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
#' @param scale_range Numeric vector of length 2 defining the minimum and maximum sizes of the bubbles. 
#' Defaults to \code{c(5, 30)}.
#' @param number_color Character. Color of the text labels inside bubbles. Default is \code{"white"}.
#' @param bubble_colors Optional vector or palette of colors to use for bubbles by module. If \code{NA}, uses \code{module_colors}.
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
#' @return A \code{ggplot2} object is printed showing a bubble chart of 
#' variable activity distributions across modules.
#'
#' @export

plot_activity_bubble = function(
    spca_data,
    labels,
    modules = "all",
    threshold = 0.8,
    return_information = TRUE,
    cex = 20,
    scale_range = c(5, 30),
    number_color = "white",
    bubble_colors = NA){
  
  if(!is.list(spca_data)){
    stop("Invalid spca_data argument, must be a list.")
  }
  
  if(length(scale_range) != 2){
    warning("Argument scale_range should be a vector of length 2.")
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

  if(sum(modules == "all") == 1){
    module_colors = spca_data$module_colors
  } else {
    module_colors = modules
  }
  
  activity_table = compute_activity(spca_data = spca_data, labels = labels, threshold = threshold, return_information = return_information)
  
  all_combos = expand.grid(
    module = module_colors,
    activity = unique(activity_table$activity),
    stringsAsFactors = FALSE
  )
  
  data = activity_table |>
    mutate(activity = factor(.data$activity)) |>
    group_by(.data$module, .data$activity) |>
    summarize(counts = n()) |>
    right_join(all_combos, by = c("module", "activity")) |>
    mutate(counts = ifelse(is.na(.data$counts), 0, .data$counts))
  
  if(any(is.na(spca_data$module_colors))) {
    x_tick_labels = element_blank()
    x_axis_label = element_blank()
  } else {
    x_tick_labels = element_text(size = cex)
    x_axis_label = element_text(size = cex)
  }

  print(
    ggplot(data,
         aes(x = .data$module, 
             y = .data$activity,
             colour = .data$module,
             size = .data$counts)) +
    geom_point() +
    geom_text(aes(label = .data$counts), 
              colour = number_color, 
              size = 4) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = cex) + 
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.ticks = element_blank()) +
    scale_size(range = scale_range) + 
    scale_color_manual(values = if(is.na(bubble_colors)){sort(module_colors)}else{bubble_colors}) +
    theme(axis.title.x = x_tick_labels,
          axis.text.x = x_axis_label)
  )

}