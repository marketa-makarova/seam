#' Visualize variable–trait correlations as lollipop plot
#'
#' This function creates lollipop-style visualizations of 
#' correlations between variables and external clinical traits. 
#' If groups are present, it highlights 
#' differences between them using SPCA as feature selection method. Optionally, 
#' weighted correlations can be used to gain insight into how variables contribute
#' to sparse eigengene-trait correlations. If modules are present, 
#' each module is plotted separately.
#'
#' @param module_trait_correlations A list of module trait correlations
#' computed using the \code{compute_module_trait_correlations} function.
#' @param spca_data A list containing SPCA decomposition results from the 
#' function \code{spca_all_modules} or \code{spca_all_modules_grid}.
#' @param type Character. Specifies which type of plot to create, options are
#' \code{"weighted"}, \code{"classic"} or \code{"both"}. Defaults to \code{"both"}.
#' @param interest_threshold Numeric. Only module–trait correlations greater
#' than this value (in absolute terms) are considered interesting and plotted. 
#' Default is 0.6.
#' @param cex Numeric. Font size for plot elements. 
#' Passed to \code{theme_classic()} argument. 
#' Defaults to 20.
#' @param cex_point Numeric. Size of points visualising correlation. 
#' Defaults to 20.
#' @param labels Character vector. Provides group labels. 
#' These labels are used for activity assignment.
#' @param n Integer. The number of top variables 
#' (by absolute weighted correlation) to label on the y-axis. 
#' Defaults to 5.
#' @param xlab Character. X-axis label.
#' @param ylab Character. Y-axis label.
#' @param col Color or vector of colors if groups are present. 
#' Used for plotting points. 
#' @param col_connective Color of connective lines between correlation values. 
#' Defaults to 8 (grey).
#' @param col_active Color for variables that are active in both groups. 
#' Default is 3 (green).
#' @param y_dodge Integer. Used to spread out y-axis labels to avoid label overlap. 
#' Defaults to 2.
#'
#' @details
#' The behavior of the plot depends on the value of the parameter \code{type}:
#' \itemize{
#'   \item \code{"classic"} — plots standard variable–trait correlations as a lollipop plot.
#'   Correlations are computed directly between each variable and the trait.
#' 
#'   \item \code{"weighted"} — plots weighted variable–trait correlations, defined as the
#'   variable–trait correlation multiplied by the SPCA weight assigned to the variable.
#'   This emphasizes variables that are both strongly associated with the trait and
#'   strongly weighted by SPCA.
#' 
#'   \item \code{"both"} — produces a two-panel plot showing \code{"classic"} (unweighted)
#'   and \code{"weighted"} correlations side by side for direct comparison.
#' }
#'
#' If multiple groups are present, variable-trait correlations are computed separately for each group.
#' Group-specific correlations are displayed within the same plot and are differentiated
#' by color to facilitate comparison of group-wise association patterns.
#'
#' The parameter \code{n} controls how many variables are labeled on the y-axis.
#' Specifically, only the top \code{n} variables with the strongest absolute (weighted or
#' unweighted, depending on \code{type}) correlations are labeled.
#' If multiple labeled variables would overlap on the y-axis, the \code{y_dodge} parameter
#' is used to vertically offset labels to improve readability.
#'
#' This visualization facilitates identification of biologically relevant drivers of
#' module–trait associations. In particular, weighted correlations can reveal key driver
#' variables that differentiate between groups.
#' 
#' @export
insight_correlations_lollipop = function(
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
    col,
    col_connective = 8,
    col_active = 3,
    y_dodge = 2){
  
  if(!(type %in% c("classic", "weighted", "both"))){
    stop("Type can be either classic, weighted, or both.")
  }
  
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
        variable_id = factor(.data$variable_id)
      ) |>
      mutate(
        variable_id = reorder(.data$variable_id, .data$num, mean)
      ) |>
      filter(.data$activity != "dead")
    
    new_cols = setNames(
      lapply(labels, function(lbl) {
        cors[[paste0("w_", lbl)]] * cors[[paste0("cor_", lbl)]]
      }),
      paste0("w_cor_", labels)
    )
    
    cors = bind_cols(cors, as.data.frame(new_cols))
    
    # find top correlations
    top_n_cors = cors |>
      select(all_of(paste0("w_cor_", labels))) |>
      unlist(use.names = FALSE) |> abs() |> sort() |> tail(n)
    
    top_n_labels = cors |>
      rowwise() |>
      mutate(top_n_id = if_else(
        any(abs(c_across(all_of(paste0("w_cor_", labels)))) %in% top_n_cors),
        as.character(.data$variable_id),
        "")) |>
      ungroup() |>
      pull(.data$top_n_id)
    
    # title
    if (sum(is.na(spca_data$module_colors)) == 0) {
      
      cor_strings <- vapply(
        seq_len(nsets),
        function(set) {
          paste0(
            labels[set],
            " = ",
            round(
              module_trait_correlations[[1]][[set]][positions[i, 1], positions[i, 2]],
              2
            )
          )
        },
        character(1)
      )
      
      classic_title <- paste0(
        positions$module[i],
        " module correlated with ",
        positions$trait[i],
        "\nSPCA correlations: ",
        paste(cor_strings, collapse = ", ")
      )
    } else {
      
      cor_strings <- vapply(
        seq_len(nsets),
        function(set) {
          paste0(
            labels[set],
            " = ",
            round(
              module_trait_correlations[[1]][[set]][positions[i, 1], positions[i, 2]],
              2
            )
          )
        },
        character(1)
      )
      
      classic_title <- paste0(
        "correlation with ",
        positions$trait[i],
        "\nSPCA correlations: ",
        paste(cor_strings, collapse = ", ")
      )
    }
    
    
    ggcors = rbind(
      cors |> select(.data$variable_id, .data$activity, .data$module, all_of(paste0("cor_", labels))) |> mutate(weighted = "unweighted"),
      cors |> select(.data$variable_id, .data$activity, .data$module, all_of(paste0("w_cor_", labels))) |> mutate(weighted = "weighted") |> rename_with(~ sub("^w_cor_", "cor_", .x), starts_with("w_cor_"))
    ) |>
      tidyr::pivot_longer(
        cols = starts_with("cor_"),   # all columns that begin with "cor_"
        names_to = "label",           # new column for names
        values_to = "cor"             # new column for values
      ) |>
      mutate(label = gsub("cor_", "", .data$label))  # remove "cor_" prefix for cleaner labels
    
    # plots
    if(type == "weighted"){
      print(
        ggplot(data = ggcors |> filter(.data$weighted == "weighted")) +
          
          # connective lines 
          geom_segment(aes(x = .data$variable_id, y = cor, yend = 0, color = .data$label), na.rm = TRUE) +
          geom_segment(
            data = ggcors |> filter(.data$weighted == "weighted" & .data$activity == "active"),
            mapping = aes(x = .data$variable_id, y = cor, yend = 0),
            color = col_active,
            na.rm = TRUE
          ) +
          
          # all points coloured by activity
          geom_point(aes(x = .data$variable_id, y = cor, color = .data$label), size = cex_point, na.rm = TRUE) +
          
          coord_flip() +
          theme_classic(base_size = cex) +
          scale_x_discrete(labels = top_n_labels, guide = guide_axis(n.dodge = y_dodge)) +
          labs(x = xlab, y = ylab) +
          ggtitle(classic_title) +
          geom_hline(yintercept = 0, color = col_connective, linewidth = 1) +
          guides(size = "none") +
          scale_color_discrete(name = "Group")
        
      )
    }
    
    if(type == "classic"){
      print(
        ggplot(data = ggcors |> filter(.data$weighted == "unweighted")) +
          
          # connective lines 
          geom_segment(aes(x = .data$variable_id, y = .data$cor, yend = 0, color = .data$label), na.rm = TRUE, color = "grey") +
          
          # all points coloured by activity
          geom_point(aes(x = .data$variable_id, y = .data$cor, color = .data$label), size = cex_point, na.rm = TRUE) +
          
          coord_flip() +
          theme_classic(base_size = cex) +
          scale_x_discrete(labels = top_n_labels, guide = guide_axis(n.dodge = y_dodge)) +
          labs(x = xlab, y = ylab) +
          ggtitle(classic_title) +
          geom_hline(yintercept = 0, color = col_connective, linewidth = 1) +
          guides(size = "none") +
          scale_color_discrete(name = "Group")
        
      )
    }
    
    if(type == "both"){
      print(
        ggplot(data = ggcors) +
          
          # connective lines for weighted
          geom_segment(
            data = ggcors |> filter(.data$weighted == "weighted"),
            mapping = aes(x = .data$variable_id, y = cor, yend = 0, color = .data$label), na.rm = TRUE) +
          geom_segment(
            data = ggcors |> filter(.data$weighted == "weighted" & .data$activity == "active"),
            mapping = aes(x = .data$variable_id, y = cor, yend = 0),
            color = col_active,
            na.rm = TRUE
          ) +
          
          # connective lines for classic
          geom_segment(
            data = ggcors |> filter(.data$weighted == "unweighted"),
            mapping = aes(x = .data$variable_id, y = cor, yend = 0, color = .data$label), na.rm = TRUE, color = "grey") +
          
          # all points coloured by activity
          geom_point(data = ggcors, mapping = aes(x = .data$variable_id, y = cor, color = .data$label), size = cex_point, na.rm = TRUE) +
          
          coord_flip() +
          theme_classic(base_size = cex) +
          scale_x_discrete(labels = top_n_labels, guide = guide_axis(n.dodge = y_dodge)) +
          labs(x = xlab, y = ylab) +
          ggtitle(classic_title) +
          geom_hline(yintercept = 0, color = col_connective, linewidth = 1) +
          guides(size = "none") +
          scale_color_discrete(name = "Group") +
          facet_wrap(~ .data$weighted, scales = "free")
      )
    }
    
  }
  
}