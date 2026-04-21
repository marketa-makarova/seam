#' Compute activity statuses of variables using SPCA
#'
#' This function evaluates variable activity based on sparse principal component 
#' analysis (SPCA) weights. It initially labels each variable's activity as 
#' either **active** in the case it has been assigned none-zero weights or **dead** 
#' in case it has been assigned weight equal to zero. 
#' It optionally performs a correlation-based 
#' refinement when the parameter \code{return_information = TRUE}, 
#' to detect so called **undead** variables -- those that appear inactive 
#' (**dead**) but are highly correlated with active ones.
#'
#' @param spca_data A list containing SPCA decomposition results from the 
#' function \code{spca_all_modules} or \code{spca_all_modules_grid}.
#' @param return_information Logical. If \code{FALSE}, returns only raw activity 
#' assignments. If \code{TRUE} performs additional correlation-based checks to return
#' dead variables highly correlated with actives one to **undead** activity status. 
#' Defaults to \code{TRUE}.
#' @param threshold Numeric. The correlation threshold used to determine 
#' whether a **dead** variable highly correlated with an active variable, 
#' suggesting it should be assigned the status **undead**. 
#' Default value is 0.8.
#' @param labels Character vector. Provides group labels. 
#' These labels are used for activity assignment.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{variable_id}{Variable identifiers/names.}
#'   \item{module}{The assigned module/cluster if applicable.}
#'   \item{weight columns}{SPCA weight of the variable. If several groups are present, 
#'   there will be several of these column with the suffix \code{_label}.}
#'   \item{activity}{Activity status.}
#' }
#'
#' @details
#' In case of no groups, activity labels are assigned as follows:
#' \itemize{
#'   \item **active** — variables assigned non-zero weights by SPCA.
#'   \item **dead** — variables assigned zero weights by SPCA.
#' }
#' 
#' In case of two groups, labelled \code{c("A", "B")}, activity labels are assigned as follows:
#' \itemize{
#'   \item **both_dead** — variables assigned zero weights by SPCA performed on the expression data of *both* groups.
#'   \item **both_active** — variables assigned non-zero weights by SPCA in *both* groups.
#'   \item **active_A** — variables assigned non-zero weights only in group \code{"A"} and zero weights in group \code{"B"}.
#'   \item **active_B** — variables assigned non-zero weights only in group \code{"B"} and zero weights in group \code{"A"}.
#' }
#'
#' If \code{return_information == TRUE}, the function additionally checks whether any
#' **both_dead** variables are highly correlated (Pearson's \eqn{r > threshold})
#' with variables that are active in at least one group.
#' Such variables are relabeled as follows:
#' \itemize{
#'   \item **undead_A** — correlated with active variables in group \code{"A"} only.
#'   \item **undead_B** — correlated with active variables in group \code{"B"} only.
#'   \item **both_undead** — correlated with active variables in *both* groups.
#' }
#' 
#' For more than two groups, the function applies analogous logic. Specifically:
#' \itemize{
#'   \item **active** — varaibles assigned non-zero weights by SPCA
#'   in *at least two* groups.
#'   \item **dead** — varibales assigned zero weights by SPCA
#'   in *all* groups.
#'   \item If a variable is active **only in one group**, it is labeled accordingly
#'   (e.g., \code{active_A}, \code{active_B}, \code{active_C}).
#'   \item If \code{return_information == TRUE}, dead variables may be relabeled as
#'   **undead** variants based on high correlation with active variables in one or more groups,
#'   following the same principles as in the two-group case.
#' }
#'
#' @export

compute_activity = function(
    spca_data, 
    return_information = TRUE, 
    threshold = 0.8,
    labels){
  
  if(!is.numeric(threshold)){
    stop("Invalid threshold argument value.")
  }
  
  if(!is.list(spca_data)){
    stop("Invalid spca_data argument, must be a list.")
  }
  
  if(threshold <= 0 | threshold >= 1){
    warning("Threshold argument value should be between 0 and 1.")
  }
  
  if(length(labels) != length(spca_data$PCs)){
    stop("Length of labels does not coincide with number of groups in spca_data.")
  }
  
  if(any(grepl("-", labels))){
    stop("Do not use - in labels.")
  }

  # extract information
  nsets = length(spca_data$vs)
  number_of_modules = length(spca_data$module_matrices[[1]])
  
  # initialize activity table
  activity_table = data.frame(
    variable_id = unlist(lapply(spca_data$module_matrices[[1]], colnames)),
    module = rep(names(spca_data$module_matrices[[1]]), lapply(spca_data$module_matrices[[1]], ncol))
  )
  
  # save weights into activity table and name the columns appropriately
  for(set in 1:nsets){
    activity_table = cbind(activity_table, unlist(spca_data$vs[[set]]))
  }
  names(activity_table)[(ncol(activity_table) - nsets + 1):ncol(activity_table)] = paste("w", labels, sep = "_")
  
  # compute activity status
  activity_table = activity_table |>
    mutate(
      weights = rowSums(across(starts_with("w_")) != 0)) |> # number of non-zero weights
    mutate(
      activity = case_when(
        .data$weights == 0 ~ "dead",
        .data$weights > 1 ~ "active",
        .data$weights == 1 ~ paste0( # specific only for one group
          "active_",
          sub("^w_", "", labels[max.col(across(starts_with("w_")) != 0, ties.method = "first")]) # which group
        )
      )
    ) |>
    select(-"weights") # remove the column used for computations
  
  
  if(return_information == FALSE){ return(activity_table) }
  
  if(return_information == TRUE | return_information == "both"){
    
    # initialize the new undead group specific columns so they can be updated when iterating through modules
    for(set in 1:nsets){
      new_column_name = paste0("undead_", labels[set])
      activity_table[[new_column_name]] = activity_table$activity
    }
    
    for(module in 1:number_of_modules){
      for(set in 1:nsets){ 
        
        module_cor = abs(cor(spca_data$module_matrices[[set]][[module]]))
        diag(module_cor) = 0
        
        # find which variables have high correlations with other variables (> threshold) 
        module_cor_threshold = module_cor |> 
          data.frame() |> 
          select(where(~ any(. > threshold))) |> 
          filter(if_any(everything(), ~ . > threshold)) 
        
        # if nothing passes threshold, move on
        if(length(names(module_cor_threshold)) == 0){
          next
        } else { # there are some undead, so lets identify them
          variable_activity_table_threshold = activity_table |> 
            filter(.data$variable_id %in% rownames(module_cor_threshold)) 
          
          # check which ones are classified as dead
          dead = variable_activity_table_threshold |> 
            filter(.data$activity == "dead") |> 
            pull(.data$variable_id)
          
          if(length(dead) == 0){
            next
          }
          
          # check which ones of these are correlated with some active variable (highly correlated but both dead are of no interest)
          full_rowsums = (module_cor_threshold > threshold) |> rowSums() 
          full_rowsums = data.frame(variable_id = names(full_rowsums), rowsums = full_rowsums)
          rowsums = (module_cor_threshold |> 
                       filter(names(module_cor_threshold) %in% dead) |> 
                       select(all_of(dead)) > threshold) |> 
            rowSums() # same as colSums due to symmetry
          both_dead_rowsums = data.frame(variable_id = names(rowsums), both_dead_rowsums = rowsums)
          
          # dont delete this additional skip, it breaks without it
          if(nrow(both_dead_rowsums) == 0){ 
            next
          }
          
          undead_ids = both_dead_rowsums |> 
            left_join(by = join_by('variable_id'), full_rowsums) |> 
            filter(both_dead_rowsums != rowsums) |> 
            pull(.data$variable_id)
          
          undead_label = paste("undead", labels[set], sep = "_")
          
          # modify the activity of those returning from the dead for interpretation to "undead"
          activity_table[activity_table$variable_id %in% undead_ids, paste0("undead_", labels[set])] = undead_label  
        }
      }
    }
    
    undead_cols <- names(activity_table)[grepl("^undead_", names(activity_table))]
    
    activity_table <- activity_table |>
      mutate(
        undead_counts = rowSums(as.data.frame(across(all_of(undead_cols), ~ grepl("undead", .x))))) |>
      mutate(
        undead_label  = case_when(
          .data$undead_counts == 1 ~ sub("^undead_", "", undead_cols[
            max.col(across(all_of(undead_cols), ~ grepl("undead", .x)), ties.method = "first")
          ]),
          TRUE ~ NA_character_
        )) |>
      mutate(activity = case_when(
        .data$undead_counts == 0 ~ activity,                        # keep original activity if nothing is undead
        .data$undead_counts > 1 ~ "undead",                         # undead status for undead in multiple groups
        .data$undead_counts == 1 ~ paste0("undead_", undead_label)  # one undead → give specific label
      )
      ) |>
      select(- all_of(undead_cols)) |>
      select(- c("undead_counts", "undead_label"))
    
    return(activity_table)
  }
}