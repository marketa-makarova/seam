#' Compute correlations between clinical traits and variables
#'
#' This function computes Pearson correlations between variables and 
#' external clinical trait data. It also computes weighted correlations based
#' on SPCA assigned weights.
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
#' @param trait_data A list containing trait data.frame(s) with 
#' clinical variables in columns and samples in rows. Rows must be named
#' in correspondence to sample names in expression data.
#' Length of the list corresponds to the number of groups.
#'
#' @return A list with of data frames with the following columns:
#' \describe{
#'   \item{variable_id}{Variable identifiers/names.}
#'   \item{module}{The assigned module/cluster if applicable.}
#'   \item{activity}{Activity status.}
#'   \item{cor}{Correlation with a given clinical trait.}
#'   \item{w_cor}{Weighted correlation with a given clinical trait.}
#' }
#'
#' @details
#' Both weighted and non-weighted Pearson's correlations between variables and
#' traits are computed. Weights are taken as loadings from SPCA implemented by
#' the \code{\link[elasticnet]{spca}}. install_github("kbroman/broman")For details on activity assignment, refer 
#' to \code{\link{compute_activity}}. 
#' 
#' @export
compute_variable_trait_correlations = function(
    trait_data, 
    spca_data,
    return_information = TRUE, 
    threshold = 0.8,
    labels
  ){
  
  if(!is.list(spca_data)){
    stop("Invalid spca_data argument, must be a list.")
  }
  
  if(!is.list(trait_data)){
    stop("Invalid trait_data argument, must be a list.")
  }
  
  if(length(spca_data$PCs) != length(trait_data)){
    stop("Incompatable number of groups in trait_data and spca_data.")
  }
  
  n_sets = length(spca_data$module_matrices)
  
  for(i in 1:n_sets){
    if(nrow(trait_data[[i]]) != nrow(spca_data$PCs[[i]])){
      stop("Number of variables in spca_data does not coincide with length of module_assignment.")
    }
  }
  
  variable_trait_cors_list = list()
  
  activity = compute_activity(
    spca_data, 
    return_information = T, 
    threshold = 0.8,
    labels = c("a", "b")
  )
  
  multiexpression = lapply(spca_data$module_matrices, function(m) do.call(cbind, m))
  
  for(set in 1:n_sets){
    traits = trait_data[[set]]
    
    for(trait in 1:ncol(traits)){
      trait_name = colnames(traits)[trait]
      variable_trait_cors = as.vector(cor(traits[, trait], multiexpression[[set]]))
      activity[[paste0(trait_name, "_cor")]] = variable_trait_cors
      activity[[paste0(trait_name, "_w_cor")]] = variable_trait_cors * activity[, 2 + set]
    }
    
    variable_trait_cors_list[[set]] = activity
    
  }
    
  return(variable_trait_cors_list)
  
}
