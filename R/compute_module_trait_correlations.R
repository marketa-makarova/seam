#' Compute correlations between clinical traits and sparse or classic eigengenes
#'
#' This function computes Pearson correlations between eigengenes and 
#' external clinical trait data. If the data is clustered into modules, it
#' computes these correlations for each module. If several groups are present, 
#' calculations are performed for each group separately.
#'
#' @param spca_data A list containing SPCA decomposition results from the 
#' function \code{spca_all_modules} or \code{spca_all_modules_grid}.
#' @param trait_data A list containing trait data.frame(s) with 
#' clinical variables in columns and samples in rows. Rows must be named
#' in correspondence to sample names in expression data.
#' Length of the list corresponds to the number of groups.
#' @param module_assignment Character vector. Assigns
#'  module/cluster to each variable by color. If no clustering is present, this
#'  argument may be set to \code{NA}.
#' @param mode Character. Specifies which method to use: \code{"spca"},
#'  \code{"classic"}, or \code{"both"} for comparison.
#'
#' @return A list with the following elements, depending on the mode:
#' \describe{
#'   \item{\code{mtc_spca}}{A list of module–trait correlation matrices.}
#'   \item{\code{mtc_classic}}{A list of module–trait correlation matrices.}
#'   \item{\code{trait_data}}{The input trait data.}
#'   \item{\code{module_sizes}}{A vector or list indicating the number
#'   of variables in modules and groups if they are present.}
#'   \item{\code{reduced_module_sizes}}{A vector or list indicating the number
#'   of variables in modules selected by SPCA.}
#' }
#'
#' @details
#' In \code{"classic"} mode, the function computes correlations between
#' classic eigengenes and clinical trait data using Pearson correlation.
#' 
#' In \code{"spca"} mode, it uses SPCA-derived principal components 
#' instead of classic eigengenes.
#' 
#' In \code{"both"} mode, the function returns results for both methods.
#'
#' @export
compute_module_trait_correlations = function(
    trait_data, 
    spca_data,
    module_assignment, 
    mode = "spca"){
  
  if(!(mode %in% c("classic", "spca", "both"))){
    stop("Type can be either classic, spca, or both.")
  }
  
  if(!is.list(spca_data)){
    stop("Invalid spca_data argument, must be a list.")
  }
  
  if(!is.list(trait_data)){
    stop("Invalid trait_data argument, must be a list.")
  }
  
  if(length(spca_data$PCs) != length(trait_data)){
    stop("Incompatable number of groups in trait_data and spca_data.")
  }
  
  if(sum(unlist(lapply(spca_data$vs[[1]], length))) != length(module_assignment)){
    stop("Number of variables in spca_data does not coincide with length of module_assignment.")
  }
  
  n_sets = length(spca_data$module_matrices)
  
  for(i in 1:n_sets){
    if(nrow(trait_data[[i]]) != nrow(spca_data$PCs[[i]])){
      stop("Number of variables in spca_data does not coincide with length of module_assignment.")
    }
  }

  module_colors = spca_data$module_colors
  mtc_classic = mtc_spca = list() # mtc stands for module-trait-correlation   
  module_sizes = lapply(spca_data$vs[[1]], length)
  reduced_module_sizes = lapply(spca_data$vs, (\(k) lapply(k, (\(x) sum(x != 0)))))

  
  if(mode == "classic" | mode == "CLASSIC"){
    
    for(set in 1:n_sets){
      mtc_classic[[set]] = cor(compute_eigengenes(
        module_matrices = spca_data$module_matrices, 
        module_assignment = module_assignment, 
        module_colors = module_colors)[[set]], trait_data[[set]], use = "p")
      rownames(mtc_classic[[set]]) = module_colors
    }
    return(list(mtc_classic, trait_data = trait_data, module_sizes = module_sizes, reduced_module_sizes = reduced_module_sizes))
  }
  
  if(mode == "SPCA" | mode == "spca" ){
    
    for(set in 1:n_sets){
      mtc_spca[[set]] = cor(spca_data$PCs[[set]], trait_data[[set]], use = "p")
      rownames(mtc_spca[[set]]) = module_colors
    }
    return(list(mtc_spca, trait_data = trait_data, module_sizes = module_sizes, reduced_module_sizes = reduced_module_sizes))
  }
  
  if(mode == "both"){
    
    for(set in 1:n_sets){
      mtc_classic[[set]] = cor(compute_eigengenes(
        module_matrices = spca_data$module_matrices, 
        module_assignment = module_assignment, 
        module_colors = module_colors)[[set]], trait_data[[set]], use = "p")
      mtc_spca[[set]] = cor(spca_data$PCs[[set]], trait_data[[set]], use = "p")
      rownames(mtc_classic[[set]]) = module_colors
      rownames(mtc_spca[[set]]) = module_colors
    }
    return(list(mtc_classic = mtc_classic, mtc_spca = mtc_spca, trait_data = trait_data, module_sizes = module_sizes, reduced_module_sizes = reduced_module_sizes))
  }
}
