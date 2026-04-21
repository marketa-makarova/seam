#' Perform SPCA across all groups and modules with different sparsity
#'
#' Applies sparse PCA to each module in each dataset, 
#' storing weights, explained variance, and the first principal component.
#'
#' @param expression_data A list containing expression matrices for each group.
#' @param module_assignment Character vector. Assigns
#'  module/cluster to each variable by color. If no clustering is present, this
#'  argument may be set to \code{NA}.
#' @param modules Character or character vector. Allows for module selection. Defaults to \code{"all"}.
#' @param lambdas A numeric vector of lambda values specifying the value for each module/cluster. Sparsity penalty parameter passed to \code{\link[elasticnet]{spca}}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{PCs}}{List of matrices. Each matrix contains the first SPCA component (PC1) for each module, per dataset.}
#'   \item{\code{var_exp}}{List of matrices. Each matrix holds the variance explained by PC1 in each module.}
#'   \item{\code{vs}}{List of lists. Each inner list contains loading vectors (weights) for PC1 in each module.}
#'   \item{\code{module_assignment}}{Original vector of module assignments.}
#'   \item{\code{module_matrices}}{List of lists. Each sublist contains the expression matrix for each module and group.}
#'   \item{\code{module_colors}}{Module color order used.}
#'   \item{\code{module_sizes}}{List of vectors of module sizes (number of variables per module).}
#'   \item{\code{module_signs}}{List of vectors of sign changes for eigengenes.}
#' }
#'
#' @details
#' Each module is extracted based on \code{module_assignment}. 
#' If missing values are present, they are imputed using KNN imputation from \pkg{impute}.
#'
#' The SPCA is performed using the \code{\link[elasticnet]{spca}} 
#' function from the \pkg{elasticnet} package, with a single component (PC1) extracted per module.
#' The sign of each PC is adjusted to match the average expression direction for interpretability.
#'
#' @export

spca_all_modules_grid <- function(
  expression_data, 
  module_assignment, 
  modules = "all",
  lambdas){
  
  if(any(lambdas < 0)) { stop("Select non-negative lambda.") }
  
  if(!is.list(expression_data)){
    stop("Invalid expression_data argument, must be a list.")
  }
  
  nsets = length(expression_data)
  
  if(any(!is.na(module_assignment))){
    for(i in 1:nsets){
      if(length(module_assignment) != ncol(expression_data[[i]])){
        stop("Number of variables in expression_data does not coincide with length of module_assignment.")
      }
    }
  }

  PCs = list() # principal components
  var_exp_list = list() # explained variance
  vs = list() # weights, eigenvectors, right singular vectors
  xs = list()


  if(sum(modules == "all") == 1){
    module_colors = unique(module_assignment)
  } else {
    module_colors = modules
  }
  
  module_matrices = list(length = nsets)
  module_signs = list(length = nsets)
  for(set in 1:nsets){
    module_matrices[[set]] = list(length = length(module_colors))
    module_signs[[set]] = numeric(length = length(unique(module_assignment)))
  }
  
  if(any(is.na(module_assignment))){
    for(set in 1:nsets){
      
      var_exp_matrix = matrix(NA, nrow = 1, ncol = 1) # save the variance explained here
      x = as.matrix(expression_data[[set]])
      PC1 = matrix(NA, nrow = nrow(x), ncol = 1)
      curr_vs = list()
      i = 1
      xs[[set]] = list()
      
      # if there are NAs, fill them using KNN
      if(sum(is.na(x)) != 0){
        xt = t(x)
        xtimp = impute.knn(xt, k = min(10, nrow(xt) - 1))
        x = t(xtimp$data)
      }
      
      # perform sparse pca
      spca_list = spca(x, K = 1, type = "predictor", sparse = "penalty", para = lambdas[1])
      
      # save variables
      y = x %*% spca_list$loadings
      PC1[, 1] = if(cor(y[, 1], apply(x, 1, sum) / ncol(x)) > 0){ y[, 1] } else { - y[, 1]}
      var_exp_matrix[1, ] = spca_list$pev
      curr_vs[[1]] = spca_list$loadings |> as.vector()
      
      module_matrices[[set]][[1]] = x
      PCs[[set]] = PC1
      var_exp_list[[set]] = var_exp_matrix
      vs[[set]] = curr_vs
      
    }
  } else {
    
    for(set in 1:nsets){
      var_exp_matrix = matrix(NA, nrow = length(module_colors), ncol = 1) # save the variance explained here
      fullx = expression_data[[set]]
      PC1 = matrix(NA, nrow = nrow(fullx), ncol = length(module_colors))
      curr_vs = list()
      i = 1
      xs[[set]] = list()
      
      # iterate over modules/colors
      for(color in module_colors){
        x = as.matrix(fullx[, as.character(module_assignment) == as.character(color)])
        
        # if there are NAs, fill them using KNN
        if(sum(is.na(x)) != 0){
          xt = t(x)
          xtimp = impute.knn(xt, k = min(10, nrow(xt) - 1))
          x = t(xtimp$data)
        }
        
        # perform sparse pca
        spca_list = spca(x, K = 1, type = "predictor", sparse = "penalty", para = lambdas[which(module_colors == color)])
        
        # save variables
        y = x %*% spca_list$loadings
        PC1[, which(module_colors == color)] = if(cor(y[, 1], apply(x, 1, sum) / ncol(x)) > 0){ y[, 1] } else { - y[, 1]}
        var_exp_matrix[which(module_colors == color), ] = spca_list$pev
        curr_vs[[which(module_colors == color)]] = spca_list$loadings |> as.vector()
        
        module_matrices[[set]][[which(module_colors == color)]] = x
        module_signs[[set]][which(module_colors == color)] = if(cor(y[, 1], apply(x, 1, sum) / ncol(x)) > 0){1}else{-1}
        
      }
      
      colnames(PC1) = module_colors
      names(module_matrices[[set]]) = module_colors
      
      PCs[[set]] = PC1
      var_exp_list[[set]] = var_exp_matrix
      vs[[set]] = curr_vs
      
    }
  }
  
  module_sizes = unlist(lapply(module_matrices[[1]], ncol))
  names(module_sizes) = module_colors
  
  return(list(
    PCs = PCs, 
    var_exp = var_exp_list, 
    vs = vs, 
    module_assignment = module_assignment, 
    module_matrices = module_matrices, 
    module_colors = module_colors,
    module_sizes = module_sizes,
    module_signs = module_signs,
    lambdas = lambdas))
}