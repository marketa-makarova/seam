#' Perform SPCA across groups and modules
#'
#' Fits SPCA to the original data and optionally to bootstrap resamples.
#' Input expression matrices must contain no missing values.
#'
#' @param expression_data List of numeric expression matrices (samples in rows).
#' @param module_assignment Vector of modules assigned to variables. One module per variable.
#' @param lambda Non-negative sparsity penalty passed to the `para` argument of [elasticnet::spca()].
#' @param modules Module labels to plot, or `"all"` (default).
#' @param bootstrap Logical. Whether to compute bootstrap replicates.
#' @param n Number of bootstrap replicates.
#'
#' @return List with original PCs, var_exp, vs, module_matrices,
#' module_colors, module_sizes, module_assignment, module_signs, and
#' optional bootstrap list (PCs, var_exp, vs, module_signs, samples).
#' 
#' @details Bootstrap resampling may be used to assess the stability of 
#' variable selection under repeated sampling. 
#' 
#' #' @examples
#' data("seam_demo", package = "seam")
#' fit = spca_all_modules(
#'  seam_demo$expression_data,
#'  seam_demo$module_assignment,
#'  lambda = 0.05, 
#'  bootstrap = FALSE
#' )
#' fit$var_exp
#' 
#' @export
spca_all_modules = function(expression_data, module_assignment, lambda,
                             modules = "all", bootstrap = TRUE, n = 100){

  if(!is.list(expression_data) || length(expression_data) == 0 || any(vapply(expression_data, anyNA, logical(1)))){
    stop("expression_data must be a non-empty list of matrices without missing values.")
  }
  if(length(lambda) != 1 || is.na(lambda) || !is.finite(lambda) || lambda < 0){
    stop("Select a non-negative lambda.")
  }
  if(bootstrap && (length(n) != 1 || is.na(n) || n < 1 || n != as.integer(n))){
    stop("n must be a positive integer.")
  }
  expression_data = lapply(expression_data, as.matrix)
  n_sets = length(expression_data)
  for(set in seq_len(n_sets)){
    x = expression_data[[set]]
    if(!is.numeric(x) || anyNA(x) || any(!is.finite(x))){
      stop("Expression matrices must contain finite numeric values without missing data.")
    }
    if(nrow(x) < 2 || ncol(x) < 2){
      stop("Each expression matrix must have at least two rows and columns.")
    }
  }
  no_modules = length(module_assignment) == 1 && is.na(module_assignment)
  if(no_modules){
    module_colors = "all"
  } else {
    if(anyNA(module_assignment) || any(vapply(expression_data, ncol, integer(1)) != length(module_assignment))){
      stop("module_assignment must contain one non-missing label per variable.")
    }
    module_colors = if(length(modules) == 1 && modules == "all"){
      unique(as.character(module_assignment))
    } else {
      as.character(modules)
    }
    if(anyDuplicated(module_colors) || !all(module_colors %in% module_assignment)){
      stop("Selected modules must be unique and present in module_assignment.")
    }
  }

  # prepare variables
  n_modules = length(module_colors)
  PCs = vector("list", n_sets)
  var_exp = vector("list", n_sets)
  vs = vector("list", n_sets)
  module_matrices = vector("list", n_sets)
  module_signs = vector("list", n_sets)
  boot = NULL
  if(bootstrap){
    boot = list(
      PCs = vector("list", n_sets),
      var_exp = vector("list", n_sets),
      vs = vector("list", n_sets),
      module_signs = vector("list", n_sets),
      samples = vector("list", n_sets)
    )
  }

  for(set in seq_len(n_sets)){
    full_x = expression_data[[set]]
    n_samples = nrow(full_x)
    module_matrices[[set]] = vector("list", n_modules)
    vs[[set]] = vector("list", n_modules)
    names(module_matrices[[set]]) = names(vs[[set]]) = module_colors
    PCs[[set]] = matrix(NA_real_, n_samples, n_modules, dimnames = list(rownames(full_x), module_colors))
    var_exp[[set]] = matrix(NA_real_, n_modules, 1, dimnames = list(module_colors, NULL))
    module_signs[[set]] = setNames(numeric(n_modules), module_colors)

    for(i in seq_len(n_modules)){
      x = if(no_modules) full_x else full_x[, as.character(module_assignment) == module_colors[i], drop = FALSE]
      module_matrices[[set]][[i]] = x

      fit = elasticnet::spca(x, K = 1, type = "predictor", sparse = "penalty", para = lambda)
      raw_weights = as.vector(fit$loadings)
      raw_scores = as.vector(x %*% raw_weights)
      direction = suppressWarnings(cor(raw_scores, rowMeans(x)))
      sign = if(is.finite(direction) && direction > 0) 1 else -1

      weights = sign * raw_weights
      names(weights) = colnames(x)
      vs[[set]][[i]] = weights
      PCs[[set]][, i] = sign * raw_scores
      var_exp[[set]][i, 1] = fit$pev
      module_signs[[set]][i] = sign
    }

    if(bootstrap){
      boot$PCs[[set]] = array(NA_real_, c(n_samples, n_modules, n), dimnames = list(NULL, module_colors, NULL))
      boot$var_exp[[set]] = matrix(NA_real_, n, n_modules, dimnames = list(NULL, module_colors))
      boot$module_signs[[set]] = matrix(NA_real_, n, n_modules, dimnames = list(NULL, module_colors))
      boot$samples[[set]] = matrix(NA_integer_, n, n_samples)
      boot$vs[[set]] = vector("list", n_modules)
      names(boot$vs[[set]]) = module_colors

      for(i in seq_len(n_modules)){
        x = module_matrices[[set]][[i]]
        boot$vs[[set]][[i]] = matrix(NA_real_, n, ncol(x), dimnames = list(NULL, colnames(x)))
      }

      for(b in seq_len(n)){
        sampled_rows = sample.int(n_samples, n_samples, replace = TRUE)
        boot$samples[[set]][b, ] = sampled_rows

        for(i in seq_len(n_modules)){
          x = module_matrices[[set]][[i]][sampled_rows, , drop = FALSE]
          fit = elasticnet::spca(x, K = 1, type = "predictor", sparse = "penalty", para = lambda)
          raw_weights = as.vector(fit$loadings)

          # align each bootstrap component to the original loading vector
          sign = if(sum(raw_weights * vs[[set]][[i]]) >= 0) 1 else -1
          boot$vs[[set]][[i]][b, ] = sign * raw_weights
          boot$PCs[[set]][, i, b] = sign * as.vector(x %*% raw_weights)
          boot$var_exp[[set]][b, i] = fit$pev
          boot$module_signs[[set]][b, i] = sign
        }
      }
    }
  }

  module_sizes = vapply(module_matrices[[1]], ncol, integer(1))

  return(list(
    PCs = PCs,
    var_exp = var_exp,
    vs = vs,
    module_assignment = module_assignment,
    module_matrices = module_matrices,
    module_colors = module_colors,
    module_sizes = module_sizes,
    module_signs = module_signs,
    bootstrap = boot
  ))
}
