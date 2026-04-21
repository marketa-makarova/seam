compute_eigengenes = function(
    module_matrices, 
    module_assignment, 
    module_colors){
  
  multiexpression = lapply(module_matrices, function(m) do.call(cbind, m))
  
  nsets = length(multiexpression)
  MEs = list()
  
  for(set in 1:nsets){
    MEs[[set]] = matrix(NA, nrow = nrow(multiexpression[[set]]), ncol = length(module_colors))
  }
  
  if((sum(is.na(module_assignment)) != 0) | (sum(is.na(module_colors)) != 0)){
    for(set in 1:nsets){
      expr = multiexpression[[set]]
      ME = svd(t(scale(expr)))$v[, 1]
      MEs[[set]][, 1] = if(cor(ME, apply(expr, 1, sum) / ncol(expr)) > 0){ ME } else { - ME}
      }
    } else {
    for(set in 1:nsets){
      for(i in 1:length(module_colors)){
        module_expr = module_matrices[[set]][[i]]
        ME = svd(t(scale(module_expr)))$v[, 1]
        MEs[[set]][, i] = if(cor(ME, apply(module_expr, 1, sum) / ncol(module_expr)) > 0){ ME } else { - ME}
      }
    }
  }

  return(MEs)
}