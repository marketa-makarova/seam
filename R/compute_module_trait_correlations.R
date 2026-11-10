#' Compute Pearson's correlations between clinical traits and sparse or classic eigengenes
#'
#' The original correlations are always returned. If spca_data contains
#' bootstrap replicates, bootstrap correlations are returned separately.
#' Each bootstrap replicate uses the same sampled rows for eigengenes and traits.
#'
#' @param trait_data List of numeric trait data frames, one per group, with rows
#'   in the same order as the corresponding expression matrix.
#' @param spca_data Result of spca_all_modules() or spca_all_modules_grid().
#' @param module_assignment Vector of modules assigned to variables. One module per variable.
#' @param mode One of "spca", "classic", or "both".
#' 
#' @return A list containing mtc_spca and/or mtc_classic (lists of original
#'   module-by-trait correlation matrices), trait_data, module_sizes,
#'   reduced_module_sizes, and bootstrap. The bootstrap element is NULL if
#'   no resampling was performed.
#' 
#' @examples
#' data("seam_demo", package = "seam")
#' fit = spca_all_modules(
#'  seam_demo$expression_data,
#'  seam_demo$module_assignment,
#'  lambda = 0.05, 
#'  bootstrap = FALSE
#' )
#' mtc = compute_module_trait_correlations(
#'  seam_demo$trait_data, 
#'  fit, seam_demo$module_assignment,
#'  mode = "both"
#' )
#' mtc$mtc_spca[[1]]
#' 
#' @export
#' 
compute_module_trait_correlations = function(
    trait_data,
    spca_data,
    module_assignment,
    mode = "spca"){

  mode = tolower(mode)
  if(!mode %in% c("classic", "spca", "both")){
    stop("Mode must be 'classic', 'spca', or 'both'.")
  }
  if(!is.list(spca_data) || !is.list(trait_data)){
    stop("spca_data and trait_data must be lists.")
  }
  n_sets = length(spca_data$PCs)
  if(length(trait_data) != n_sets){
    stop("Incompatible number of groups in trait_data and spca_data.")
  }

  # check and prepare trait data
  traits = lapply(seq_len(n_sets), function(set){
    x = as.matrix(trait_data[[set]])
    if(!is.numeric(x)) stop("Trait columns must be numeric.")
    if(nrow(x) != nrow(spca_data$PCs[[set]])) stop("Trait and eigengene sample counts do not match.")
    if(!is.null(rownames(x)) && !is.null(rownames(spca_data$PCs[[set]])) && !identical(rownames(x), rownames(spca_data$PCs[[set]]))){
      stop("Trait rows must be in the same order as expression samples.")
    }
    return(x)
  })

  # prepare necessary variables
  module_colors = spca_data$module_colors
  n_modules = length(module_colors)
  if(length(module_assignment) != 1 || !is.na(module_assignment)){
    if(!all(module_colors %in% module_assignment)){
      stop("module_assignment does not match the fitted modules.")
    }
  }
  module_sizes = as.list(spca_data$module_sizes)
  reduced_module_sizes = lapply(spca_data$vs, function(group){
    lapply(group, function(weights) sum(weights != 0))
  })
  use_spca = mode %in% c("spca", "both")
  use_classic = mode %in% c("classic", "both")
  mtc_spca = if(use_spca) vector("list", n_sets) else NULL
  mtc_classic = if(use_classic) vector("list", n_sets) else NULL

  # calculate classic eigengenes
  if(use_classic){
    classic_PCs = compute_eigengenes(
      module_matrices = spca_data$module_matrices,
      module_assignment = module_assignment,
      module_colors = module_colors
    )
  }

  # compute correlations
  for(set in seq_len(n_sets)){
    if(use_spca){
      mtc_spca[[set]] = cor(spca_data$PCs[[set]], traits[[set]], use = "pairwise.complete.obs")
      rownames(mtc_spca[[set]]) = module_colors
    }
    if(use_classic){
      mtc_classic[[set]] = cor(classic_PCs[[set]], traits[[set]], use = "pairwise.complete.obs")
      rownames(mtc_classic[[set]]) = module_colors
    }
  }

  # bootstrap correlations, if bootstrap is present
  boot = NULL
  if(!is.null(spca_data$bootstrap)){
    samples = spca_data$bootstrap$samples
    n_boot = nrow(samples[[1]])
    boot = list(
      mtc_spca = if(use_spca) vector("list", n_sets) else NULL,
      mtc_classic = if(use_classic) vector("list", n_sets) else NULL
    )

    for(set in seq_len(n_sets)){
      n_traits = ncol(traits[[set]])
      dims = list(module_colors, colnames(traits[[set]]), NULL)
      if(use_spca){
        boot$mtc_spca[[set]] = array(NA_real_, c(n_modules, n_traits, n_boot), dimnames = dims)
      }
      if(use_classic){
        boot$mtc_classic[[set]] = array(NA_real_, c(n_modules, n_traits, n_boot), dimnames = dims)
      }

      for(b in seq_len(n_boot)){
        rows = samples[[set]][b, ]
        sampled_traits = traits[[set]][rows, , drop = FALSE]

        if(use_spca){
          scores = matrix(spca_data$bootstrap$PCs[[set]][, , b], nrow = length(rows), ncol = n_modules)
          boot$mtc_spca[[set]][, , b] = cor(scores, sampled_traits, use = "pairwise.complete.obs")
        }

        if(use_classic){
          sampled_modules = lapply(spca_data$module_matrices[[set]], function(x){ x[rows, , drop = FALSE] })
          sampled_classic = compute_eigengenes(
            module_matrices = list(sampled_modules),
            module_assignment = module_assignment,
            module_colors = module_colors
          )[[1]]
          # orient resampled classic scores to the full-data eigengenes so they are comparable
          for(i in seq_len(n_modules)){
            ref = classic_PCs[[set]][rows, i]
            similarity = suppressWarnings(cor(sampled_classic[, i], ref))
            if(is.finite(similarity) && similarity < 0){
              sampled_classic[, i] = - sampled_classic[, i]
            }
          }
          boot$mtc_classic[[set]][, , b] = cor(sampled_classic, sampled_traits, use = "pairwise.complete.obs")
        }
      }
    }
  }

  result = list(
    trait_data = trait_data,
    module_sizes = module_sizes,
    reduced_module_sizes = reduced_module_sizes,
    bootstrap = boot
  )
  if(use_spca) result$mtc_spca = mtc_spca
  if(use_classic) result$mtc_classic = mtc_classic
  return(result)
}
