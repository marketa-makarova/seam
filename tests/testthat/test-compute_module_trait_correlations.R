spca_data <- list( # one group
  PCs = list( # one group
    matrix( # two modules
      1:4,
      nrow = 2,
      dimnames = list(c("s1", "s2"), c("PC1", "PC2"))
    )
  ),
  vs = list(
    list(blue = c(gene1 = 1, gene2 = 0), red = c(gene1 = 1, gene2 = 1)) # weights for 2 modules
  ),
  module_matrices = list(
    list(
      blue = matrix(c(1, 2, 3, 4), nrow = 2, dimnames = list(c("id1", "id2"), c("gene1", "gene2"))),
      red = matrix(c(2, 1, 4, 3), nrow = 2, dimnames = list(c("id1", "id2"), c("gene1", "gene2")))
    )
  ),
  module_colors = c("blue", "red"), # two modules
  module_sizes = c(4, 4),
  var_expr = list(c(blue = 0.7, red = 0.6))
)

trait_data <- list(
    data.frame( # two traits
      trait1 = c(10, 20),
      trait2 = c(5, 15),
      row.names = c("s1", "s2")
    )
)

module_assignment = c("blue", "blue", "red", "red")

test_that("compute_module_trait_correlations works in spca mode", {
  
  result <- compute_module_trait_correlations(
    trait_data = trait_data,
    spca_data = spca_data,
    module_assignment = module_assignment,
    mode = "spca"
  )
  
  expect_type(result, "list")
  expect_named(result, c("", "trait_data", "module_sizes", "reduced_module_sizes"))
  
  mtc <- result[[1]][[1]]
  
  expect_true(is.matrix(mtc))
  expect_equal(dim(mtc), c(2, 2))  # modules / traits
  expect_equal(rownames(mtc), c("blue", "red"))
  expect_equal(colnames(mtc), c("trait1", "trait2"))
  
  expect_true(all(is.finite(mtc)))
})

test_that("compute_module_trait_correlations works in classic mode", {
  
  result <- compute_module_trait_correlations(
    trait_data = trait_data,
    spca_data = spca_data,
    module_assignment = module_assignment,
    mode = "classic"
  )
  
  mtc <- result[[1]][[1]]
  
  expect_true(is.matrix(mtc))
  expect_equal(dim(mtc), c(2, 2))
  expect_equal(rownames(mtc), c("blue", "red"))
  expect_equal(colnames(mtc), c("trait1", "trait2"))
})

test_that("compute_module_trait_correlations works in both mode", {
  
  result <- compute_module_trait_correlations(
    trait_data = trait_data,
    spca_data = spca_data,
    module_assignment = module_assignment,
    mode = "both"
  )
  
  expect_named(
    result,
    c("mtc_classic", "mtc_spca", "trait_data", "module_sizes", "reduced_module_sizes")
  )
  
  expect_true(is.matrix(result$mtc_classic[[1]]))
  expect_true(is.matrix(result$mtc_spca[[1]]))
  
  expect_equal(
    dim(result$mtc_classic[[1]]),
    dim(result$mtc_spca[[1]])
  )
})

