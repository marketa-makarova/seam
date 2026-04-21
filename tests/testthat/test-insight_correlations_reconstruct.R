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


trait_data <- list( # one group
    data.frame( # two traits
      trait1 = c(10, 7),
      trait2 = c(5, 4),
      row.names = c("s1", "s2")
  )
)

module_trait_correlations <- list(
  mtc_spca = list(matrix(c(0.1, 0.1), nrow = 2)),
  trait_data = list(
    data.frame( # two traits
      trait1 = c(10, 7),
      trait2 = c(5, 4),
      row.names = c("s1", "s2")
    ))
)

test_that("insight_correlations_reconstruct returns message when no correlations exceed threshold", {
  
  result <- insight_correlations_reconstruct(
    module_trait_correlations = module_trait_correlations,
    spca_data = spca_data,
    interest_threshold = 0.99,  # too high interest threshold 
    labels = "A",
    col = c("red", "blue")
  )
  
  expect_type(result, "character")
  expect_match(result, "No interesting correlations")
})

test_that("insight_correlations_reconstruct validates interest_threshold", {

  expect_error(
    insight_correlations_reconstruct(
      module_trait_correlations = module_trait_correlations,
      spca_data = spca_data,
      interest_threshold = "high",
      labels = c("A", "B"),
      col = c("red", "blue")
    ),
    "Invalid interest_threshold"
  )
})

test_that("insight_correlations_reconstruct checks label length", {
 
  expect_error(
    insight_correlations_reconstruct(
      module_trait_correlations = module_trait_correlations,
      spca_data = spca_data,
      labels = c("A", "B", "C"),  # wrong length
      col = c("red", "blue")
    ),
    "Length of labels does not coincide"
  )
})

test_that("insight_correlations_reconstruct produces ggplot output", {
  
  expect_silent({
    insight_correlations_reconstruct(
      module_trait_correlations = module_trait_correlations,
      spca_data = spca_data,
      interest_threshold = 0.2,
      top_n_upper = 3, 
      top_n_lower = 2, 
      labels = "A",
      col = c("red", "blue")
    )
  })
})
