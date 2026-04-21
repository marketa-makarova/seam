mat <- matrix(c(0.3, 0.6, 0.8, 0.2), nrow = 2)
rownames(mat) <- c("red", "blue")
colnames(mat) <- c("Trait1", "Trait2")

module_trait_correlations <- list(
  list(mat),
  list(mat),
  mtc_spca = list(mat),
  module_sizes = c(10, 12),
  reduced_module_sizes = list(c(8, 9))
)

test_that("plot_module_trait_correlations rejects invalid mode", {
  expect_error(
    plot_module_trait_correlations(
      module_trait_correlations = list(),
      labels = "G1",
      mode = "invalid"
    ),
    "Type can be either classic, spca, or both"
  )
})

test_that("plot_module_trait_correlations warns for invalid mfrow length", {
  
  expect_error(
    plot_module_trait_correlations(
      module_trait_correlations = module_trait_correlations,
      labels = "G1",
      mfrow = c(1, 2, 3)
    ),
    "Argument mfrow should be a vector of length 2"
  )
})

test_that("plot_module_trait_correlations warns on labels length mismatch", {
  
  expect_error(
    plot_module_trait_correlations(
      module_trait_correlations = module_trait_correlations,
      labels = c("G1", "G2")
    ),
    "Length of labels does not coincide"
  )
})

test_that("plot_module_trait_correlations rejects single module selection", {

  expect_error(
    plot_module_trait_correlations(
      module_trait_correlations = module_trait_correlations,
      labels = "G1",
      modules = "red"
    ),
    "Select at least two modules"
  )
})

test_that("plot_module_trait_correlations runs silently in both mode", {

  expect_silent(
    plot_module_trait_correlations(
      module_trait_correlations = module_trait_correlations,
      labels = "G1",
      mode = "both",
      plot_avg_differences = TRUE
    )
  )
})

test_that("plot_module_trait_correlations runs without RMSE when plot_avg_differences = FALSE", {
  
  expect_silent(
    plot_module_trait_correlations(
      module_trait_correlations = module_trait_correlations,
      labels = "G1",
      mode = "both",
      plot_avg_differences = FALSE
    )
  )
})
