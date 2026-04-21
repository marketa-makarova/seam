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
      red = matrix(c(2, 1, 4, 3), nrow = 2, dimnames = list(c("id1", "id2"), c("gene3", "gene4")))
    )
  ),
  module_colors = c("blue", "red"), # two modules
  module_sizes = c(4, 4),
  var_expr = list(c(blue = 0.7, red = 0.6))
)

labels <- "A"

test_that("plot_module_heatmaps validates spca_data input", {
  expect_error(
    plot_module_heatmaps(
      spca_data = matrix(1:4, nrow = 2),
      labels = "G1"
    ),
    "Invalid spca_data argument"
  )
})

test_that("plot_module_heatmaps checks labels length", {

  expect_error(
    plot_module_heatmaps(
      spca_data = spca_data,
      labels = c("G1", "G2")
    ),
    "Length of labels does not coincide"
  )
})

test_that("plot_module_heatmaps rejects non-numeric threshold", {

  expect_error(
    plot_module_heatmaps(
      spca_data = spca_data,
      labels = "G1",
      threshold = "high"
    ),
    "Invalid threshold argument value"
  )
})

test_that("plot_module_heatmaps warns twice for threshold outside (0, 1)", {
  
  warnings <- capture_warnings(
    plot_module_heatmaps(
      spca_data = spca_data,
      labels = "G1",
      threshold = 1.5
    )
  )
  
  expect_length(warnings, 2)
  expect_true(all(grepl(
    "Threshold argument value should be between 0 and 1.",
    warnings
  )))
})

test_that("plot_module_heatmaps runs silently with valid input", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  
  expect_silent(
    plot_module_heatmaps(
      spca_data = spca_data,
      labels = "G1",
      return_information = TRUE
    )
  )
})

test_that("plot_module_heatmaps works with return_information = FALSE", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  
  expect_silent(
    plot_module_heatmaps(
      spca_data = spca_data,
      labels = "G1",
      return_information = FALSE
    )
  )
})

test_that("plot_module_heatmaps supports module selection", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  expect_silent(
    plot_module_heatmaps(
      spca_data = spca_data,
      labels = "G1",
      modules = "blue"
    )
  )
})
