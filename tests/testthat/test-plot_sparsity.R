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

labels <- "A"

test_that("plot_sparsity errors on invalid spca_data", {
  expect_error(
    plot_sparsity(
      spca_data = "not_a_list",
      labels = labels
    ),
    "Invalid spca_data argument"
  )
})

test_that("plot_sparsity checks label length", {
  expect_error(
    plot_sparsity(
      spca_data = spca_data,
      labels = c("A", "B", "C")
    ),
    "Length of labels does not coincide"
  )
})

test_that("plot_sparsity runs silently with all modules", {
  expect_silent(
    plot_sparsity(
      spca_data = spca_data,
      labels = labels
    )
  )
})

test_that("plot_sparsity supports module selection by color", {
  expect_silent(
    plot_sparsity(
      spca_data = spca_data,
      labels = labels,
      modules = "blue"
    )
  )
})

test_that("plot_sparsity runs without legend and module size labels", {
  expect_silent(
    plot_sparsity(
      spca_data = spca_data,
      labels = labels,
      show_legend = FALSE,
      show_size = FALSE
    )
  )
})
