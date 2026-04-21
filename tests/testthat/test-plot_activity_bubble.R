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

test_that("plot_activity_bubble errors on invalid spca_data", {
  expect_error(
    plot_activity_bubble(
      spca_data = "not a list",
      labels = labels
    ),
    "Invalid spca_data argument"
  )
})

test_that("plot_activity_bubble errors on non-numeric threshold", {
  expect_error(
    plot_activity_bubble(
      spca_data = spca_data,
      labels = labels,
      threshold = "high"
    ),
    "Invalid threshold argument"
  )
})


test_that("plot_activity_bubble warns on threshold outside (0,1)", {
  expect_warning(
    expect_warning(
      plot_activity_bubble(
        spca_data = spca_data,
        labels = labels,
        threshold = 1.5
      ),
      "Threshold argument value should be between 0 and 1"
    ),
    "Threshold argument value should be between 0 and 1"
  )
})

test_that("plot_activity_bubble errors when label length does not match groups", {
  expect_error(
    plot_activity_bubble(
      spca_data = spca_data,
      labels = c("A", "B")
    ),
    "Length of labels does not coincide"
  )
})

test_that("plot_activity_bubble runs silently with valid input", {
  expect_silent(
    plot_activity_bubble(
      spca_data = spca_data,
      labels = labels
    )
  )
})

test_that("plot_activity_bubble works with return_information = FALSE", {
  expect_silent(
    plot_activity_bubble(
      spca_data = spca_data,
      labels = labels,
      return_information = FALSE
    )
  )
})
