expr <- list(matrix(rnorm(40), nrow = 10,
                     dimnames = list(NULL, paste0("V", 1:4)))
)

test_that("spca_all_modules rejects negative lambda", {
  expect_error(
    spca_all_modules(
      expression_data = expr,
      module_assignment = "A",
      lambda = -1
    ),
    "Select non-negative lambda"
  )
})

test_that("spca_all_modules validates expression_data", {
  expect_error(
    spca_all_modules(
      expression_data = matrix(1:4, nrow = 2),
      module_assignment = c("A", "B"),
      lambda = 0
    ),
    "Invalid expression_data argument"
  )
})

test_that("spca_all_modules checks module_assignment length", {
  
  expect_error(
    spca_all_modules(
      expression_data = expr,
      module_assignment = c("A", "B", "C"),
      lambda = 0
    ),
    "Number of variables in expression_data does not coincide"
  )
})

test_that("spca_all_modules runs with multiple modules", {
  skip_if_not_installed("elasticnet")
  
  res <- spca_all_modules(
    expression_data = expr,
    module_assignment = c("red", "red", "blue", "blue"),
    lambda = 0.1
  )
  
  expect_type(res, "list")
  expect_named(res, c(
    "PCs", "var_exp", "vs", "module_assignment",
    "module_matrices", "module_colors", "module_sizes", "module_signs"
  ))
  
  expect_equal(res$module_colors, c("red", "blue"))
  expect_equal(res$module_sizes, c(red = 2, blue = 2))
  
  expect_equal(ncol(res$PCs[[1]]), 2)
  expect_equal(length(res$vs[[1]]), 2)
})

test_that("spca_all_modules works with NA module_assignment", {
  skip_if_not_installed("elasticnet")
  
  res <- spca_all_modules(
    expression_data = expr,
    module_assignment = NA,
    lambda = 0.2
  )
  
  expect_equal(res$module_colors, NA)
  
  expect_equal(ncol(res$PCs[[1]]), 1)
  expect_length(res$vs[[1]], 1)
})