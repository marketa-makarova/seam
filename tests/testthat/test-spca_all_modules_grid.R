expr <- list(matrix(rnorm(40), nrow = 10,
                     dimnames = list(NULL, paste0("V", 1:4)))
)

test_that("spca_all_modules_grid rejects negative lambdas", {
  expect_error(
    spca_all_modules_grid(
      expression_data = expr,
      module_assignment = "A",
      lambdas = -1
    ),
    "Select non-negative lambda"
  )
})

test_that("spca_all_modules_grid validates expression_data", {
  expect_error(
    spca_all_modules_grid(
      expression_data = matrix(1:4, nrow = 2),
      module_assignment = c("A", "B"),
      lambdas = c(0, 0)
    ),
    "Invalid expression_data argument"
  )
})

test_that("spca_all_modules_grid checks module_assignment length", {
  
  expect_error(
    spca_all_modules_grid(
      expression_data = expr,
      module_assignment = c("A", "B", "C"),
      lambdas = c(0, 0, 0)
    ),
    "Number of variables in expression_data does not coincide"
  )
})

test_that("spca_all_modules_grid runs with multiple modules and lambdas", {
  skip_if_not_installed("elasticnet")
  
  module_assignment <- c("red", "red", "blue", "blue")
  lambdas <- c(red = 0.1, blue = 0.5)
  
  res <- spca_all_modules_grid(
    expression_data = expr,
    module_assignment = module_assignment,
    lambdas = lambdas
  )
  
  expect_equal(res$module_colors, c("red", "blue"))
  expect_equal(res$module_sizes, c(red = 2, blue = 2))
  expect_equal(res$lambdas, lambdas)
  
  expect_equal(ncol(res$PCs[[1]]), 2)
  expect_equal(length(res$vs[[1]]), 2)
})

test_that("spca_all_modules_grid works with NA module_assignment", {
  skip_if_not_installed("elasticnet")
  
  res <- spca_all_modules_grid(
    expression_data = expr,
    module_assignment = NA,
    lambdas = 0.3
  )
  
  expect_equal(res$module_colors, NA)
  expect_equal(res$lambdas, 0.3)
  
  expect_equal(ncol(res$PCs[[1]]), 1)
  expect_length(res$vs[[1]], 1)
})