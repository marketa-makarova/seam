expr <- list(
  matrix(rnorm(20), nrow = 10),
  matrix(rnorm(20), nrow = 10)
)

module_assignment <- c("red", "blue")
lambda_grid <- c(0, 0.5)


test_that("plot_lambdas validates expression_data input", {
  expect_error(
    plot_lambdas(
      expression_data = matrix(1:4, nrow = 2),
      module_assignment = c("A", "B"),
      lambda_grid = c(0, 1),
      labels = "G1"
    ),
    "Invalid expression_data argument"
  )
})

test_that("plot_lambdas rejects negative lambda values", {
  expr <- list(
    matrix(rnorm(20), nrow = 10)
  )
  
  expect_error(
    plot_lambdas(
      expression_data = expr,
      module_assignment = rep("A", 2),
      lambda_grid = c(-0.1, 0.2),
      labels = "G1"
    ),
    "Select non-negative lambdas"
  )
})

test_that("plot_lambdas checks module_assignment length", {
  expr <- list(
    matrix(rnorm(30), ncol = 10)
  )
  
  expect_error(
    plot_lambdas(
      expression_data = expr,
      module_assignment = c("A", "B", "C"), # wrong length
      lambda_grid = c(0, 1),
      labels = "G1"
    ),
    "Number of variables in expression_data does not coincide with length of module_assignment."
  )
})

test_that("plot_lambdas checks labels length", {
  expr <- list(
    matrix(rnorm(20), nrow = 10),
    matrix(rnorm(20), nrow = 10)
  )
  
  expect_error(
    plot_lambdas(
      expression_data = expr,
      module_assignment = c("A", "B"),
      lambda_grid = c(0, 1),
      labels = "OnlyOneLabel"
    ),
    "Length of labels does not coincide"
  )
})

test_that("plot_lambdas runs successfully with valid input", {
  skip_if_not_installed("ggplot2")
  
  expect_silent(
    plot_lambdas(
      expression_data = expr,
      module_assignment = module_assignment,
      lambda_grid = lambda_grid,
      labels = c("G1", "G2")
    )
  )
})